module SinglePendulum

using DifferentialEquations, GLMakie, Elliptic, Logging
using LinearAlgebra, Dates
using RigidBodyDynamics
using StaticArrays
using Rotations
using Colors
using GLMakie.Makie: Mouse

# Global reference to the current screen for GUI management
const CURRENT_SCREEN = Ref{Union{Nothing, GLMakie.Screen}}(nothing)
const PHASE_MAX_OMEGA = 3π/2  #
const control_label_fontsize = 16
# Struct to hold all time series data for atomic observable updates
struct TimeSeriesState
    ts::Vector{Float32}
    thetas::Vector{Float32}
    omegas::Vector{Float32}
    thetas_sho::Vector{Float32}
    omegas_sho::Vector{Float32}
end

# Mutable wrapper for mechanism components (allows replacement in closures)
mutable struct MechanismWrapper
    mech::Any
    state::Any
    result::Any
    joint::Any
end

# Mutable struct to expose all GUI controls for programmatic access
mutable struct GUIControls
    # Figure and screen
    fig::Any
    screen::Any

    # Buttons
    run_btn::Any
    rand_btn::Any
    reset_btn::Any
    clear_phase_btn::Any
    # normalize_toggle::Any

    # Sliders
    damp_sld::Any
    L_sld::Any
    m_sld::Any
    speed_sld::Any

    # Observables
    is_running::Any
    is_normalized::Any
    L_obs::Any
    θ_obs::Any
    ω_obs::Any
    t_now::Any
    ke_proportion::Any
    pe_proportion::Any

    # Functions
    set_ic!::Any
    rebuild_mechanism!::Any
    set_phase_ic!::Any
end

# Global reference to GUI controls for programmatic access (defined after GUIControls struct)
const CURRENT_GUI = Ref{Union{Nothing, GUIControls}}(nothing)

function run_gui()

    # ==========================
    # Logging setup
    # ==========================
    global screenshot_dir = abspath(joinpath(@__DIR__, "..", "screenshots"))
    mkpath(screenshot_dir)
    log_dir = abspath(joinpath(@__DIR__, "..", "logs"))
    mkpath(log_dir)
    log_file = joinpath(log_dir, "pendulum_debug.log")

    # Create file logger that appends to log file
    file_logger = SimpleLogger(open(log_file, "a"))

    # Set global logger - debug messages go to file only
    global_logger(file_logger)

    @info "========== Pendulum Simulation Started =========="
    @info "Log file: $log_file"
    @info "Timestamp: $(now())"

    # ==========================
    # Physical model (RigidBodyDynamics)
    # ==========================
    g = 9.81
    ρ_rod = 0.5  # Linear density of rod (kg/m) - allows exploring simple vs physical pendulum regimes

    # Function to create mechanism with given length and bob mass
    function create_mechanism(L::Float64, m_bob::Float64)
        frame = CartesianFrame3D("link")

        # Rod mass scales with length
        m_rod = ρ_rod * L

        # Rod inertia (uniform rod about one end)
        Irod_com = Diagonal(SVector((1 / 12) * m_rod * L^2, 0.0, (1 / 12) * m_rod * L^2))
        I_rod = SpatialInertia(frame;
            mass=m_rod,
            com=SVector(0.0, -L / 2, 0.0),  # Center of mass at midpoint
            moment_about_com=Matrix(Irod_com))

        # Bob as point mass at end
        I_bob = SpatialInertia(frame;
            mass=m_bob,
            com=SVector(0.0, -L, 0.0),  # At the end of the rod
            moment_about_com=zeros(3, 3))

        # Combined inertia
        I_link = I_rod + I_bob
        link = RigidBody("link", I_link)

        mech = Mechanism(RigidBody{Float64}("world"); gravity=SVector(0.0, -g, 0.0))
        ground = root_body(mech)

        axis = SVector(0.0, 0.0, 1.0)
        joint = Joint("hinge", Revolute{Float64}(axis))

        Reye = one(RotMatrix{3,Float64})
        pzero = SVector(0.0, 0.0, 0.0)

        attach!(mech, ground, link, joint;
            joint_pose=Transform3D(frame_before(joint), default_frame(ground), Reye, pzero),
            successor_pose=Transform3D(default_frame(link), frame_after(joint), Reye, pzero)
        )

        return mech, MechanismState(mech), DynamicsResult(mech), joint
    end

    # Initial rod length and bob mass
    L_obs = Observable(1.0)
    m_bob_initial = 0.8  # Initial bob mass in kg

    # Create initial mechanism (wrapped in mutable struct for mutability in closures)
    mech_init, state_init, result_init, joint_init = create_mechanism(L_obs[], m_bob_initial)
    mw = MechanismWrapper(mech_init, state_init, result_init, joint_init)

    # ==========================
    # Visualization (Makie)
    # ==========================
    GLMakie.activate!(;float=false, focus_on_show=true)
    set_theme!(theme_dark())

    fig = Figure(size=(1000, 1000), figure_padding=20)

    # Phase space plot (LEFT, 70% width) - 1:1 aspect ratio
    newAxPhase() = begin
        ax = Axis(fig[1, 1],
            title="Phase space (θ, ω)",
            xlabel="θ [rad]",
            ylabel="ω/ω₀",
            aspect=1,  # 1:1 aspect ratio for circular SHO orbits
            xticks=([-π, -π/2, 0, π/2, π], ["-π", "-π/2", "0", "π/2", "π"]),
            yticks=([-3π, -2π, -π, 0, π, 2π, 3π], ["-3π", "-2π", "-π", "0", "π", "2π", "3π"]))
        # Disable default mouse interactions (scaling, panning, etc.)
        deregister_interacis_oscillating[]tion!(ax, :rectanglezoom)
        ax
    end
    axPhase = Observable(newAxPhase())
    xlims!(axPhase[], -π, π)  # Exactly -π to π horizontally
    ylims!(axPhase[], -PHASE_MAX_OMEGA, PHASE_MAX_OMEGA)  # Larger limits for rotating orbits

    # Observable for normalized coordinates (defined early for use in separatrix)
    is_normalized = Observable(true)  # Start with normalized phase space

    # Plot separatrix (curve separating oscillation from rotation)
    # For a physical pendulum: E = ½Iω² + mgL_cm(1 - cos(θ))
    # At separatrix: E_crit = 2mgL_cm (energy at θ=π, ω=0)
    # Therefore: ω² = 2mgL_cm(1 + cos(θ))/I
    # Need to compute I and L_cm from the mechanism
    θ_sep = range(-π, π, length=200)
    
    # Function to compute separatrix given L and m_bob
    function compute_separatrix(L::Float64, m_bob::Float64)
        m_rod = ρ_rod * L
        
        # Moment of inertia about pivot (parallel axis theorem)
        # Rod: I_rod = (1/3)m_rod*L² (uniform rod about end)
        I_rod = (1/3) * m_rod * L^2
        # Bob: I_bob = m_bob*L² (point mass at distance L)
        I_bob = m_bob * L^2
        I_total = I_rod + I_bob
        
        # Center of mass distance from pivot
        L_cm = (m_rod * L/2 + m_bob * L) / (m_rod + m_bob)
        
        # Total mass
        m_total = m_rod + m_bob
        
        # Separatrix: ω² = 2*m_total*g*L_cm*(1 + cos(θ))/I_total
        ω_sep_upper = [sqrt(max(2 * m_total * g * L_cm * (1 + cos(θ)) / I_total, 0)) for θ in θ_sep]
        ω_sep_lower = [-sqrt(max(2 * m_total * g * L_cm * (1 + cos(θ)) / I_total, 0)) for θ in θ_sep]
        
        # Normalized version (divide by ω₀ = √(m_total*g*L_cm/I_total))
        ω₀ = sqrt(m_total * g * L_cm / I_total)
        ω_sep_norm_upper = ω_sep_upper ./ ω₀
        ω_sep_norm_lower = ω_sep_lower ./ ω₀
        
        return (ω_sep_upper, ω_sep_lower, ω_sep_norm_upper, ω_sep_norm_lower)
    end
    
    # Compute initial separatrix (use same initial mass as mechanism)
    L_initial = L_obs[]
    m_initial = m_bob_initial  # Use the actual initial bob mass
    sep_initial = compute_separatrix(L_initial, m_initial)
    
    # Make observables that update when L or m change
    ω_sep_physical_upper = Observable(sep_initial[1])
    ω_sep_physical_lower = Observable(sep_initial[2])
    ω_sep_normalized_upper = Observable(sep_initial[3])
    ω_sep_normalized_lower = Observable(sep_initial[4])

    # Plot separatrix (will update based on is_normalized and L_obs)
    lines!(axPhase[], θ_sep, @lift($is_normalized ? ω_sep_normalized_upper : $ω_sep_physical_upper),
           color=(:white, 0.3), linestyle=:dash, linewidth=2)
    lines!(axPhase[], θ_sep, @lift($is_normalized ? ω_sep_normalized_lower : $ω_sep_physical_lower),
           color=(:white, 0.3), linestyle=:dash, linewidth=2)

    # Time series plots (left side only, narrower)
    axTheta = Axis(fig[2, 1],
        title="θ(t)",
        xlabel="",
        ylabel="θ [rad]",
        xticklabelsvisible=false)  # Hide x-tick labels on top plot

    axOmega = Axis(fig[3, 1],
        title="ω(t)",
        xlabel="t [s]",
        ylabel="ω [rad/s]")

    # Link x-axes so they zoom/pan together
    linkxaxes!(axTheta, axOmega)

    # Calculate period for initial condition and set xlim to show 2 periods
    # θ0_initial = deg2rad(120.0)
    # T_sho = 2π * sqrt(L_obs[] / g)  # Small angle approximation
    # T_exact = T_sho * (2 / π) * Elliptic.K(sin(θ0_initial / 2)^2)  # Large angle period
    # two_periods = 2 * T_exact
    # set_ic!(50.0,0.0)

    # Set x-ticks at even divisions of the period (0, T/2, T, 3T/2, 2T)
    # xtick_values = [0, T_exact/2, T_exact, 3*T_exact/2, 2*T_exact]
    # xtick_labels = ["0", "T/2", "T", "3T/2", "2T"]
    # axTheta.xticks = (xtick_values, xtick_labels)
    # axOmega.xticks = (xtick_values, xtick_labels)

    # Shared x-axis limits (set once, don't update dynamically)
    # xlims!(axTheta, 0, two_periods)
    # xlims!(axOmega, 0, two_periods)
    # ylims!(axTheta, -π-π/8, π+π/8)  # Initial y-limits with margin

    # Enable y-axis autoscaling to data range
    # limits!(axTheta, 0, two_periods, nothing, nothing)  # Auto y-limits
    # limits!(axOmega, 0, two_periods, nothing, nothing)  # Auto y-limits

    # Right column top: Subgrid for pendulum + energy bars
    pend_energy_grid = GridLayout(fig[1, 2], tellwidth=true, tellheight=false, alignmode = Inside())

    # Pendulum visualization
    axPend = Axis(pend_energy_grid[1, 1:5],
        title="Pendulum", aspect=1
    )
        # width=Auto(),
        # height=Auto(),
        # alignmode=Mixed(left=0, right=0, bottom=0, top=0)
    # )
    # Pendulum limits scale with rod length (1.3x for margin)
    pend_lim = 2.0 * L_obs[]  # Reduced from 2.3x for tighter view
    xlims!(axPend, -pend_lim * 1.2, pend_lim * 1.2)  # Moderate horizontal margin
    ylims!(axPend, -pend_lim * 1.2, pend_lim * 1.2)  # Match vertical to horizontal

    # Update pendulum limits when L changes
    on(L_obs) do L
        lim = 1.3 * L
        # xlims!(axPend, -lim, lim)
        # ylims!(axPend, -lim, lim)
    end

    # hidedecorations!(axPend, grid=false, ticks=false)
    # Energy display as proportion bars (below pendulum)
    ke_proportion = Observable(0.5)  # Fraction of total energy
    pe_proportion = Observable(0.5)

    # KE bar (blue)
    Label(pend_energy_grid[2, 1], "KE", halign=:left, fontsize=control_label_fontsize, font=:bold, padding=(5, 0, 0, 0))
    ax_ke = Axis(pend_energy_grid[3, 1:5], height=20, limits=(0,1,0,1))
    hidedecorations!(ax_ke)
    hidespines!(ax_ke)
    ke_rect = @lift(Rect(0, 0, $ke_proportion,1))
    poly!(ax_ke, ke_rect, color=:dodgerblue, strokewidth=1, strokecolor=:white)

    # PE bar (green)
    Label(pend_energy_grid[4, 1], "PE", halign=:left, fontsize=control_label_fontsize, font=:bold, padding=(5, 0, 0, 0))
    ax_pe = Axis(pend_energy_grid[5, 1:5], height=20, limits=(0,1,0,1))
    hidedecorations!(ax_pe)
    hidespines!(ax_pe)
    pe_rect = @lift(Rect(0, 0, $pe_proportion, 1))
    poly!(ax_pe, pe_rect, color=:green, strokewidth=1, strokecolor=:white)

    # Spacing within pendulum+energy grid
    rowgap!(pend_energy_grid, 2)
    # rowsize!(pend_energy_grid, 1, Auto(true))  # Pendulum gets remaining space
    # colsize!(pend_energy_grid, 1, Auto(true))  # Column expands to fill width
    trim!(pend_energy_grid)
    # Right column bottom: Control panel
    controls_grid = GridLayout(fig[2:3, 2], tellwidth=false, tellheight=false)
    
    # Buttons
    button_grid = GridLayout(controls_grid[1, 1:2], tellwidth=false, tellheight=false)
    run_btn = Button(button_grid[1, 1], label="Run", width=60)
    rand_btn = Button(button_grid[1,2], label="Random", width=60)
    reset_btn = Button(button_grid[1,3], label="Reset", width=60)
    clear_phase_btn = Button(button_grid[1,4], label="Clear", width=60)

    # Normalize toggle with label
    toggle_layout = GridLayout(controls_grid[2, :], halign=:center, tellwidth=false, tellheight=false)
    normalize_toggle = Toggle(toggle_layout[1, 1], active=false)
    Label(toggle_layout[1, 2], "ω/ω₀", halign=:center, fontsize=control_label_fontsize)
    is_running = Observable(false)  # Start paused
    row = 2
    # Sliders with labels (right side of buttons)
    function make_control(label, ctrl_fn)
        row += 1
        Label(controls_grid[row, 1], label, halign=:right, fontsize=control_label_fontsize)
        ctrl_fn(row)
    end

    damp_sld = make_control("damping", row -> begin
        Slider(controls_grid[row,2], range=0:0.01:0.35, startvalue=0.0)
    end)
    L_sld = make_control("L", row -> begin
        Slider(controls_grid[3, 2], range=0.2:0.1:2.0, startvalue=1.0)
    end)
    m_sld = make_control("m", row -> begin
        Slider(controls_grid[4, 2], range=0.1:0.1:2.0, startvalue=0.8)
    end)
    speed_sld = make_control("speed", row -> begin
        Slider(controls_grid[5, 2], range=0.01:0.1:2.0, startvalue=0.5)
    end)
    # Spacing within controls
    rowgap!(controls_grid, 15)
    colgap!(controls_grid, 5)

    # Add spacing between plots and controls
    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 10)

    # Row sizes - emphasize top row with phase plot and pendulum
    rowsize!(fig.layout, 1, Relative(0.70))  # Top row: phase + pendulum (70%)
    rowsize!(fig.layout, 2, Relative(0.15))  # Time series + controls top (15%)
    rowsize!(fig.layout, 3, Relative(0.15))  # Time series + controls bottom (15%)

    # Column sizes - phase plot gets 60%, pendulum gets 40%
    colsize!(fig.layout, 1, Relative(0.60))  # Phase plot (left)
    colsize!(fig.layout, 2, Relative(0.40))  # Pendulum (right)

    # Display and set window to floating mode at upper left
    screen = display(GLMakie.Screen(), fig)

    # Small delay to ensure window is created before positioning
    sleep(0.1)

    # Position window at upper left
    GLMakie.GLFW.SetWindowPos(screen.glscreen, 50, 50)  # Position at (50, 50) pixels from top-left

    # Store global reference for GUI management
    CURRENT_SCREEN[] = screen

    # ==========================
    # Screenshot task (DISABLED - was causing input focus stealing)
    # Use take_screenshot!() to manually capture screenshots
    # ==========================
    # @async begin
    #     while isopen(fig.scene)
    #         try
    #             # Create timestamped filename
    #             timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    #             screenshot_path = joinpath(screenshot_dir, "pendulum_$timestamp.png")

    #             # Save screenshot
    #             save(screenshot_path, fig)

    #             # Clean up old screenshots (keep only last 50)
    #             all_screenshots = filter(f -> endswith(f, ".png") && startswith(f, "pendulum_"), readdir(screenshot_dir))
    #             if length(all_screenshots) > 50
    #                 # Sort by filename (which includes timestamp, so chronological)
    #                 sort!(all_screenshots)
    #                 # Delete oldest ones
    #                 num_to_delete = length(all_screenshots) - 50
    #                 for i in 1:num_to_delete
    #                     old_file = joinpath(screenshot_dir, all_screenshots[i])
    #                     rm(old_file)
    #                 end
    #             end
    #         catch e
    #             @warn "Screenshot failed: $e"
    #         end
    #         sleep(10.0)  # Save every 10 seconds
    #     end
    # end

    # ==========================
    # Observables / data buffers
    # ==========================
    θ_obs = Observable(0.0)
    ω_obs = Observable(0.0)
    t_now = Observable(0.0)

    rod_pts = @lift(Point2f[(0, 0), (Float32($L_obs * sin($θ_obs)), Float32(-$L_obs * cos($θ_obs)))])
    bob_pt = @lift(Point2f(Float32($L_obs * sin($θ_obs)), Float32(-$L_obs * cos($θ_obs))))
    lines!(axPend, rod_pts, color=:gray80, linewidth=5)
    scatter!(axPend, bob_pt, color=:orange, markersize=25)

    # Mutable buffers for time series data
    ts = Float32[]
    thetas = Float32[]
    omegas = Float32[]

    # Single observable holding all time series data (for atomic updates)
    timeseries_state = Observable(TimeSeriesState(Float32[], Float32[], Float32[], Float32[], Float32[]))

    # SHO initial conditions
    θ0_sho = Observable(0.0)
    ω0_sho = Observable(0.0)

    # Flag to track if trajectory is oscillating or rotating
    is_oscillating = Observable(true)  # true if inside separatrix, false if rotating

    # Current position markers on time series (wrap within 2-period window)
    current_t = Observable([0.0f0])
    current_theta = Observable([0.0f0])
    current_omega = Observable([0.0f0])

    # Plot actual theta solution
    lines!(axTheta, @lift($timeseries_state.ts), @lift($timeseries_state.thetas), color=:cyan, linewidth=3, label="Actual")

    # Plot SHO theta solution (semi-transparent) - only visible for oscillating motion
    sho_theta_line = lines!(axTheta, @lift($timeseries_state.ts), @lift($timeseries_state.thetas_sho),
                            color=(:dodgerblue, 0.3), linewidth=2, label="SHO", linestyle=:dash,
                            visible=is_oscillating)

    # Fill between theta plots - only visible for oscillating motion
    theta_band = band!(axTheta, @lift($timeseries_state.ts), @lift($timeseries_state.thetas),
                       @lift($timeseries_state.thetas_sho), color=(:skyblue, 0.15),
                       visible=is_oscillating)

    axislegend(axTheta, position=:rt)

    # Plot actual omega solution
    lines!(axOmega, @lift($timeseries_state.ts), @lift($timeseries_state.omegas), color=:cyan, linewidth=3, label="Actual")

    # Plot SHO omega solution (semi-transparent) - only visible for oscillating motion
    sho_omega_line = lines!(axOmega, @lift($timeseries_state.ts), @lift($timeseries_state.omegas_sho),
                            color=(:dodgerblue, 0.3), linewidth=2, label="SHO", linestyle=:dash,
                            visible=is_oscillating)

    # Fill between omega plots - only visible for oscillating motion
    omega_band = band!(axOmega, @lift($timeseries_state.ts), @lift($timeseries_state.omegas),
                       @lift($timeseries_state.omegas_sho), color=(:skyblue, 0.15),
                       visible=is_oscillating)

    axislegend(axOmega, position=:rt)

    # Add current position markers (orange scatter points that wrap around)
    # scatter!(axTheta, current_t, current_theta, color=:orange, markersize=5, marker=:circle, strokewidth=2, strokecolor=:black)
    # scatter!(axOmega, current_t, current_omega, color=:orange, markersize=5, marker=:circle, strokewidth=2, strokecolor=:black)
    # --- Phase traces ---
    palette = [RGB(0.90, 0.30, 0.25), RGB(0.20, 0.70, 0.90), RGB(0.30, 0.85, 0.40),
        RGB(0.55, 0.45, 0.85), RGB(0.95, 0.70, 0.20), RGB(0.40, 0.80, 0.70)]

    traj_idx = Ref(0)
    θ_phase_cur = Ref(Observable(Float32[]))
    ω_phase_cur = Ref(Observable(Float32[]))

    # Arrow showing current position - simple arrowhead following the orbit
    arrow_θ = Observable([0.0f0])
    arrow_ω = Observable([0.0f0])
    arrow_u = Observable([0.1f0])  # Small direction vector
    arrow_v = Observable([0.1f0])

    # Just an arrowhead (no visible shaft) - magenta to stand out, much larger
    arrows2d!(axPhase[], arrow_θ, arrow_ω, arrow_u, arrow_v,
              lengthscale=5.0, color=:magenta,
              shaftwidth=0, tipwidth=10, tiplength=20)

    function new_phase_trace!(axPhase, palette, traj_idx, θref::Ref{Observable{Vector{Float32}}},
        ωref::Ref{Observable{Vector{Float32}}})
        traj_idx[] += 1
        θref[] = Observable(Float32[])
        ωref[] = Observable(Float32[])
        col = palette[mod1(traj_idx[], length(palette))]
        lines!(axPhase[], θref[], ωref[], color=col, linewidth=2)
    end

    function clear_time_plots!(ts, thetas, omegas, timeseries_state, t_now)
        empty!(ts)
        empty!(thetas)
        empty!(omegas)
        timeseries_state[] = TimeSeriesState(Float32[], Float32[], Float32[], Float32[], Float32[])
        t_now[] = 0.0
    end

    # Observable for current period (updates when IC changes)
    current_period = Observable(0.0)

    # Track previous angle and unwrapped angle for continuous time series
    θ_prev = Ref(0.0)
    θ_unwrapped = Ref(0.0)

    # ==========================
    # IC handling
    # ==========================
    function set_ic!(θ0::Float64, ω0::Float64)
        set_configuration!(mw.state, mw.joint, θ0)
        set_velocity!(mw.state, mw.joint, ω0)
        θ_obs[] = θ0
        ω_obs[] = ω0

        # Initialize unwrapped angle tracking
        θ_prev[] = θ0
        θ_unwrapped[] = 0.0

        # Initialize phase space arrow (normalized by ω₀ = √(g/L))
        θ_wrapped_init = mod(θ0 + π, 2π) - π
        ω₀ = sqrt(g / L_obs[])
        arrow_θ[] = [Float32(θ_wrapped_init)]
        arrow_ω[] = [Float32(is_normalized[] ? ω0 / ω₀ : ω0)]

        # Normalize initial direction
        ω0_phase = is_normalized[] ? ω0 / ω₀ : ω0
        init_dir_mag = sqrt(Float32(ω0_phase)^2) + Float32(1e-6)
        arrow_u[] = [Float32(ω0_phase) / init_dir_mag]
        arrow_v[] = [Float32(0.0)]  # Initial acceleration is zero

        clear_time_plots!(ts, thetas, omegas, timeseries_state, t_now)
        new_phase_trace!(axPhase, palette, traj_idx, θ_phase_cur, ω_phase_cur)

        # Initialize current position markers at t=0
        current_t[] = [0.0f0]
        current_theta[] = [Float32(θ0)]
        current_omega[] = [Float32(ω0)]

        # Calculate period for this IC
        T_sho_local = 2π * sqrt(L_obs[] / g)
        T_exact_local = T_sho_local * (2 / π) * Elliptic.K(sin(abs(θ0) / 2)^2)
        current_period[] = T_exact_local

        # Calculate maximum amplitudes from energy conservation
        # Total energy (normalized by mgL): E = (1/2)*(L/g)*ω^2 + (1 - cos(θ))
        E_normalized = 0.5 * (L_obs[] / g) * ω0^2 + (1.0 - cos(θ0))

        # Check if trajectory is oscillating (E < 2) or rotating (E >= 2)
        # The separatrix energy is E = 2 (the energy at the unstable equilibrium θ=π, ω=0)
        is_oscillating[] = E_normalized < 2.0

        # Maximum angle (when ω=0): E = 1 - cos(θ_max)
        cos_θ_max = 1.0 - E_normalized
        θ_max = if cos_θ_max >= -1.0
            acos(max(-1.0, cos_θ_max))  # Clamp to valid range
        else
            π  # Over-the-top motion
        end

        # Maximum velocity (when θ=0): E = (1/2)*(L/g)*ω_max^2
        ω_max = sqrt(2.0 * g * E_normalized / L_obs[])

        # SHO maximum amplitudes with arbitrary initial conditions
        # Amplitude: A = sqrt(θ0^2 + (ω0/ω_sho)^2)
        ω_sho = 2π / T_sho_local
        θ_max_sho = sqrt(θ0^2 + (ω0 / ω_sho)^2)
        ω_max_sho = sqrt(θ0^2 * ω_sho^2 + ω0^2)

        # Determine y-limits based on whether trajectory is oscillating or rotating
        if is_oscillating[]
            # Oscillating: use bounded limits based on energy
            θ_limit = max(θ_max, θ_max_sho) * 1.1  # 10% margin
            ω_limit = max(ω_max, ω_max_sho) * 1.1  # 10% margin
        else
            # Rotating: θ increases/decreases continuously
            # Set θ range to cover ~4 complete rotations (8π total range)
            # Direction determined by sign of ω0
            if ω0 > 0
                θ_limit = 4π * 1.1  # Positive rotation
            else
                θ_limit = 4π * 1.1  # Negative rotation (magnitude)
            end
            # ω stays roughly constant for rotation, use a reasonable range
            ω_limit = max(abs(ω0), ω_max) * 1.2
        end

        # Update axes with calculated limits (4 periods now)
        limits!(axTheta, 0, 4 * T_exact_local, -θ_limit, θ_limit)
        limits!(axOmega, 0, 4 * T_exact_local, -ω_limit, ω_limit)

        # Update x-ticks for 4 periods
        xtick_values = [0, T_exact_local, 2*T_exact_local, 3*T_exact_local, 4*T_exact_local]
        xtick_labels = ["0", "T", "2T", "3T", "4T"]
        axTheta.xticks = (xtick_values, xtick_labels)
        axOmega.xticks = (xtick_values, xtick_labels)

        # Store SHO initial conditions
        θ0_sho[] = θ0
        ω0_sho[] = ω0
        
        # Add initial data point so pendulum is visible even when paused
        push!(ts, Float32(0.0))
        push!(thetas, Float32(θ0))
        push!(omegas, Float32(ω0))

        # Update phase space trace
        θvec = θ_phase_cur[][]
        ωvec = ω_phase_cur[][]
        ω₀_phase = sqrt(g / L_obs[])
        push!(θvec, Float32(θ0))
        push!(ωvec, Float32(is_normalized[] ? ω0 / ω₀_phase : ω0))
        θ_phase_cur[][] = θvec
        ω_phase_cur[][] = ωvec

        # Update time series state atomically
        timeseries_state[] = TimeSeriesState(
            Float32[0.0],
            Float32[θ0],
            Float32[ω0],
            Float32[θ0],  # Initialize SHO with same IC
            Float32[ω0]   # Initialize SHO omega with same IC
        )

        θ_phase_cur[][] = θvec
        ω_phase_cur[][] = ωvec
    end

    """
        set_phase_ic!(θ::Float64, ω::Float64)

    Sets the initial conditions (angle θ and angular velocity ω) for the pendulum.
    This function can be used programmatically to control the pendulum's starting state.
    """
    function set_phase_ic!(θ::Float64, ω::Float64)
        set_ic!(θ, ω)
        # Ensure the simulation starts running when ICs are set programmatically
        is_running[] = true
        run_btn.label[] = "Pause"
        @info "Initial conditions set programmatically" θ=θ ω=ω
    end

    # Button callbacks (now that all variables are defined)
    on(run_btn.clicks) do _
        is_running[] = !is_running[]
        run_btn.label[] = is_running[] ? "Pause" : "Run"
        @debug "Run button clicked" is_running=is_running[] data_points=length(ts)
    end
    
    on(rand_btn.clicks) do _
        # Random IC over full range: θ ∈ [-π, π], ω ∈ [-π, π]
        θ0 = rand() * 2π - π
        ω0 = rand() * 2π - π
        set_ic!(θ0, ω0)
        is_running[] = true
        run_btn.label[] = "Pause"
    end

    on(reset_btn.clicks) do _
        set_ic!(deg2rad(120.0), 0.0)
    end

    on(clear_phase_btn.clicks) do _
        # Clear all plot content
        empty!(axPhase[])

        # Re-add separatrix
        lines!(axPhase[], θ_sep, @lift($is_normalized ? ω_sep_normalized_upper : $ω_sep_physical_upper),
               color=(:white, 0.3), linestyle=:dash, linewidth=2)
        lines!(axPhase[], θ_sep, @lift($is_normalized ? ω_sep_normalized_lower : $ω_sep_physical_lower),
               color=(:white, 0.3), linestyle=:dash, linewidth=2)

        # Re-add arrow
        arrows2d!(axPhase[], arrow_θ, arrow_ω, arrow_u, arrow_v,
                  lengthscale=5.0, color=:magenta,
                  shaftwidth=0, tipwidth=10, tiplength=20)

        # Clear trajectory data and start new trace
        θ_phase_cur[][] = Float32[]
        ω_phase_cur[][] = Float32[]
        traj_idx[] = 0
        new_phase_trace!(axPhase, palette, traj_idx, θ_phase_cur, ω_phase_cur)
    end

    on(normalize_toggle.active) do active
        is_normalized[] = active
        # Update axis ylabel and limits without recreating the axis
        if is_normalized[]
            axPhase[].ylabel = "ω/ω₀"
            ylims!(axPhase[], -PHASE_MAX_OMEGA, PHASE_MAX_OMEGA)
        else
            axPhase[].ylabel = "ω [rad/s]"
            ylims!(axPhase[], -PHASE_MAX_OMEGA*sqrt(4), PHASE_MAX_OMEGA*sqrt(4))
        end
        # Clear all trajectory data to avoid mixing normalized/unnormalized
        clear_time_plots!(ts, thetas, omegas, timeseries_state, t_now)
        θ_phase_cur[][] = Float32[]
        ω_phase_cur[][] = Float32[]
        traj_idx[] = 0
        new_phase_trace!(axPhase, palette, traj_idx, θ_phase_cur, ω_phase_cur)
    end

    # Helper function to rebuild mechanism with new parameters
    function rebuild_mechanism!(new_L, new_m)
        # Get current state
        θ_current = configuration(mw.state)[1]
        ω_current = velocity(mw.state)[1]
        was_running = is_running[]

        # Pause during rebuild
        is_running[] = false
        run_btn.label[] = "Run"

        # Rebuild mechanism
        new_mech, new_state, new_result, new_joint = create_mechanism(new_L, new_m)

        # Update mutable wrapper
        mw.mech = new_mech
        mw.state = new_state
        mw.result = new_result
        mw.joint = new_joint

        # Update separatrix for new L and m
        sep_new = compute_separatrix(new_L, new_m)
        ω_sep_physical_upper[] = sep_new[1]
        ω_sep_physical_lower[] = sep_new[2]
        ω_sep_normalized_upper[] = sep_new[3]
        ω_sep_normalized_lower[] = sep_new[4]

        # Restore state
        set_ic!(θ_current, ω_current)

        if was_running
            is_running[] = true
            run_btn.label[] = "Pause"
        end
    end

    # L slider callback - rebuild mechanism with new length
    on(L_sld.value) do new_L
        L_obs[] = new_L
        rebuild_mechanism!(new_L, m_sld.value[])
    end

    # Mass slider callback - rebuild mechanism with new bob mass
    on(m_sld.value) do new_m
        rebuild_mechanism!(L_obs[], new_m)
    end

    new_phase_trace!(axPhase, palette, traj_idx, θ_phase_cur, ω_phase_cur)
    set_ic!(deg2rad(120.0), 0.0)

    # Spacebar to toggle run/pause
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.space
            is_running[] = !is_running[]
            run_btn.label[] = is_running[] ? "Pause" : "Run"
            @debug "Spacebar pressed" is_running=is_running[]
        end
    end

    # Click in phase plot to set initial conditions
    on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            # Check if mouse is over phase plot
            if Makie.is_mouseinside(axPhase[])
                # Get mouse position in data coordinates
                mp = Makie.mouseposition(axPhase[])
                θ_clicked = mp[1]
                ω_clicked_display = mp[2]

                # Clamp θ to phase plot limits
                θ_clicked = clamp(θ_clicked, -π, π)

                # Convert ω from display coordinates to physical coordinates
                if is_normalized[]
                    # In normalized mode, ω_clicked_display is ω/ω₀
                    ω₀ = sqrt(g / L_obs[])
                    ω_clicked = ω_clicked_display * ω₀  # Convert back to physical
                    ω_clicked = clamp(ω_clicked, -12 * ω₀, 12 * ω₀)
                else
                    # In physical mode, use directly
                    ω_clicked = clamp(ω_clicked_display, -20, 20)
                end

                # Set new initial conditions (always in physical units)
                set_ic!(θ_clicked, ω_clicked)
                is_running[] = true
                run_btn.label[] = "Pause"

                @debug "Phase plot clicked" θ=θ_clicked ω=ω_clicked
            end
        end
    end

    # ==========================
    # Simulation loop
    # ==========================
    dt = 0.002
    τ = zeros(eltype(configuration(mw.state)), num_velocities(mw.mech))

    sim_task = @async begin
        sleep(0.1)  # Let the GUI initialize first
        @info "RigidBodyDynamics + Makie running" mode="Run/Pause / Random IC"
        @debug "Initial state" is_running=is_running[] ts_length=length(ts) thetas_length=length(thetas)

        # Frame rate control for smooth 120fps with adjustable speed
        target_fps = 120
        frame_time = 1.0 / target_fps

        update_every = 4  # Update plots every 4 frames for better performance
        k = 0
        frame_count = 0

        while isopen(fig.scene)
            if is_running[]
                # Do multiple physics steps per frame for smooth motion
                c = damp_sld.value[]
                speed_multiplier = speed_sld.value[]
                steps_per_frame = max(1, round(Int, frame_time * speed_multiplier / dt))

                for step in 1:steps_per_frame
                    q = configuration(mw.state)
                    v = velocity(mw.state)

                    # Velocity Verlet integration (2nd order symplectic integrator)
                    # Step 1: Compute acceleration at current state
                    τ .= 0
                    τ[1] = -c * v[1]
                    dynamics!(mw.result, mw.state, τ)
                    a = mw.result.v̇

                    # Step 2: Half-step velocity update
                    v_half = v .+ 0.5 .* dt .* a

                    # Step 3: Full-step position update using half-step velocity
                    q .= q .+ dt .* v_half
                    set_configuration!(mw.state, q)

                    # Step 4: Compute acceleration at new position
                    τ .= 0
                    τ[1] = -c * v_half[1]
                    dynamics!(mw.result, mw.state, τ)
                    a_new = mw.result.v̇

                    # Step 5: Complete velocity update
                    v .= v_half .+ 0.5 .* dt .* a_new
                    set_velocity!(mw.state, v)
                    t_now[] = t_now[] + dt

                    # Track wrapped angle for phase space (always in [-π, π])
                    newq = q[1]
                    while newq> π || newq < -π
                        if newq > π
                            newq -= 2π
                        elseif newq < -π
                            newq += 2π
                        end
                    end
                    newv = v[1]

                    # Track unwrapped angle for continuous time series
                    # Detect wraps and accumulate offset
                    θ_diff = newq - θ_prev[]
                    if θ_diff > π
                        θ_unwrapped[] -= 2π
                    elseif θ_diff < -π
                        θ_unwrapped[] += 2π
                    end
                    θ_prev[] = newq
                    θ_continuous = θ_unwrapped[] + newq

                    # Only add time series data if we're within first 4 periods
                    four_periods_duration = 4 * current_period[]
                    if t_now[] <= four_periods_duration
                        push!(ts, Float32(t_now[]))
                        push!(thetas, Float32(θ_continuous))  # Use unwrapped angle
                        push!(omegas, Float32(newv))
                    end
                    θvec = θ_phase_cur[][]
                    ωvec = ω_phase_cur[][]

                    # Normalize ω for phase space if needed
                    ω₀ = sqrt(g / L_obs[])
                    newv_phase = is_normalized[] ? newv / ω₀ : newv

                    # Detect wrapping discontinuity and insert NaN to break the line
                    if !isempty(θvec)
                        θ_last = θvec[end]
                        if abs(newq - θ_last) > π
                            # Wrapping detected - insert NaN to break the line
                            push!(θvec, NaN32)
                            push!(ωvec, NaN32)
                        end
                    end

                    push!(θvec, Float32(newq))
                    push!(ωvec, Float32(newv_phase))
                    θ_phase_cur[][] = θvec
                    ω_phase_cur[][] = ωvec
                end  # End of for loop over steps_per_frame

                # Update observables after all physics steps for smooth pendulum motion
                θ_obs[] = configuration(mw.state)[1]
                ω_obs[] = velocity(mw.state)[1]

                # Update phase space arrow (current position and direction)
                # Wrap θ to [-π, π] using same method as time series data
                θ_wrapped = configuration(mw.state)[1]
                while θ_wrapped > π || θ_wrapped < -π
                    if θ_wrapped > π
                        θ_wrapped -= 2π
                    elseif θ_wrapped < -π
                        θ_wrapped += 2π
                    end
                end
                ω_val = Float32(velocity(mw.state)[1])
                α_val = Float32(mw.result.v̇[1])

                # Normalize position for phase space display
                ω₀ = sqrt(g / L_obs[])
                ω_val_phase = is_normalized[] ? ω_val / ω₀ : ω_val

                # Update arrow position
                arrow_θ[] = [θ_wrapped]
                arrow_ω[] = [ω_val_phase]

                # Direction in phase space: (dθ/dt, d(ω or ω/ω₀)/dt)
                # dθ/dt = ω (always)
                # d(ω/ω₀)/dt = α/ω₀ (in normalized coords) or dω/dt = α (in physical coords)
                dir_u = ω_val
                dir_v = is_normalized[] ? α_val / ω₀ : α_val
                dir_mag = sqrt(dir_u^2 + dir_v^2) + 1e-6  # Avoid division by zero
                arrow_u[] = [dir_u / dir_mag]  # Normalized horizontal component
                arrow_v[] = [dir_v / dir_mag]  # Normalized vertical component

                # Calculate energies in Joules using mass slider
                # KE = (1/2) * m * L^2 * ω^2
                # PE = m * g * L * (1 - cos(θ))
                θ_current = θ_obs[]
                ω_current = ω_obs[]
                m_current = m_sld.value[]

                ke_joules = 0.5 * m_current * L_obs[]^2 * ω_current^2
                pe_joules = m_current * g * L_obs[] * (1.0 - cos(θ_current))

                # Calculate energy proportions for bar display
                total_energy = ke_joules + pe_joules
                if total_energy > 1e-10  # Avoid division by zero
                    ke_proportion[] = ke_joules / total_energy
                    pe_proportion[] = pe_joules / total_energy
                else
                    ke_proportion[] = 0.5
                    pe_proportion[] = 0.5
                end

                # Update current position markers with time wrapping
                # Interpolate from recorded time series to show position on the original trajectory
                four_periods_duration = 4 * current_period[]
                t_wrapped = mod(t_now[], four_periods_duration)

                # Find the marker values by interpolating from the recorded time series
                if !isempty(ts) && length(ts) > 1
                    # Find where t_wrapped falls in the time series
                    idx = searchsortedlast(ts, Float32(t_wrapped))
                    if idx >= 1 && idx < length(ts)
                        # Linear interpolation between ts[idx] and ts[idx+1]
                        t1, t2 = ts[idx], ts[idx+1]
                        if t2 > t1  # Avoid division by zero
                            α = (Float32(t_wrapped) - t1) / (t2 - t1)
                            θ_marker = thetas[idx] * (1 - α) + thetas[idx+1] * α
                            ω_marker = omegas[idx] * (1 - α) + omegas[idx+1] * α
                            current_t[] = [Float32(t_wrapped)]
                            current_theta[] = [θ_marker]
                            current_omega[] = [ω_marker]
                        end
                    elseif idx == length(ts)
                        # Use last data point
                        current_t[] = [ts[end]]
                        current_theta[] = [thetas[end]]
                        current_omega[] = [omegas[end]]
                    elseif idx == 0 && !isempty(ts)
                        # Use first data point
                        current_t[] = [ts[1]]
                        current_theta[] = [thetas[1]]
                        current_omega[] = [omegas[1]]
                    end
                end

                frame_count += 1
                k += 1

                # Update observable every update_every frames
                if k % update_every == 0
                    # Create synchronized copies for atomic update
                    n_points = length(ts)
                    ts_copy = Vector{Float32}(undef, n_points)
                    thetas_copy = Vector{Float32}(undef, n_points)
                    omegas_copy = Vector{Float32}(undef, n_points)
                    thetas_sho = Vector{Float32}(undef, n_points)
                    omegas_sho = Vector{Float32}(undef, n_points)

                    # Calculate SHO solution and copy data in one pass
                    # General SHO solution with initial conditions θ(0)=θ0, ω(0)=ω0:
                    # θ(t) = θ0*cos(ω_sho*t) + (ω0/ω_sho)*sin(ω_sho*t)
                    # ω(t) = -θ0*ω_sho*sin(ω_sho*t) + ω0*cos(ω_sho*t)
                    # Use idealized simple pendulum frequency (compare to physical pendulum simulation)
                    L_current = L_obs[]
                    ω_sho = sqrt(g / L_current)
                    for i in 1:n_points
                        ts_copy[i] = ts[i]
                        thetas_copy[i] = thetas[i]
                        omegas_copy[i] = omegas[i]

                        t = ts[i]
                        θ_raw = θ0_sho[] * cos(ω_sho * t) + (ω0_sho[] / ω_sho) * sin(ω_sho * t)
                        thetas_sho[i] = θ_raw  # Keep unwrapped for continuous display
                        omegas_sho[i] = -θ0_sho[] * ω_sho * sin(ω_sho * t) + ω0_sho[] * cos(ω_sho * t)
                    end

                    # Single atomic update - all arrays guaranteed to have same length
                    timeseries_state[] = TimeSeriesState(ts_copy, thetas_copy, omegas_copy, thetas_sho, omegas_sho)
                end
            end
            sleep(frame_time)  # 60fps frame rate for smooth motion
        end
    end

    # Return all GUI controls for programmatic access
    gui = GUIControls(
        fig, screen,
        run_btn, rand_btn, reset_btn, clear_phase_btn, #
        damp_sld, L_sld, m_sld, speed_sld,
        is_running, is_normalized, L_obs, θ_obs, ω_obs, t_now, ke_proportion, pe_proportion,
        set_ic!, rebuild_mechanism!, set_phase_ic!
    )
    CURRENT_GUI[] = gui
    return gui
end

"""
    get_screenshot_path() -> String

Returns the path to the most recent screenshot saved by the simulation.
Screenshots are saved to ~/.pendulum_screenshots/ with timestamps.
The circular buffer keeps the last 50 screenshots.
"""
function get_screenshot_path()
    global screenshot_dir
    if !isdir(screenshot_dir)
        @warn "Screenshot directory not found: $screenshot_dir"
        return ""
    end

    # Get all screenshots and return the most recent one
    all_screenshots = filter(f -> endswith(f, ".png") && startswith(f, "pendulum_"), readdir(screenshot_dir))
    if isempty(all_screenshots)
        @warn "No screenshots found in $screenshot_dir"
        return ""
    end

    # Sort and return the latest (last in sorted order)
    sort!(all_screenshots)
    joinpath(screenshot_dir, last(all_screenshots))
end

"""
    get_all_screenshot_paths(; n::Int=50) -> Vector{String}

Returns paths to the n most recent screenshots (default: 50).
Returns them in chronological order (oldest to newest).
Useful for reviewing simulation history.
"""
function get_all_screenshot_paths(; n::Int=50)
    global screenshot_dir
    if !isdir(screenshot_dir)
        @warn "Screenshot directory not found: $screenshot_dir"
        return String[]
    end

    # Get all screenshots
    all_screenshots = filter(f -> endswith(f, ".png") && startswith(f, "pendulum_"), readdir(screenshot_dir))
    if isempty(all_screenshots)
        return String[]
    end

    # Sort chronologically and take last n
    sort!(all_screenshots)
    recent = all_screenshots[max(1, end-n+1):end]

    # Return full paths
    [joinpath(screenshot_dir, f) for f in recent]
end

"""
    view_latest_screenshot()

Opens the latest screenshot in your default image viewer.
"""
function view_latest_screenshot()
    path = get_screenshot_path()
    if isfile(path)
        run(`open $path`)  # macOS command to open file
        @info "Opening screenshot" path=path
    else
        @warn "No screenshot found at $path. Make sure run_gl2() is running."
    end
end

"""
    restart_gui!()

Closes the current GUI window (if any) and starts a fresh one.
This allows agents to reload the GUI after code changes without manual intervention.

# Example
```julia
using PendulumSim
PendulumSim.SinglePendulum.run_gui()

# ... make some code changes ...

PendulumSim.SinglePendulum.restart_gui!()  # Close old window and open new one
```
"""
function restart_gui!()
    # Close existing window if present
    if CURRENT_SCREEN[] !== nothing
        try
            # Check if window is still valid before closing
            if CURRENT_SCREEN[].glscreen.handle != C_NULL
                GLMakie.GLFW.SetWindowShouldClose(CURRENT_SCREEN[].glscreen, true)
                @debug "Closed existing GUI window"
            end
        catch e
            @debug "Window already closed or invalid" exception=e
        finally
            CURRENT_SCREEN[] = nothing
        end
    end

    # Small delay to let window close cleanly
    sleep(0.2)

    # Start fresh GUI
    @info "Starting new GUI..."
    gui = try
        run_gui()
    catch e
        @warn "Failed to restart GUI with run_gui(), trying direct call" exception=e
        # If Revise changed the function signature, call it directly from the module
        Base.invokelatest(run_gui)
    end

    return gui
end

"""
    take_screenshot!() -> String

Takes a screenshot of the current GUI and returns the file path.
Useful for programmatically inspecting the GUI state.

# Example
```julia
gui = PendulumSim.CURRENT_GUI[]
gui.damp_sld.value[] = 0.02  # Set damping
path = PendulumSim.take_screenshot!()
# Now inspect the screenshot at path
```
"""
function take_screenshot!()
    if CURRENT_GUI[] === nothing || CURRENT_GUI[].fig === nothing
        error("No GUI is currently running. Call run_gui() or restart_gui!() first.")
    end

    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS-sss")
    screenshot_path = joinpath(screenshot_dir, "manual_$timestamp.png")
    save(screenshot_path, CURRENT_GUI[].fig)
    return screenshot_path
end

end # module SinglePendulum
