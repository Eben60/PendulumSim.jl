module PendulumSim

using DifferentialEquations, GLMakie, Elliptic, Logging
using LinearAlgebra, Dates
using RigidBodyDynamics
using StaticArrays
using Rotations
using Colors

# Global reference to the current screen for GUI management
const CURRENT_SCREEN = Ref{Union{Nothing, GLMakie.Screen}}(nothing)

function run_gl2()
    # ==========================
    # Logging setup
    # ==========================
    log_dir = joinpath(homedir(), ".pendulum_screenshots")
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
    L = 1.0
    m_rod = 0.10
    m_bob = 0.80

    frame = CartesianFrame3D("link")

    Irod_com = Diagonal(SVector((1 / 12) * m_rod * L^2, 0.0, (1 / 12) * m_rod * L^2))
    I_rod = SpatialInertia(frame;
        mass=m_rod,
        com=SVector(0.0, -L / 2, 0.0),
        moment_about_com=Matrix(Irod_com))

    I_bob = SpatialInertia(frame;
        mass=m_bob,
        com=SVector(0.0, -L, 0.0),
        moment_about_com=zeros(3, 3))

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

    state = MechanismState(mech)
    result = DynamicsResult(mech)

    # ==========================
    # Visualization (Makie)
    # ==========================
    GLMakie.activate!()
    set_theme!(theme_dark())

    fig = Figure(resolution=(1400, 800), figure_padding=20)
    
    # Pendulum visualization (left, top)
    axPend = Axis(fig[1, 1], 
        title="Pendulum (plan view)", 
        xlabel="x [m]", 
        ylabel="y [m]", 
        aspect=DataAspect())
    xlims!(axPend, -1.3, 1.3)
    ylims!(axPend, -1.3, 0.3)
    hidedecorations!(axPend, grid=false, ticks=false)

    # Phase space plot (right, top) - aspect ratio 1:1
    newAxPhase() = begin
        ax = Axis(fig[1, 2], 
            title="Phase space (θ, ω)", 
            xlabel="θ [rad]", 
            ylabel="ω [rad/s]", 
            aspect=1)  # Force 1:1 aspect ratio
        ax
    end
    axPhase = Observable(newAxPhase())
    xlims!(axPhase[], -π, π)
    ylims!(axPhase[], -6.5, 6.5)

    # Time series plots (left column, rows 2 and 3)
    axTheta = Axis(fig[2, 1], 
        title="θ(t)", 
        xlabel="", 
        ylabel="θ [rad]")
    
    axOmega = Axis(fig[3, 1],
        title="ω(t)",
        xlabel="t [s]",
        ylabel="ω [rad/s]")
    
    # Calculate period for initial condition and set xlim to show 2 periods
    θ0_initial = deg2rad(120.0)
    T_sho = 2π * sqrt(L / g)  # Small angle approximation
    T_exact = T_sho * (2 / π) * Elliptic.K(sin(θ0_initial / 2)^2)  # Large angle period
    two_periods = 2 * T_exact
    
    xlims!(axTheta, 0, two_periods)
    ylims!(axTheta, -π, π)
    xlims!(axOmega, 0, two_periods)
    ylims!(axOmega, -7, 7)

    # Control panel (right column, rows 2-3) with better spacing
    controls = GridLayout(fig[2:3, 2], tellwidth=false, tellheight=false)
    
    # Buttons row (callbacks will be set up after data arrays are defined)
    run_btn = Button(controls[1, 1:3], label="Run", width=Auto())
    is_running = Observable(false)  # Start paused

    rand_btn = Button(controls[2, 1], label="Random IC")
    reset_btn = Button(controls[2, 2], label="Reset")
    clear_phase_btn = Button(controls[2, 3], label="Clear phase")

    # Spacer
    Label(controls[3, 1:3], " ")
    
    # Sliders with labels
    lbl1 = Label(controls[4, 1:3], "Max |θ₀| (deg)", halign=:left)
    θmax_sld = Slider(controls[5, 1:3], range=0:1:170, startvalue=120)
    
    lbl2 = Label(controls[6, 1:3], "Max |ω₀| (rad/s)", halign=:left)
    ωmax_sld = Slider(controls[7, 1:3], range=0:0.1:6.0, startvalue=0.5)
    
    lbl3 = Label(controls[8, 1:3], "Damping c (N·m·s)", halign=:left)
    c_sld = Slider(controls[9, 1:3], range=0:0.001:0.05, startvalue=0.0)
    
    # Add spacing between plots and controls
    rowgap!(fig.layout, 20)
    colgap!(fig.layout, 20)

    # Display and set window to floating mode at upper left
    screen = display(GLMakie.Screen(), fig)
    GLMakie.GLFW.SetWindowAttrib(screen.glscreen, GLMakie.GLFW.FLOATING, true)
    GLMakie.GLFW.SetWindowPos(screen.glscreen, 50, 50)  # Position at (50, 50) pixels from top-left
    
    # Store global reference for GUI management
    CURRENT_SCREEN[] = screen

    # ==========================
    # Screenshot task (async, every 1 second, circular buffer of 50)
    # ==========================
    screenshot_dir = joinpath(homedir(), ".pendulum_screenshots")
    mkpath(screenshot_dir)

    @async begin
        while isopen(fig.scene)
            try
                # Create timestamped filename
                timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
                screenshot_path = joinpath(screenshot_dir, "pendulum_$timestamp.png")

                # Save screenshot
                save(screenshot_path, fig)

                # Clean up old screenshots (keep only last 50)
                all_screenshots = filter(f -> endswith(f, ".png") && startswith(f, "pendulum_"), readdir(screenshot_dir))
                if length(all_screenshots) > 50
                    # Sort by filename (which includes timestamp, so chronological)
                    sort!(all_screenshots)
                    # Delete oldest ones
                    num_to_delete = length(all_screenshots) - 50
                    for i in 1:num_to_delete
                        old_file = joinpath(screenshot_dir, all_screenshots[i])
                        rm(old_file)
                    end
                end
            catch e
                @warn "Screenshot failed: $e"
            end
            sleep(1.0)  # Save every 1 second
        end
    end

    # ==========================
    # Observables / data buffers
    # ==========================
    θ_obs = Observable(0.0)
    ω_obs = Observable(0.0)
    t_now = Observable(0.0)

    rod_pts = @lift(Point2f[(0, 0), (Float32(L * sin($θ_obs)), Float32(-L * cos($θ_obs)))])
    bob_pt = @lift(Point2f(Float32(L * sin($θ_obs)), Float32(-L * cos($θ_obs))))
    lines!(axPend, rod_pts, color=:gray80, linewidth=5)
    scatter!(axPend, bob_pt, color=:tomato, markersize=20)

    ts = Float32[]
    thetas = Float32[]
    omegas = Float32[]
    ts_obs = Observable(copy(ts))
    thetas_obs = Observable(copy(thetas))
    omegas_obs = Observable(copy(omegas))
    
    # SHO solution observables
    θ0_sho = Observable(θ0_initial)
    ω0_sho = Observable(0.0)
    thetas_sho_obs = Observable(Float32[])
    omegas_sho_obs = Observable(Float32[])
    
    # Plot actual theta solution
    lines!(axTheta, ts_obs, thetas_obs, color=:cyan, linewidth=3, label="Actual")
    
    # Plot SHO theta solution (semi-transparent)
    lines!(axTheta, ts_obs, thetas_sho_obs, color=(:orange, 0.5), linewidth=2, label="SHO", linestyle=:dash)
    
    # Fill between theta plots
    band!(axTheta, ts_obs, thetas_obs, thetas_sho_obs, color=(:yellow, 0.2))
    
    axislegend(axTheta, position=:rt)
    
    # Plot actual omega solution
    lines!(axOmega, ts_obs, omegas_obs, color=:cyan, linewidth=3, label="Actual")
    
    # Plot SHO omega solution (semi-transparent)
    lines!(axOmega, ts_obs, omegas_sho_obs, color=(:orange, 0.5), linewidth=2, label="SHO", linestyle=:dash)
    
    # Fill between omega plots
    band!(axOmega, ts_obs, omegas_obs, omegas_sho_obs, color=(:yellow, 0.2))
    
    axislegend(axOmega, position=:rt)
    # --- Phase traces ---
    palette = [RGB(0.90, 0.30, 0.25), RGB(0.20, 0.70, 0.90), RGB(0.30, 0.85, 0.40),
        RGB(0.55, 0.45, 0.85), RGB(0.95, 0.70, 0.20), RGB(0.40, 0.80, 0.70)]

    traj_idx = Ref(0)
    θ_phase_cur = Ref(Observable(Float32[]))
    ω_phase_cur = Ref(Observable(Float32[]))

    function new_phase_trace!(axPhase, palette, traj_idx, θref::Ref{Observable{Vector{Float32}}},
        ωref::Ref{Observable{Vector{Float32}}})
        traj_idx[] += 1
        θref[] = Observable(Float32[])
        ωref[] = Observable(Float32[])
        col = palette[mod1(traj_idx[], length(palette))]
        lines!(axPhase[], θref[], ωref[], color=col, linewidth=2)
    end

    function clear_time_plots!(ts, thetas, omegas, ts_obs, thetas_obs, omegas_obs, t_now)
        empty!(ts)
        empty!(thetas)
        empty!(omegas)
        ts_obs[] = ts
        thetas_obs[] = thetas
        omegas_obs[] = omegas
        t_now[] = 0.0
    end

    # Observable for current period (updates when IC changes)
    current_period = Observable(T_exact)
    
    # ==========================
    # IC handling
    # ==========================
    function set_ic!(θ0::Float64, ω0::Float64)
        set_configuration!(state, joint, θ0)
        set_velocity!(state, joint, ω0)
        θ_obs[] = θ0
        ω_obs[] = ω0
        clear_time_plots!(ts, thetas, omegas, ts_obs, thetas_obs, omegas_obs, t_now)
        new_phase_trace!(axPhase, palette, traj_idx, θ_phase_cur, ω_phase_cur)
        
        # Calculate period for this IC
        T_sho_local = 2π * sqrt(L / g)
        T_exact_local = T_sho_local * (2 / π) * Elliptic.K(sin(abs(θ0) / 2)^2)
        current_period[] = T_exact_local
        
        # Update x-axis to show 2 periods
        xlims!(axTheta, 0, 2 * T_exact_local)
        
        # Store SHO initial conditions
        θ0_sho[] = θ0
        ω0_sho[] = ω0
        
        # Add initial data point so pendulum is visible even when paused
        push!(ts, Float32(0.0))
        push!(thetas, Float32(θ0))
        push!(omegas, Float32(ω0))
        θvec = θ_phase_cur[][]
        ωvec = ω_phase_cur[][]
        push!(θvec, Float32(θ0))
        push!(ωvec, Float32(ω0))
        ts_obs[] = ts
        thetas_obs[] = thetas
        omegas_obs[] = omegas
        thetas_sho_obs[] = Float32[θ0]  # Initialize SHO with same IC
        omegas_sho_obs[] = Float32[ω0]  # Initialize SHO omega with same IC
        θ_phase_cur[][] = θvec
        ω_phase_cur[][] = ωvec
    end

    # Button callbacks (now that all variables are defined)
    on(run_btn.clicks) do _
        is_running[] = !is_running[]
        run_btn.label[] = is_running[] ? "Pause" : "Run"
        println("Run button clicked: is_running = $(is_running[]), data points = $(length(ts))")
    end
    
    on(rand_btn.clicks) do _
        θmax = deg2rad(θmax_sld.value[])
        ωmax = ωmax_sld.value[]
        θ0 = rand() * 2θmax - θmax
        ω0 = rand() * 2ωmax - ωmax
        set_ic!(θ0, ω0)
        is_running[] = true
        run_btn.label[] = "Pause"
    end

    on(reset_btn.clicks) do _
        set_ic!(deg2rad(120.0), 0.0)
    end

    on(clear_phase_btn.clicks) do _
        empty!(axPhase[])
        axPhase = Observable(newAxPhase())
        
        xlims!(axPhase[], -π, π)
        ylims!(axPhase[], -6.5, 6.5)
        θ_phase_cur[][] = Float32[]
        ω_phase_cur[][] = Float32[]
        traj_idx[] = 0
        new_phase_trace!(axPhase, palette, traj_idx, θ_phase_cur, ω_phase_cur)
    end

    new_phase_trace!(axPhase, palette, traj_idx, θ_phase_cur, ω_phase_cur)
    set_ic!(deg2rad(120.0), 0.0)

    # ==========================
    # Simulation loop
    # ==========================
    dt = 0.002
    τ = zeros(eltype(configuration(state)), num_velocities(mech))

    sim_task = @async begin
        sleep(0.1)  # Let the GUI initialize first
        println("RigidBodyDynamics + Makie running... (Run/Pause / Random IC).")
        println("Initial is_running state: $(is_running[])")
        println("Initial data: ts=$(length(ts)), thetas=$(length(thetas))")
        update_every = 2
        k = 0
        loop_count = 0
        iter = 0
        while isopen(fig.scene)
            iter += 1
            if iter <= 3 || (iter % 1000 == 0)
                println("Loop iteration $iter: is_running=$(is_running[]), data points=$(length(ts))")
            end
            if is_running[]
                loop_count += 1
                if loop_count == 1 || loop_count % 500 == 0
                    println("Simulation running: loop_count=$loop_count, data points=$(length(ts))")
                end
                c = c_sld.value[]
                q = configuration(state)
                v = velocity(state)
                τ .= 0
                τ[1] = -c * v[1]
                dynamics!(result, state, τ)
                v .= v .+ dt .* result.v̇
                q .= q .+ dt .* v
                set_velocity!(state, v)
                set_configuration!(state, q)
                θ_obs[] = q[1]
                ω_obs[] = v[1]
                t_now[] = t_now[] + dt
                push!(ts, Float32(t_now[]))
                newq = q[1]
                while newq> π || newq < -π
                    if newq > π
                        newq -= 2π
                    elseif newq < -π
                        newq += 2π
                    end
                end
                newv = v[1]
                push!(thetas, Float32(newq))
                push!(omegas, Float32(newv))
                θvec = θ_phase_cur[][]
                ωvec = ω_phase_cur[][]
                push!(θvec, Float32(newq))
                push!(ωvec, Float32(newv))
                θ_phase_cur[][] = θvec
                ω_phase_cur[][] = ωvec
                k += 1
                if k % update_every == 0
                    # Implement rolling 2-period window
                    # Calculate how many points correspond to 2 periods
                    two_periods_duration = 2 * current_period[]
                    max_points = ceil(Int, two_periods_duration / dt) + 100  # Add buffer

                    # Trim arrays if they exceed the window size
                    if length(ts) > max_points
                        n_to_remove = length(ts) - max_points
                        deleteat!(ts, 1:n_to_remove)
                        deleteat!(thetas, 1:n_to_remove)
                        deleteat!(omegas, 1:n_to_remove)
                    end

                    # Fix DimensionMismatch: capture length once to ensure all arrays match
                    n_points = length(ts)
                    ts_obs[] = ts[1:n_points]
                    thetas_obs[] = thetas[1:n_points]
                    omegas_obs[] = omegas[1:n_points]

                    # Update x-axis to show rolling window
                    if n_points > 0
                        t_start = ts[1]
                        t_end = ts[end]
                        xlims!(axTheta, t_start, max(t_end, t_start + two_periods_duration))
                        xlims!(axOmega, t_start, max(t_end, t_start + two_periods_duration))
                    end

                    # Calculate SHO solution at current time points
                    # Use period-matched frequency so SHO has same period as actual motion
                    ω_sho = 2π / current_period[]
                    ts_slice = ts[1:n_points]

                    # SHO uses absolute time from start (t=0) with initial conditions
                    # This ensures periods stay aligned with actual motion
                    thetas_sho_raw = [θ0_sho[] * cos(ω_sho * t) for t in ts_slice]

                    # Wrap SHO theta values to [-π, π] to match actual pendulum wrapping
                    thetas_sho = map(thetas_sho_raw) do θ
                        θ_wrapped = θ
                        while θ_wrapped > π
                            θ_wrapped -= 2π
                        end
                        while θ_wrapped < -π
                            θ_wrapped += 2π
                        end
                        θ_wrapped
                    end

                    omegas_sho = [-θ0_sho[] * ω_sho * sin(ω_sho * t) for t in ts_slice]
                    thetas_sho_obs[] = Float32.(thetas_sho)
                    omegas_sho_obs[] = Float32.(omegas_sho)
                end
            end
            sleep(dt * 0.75)  # Run slightly faster than real-time (25% speedup)
        end
    end
end

"""
    get_screenshot_path() -> String

Returns the path to the most recent screenshot saved by the simulation.
Screenshots are saved to ~/.pendulum_screenshots/ with timestamps.
The circular buffer keeps the last 50 screenshots.
"""
function get_screenshot_path()
    screenshot_dir = joinpath(homedir(), ".pendulum_screenshots")
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
    screenshot_dir = joinpath(homedir(), ".pendulum_screenshots")
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
        println("Opening: $path")
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
PendulumSim.run_gl2()

# ... make some code changes ...

PendulumSim.restart_gui!()  # Close old window and open new one
```
"""
function restart_gui!()
    # Close existing window if present
    if CURRENT_SCREEN[] !== nothing
        try
            close(CURRENT_SCREEN[].glscreen)
            println("Closed existing GUI window")
        catch e
            @warn "Failed to close existing window: $e"
        end
        CURRENT_SCREEN[] = nothing
    end
    
    # Small delay to let window close cleanly
    sleep(0.2)
    
    # Start fresh GUI
    println("Starting new GUI...")
    run_gl2()
    
    return nothing
end

end # module PendulumSim
