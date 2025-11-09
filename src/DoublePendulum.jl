module DoublePendulum

using DifferentialEquations, GLMakie, LinearAlgebra, Dates, Colors
using StaticArrays
using Logging

"""
    run_double_pendulum_gui()

Launch the full GUI double pendulum simulation with logging and dynamic 3D phase space projection.
This maintains screenshot tasks, observables, and energy visualization.
"""
function run_double_pendulum_gui()
    @info "Launching full Double Pendulum GUI..."
    g = 9.81
    m₁_obs = Observable(1.0)
    m₂_obs = Observable(1.0)
    L₁_obs = Observable(1.0)
    L₂_obs = Observable(1.0)

    # Projection selection
    projection_options = ["drop θ₁", "drop ω₁", "drop θ₂", "drop ω₂"]
    projection_choice = Observable(projection_options[1])

    # Logging and screenshots
    screenshot_dir = abspath(joinpath(@__DIR__, "..", "screenshots"))
    mkpath(screenshot_dir)
    log_dir = abspath(joinpath(@__DIR__, "..", "logs"))
    mkpath(log_dir)
    logger_path = joinpath(log_dir, "double_pendulum.log")
    global_logger(SimpleLogger(open(logger_path, "a")))
    @info "Double Pendulum Simulation started" log=logger_path

    # Default initial conditions (with some interesting motion)
    θ₁₀, ω₁₀ = π / 3, 0.0  # Start at 60 degrees from vertical
    θ₂₀, ω₂₀ = π / 2, 0.0  # Second pendulum at 90 degrees
    u₀ = [θ₁₀, ω₁₀, θ₂₀, ω₂₀]

    # Simulation parameters
    g_obs = Observable(9.81)
    damping_obs = Observable(0.0)
    time_speed_obs = Observable(1.0)

    function f!(du, u, p, t)
        θ₁, ω₁, θ₂, ω₂ = u
        g, m₁, m₂, L₁, L₂, damping = p
        Δ = θ₁ - θ₂
        denom₁ = (m₁ + m₂) * L₁ - m₂ * L₁ * cos(Δ)^2
        denom₂ = (L₂ / L₁) * denom₁
        du[1] = ω₁
        du[2] = (m₂ * L₁ * ω₁^2 * sin(Δ) * cos(Δ) +
                 m₂ * g * sin(θ₂) * cos(Δ) +
                 m₂ * L₂ * ω₂^2 * sin(Δ) -
                 (m₁ + m₂) * g * sin(θ₁) -
                 damping * ω₁) / denom₁
        du[3] = ω₂
        du[4] = (-m₂ * L₂ * ω₂^2 * sin(Δ) * cos(Δ) +
                 (m₁ + m₂) * (g * sin(θ₁) * cos(Δ) -
                 L₁ * ω₁^2 * sin(Δ) - g * sin(θ₂)) -
                 damping * ω₂) / denom₂
    end

    # Energy calculation function
    function calculate_energy(θ₁, ω₁, θ₂, ω₂, m₁, m₂, L₁, L₂, g)
        # Positions of the bobs
        x₁ = L₁ * sin(θ₁)
        y₁ = -L₁ * cos(θ₁)
        x₂ = x₁ + L₂ * sin(θ₂)
        y₂ = y₁ - L₂ * cos(θ₂)

        # Velocities of the bobs
        vx₁ = L₁ * ω₁ * cos(θ₁)
        vy₁ = L₁ * ω₁ * sin(θ₁)
        vx₂ = vx₁ + L₂ * ω₂ * cos(θ₂)
        vy₂ = vy₁ + L₂ * ω₂ * sin(θ₂)

        # Kinetic energy
        KE = 0.5 * m₁ * (vx₁^2 + vy₁^2) + 0.5 * m₂ * (vx₂^2 + vy₂^2)

        # Potential energy (using y=0 as reference at top)
        PE = m₁ * g * y₁ + m₂ * g * y₂

        return KE, PE
    end

    # GUI setup
    GLMakie.activate!(; float=true, focus_on_show=true)
    set_theme!(theme_dark())
    fig = Figure(resolution=(1400, 700))  # Wider to accommodate 3 columns

    # Layout grid for structured UI
    grid = fig[1, 1] = GridLayout()

    # Top: Title section and controls
    Label(grid[1, 1], "Double Pendulum Simulation", halign=:center, fontsize=18, font=:bold)

    # Status displays
    fps_obs = Observable(0.0)
    energy_drift_obs = Observable(0.0)
    status_label = Label(grid[1, 2], @lift("FPS: $(round($fps_obs, digits=1)) | ΔE: $(round($energy_drift_obs*100, digits=3))%"),
                        halign=:center, fontsize=12)

    projection_dropdown = Menu(grid[1, 3], options=projection_options, width=140, default=projection_options[1])

    # Middle: visualization panels (3 columns now)
    ax3d = LScene(grid[2, 1], scenekw=(clear=true,))
    Label(grid[2, 1, Top()], "3D Phase Space", halign=:center, fontsize=12, font=:bold)

    axPend = Axis(grid[2, 2], aspect=DataAspect(), title="Double Pendulum Motion")

    # Phase space plot for first pendulum (θ₁ vs ω₁) - click to set IC
    axPhase = Axis(grid[2, 3], xlabel="θ₁ (rad)", ylabel="ω₁ (rad/s)", title="Phase Space (Click to Set IC)")
    xlims!(axPhase, -π, π)
    ylims!(axPhase, -10, 10)

    # Phase space trajectory marker
    phase_pt = Observable(Point2f(π/3, 0.0))
    phase_trail = Observable(Point2f[])

    # Draw phase space trajectory and current point
    lines!(axPhase, phase_trail, color=(:cyan, 0.4), linewidth=1)
    scatter!(axPhase, phase_pt, color=:magenta, markersize=12)

    idx = Observable(1) # Make idx an Observable

    # Trail for second bob (to show chaotic motion)
    trail_length = 100
    trail_pts = Observable(Point2f[])

    # 3D phase space trail
    phase3d_trail = Observable(Point3f[])

    xlims!(axPend, -2.2, 2.2)
    ylims!(axPend, -2.2, 2.2)
    hidedecorations!(axPend)

    onany(L₁_obs, L₂_obs) do L1, L2
        max_len = (L1 + L2) * 1.1
        xlims!(axPend, -max_len, max_len)
        ylims!(axPend, -max_len, max_len)
    end

    # Bottom: Energy bar and display
    energy_grid = GridLayout(grid[3, 1:3])

    # Energy values and trends
    ke_val_obs = Observable(0.0)
    pe_val_obs = Observable(0.0)
    total_e_obs = Observable(0.0)
    ke_trend_obs = Observable("")  # "↑", "↓", or ""
    pe_trend_obs = Observable("")

    # Energy labels
    ke_label = Label(energy_grid[1, 1], @lift("KE: $(round($ke_val_obs, digits=2)) J $($ke_trend_obs)"),
                     halign=:left, color=:deepskyblue, fontsize=12)
    pe_label = Label(energy_grid[1, 2], @lift("PE: $(round($pe_val_obs, digits=2)) J $($pe_trend_obs)"),
                     halign=:center, color=:limegreen, fontsize=12)
    total_label = Label(energy_grid[1, 3], @lift("Total: $(round($total_e_obs, digits=2)) J"),
                       halign=:right, fontsize=12)

    # Energy bar
    axKE = Axis(energy_grid[2, 1:3], height=20, title="", limits=(0, 1, 0, 1))
    ke_fraction_obs = Observable(0.5)
    ke_bar = @lift(Rect(0, 0, $ke_fraction_obs, 1))
    pe_bar = @lift(Rect($ke_fraction_obs, 0, 1 - $ke_fraction_obs, 1))
    poly!(axKE, ke_bar, color=:deepskyblue)
    poly!(axKE, pe_bar, color=:limegreen)
    hidespines!(axKE)
    hidedecorations!(axKE)

    # Control buttons layout
    controls = GridLayout(grid[4, 1:3])
    run_btn = Button(controls[1, 1], label="Run", width=60)
    pause_btn = Button(controls[1, 2], label="Pause", width=60)
    reset_btn = Button(controls[1, 3], label="Reset", width=60)
    rand_btn = Button(controls[1, 4], label="Random", width=60)
    clear_btn = Button(controls[1, 5], label="Clear", width=60)

   # Sliders for physical parameters
   Label(controls[2, 1], "m₁:", halign=:left)
   slider_m1 = Slider(controls[2, 2:5], range=0.1:0.1:5.0, startvalue=m₁_obs[])
   connect!(m₁_obs, slider_m1.value)

   Label(controls[3, 1], "m₂:", halign=:left)
   slider_m2 = Slider(controls[3, 2:5], range=0.1:0.1:5.0, startvalue=m₂_obs[])
   connect!(m₂_obs, slider_m2.value)

   Label(controls[4, 1], "L₁:", halign=:left)
   slider_L1 = Slider(controls[4, 2:5], range=0.1:0.1:5.0, startvalue=L₁_obs[])
   connect!(L₁_obs, slider_L1.value)

   Label(controls[5, 1], "L₂:", halign=:left)
   slider_L2 = Slider(controls[5, 2:5], range=0.1:0.1:5.0, startvalue=L₂_obs[])
   connect!(L₂_obs, slider_L2.value)

   Label(controls[6, 1], "g:", halign=:left)
   slider_g = Slider(controls[6, 2:5], range=0.0:0.5:20.0, startvalue=g_obs[])
   connect!(g_obs, slider_g.value)

   Label(controls[7, 1], "Damping:", halign=:left)
   slider_damping = Slider(controls[7, 2:5], range=0.0:0.01:0.5, startvalue=damping_obs[])
   connect!(damping_obs, slider_damping.value)

   Label(controls[8, 1], "Speed:", halign=:left)
   slider_speed = Slider(controls[8, 2:5], range=0.1:0.1:3.0, startvalue=time_speed_obs[])
   connect!(time_speed_obs, slider_speed.value)

    # Connect projection dropdown to observable
    connect!(projection_choice, projection_dropdown.selection)

    is_running = Observable(false)

    # Function to create and solve the ODEProblem
    function solve_pendulum(m₁, m₂, L₁, L₂, g_val, damp_val)
        p = [g_val, m₁, m₂, L₁, L₂, damp_val]
        prob = ODEProblem(f!, copy(u₀), (0, 60), p)
        solve(prob, Tsit5(), abstol=1e-9, reltol=1e-9, saveat=0.02)
    end

    # Create an observable for the solution - MUST be defined before @lift uses it
    sol_obs = Observable(solve_pendulum(m₁_obs[], m₂_obs[], L₁_obs[], L₂_obs[], g_obs[], damping_obs[]))

    # Now we can create the pendulum line visualization that depends on sol_obs
    pend_line = @lift begin
        sol = $sol_obs
        current_idx = $idx
        L1 = $L₁_obs
        L2 = $L₂_obs

        θ₁_val = sol[1, current_idx]
        θ₂_val = sol[3, current_idx]

        x₁ = L1 * sin(θ₁_val)
        y₁ = -L1 * cos(θ₁_val)
        x₂ = x₁ + L2 * sin(θ₂_val)
        y₂ = y₁ - L2 * cos(θ₂_val)

        Point2f[(0, 0), (x₁, y₁), (x₂, y₂)]
    end

    # Draw trail first (so it's behind the pendulum)
    lines!(axPend, trail_pts, color=(:orange, 0.3), linewidth=1)

    # Draw pendulum rods and bobs
    lines!(axPend, pend_line, color=:gray80, linewidth=4)
    scatter!(axPend, @lift([$pend_line[2], $pend_line[3]]), color=:orange, markersize=15)

    # React to changes in parameters
    onany(m₁_obs, m₂_obs, L₁_obs, L₂_obs, g_obs, damping_obs) do m₁, m₂, L₁, L₂, g_val, damp_val
        sol_obs[] = solve_pendulum(m₁, m₂, L₁, L₂, g_val, damp_val)
        idx[] = 1  # Reset animation index when parameters change
        empty!(phase3d_trail[])  # Clear 3D trail when parameters change
    end

    on(run_btn.clicks) do _
        is_running[] = true
    end
    on(pause_btn.clicks) do _
        is_running[] = false
    end
    on(reset_btn.clicks) do _
        θ₁₀, ω₁₀, θ₂₀, ω₂₀ = π / 3, 0.0, π / 2, 0.0
        u₀ .= [θ₁₀, ω₁₀, θ₂₀, ω₂₀]
        sol_obs[] = solve_pendulum(m₁_obs[], m₂_obs[], L₁_obs[], L₂_obs[], g_obs[], damping_obs[])
        idx[] = 1
        empty!(trail_pts[])
        empty!(phase3d_trail[])
        @info "Simulation reset"
    end
    on(rand_btn.clicks) do _
        θ₁₀, ω₁₀, θ₂₀, ω₂₀ = (rand() * 2π - π), (rand() * 2π - π), (rand() * 2π - π), (rand() * 2π - π)
        u₀ .= [θ₁₀, ω₁₀, θ₂₀, ω₂₀]
        sol_obs[] = solve_pendulum(m₁_obs[], m₂_obs[], L₁_obs[], L₂_obs[], g_obs[], damping_obs[])
        idx[] = 1
        empty!(trail_pts[])
        empty!(phase3d_trail[])
        @info "Randomized IC" θ₁₀ θ₂₀ ω₁₀ ω₂₀
    end
    on(clear_btn.clicks) do _
        empty!(ax3d.scene.plots)
        empty!(trail_pts[])
        empty!(phase_trail[])
        empty!(phase3d_trail[])
        @info "Cleared 3D plot and trails"
    end

    # Click interaction to set initial conditions
    # Register on the specific axis instead of the whole figure for cleaner interaction
    register_interaction!(axPhase, :click_to_set_ic) do event::MouseEvent, axis
        if event.type === MouseEventTypes.leftclick
            # Get click position in data coordinates
            θ₁_new = event.data[1]
            ω₁_new = event.data[2]

            # Clamp to reasonable bounds
            θ₁_new = clamp(θ₁_new, -π, π)
            ω₁_new = clamp(ω₁_new, -10, 10)

            # Set new initial conditions (keep θ₂ and ω₂ the same)
            u₀[1] = θ₁_new
            u₀[2] = ω₁_new

            # Resolve and reset
            sol_obs[] = solve_pendulum(m₁_obs[], m₂_obs[], L₁_obs[], L₂_obs[], g_obs[], damping_obs[])
            idx[] = 1
            empty!(trail_pts[])
            empty!(phase_trail[])
            empty!(phase3d_trail[])

            @info "Set IC from phase space click" θ₁=θ₁_new ω₁=ω₁_new
            return Consume(true)
        end
        return Consume(false)
    end

    # Track previous energy for trend calculation
    prev_ke = Ref(0.0)
    prev_pe = Ref(0.0)
    initial_energy = Ref(0.0)
    last_frame_time = Ref(time())
    frame_count = Ref(0)
    fps_update_interval = 0.5  # Update FPS every 0.5 seconds

    # Screenshot task (every 10 seconds, keep last 50)
    screenshot_task = @async begin
        max_screenshots = 50
        screenshot_interval = 10.0
        last_screenshot_time = Ref(time())
        sleep(1.0)  # Give the window time to open

        while true  # Run indefinitely
            current_time = time()
            if current_time - last_screenshot_time[] >= screenshot_interval
                try
                    # Get all screenshot files
                    files = filter(f -> endswith(f, ".png"), readdir(screenshot_dir, join=true))
                    sort!(files, by=mtime)

                    # Remove oldest if we have too many
                    while length(files) >= max_screenshots
                        rm(files[1])
                        deleteat!(files, 1)
                    end

                    # Save new screenshot
                    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
                    filename = joinpath(screenshot_dir, "double_pendulum_$(timestamp).png")
                    save(filename, fig)
                    @info "Screenshot saved" filename
                    last_screenshot_time[] = current_time
                catch e
                    @warn "Screenshot failed" exception=e
                end
            end
            sleep(1.0)
        end
    end

    # Integration and animation loop
    @async begin
        sleep(0.5)  # Give the window time to open
        while true  # Run indefinitely - user can close window to stop
            if is_running[]
                sol = sol_obs[]
                N = length(sol.t)
                current_idx = idx[]
                idx[] = (current_idx % N) + 1

                # Update energy calculations
                θ₁, ω₁, θ₂, ω₂ = sol[:, current_idx]
                KE, PE = calculate_energy(θ₁, ω₁, θ₂, ω₂, m₁_obs[], m₂_obs[], L₁_obs[], L₂_obs[], g_obs[])
                total_E = KE + PE

                # Track initial energy (when starting or resetting)
                if current_idx == 1
                    initial_energy[] = total_E
                end

                # Calculate energy drift percentage
                if initial_energy[] > 0
                    energy_drift_obs[] = (total_E - initial_energy[]) / initial_energy[]
                end

                # Calculate trends
                ke_trend = KE > prev_ke[] + 0.01 ? "↑" : (KE < prev_ke[] - 0.01 ? "↓" : "")
                pe_trend = PE > prev_pe[] + 0.01 ? "↑" : (PE < prev_pe[] - 0.01 ? "↓" : "")

                # Update FPS
                frame_count[] += 1
                current_time = time()
                if current_time - last_frame_time[] >= fps_update_interval
                    fps_obs[] = frame_count[] / (current_time - last_frame_time[])
                    frame_count[] = 0
                    last_frame_time[] = current_time
                end

                # Update observables
                ke_val_obs[] = KE
                pe_val_obs[] = PE
                total_e_obs[] = total_E
                ke_fraction_obs[] = total_E > 0 ? KE / total_E : 0.5
                ke_trend_obs[] = ke_trend
                pe_trend_obs[] = pe_trend

                prev_ke[] = KE
                prev_pe[] = PE

                # Update trails and phase space
                L1 = L₁_obs[]
                L2 = L₂_obs[]
                x₁ = L1 * sin(θ₁)
                y₁ = -L1 * cos(θ₁)
                x₂ = x₁ + L2 * sin(θ₂)
                y₂ = y₁ - L2 * cos(θ₂)

                # Update second bob trail
                current_trail = trail_pts[]
                push!(current_trail, Point2f(x₂, y₂))
                if length(current_trail) > trail_length
                    popfirst!(current_trail)
                end
                trail_pts[] = current_trail

                # Update phase space current point
                phase_pt[] = Point2f(θ₁, ω₁)

                # Update phase space trail
                current_phase_trail = phase_trail[]
                push!(current_phase_trail, Point2f(θ₁, ω₁))
                if length(current_phase_trail) > 500  # Keep last 500 points
                    popfirst!(current_phase_trail)
                end
                phase_trail[] = current_phase_trail

                # Update 3D phase space with growing trail
                choice = projection_choice[]

                # Add current point to 3D trail
                current_3d_pt = if choice == "drop θ₁"
                    Point3f(ω₁, θ₂, ω₂)
                elseif choice == "drop ω₁"
                    Point3f(θ₁, θ₂, ω₂)
                elseif choice == "drop θ₂"
                    Point3f(θ₁, ω₁, ω₂)
                else
                    Point3f(θ₁, ω₁, θ₂)
                end

                current_3d_trail = phase3d_trail[]
                push!(current_3d_trail, current_3d_pt)
                if length(current_3d_trail) > 1000  # Keep last 1000 points
                    popfirst!(current_3d_trail)
                end
                phase3d_trail[] = current_3d_trail

                # Draw the growing trail with color gradient
                if length(current_3d_trail) > 1
                    n_pts = length(current_3d_trail)
                    colors = [RGBAf(i/n_pts, 0.3, 1.0 - i/n_pts, 0.6) for i in 1:n_pts]

                    empty!(ax3d.scene.plots)
                    lines!(ax3d, current_3d_trail, color=colors, linewidth=2)
                    scatter!(ax3d, [current_3d_trail[end]], color=:red, markersize=12)
                end
                sleep(0.016 / time_speed_obs[])  # Adjust sleep based on speed slider
            else
                sleep(0.05)
            end
        end
    end

    display(fig)
    return (
        fig = fig,
        ax3d = ax3d,
        axPend = axPend,
        pend_line = pend_line,
        projection_choice = projection_choice,
        controls = (
            run_btn = run_btn,
            pause_btn = pause_btn,
            reset_btn = reset_btn,
            rand_btn = rand_btn,
            clear_btn = clear_btn
        ),
        state = (
            is_running = is_running,
            u₀ = u₀
        )
    )
end
end # module DoublePendulum