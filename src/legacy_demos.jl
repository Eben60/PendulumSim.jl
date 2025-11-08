function run()
    # pendulum params
    g = 9.81
    L = 1.0
    β = 0.02          # small viscous damping (s⁻¹)
    θ0 = deg2rad(120) # initial angle
    ω0 = 0.0

    function pend!(du, u, p, t)
        θ, ω = u
        du[1] = ω
        du[2] = -(g / L) * sin(θ) - β * ω
    end

    u0 = [θ0, ω0]
    tspan = (0.0, 10.0)
    sol = solve(ODEProblem(pend!, u0, tspan), Tsit5(); reltol=1e-9, abstol=1e-9)

    # exact large-angle period for comparison (β≈0)
    T0 = 2π * sqrt(L / g)
    Tθ = T0 * (2 / π) * elliptic_k(sin(θ0 / 2)^2)  # Elliptic.K(m) with m = sin^2(θ0/2)

    # quick animation
    f = Figure(backgroundcolor=:black)
    ax = Axis(f[1, 1], aspect=1, xlabel="x", ylabel="y",
        xgridvisible=false, ygridvisible=false, xticklabelsvisible=false, yticklabelsvisible=false,
        xticksvisible=false, yticksvisible=false, limits=(-1.2, 1.2, -1.2, 0.2))

    rod = lines!(ax, [0, 0], [0, 0], color=:gray80, linewidth=4)
    bob = scatter!(ax, [0.0], [0.0], color=:tomato, markersize=20)

    xs = @lift @. L * sin(sol($t, idxs=1))
    ys = @lift @. -L * cos(sol($t, idxs=1))

    record(f, "pendulum.mp4", range(0, stop=sol.t[end], length=600)) do tt
        θ = sol(tt, idxs=1)
        x = L * sin(θ)
        y = -L * cos(θ)
        rod[1] = ([0, x], [0, y])
        bob[1] = ([x], [y])
    end

    @info("Predicted large-angle period T(θ0) ≈ $(round(Tθ,digits=3)) s")
end

using RigidBodyDynamics, MeshCat, MeshCatMechanisms, StaticArrays
using GeometryBasics: Point3f, Cylinder  # <-- add this import

using RigidBodyDynamics.Spatial: ⊕
using Rotations: RotMatrix
using Colors: RGB
function run_rbd()

    # ----- build the link inertia (same as you had) -----
    L = 1.0
    m_rod, m_bob = 0.10, 0.80
    frame = CartesianFrame3D("link")

    Irod_com = Diagonal(SVector((1 / 12) * m_rod * L^2, 0.0, (1 / 12) * m_rod * L^2))
    I_rod = SpatialInertia(frame; mass=m_rod,
        com=SVector(0.0, -L / 2, 0.0),
        moment_about_com=Matrix(Irod_com))

    I_bob = SpatialInertia(frame; mass=m_bob,
        com=SVector(0.0, -L, 0.0),
        moment_about_com=zeros(3, 3))

    I_link = I_rod + I_bob
    link = RigidBody("link", I_link)

    # ----- mechanism & joint -----
    gvec = SVector(0.0, -9.81, 0.0)   # gravity vector in world frame

    mech = Mechanism(RigidBody{Float64}("world"); gravity=gvec)
    ground = root_body(mech)

    axis = SVector(0.0, 0.0, 1.0)
    joint = Joint("hinge", Revolute{Float64}(axis))

    # link already constructed with I_link:
    # link = RigidBody("link", I_link)

    Reye = one(RotMatrix{3,Float64})
    pzero = SVector(0.0, 0.0, 0.0)

    attach!(mech, ground, link, joint;
        # from = joint's BEFORE frame, to = a frame fixed on the parent (ground)
        joint_pose=Transform3D(frame_before(joint), default_frame(ground), Reye, pzero),

        # from = a frame fixed on the child (link), to = joint's AFTER frame
        successor_pose=Transform3D(default_frame(link), frame_after(joint), Reye, pzero)
    )

    # Sanity: set initial state and gravity
    state = MechanismState(mech)
    set_configuration!(state, joint, deg2rad(120.0))
    set_velocity!(state, joint, 0.0)
    # visualize
    # ----- MeshCat visualizer (NOTE: visualize the MECHANISM, not `world`) -----
    vis = Visualizer()
    open(vis)  # opens in browser
    mvis = MechanismVisualizer(mech, vis)
    # geometry (rod + bob) in the link's body frame:
    r = 0.006f0
    rod = Cylinder(Point3f(0, 0, 0), Point3f(0, -Float32(L), 0), r)
    bob = Sphere(Point3f(0, -Float32(L), 0), Float32(0.03))

    # materials come from MeshCat:
    rodmat = MeshCat.MeshPhongMaterial(color=RGB(0.7, 0.7, 0.75))
    bobmat = MeshCat.MeshPhongMaterial(color=RGB(0.9, 0.3, 0.2))

    # ✅ attach to the *frame* of the body, not the body object; no `color` keyword
    setelement!(mvis, default_frame(link), rod, rodmat, "rod")
    setelement!(mvis, default_frame(link), bob, bobmat, "bob")

    GLMakie.activate!()
    # ---------- plotting UI ----------
    set_theme!(theme_dark())

    fig = Figure(resolution=(1000, 560))
    axθ = Axis(fig[1, 1], title="θ(t)", xlabel="t [s]", ylabel="θ [rad]")
    axph = Axis(fig[1, 2], title="Phase space", xlabel="θ [rad]", ylabel="ω [rad/s]")

    # Controls
    run_toggle = Toggle(fig[2, 1], active=false)
    run_label = Label(fig[2, 1], lift(x -> x ? "Pause" : "Run", run_toggle.active))
    rand_btn = Button(fig[2, 2], label="Random IC")
    θmax_sld = Slider(fig[3, 1], range=0:1:170, startvalue=120)     # degrees
    ωmax_sld = Slider(fig[3, 2], range=0:0.1:5.0, startvalue=0.5)   # rad/s

    # Data buffers
    ts = Float32[]
    thetas = Float32[]
    omegas = Float32[]
    tobs = Observable(copy(ts))
    thob = Observable(copy(thetas))
    omob = Observable(copy(omegas))

    l1 = lines!(axθ, tobs, thob, color=:cyan, linewidth=2)
    l2 = lines!(axph, thob, omob, color=:tomato, linewidth=2)

    display(fig)

    # ---------- helpers ----------
    Δt = 0.002f0
    τ = zeros(eltype(configuration(state)), num_velocities(mech))
    result = DynamicsResult(mech)
    t = 0.0f0

    # (re)set state from given (θ0, ω0)
    function set_ic!(θ0::Float64, ω0::Float64)
        set_configuration!(state, joint, θ0)
        set_velocity!(state, joint, ω0)
        empty!(ts)
        empty!(thetas)
        empty!(omegas)
        t = 0.0f0
        tobs[] = ts
        thob[] = thetas
        omob[] = omegas
        return nothing
    end

    # randomize on button click
    on(rand_btn.clicks) do _
        θmax = deg2rad(θmax_sld.value[])      # radians
        ωmax = ωmax_sld.value[]
        θ0 = rand() * 2θmax - θmax              # uniform in [-θmax, θmax]
        ω0 = rand() * 2ωmax - ωmax              # uniform in [-ωmax, ωmax]
        set_ic!(θ0, ω0)
    end

    # ---------- simulation task (runs while toggle is on) ----------
    @async begin
        while isopen(fig.scene)
            if run_toggle.active[]
                dynamics!(result, state, τ)
                v = velocity(state)
                v .= v .+ Δt .* result.v̇
                q = configuration(state)
                q .= q .+ Δt .* v
                set_velocity!(state, v)
                set_configuration!(state, q)

                # log a point
                t += Δt
                push!(ts, t)
                push!(thetas, Float32(q[1]))
                push!(omegas, Float32(v[1]))

                # update plots every few frames
                if length(ts) % 5 == 0
                    tobs[] = ts
                    thob[] = thetas
                    omob[] = omegas
                end

                # (optional) keep MeshCat in sync if you created mvis
                # MeshCatMechanisms.set_configuration!(mvis, q)
            end
            sleep(Δt)  # comment out to run faster than real time
        end
    end

    # start with a random IC once:
    notify(rand_btn.clicks)
end

function set_meshcat_camera!(vis; eye=SVector(1.8, 0.8, 1.2), target=SVector(0, -0.5, 0), up=SVector(0, 0, 1))
    # Build a look-at pose: camera at `eye`, looking at `target`, with `up` axis
    z = normalize(target - eye)                # forward
    x = normalize(cross(z, up))                # right
    y = cross(x, z)                            # true up
    R = RotMatrix(SMatrix{3,3,Float64}([x -y z]))  # MeshCat expects camera -Z forward
    # Position node (translation)
    MeshCat.settransform!(vis["/Cameras/default"], MeshCat.Translation(eye[1], eye[2], eye[3]))
    # Orientation node (rotation)
    MeshCat.settransform!(vis["/Cameras/default/rotated"], MeshCat.LinearMap(R))
    # Optional: tweak camera props
    MeshCat.setprop!(vis["/Cameras/default"], "near", 0.01)
    MeshCat.setprop!(vis["/Cameras/default"], "far", 100.0)
    return nothing
end