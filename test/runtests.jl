using PendulumSim
using ReTest
using DifferentialEquations

@testset "PendulumSim.jl" begin
    @testset "Module structure" begin
        # Test that submodules are accessible
        @test isdefined(PendulumSim, :SinglePendulum)
        @test isdefined(PendulumSim, :DoublePendulum)

        # Test that main functions are exported
        @test isdefined(PendulumSim, :run_single_pendulum_gui)
        @test isdefined(PendulumSim, :run_double_pendulum_gui)

        # Test that submodule functions exist
        @test isdefined(PendulumSim.SinglePendulum, :run_gui)
        @test isdefined(PendulumSim.DoublePendulum, :run_double_pendulum_gui)
    end

    @testset "Double Pendulum Physics" begin
        # Test energy calculation
        θ₁, ω₁, θ₂, ω₂ = π/4, 0.0, π/3, 0.0
        m₁, m₂ = 1.0, 1.0
        L₁, L₂ = 1.0, 1.0
        g = 9.81

        KE, PE = PendulumSim.DoublePendulum.calculate_energy(θ₁, ω₁, θ₂, ω₂, m₁, m₂, L₁, L₂, g)

        # At rest (ω₁ = ω₂ = 0), kinetic energy should be zero
        @test KE ≈ 0.0 atol=1e-10

        # Potential energy should be negative (below reference)
        @test PE < 0.0

        # Total energy should be conserved (equal to PE at rest)
        E_total = KE + PE
        @test E_total ≈ PE atol=1e-10

        # Test with motion
        θ₁, ω₁, θ₂, ω₂ = π/4, 1.0, π/3, 0.5
        KE2, PE2 = PendulumSim.DoublePendulum.calculate_energy(θ₁, ω₁, θ₂, ω₂, m₁, m₂, L₁, L₂, g)

        # With motion, kinetic energy should be positive
        @test KE2 > 0.0

        # Total energy should be higher with motion
        @test (KE2 + PE2) > E_total
    end

    @testset "Double Pendulum Equations of Motion" begin
        # Test that the ODE system can be constructed and solved
        m₁, m₂ = 1.0, 1.0
        L₁, L₂ = 1.0, 1.0
        g = 9.81
        damping = 0.0

        # Initial conditions: small angles
        u₀ = [0.1, 0.0, 0.1, 0.0]  # θ₁, ω₁, θ₂, ω₂
        tspan = (0.0, 1.0)
        p = [g, m₁, m₂, L₁, L₂, damping]

        # Create and solve the problem
        f! = PendulumSim.DoublePendulum.double_pendulum_ode!
        prob = ODEProblem(f!, u₀, tspan, p)
        sol = solve(prob, Tsit5(), abstol=1e-9, reltol=1e-9)

        # Check that solution was successful
        @test sol.retcode == :Success

        # Check that solution has reasonable values (angles don't explode)
        @test all(abs.(sol[1, :]) .< 10.0)  # θ₁ stays bounded
        @test all(abs.(sol[3, :]) .< 10.0)  # θ₂ stays bounded

        # For small angles with no damping, energy should be approximately conserved
        KE₀, PE₀ = PendulumSim.DoublePendulum.calculate_energy(u₀[1], u₀[2], u₀[3], u₀[4], m₁, m₂, L₁, L₂, g)
        E₀ = KE₀ + PE₀

        u_final = sol.u[end]
        KE_f, PE_f = PendulumSim.DoublePendulum.calculate_energy(u_final[1], u_final[2], u_final[3], u_final[4], m₁, m₂, L₁, L₂, g)
        E_f = KE_f + PE_f

        # Energy drift should be small (< 1% for short integration with tight tolerances)
        energy_drift = abs((E_f - E₀) / E₀)
        @test energy_drift < 0.01
    end

    @testset "Double Pendulum Damping" begin
        # Test that damping reduces energy
        m₁, m₂ = 1.0, 1.0
        L₁, L₂ = 1.0, 1.0
        g = 9.81
        damping = 0.1

        u₀ = [π/4, 0.5, π/3, 0.3]  # Start with some motion
        tspan = (0.0, 5.0)
        p = [g, m₁, m₂, L₁, L₂, damping]

        f! = PendulumSim.DoublePendulum.double_pendulum_ode!
        prob = ODEProblem(f!, u₀, tspan, p)
        sol = solve(prob, Tsit5(), abstol=1e-9, reltol=1e-9)

        # Calculate energy at start and end
        KE₀, PE₀ = PendulumSim.DoublePendulum.calculate_energy(u₀[1], u₀[2], u₀[3], u₀[4], m₁, m₂, L₁, L₂, g)
        E₀ = KE₀ + PE₀

        u_final = sol.u[end]
        KE_f, PE_f = PendulumSim.DoublePendulum.calculate_energy(u_final[1], u_final[2], u_final[3], u_final[4], m₁, m₂, L₁, L₂, g)
        E_f = KE_f + PE_f

        # With damping, final energy should be less than initial energy
        @test E_f < E₀

        # Angular velocities should decrease
        @test abs(u_final[2]) < abs(u₀[2])  # ω₁ decreased
        @test abs(u_final[4]) < abs(u₀[4])  # ω₂ decreased
    end

    @testset "Energy calculation edge cases" begin
        m₁, m₂ = 1.0, 1.0
        L₁, L₂ = 1.0, 1.0
        g = 9.81

        # Test vertical down position (θ = 0)
        KE, PE = PendulumSim.DoublePendulum.calculate_energy(0.0, 0.0, 0.0, 0.0, m₁, m₂, L₁, L₂, g)
        @test KE ≈ 0.0 atol=1e-10
        # PE should be at minimum (most negative)
        expected_PE = -m₁ * g * L₁ - m₂ * g * (L₁ + L₂)
        @test PE ≈ expected_PE atol=1e-10

        # Test horizontal position (θ = π/2)
        KE, PE = PendulumSim.DoublePendulum.calculate_energy(π/2, 0.0, π/2, 0.0, m₁, m₂, L₁, L₂, g)
        @test KE ≈ 0.0 atol=1e-10
        # At horizontal, PE should be higher than vertical
        @test PE > expected_PE
    end
end
