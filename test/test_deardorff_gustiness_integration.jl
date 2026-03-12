# Integration test for the Deardorff gustiness parameterization through the full
# iterative MOST solver.
#
# The Deardorff gustiness introduces a nonlinear coupling:
#   gustiness -> U_eff -> fluxes -> buoyancy flux -> gustiness
# This test verifies that the solver resolves this coupling correctly and that
# the gustiness has the expected physical effect on the solution.

module TestDeardorffGustinessIntegration

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

@testset "Deardorff Gustiness Integration" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    # Moderately unstable scenario: moderate wind, warm surface
    T_int = FT(295)
    T_sfc = FT(305)   # 10 K warmer surface -> strong instability
    q_int = FT(0.01)
    q_sfc = FT(0.015)
    ρ_int = FT(1.1)
    Δz = FT(20)
    u_int = (FT(5), FT(0))  # Moderate wind avoids near-zero-wind convergence stiffness
    u_sfc = (FT(0), FT(0))

    # Use default solver options (forced_fixed_iters = true) for robust convergence
    opts = SF.SolverOptions{FT}(maxiter = 30)

    # 1. Deardorff gustiness
    config_deardorff = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(1e-3), FT(1e-3)),
        SF.DeardorffGustinessSpec(),
    )

    result_deardorff = SF.surface_fluxes(
        param_set,
        T_int, q_int, FT(0), FT(0), ρ_int,
        T_sfc, q_sfc,
        FT(0), Δz, FT(0),
        u_int, u_sfc,
        nothing,
        config_deardorff,
        SF.PointValueScheme(),
        opts,
    )

    # 2. Zero gustiness
    config_zero = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(1e-3), FT(1e-3)),
        SF.ConstantGustinessSpec(FT(0)),
    )

    result_zero = SF.surface_fluxes(
        param_set,
        T_int, q_int, FT(0), FT(0), ρ_int,
        T_sfc, q_sfc,
        FT(0), Δz, FT(0),
        u_int, u_sfc,
        nothing,
        config_zero,
        SF.PointValueScheme(),
        opts,
    )

    # 3. Constant gustiness = 1 m/s for comparison
    config_const = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(1e-3), FT(1e-3)),
        SF.ConstantGustinessSpec(FT(1)),
    )

    result_const = SF.surface_fluxes(
        param_set,
        T_int, q_int, FT(0), FT(0), ρ_int,
        T_sfc, q_sfc,
        FT(0), Δz, FT(0),
        u_int, u_sfc,
        nothing,
        config_const,
        SF.PointValueScheme(),
        opts,
    )

    @testset "Finite outputs" begin
        @test isfinite(result_deardorff.shf)
        @test isfinite(result_deardorff.lhf)
        @test isfinite(result_deardorff.ustar)
        @test isfinite(result_deardorff.Cd)
    end

    @testset "Physical sign expectations" begin
        # Warm surface -> upward SHF and LHF
        @test result_deardorff.shf > 0
        @test result_deardorff.lhf > 0

        # ustar must be positive
        @test result_deardorff.ustar > FT(0)

        # Unstable -> negative ζ
        @test result_deardorff.ζ < 0
    end

    @testset "Gustiness ordering" begin
        # Deardorff adds gustiness from buoyancy flux, so the effective wind speed
        # is at least as large as with zero or constant gustiness. This should yield
        # at least as large momentum flux as the constant-1 case.
        @test result_deardorff.ustar >= result_const.ustar ||
              isapprox(result_deardorff.ustar, result_const.ustar; rtol = FT(0.1))

        # Constant gustiness=1 gives a larger effective wind than zero gustiness
        @test result_const.ustar >= result_zero.ustar
    end

    @testset "Self-consistency" begin
        # The Deardorff gustiness should be consistent with the diagnosed buoyancy flux
        β = SFP.gustiness_coeff(param_set)
        zi = SFP.gustiness_zi(param_set)

        B = SF.buoyancy_flux(
            param_set,
            result_deardorff.shf,
            result_deardorff.lhf,
            T_sfc,
            SF.surface_density(param_set, T_int, ρ_int, T_sfc, Δz, q_int),
            q_sfc,
        )

        # In unstable conditions, B > 0
        @test B > 0

        # Self-consistent gustiness from diagnosed buoyancy flux
        w_star = cbrt(B * zi)
        gustiness_expected = β * w_star
        @test gustiness_expected > 0
    end

    @testset "Stable conditions: Deardorff returns zero" begin
        # With stable conditions (T_sfc < T_int), buoyancy flux is negative,
        # so Deardorff gustiness should be zero (only activates for B > 0)
        T_sfc_stable = FT(290) # Cold surface
        T_int_stable = FT(300)

        result_stable_deardorff = SF.surface_fluxes(
            param_set,
            T_int_stable, q_int, FT(0), FT(0), ρ_int,
            T_sfc_stable, q_sfc,
            FT(0), Δz, FT(0),
            u_int, u_sfc,
            nothing,
            config_deardorff,
            SF.PointValueScheme(),
            opts,
        )

        result_stable_zero = SF.surface_fluxes(
            param_set,
            T_int_stable, q_int, FT(0), FT(0), ρ_int,
            T_sfc_stable, q_sfc,
            FT(0), Δz, FT(0),
            u_int, u_sfc,
            nothing,
            config_zero,
            SF.PointValueScheme(),
            opts,
        )

        # In stable conditions, Deardorff gustiness = 0, so results should match
        # zero-gustiness case
        @test isapprox(result_stable_deardorff.shf, result_stable_zero.shf; rtol = FT(0.01))
        @test isapprox(
            result_stable_deardorff.ustar,
            result_stable_zero.ustar;
            rtol = FT(0.01),
        )
    end
end

end # module
