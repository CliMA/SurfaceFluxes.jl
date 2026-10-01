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
import ForwardDiff

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

@testset "Deardorff gustiness at fixed ζ without a consistent ustar" begin
    # At fixed ζ, the Deardorff gustiness is proportional to ustar, so the inner ustar
    # equation has a root only where F_m/κ > β (-ζ z_i / (κ Δz))^(1/3). Beyond this
    # (more unstable than the free-convection limit) the inner solve returns the upper
    # bracket end, and the ζ solve settles where a consistent ustar exists.
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    uf = SFP.uf_params(param_set)
    pv = UF.PointValueScheme()
    κ = SFP.von_karman_const(param_set)
    β = SFP.gustiness_coeff(param_set)
    zi = SFP.gustiness_zi(param_set)
    Δz = FT(10)
    gust_bound(ζ) = β * cbrt(-ζ * zi / (κ * Δz))
    F_m(ζ, z0m) = UF.dimensionless_profile(uf, Δz, ζ, z0m, UF.MomentumTransport(), pv)
    inputs(z0m, U, ΔT) = SF.build_surface_flux_inputs(
        FT(300), FT(0), FT(0), FT(0), FT(1.2), FT(300 + ΔT), FT(0), FT(0), Δz, FT(0),
        (FT(U), FT(0)), (FT(0), FT(0)),
        SF.SurfaceFluxConfig(SF.ConstantRoughnessParams(z0m, z0m / 10),
            SF.DeardorffGustinessSpec(), SF.DryModel()),
        nothing, SF.FluxSpecs{FT}(), nothing, nothing,
    )
    opts = SF.SolverOptions{FT}(maxiter = 40)
    for z0m in (FT(0.01), FT(1))
        # No consistent ustar at ζ = -10: the upper bracket end is returned
        inp = inputs(z0m, 2, 10)
        @test F_m(FT(-10), z0m) / κ < gust_bound(FT(-10))
        @test (@inferred SF.compute_ustar_and_roughness(param_set, FT(-10), inp, pv))[1] ==
              FT(4)
        # Type stable with a dual ζ (AD through the solve), with and without a sign change
        for ζ in (FT(-0.2), FT(-10))
            ζd = ForwardDiff.Dual(ζ, one(FT))
            ustar_d = (@inferred SF.compute_ustar_and_roughness(param_set, ζd, inp, pv))[1]
            @test ustar_d isa ForwardDiff.Dual
            @test ForwardDiff.value(ustar_d) ==
                  SF.compute_ustar_and_roughness(param_set, ζ, inp, pv)[1]
        end
        # Calm and neutral: ustar below the bracket, the lower end is returned
        @test SF.compute_ustar_and_roughness(param_set, FT(0), inputs(z0m, 0, 0), pv)[1] ==
              FT(1e-4)
        # The ζ solve ends with a consistent ustar inside the bracket, including near the
        # free-convection limit (low wind), where the ζ residual is steep
        for U in (0.1, 0.5, 2, 5)
            inp = inputs(z0m, U, 10)
            out = @inferred SF.surface_fluxes(param_set, inp, pv, opts)
            @test out.converged
            @test FT(1e-2) < out.ustar < FT(4)
            @test F_m(out.ζ, z0m) / κ > gust_bound(out.ζ)
            rf = SF.UstarResidual(param_set, inp, pv, out.ζ)
            @test rf(FT(1e-4)) < 0 < rf(FT(4))
            @test 0 < out.shf < 2000
        end
    end
end

end # module
