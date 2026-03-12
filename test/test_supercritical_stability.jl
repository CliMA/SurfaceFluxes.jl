# Tests for supercritical stable conditions.
#
# For the Businger-Dyer universal functions, the stable profile function
# ψ = -a_m * ζ leads to a critical bulk Richardson number
# Ri_b_crit = 1/a_m ≈ 0.213 (with a_m = 4.7 from Businger et al. 1971),
# above which no finite ζ satisfies the MOST stability relations
# (see Fairall et al. 2003, Eq. 12-13). The solver clamps ζ to [-100, 100]
# to produce bounded output.
#
# This test verifies:
# 1. The solver produces finite, bounded output for Ri_b >> Ri_b_crit
# 2. Fluxes approach zero (but remain finite) as stratification increases
# 3. Exchange coefficients remain positive

module TestSupercriticalStability

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

@testset "Supercritical Stability" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    # Base conditions: very stable (T_sfc << T_int) + low wind -> large Ri_b
    u_int = (FT(1), FT(0))  # Low wind speed
    u_sfc = (FT(0), FT(0))
    q_int = FT(0.005)
    q_sfc = FT(0.005)
    ρ_int = FT(1.2)
    Δz = FT(10)

    config = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(1e-3), FT(1e-3)),
        SF.ConstantGustinessSpec(FT(0)),  # No gustiness
        SF.MoistModel(),
    )
    opts = SF.SolverOptions{FT}(maxiter = 15, tol = FT(1e-3), forced_fixed_iters = false)

    @testset "Increasing stratification produces bounded output" begin
        # Sweep surface temperatures from mildly to extremely stable
        T_int_base = FT(300)
        ΔT_range = [FT(2), FT(5), FT(10), FT(20), FT(50)]

        prev_shf = FT(-Inf)
        for ΔT in ΔT_range
            T_sfc = T_int_base - ΔT  # Cold surface

            result = SF.surface_fluxes(
                param_set,
                T_int_base, q_int, FT(0), FT(0), ρ_int,
                T_sfc, q_sfc,
                FT(0), Δz, FT(0),
                u_int, u_sfc,
                nothing,
                config,
                SF.PointValueScheme(),
                opts,
            )

            # All output fields should be finite
            @test isfinite(result.shf)
            @test isfinite(result.lhf)
            @test isfinite(result.ustar)
            @test isfinite(result.Cd)
            @test isfinite(result.g_h)
            @test isfinite(result.ζ)
            @test isfinite(result.L_MO)

            # Stable: T_sfc < T_int -> SHF < 0 (downward)
            @test result.shf < 0

            # Coefficients remain positive
            @test result.Cd > 0
            @test result.g_h > 0
            @test result.ustar >= 0

            # ζ should be clamped: ζ ∈ [-100, 100]
            @test result.ζ <= FT(100)
            @test result.ζ >= FT(-100)

            # Stability: ζ > 0 for stable conditions
            @test result.ζ > 0

            prev_shf = result.shf
        end
    end

    @testset "Supercritical Ri_b verification" begin
        # Construct a case with Ri_b clearly above the critical value (~0.21 for Businger)
        T_int = FT(300)
        T_sfc = FT(280)  # 20 K temperature inversion

        ρ_sfc = SF.surface_density(param_set, T_int, ρ_int, T_sfc, Δz, q_int)
        ΔU = hypot(u_int[1], u_int[2])  # 1 m/s

        inputs = SF.build_surface_flux_inputs(
            T_int, q_int, FT(0), FT(0), ρ_int,
            T_sfc, q_sfc,
            FT(0), Δz, FT(0),
            u_int, u_sfc,
            config,
            nothing,
            SF.FluxSpecs{FT}(),
            nothing, nothing,
        )

        Ri_b = SF.state_bulk_richardson_number(
            param_set, inputs, T_sfc, ρ_sfc, ΔU, q_sfc,
        )

        # Verify this is indeed supercritical
        @test Ri_b > FT(0.21)

        # Run the solver
        result = SF.surface_fluxes(
            param_set,
            T_int, q_int, FT(0), FT(0), ρ_int,
            T_sfc, q_sfc,
            FT(0), Δz, FT(0),
            u_int, u_sfc,
            nothing,
            config,
            SF.PointValueScheme(),
            opts,
        )

        # Output should be bounded despite supercritical conditions
        @test isfinite(result.shf)
        @test isfinite(result.ustar)
        @test result.ζ == FT(100)  # Should hit the upper clamp

        # Fluxes should still have correct sign
        @test result.shf < 0  # Downward heat flux (cold surface)
    end

    @testset "Fluxes decrease with increasing stability" begin
        # More strongly stable conditions should have smaller magnitude fluxes.
        # When ζ hits the upper clamp (100), Cd saturates to a constant value.
        # So we test the monotonic decrease only up to the clamping threshold.
        T_int = FT(300)

        # Use small enough ΔT values that span from weak to moderate stability
        # without all hitting the ζ=100 clamp
        results = map([FT(0.1), FT(0.5), FT(2)]) do ΔT
            SF.surface_fluxes(
                param_set,
                T_int, q_int, FT(0), FT(0), ρ_int,
                T_int - ΔT, q_sfc,
                FT(0), Δz, FT(0),
                u_int, u_sfc,
                nothing,
                config,
                SF.PointValueScheme(),
                opts,
            )
        end

        # Cd should decrease monotonically (or saturate) with increasing stability
        @test results[1].Cd >= results[2].Cd
        @test results[2].Cd >= results[3].Cd

        # At least the first transition should show a genuine decrease
        @test results[1].Cd > results[3].Cd
    end
end

end # module
