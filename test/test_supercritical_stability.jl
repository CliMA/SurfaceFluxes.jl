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

    @testset "Supercritical stability with coupled canopy callbacks" begin
        # Regression test for solver robustness with T_sfc/q_vap_sfc callbacks
        # (ClimaLand canopy coupling). With a strong inversion (~10 K over 10 m)
        # and weak wind, Ri_b is supercritical. The callbacks additionally make
        # the state Ri_b *increase* with ζ (colder canopy-weighted T_sfc as
        # u_star drops), so the residual decreases toward a negative plateau.
        T_int_c = FT(294.673095703125)
        q_tot_int_c = FT(0.009545918211858238)
        ρ_int_c = FT(1.157759649975361)
        T_canopy = FT(284.6341213016535)
        q_canopy = FT(0.008910620278696149)
        Δz_c = FT(10)
        d_c = FT(0.0573)
        u_int_c = (FT(1.3673065900802612), FT(0))
        AI = FT(3.1223678081630233)
        leaf_Cd = FT(0.0726)
        g_stomata = FT(5.848035071201371e-6)

        config_c = SF.SurfaceFluxConfig(
            SF.ConstantRoughnessParams(FT(0.359), FT(0.0544)),
            SF.ConstantGustinessSpec(FT(1)),
        )

        # Canopy energy balance: T_sfc between T_int and T_canopy, weighted by
        # the ratio of canopy conductance g_land to aerodynamic conductance g_h
        function update_T_sfc(ζ, param_set, thermo_params, inputs, scheme,
            u_star, z0m, z0h)
            Φ_sfc = SF.surface_geopotential(inputs)
            Φ_int = SF.interior_geopotential(param_set, inputs)
            g_h = SF.heat_conductance(param_set, ζ, u_star, inputs, z0m, z0h, scheme)
            g_land = leaf_Cd * u_star * AI
            cp_d = SFP.cp_d(param_set)
            r = g_land / g_h
            return (inputs.T_int + inputs.T_sfc_guess * r + (Φ_int - Φ_sfc) / cp_d) /
                   (1 + r)
        end

        # Canopy moisture balance with stomatal + leaf boundary-layer conductance
        function update_q_vap_sfc(ζ, param_set, thermo_params, inputs, scheme,
            T_sfc, u_star, z0m, z0h)
            g_leaf = leaf_Cd * u_star * AI
            g_land = g_stomata * g_leaf / (g_leaf + g_stomata)
            g_h = SF.heat_conductance(param_set, ζ, u_star, inputs, z0m, z0h, scheme)
            q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
            r = g_land / g_h
            return (r * inputs.q_vap_sfc_guess + q_vap_int) / (1 + r)
        end

        for uf_type in (UF.BusingerParams, UF.GryanikParams),
            solver_opts in (
                nothing,  # default: forced fixed iterations (GPU mode)
                SF.SolverOptions{FT}(
                    tol = FT(1e-4),
                    rtol = FT(1e-4),
                    maxiter = 30,
                    forced_fixed_iters = false,
                ),
            )

            ps = SFP.SurfaceFluxesParameters(FT, uf_type)
            result = SF.surface_fluxes(
                ps,
                T_int_c, q_tot_int_c, FT(0), FT(0), ρ_int_c,
                T_canopy, q_canopy,
                FT(0), Δz_c, d_c,
                u_int_c, (FT(0), FT(0)),
                nothing,
                config_c,
                SF.PointValueScheme(),
                solver_opts,
                nothing,
                update_T_sfc,
                update_q_vap_sfc,
            )

            @test isfinite(result.shf)
            @test isfinite(result.lhf)
            @test isfinite(result.ustar)
            @test isfinite(result.ζ)

            # The stability regime must be preserved: strongly stable state
            # must not come back as unstable (ζ < 0). This case is
            # supercritical for both UFs, so the solve must deterministically
            # saturate at the stable branch limit in both solver modes.
            @test result.ζ == FT(100)
            @test result.converged == false

            # Near-decoupled surface: weak downward heat flux, weak friction
            @test result.shf < 0
            @test abs(result.shf) < FT(5)
            @test FT(0) < result.ustar < FT(0.05)

            # Final surface state must lie between the canopy and air states
            @test T_canopy <= result.T_sfc <= T_int_c
        end

        # With the tolerance-checked solver, the no-root (supercritical) case
        # must deterministically saturate at the stable limit
        ps_b = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
        result_sat = SF.surface_fluxes(
            ps_b,
            T_int_c, q_tot_int_c, FT(0), FT(0), ρ_int_c,
            T_canopy, q_canopy,
            FT(0), Δz_c, d_c,
            u_int_c, (FT(0), FT(0)),
            nothing,
            config_c,
            SF.PointValueScheme(),
            SF.SolverOptions{FT}(
                tol = FT(1e-4),
                rtol = FT(1e-4),
                maxiter = 30,
                forced_fixed_iters = false,
            ),
            nothing,
            update_T_sfc,
            update_q_vap_sfc,
        )
        @test result_sat.ζ == FT(100)
        @test result_sat.converged == false
    end

    @testset "Opposite-branch root with surface-state callbacks" begin
        # With callbacks, Ri_b_state varies with ζ, so the sign of the
        # residual at neutral, f(0) = -Ri_b_state(0), is only a heuristic
        # for the branch containing the root. This synthetic callback makes
        # the state weakly stable when evaluated at ζ >= 0 (f(0) < 0 =>
        # stable branch indicated, and supercritical there, so no stable
        # root) but strongly unstable when evaluated at ζ < 0, putting the
        # only root on the unstable branch. The ζ = ∓1 opposite-branch probe
        # in solve_stability_param must find it; committing to the indicated
        # branch alone would saturate at ζ = +100 with near-zero fluxes.
        T_int_o = FT(295)
        q_o = FT(0.005)
        ρ_o = FT(1.15)
        Δz_o = FT(10)
        grav = SFP.grav(param_set)
        cp_d = SFP.cp_d(param_set)
        T_neutral = T_int_o + grav * Δz_o / cp_d

        function update_T_sfc_flip(ζ, param_set, thermo_params, inputs,
            scheme, u_star, z0m, z0h)
            # Weakly stable for ζ >= 0 (growing supercritically), strongly
            # unstable for ζ < 0
            return ζ >= 0 ? T_neutral - FT(0.1) - 3 * ζ / (1 + ζ) :
                   T_neutral - FT(0.1) + 5 * (-ζ) / (1 - ζ)
        end

        config_o = SF.SurfaceFluxConfig(
            SF.ConstantRoughnessParams(FT(0.01), FT(0.001)),
            SF.ConstantGustinessSpec(FT(1)),
        )

        for solver_opts in (
            nothing,  # default: forced fixed iterations (GPU mode)
            SF.SolverOptions{FT}(
                tol = FT(1e-4),
                rtol = FT(1e-4),
                maxiter = 30,
                forced_fixed_iters = false,
            ),
        )
            result = SF.surface_fluxes(
                param_set,
                T_int_o, q_o, FT(0), FT(0), ρ_o,
                T_neutral - FT(0.1), q_o,
                FT(0), Δz_o, FT(0),
                (FT(0.5), FT(0)), (FT(0), FT(0)),
                nothing,
                config_o,
                SF.PointValueScheme(),
                solver_opts,
                nothing,
                update_T_sfc_flip,
                nothing,
            )

            # Root found on the unstable branch, near neutral
            @test FT(-1) < result.ζ < FT(0)
            # Warm surface: upward sensible heat flux
            @test result.shf > 0
            @test isfinite(result.ustar) && result.ustar > 0
        end

        # In tolerance-checked mode the bracketed root must report converged
        result_conv = SF.surface_fluxes(
            param_set,
            T_int_o, q_o, FT(0), FT(0), ρ_o,
            T_neutral - FT(0.1), q_o,
            FT(0), Δz_o, FT(0),
            (FT(0.5), FT(0)), (FT(0), FT(0)),
            nothing,
            config_o,
            SF.PointValueScheme(),
            SF.SolverOptions{FT}(
                tol = FT(1e-4),
                rtol = FT(1e-4),
                maxiter = 30,
                forced_fixed_iters = false,
            ),
            nothing,
            update_T_sfc_flip,
            nothing,
        )
        @test result_conv.converged == true
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
