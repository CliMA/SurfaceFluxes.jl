# Stability caps: `max_heat_flux_stability` agrees with brute-force maximization of
# ζ / F_m(ζ)^3 and with the analytical point-value condition; beyond a cap, the exchange
# coefficients equal those at the cap, the returned ζ is the Obukhov value implied by
# the fluxes, the solve converges in supercritical conditions, and the sensible heat
# flux increases monotonically with the surface–air temperature difference.

module TestStabilityCap

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaParams as CP
import ForwardDiff

function solve(param_set, config, T_int, T_sfc, U, Δz, d; FT = eltype(param_set),
    update_T_sfc = nothing)
    SF.surface_fluxes(
        param_set,
        FT(T_int), FT(0.005), FT(0), FT(0), FT(1.2),
        FT(T_sfc), FT(0.005),
        FT(0), FT(Δz), FT(d),
        (FT(U), FT(0)), (FT(0), FT(0)),
        nothing,
        config,
        SF.PointValueScheme(),
        SF.SolverOptions{FT}(
            maxiter = 40,
            tol = FT(1e-6),
            rtol = FT(1e-6),
            forced_fixed_iters = false,
        ),
        nothing,
        update_T_sfc,
        nothing,
    )
end

@testset "Stability cap" begin
    for FT in (Float32, Float64)
        param_set = SFP.SurfaceFluxesParameters(FT, UF.GryanikParams)
        uf = SFP.uf_params(param_set)
        tm = UF.MomentumTransport()
        pv = SF.PointValueScheme()

        @testset "ζ_p maximizes the heat flux at fixed wind ($FT)" begin
            for (Δz_eff, z0m) in ((12.0, 2.0), (2.0, 0.01), (10.0, 1e-3), (30.0, 1e-4))
                ζ_p = SF.max_heat_flux_stability(param_set, FT(Δz_eff), FT(z0m))
                @test ζ_p isa FT
                # Type-stable (required for GPU kernels)
                @inferred SF.max_heat_flux_stability(param_set, FT(Δz_eff), FT(z0m), pv,
                    SF.ExponentialRSL(FT; z_RSL = 20.0))
                Fm(ζ) = UF.dimensionless_profile(uf, Δz_eff, ζ, z0m, tm, pv)
                ζs = exp.(range(log(1e-2), log(20), length = 100_000))
                ζ_b = ζs[argmax([ζ / Fm(ζ)^3 for ζ in ζs])]
                @test isapprox(ζ_p, ζ_b; rtol = 5e-3)
                # Analytical condition for point values
                G =
                    Fm(Float64(ζ_p)) -
                    3 * (UF.phi(uf, Float64(ζ_p), tm) -
                         UF.phi(uf, ζ_p * z0m / Δz_eff, tm))
                @test abs(G) < 5e-3 * Fm(Float64(ζ_p))
            end
            # ζ_p increases with Δz_eff / z0m
            @test SF.max_heat_flux_stability(param_set, FT(12), FT(2)) <
                  SF.max_heat_flux_stability(param_set, FT(2), FT(0.01)) <
                  SF.max_heat_flux_stability(param_set, FT(30), FT(1e-4))
            # Ranges stated in the MaxHeatFluxStabilityCap docstring (Gryanik functions)
            for (r_lo, r_hi, ζ_lo, ζ_hi) in ((2, 10, 0.14, 0.3), (30, 1000, 0.38, 0.85),
                (3000, 1e4, 0.99, 1.2), (3e4, 1e5, 1.37, 1.61))
                ζs = [
                    SF.max_heat_flux_stability(param_set, FT(10), FT(10 / r))
                    for r in (r_lo, r_hi)
                ]
                @test ζ_lo <= ζs[1] < ζs[2] <= ζ_hi
            end
            for r in (2, 1e5)
                ζ_p = SF.max_heat_flux_stability(param_set, FT(10), FT(10 / r))
                R_f = ζ_p / UF.phi(uf, ζ_p, tm)
                @test 0.08 <= R_f <= 0.23
                @test 0.07 <=
                      ζ_p / UF.dimensionless_profile(uf, FT(10), ζ_p, FT(10 / r),
                          tm, pv) <= 0.14
            end
        end

        rough = SF.ConstantRoughnessParams(FT(0.1), FT(0.01))
        gust = SF.ConstantGustinessSpec(FT(1))
        cfg_default = SF.SurfaceFluxConfig(rough, gust)
        cfg_none = SF.SurfaceFluxConfig(
            rough, gust, SF.MoistModel(), SF.NoRoughnessSubLayer(), SF.NoStabilityCap(),
        )
        cfg_const = SF.SurfaceFluxConfig(
            rough, gust, SF.MoistModel(), SF.NoRoughnessSubLayer(),
            SF.ConstantStabilityCap(FT(0.5)),
        )
        cfg_max = SF.SurfaceFluxConfig(
            rough, gust, SF.MoistModel(), SF.NoRoughnessSubLayer(),
            SF.MaxHeatFluxStabilityCap(),
        )
        Δz, d = FT(10), FT(0)
        ζ_p = SF.max_heat_flux_stability(param_set, Δz - d, FT(0.1))

        @testset "NoStabilityCap is the default ($FT)" begin
            @test cfg_default.stability_cap isa SF.NoStabilityCap
            for (T_sfc, U) in ((290, 3), (299, 5), (305, 2))
                a = solve(param_set, cfg_default, 300, T_sfc, U, Δz, d)
                b = solve(param_set, cfg_none, 300, T_sfc, U, Δz, d)
                @test a.shf == b.shf && a.ustar == b.ustar && a.ζ == b.ζ
            end
        end

        @testset "Unstable and weakly stable conditions unaffected ($FT)" begin
            for (T_sfc, U) in ((305, 2), (310, 5), (299.8, 8))
                a = solve(param_set, cfg_none, 300, T_sfc, U, Δz, d)
                b = solve(param_set, cfg_max, 300, T_sfc, U, Δz, d)
                @test a.ζ < ζ_p
                @test isapprox(b.shf, a.shf; rtol = 1e-3) &&
                      isapprox(b.ustar, a.ustar; rtol = 1e-3) &&
                      isapprox(b.ζ, a.ζ; rtol = 1e-3, atol = 1e-4)
            end
        end

        @testset "Capped coefficients and consistent ζ ($FT)" begin
            for (cfg, ζ_cap) in ((cfg_const, FT(0.5)), (cfg_max, ζ_p))
                out = solve(param_set, cfg, 300, 290, 2, Δz, d)
                @test out.converged
                @test out.ζ > ζ_cap
                @test out.Cd ≈ SF.drag_coefficient(param_set, ζ_cap, FT(0.1), Δz - d)
                Ch = SF.heat_exchange_coefficient(
                    param_set,
                    ζ_cap,
                    FT(0.1),
                    FT(0.01),
                    Δz - d,
                )
                @test out.g_h ≈ Ch * FT(2)
                # The returned ζ is the Obukhov stability parameter implied by the
                # fluxes: ζ = Ri_b F_m(ζ_cap)^2 / F_h(ζ_cap), i.e., Ri_b is linear in ζ
                Ri_b = SF.bulk_richardson_number(uf, SF.NoRoughnessSubLayer(),
                    Δz - d, out.ζ, FT(0.1), FT(0.01), pv, ζ_cap)
                Ri_b_cap = SF.bulk_richardson_number(uf, SF.NoRoughnessSubLayer(),
                    Δz - d, ζ_cap, FT(0.1), FT(0.01), pv)
                @test Ri_b ≈ Ri_b_cap * out.ζ / ζ_cap
                @test out.L_MO ≈ (Δz - d) / out.ζ rtol = 2e-2
            end
        end

        @testset "Profile recovery with a cap ($FT)" begin
            # Beyond the cap, profiles recovered with the returned effective length
            # L_eff = Δz_eff / min(ζ, ζ_cap) reproduce the forcing values; profiles
            # recovered with L_MO do not.
            gust0 = SF.ConstantGustinessSpec(FT(0))
            cfg =
                SF.SurfaceFluxConfig(rough, gust0, SF.DryModel(), SF.NoRoughnessSubLayer(),
                    SF.MaxHeatFluxStabilityCap())
            sfc(T_sfc) = SF.surface_fluxes(param_set, FT(300), FT(0), FT(0), FT(0),
                FT(1.2), FT(T_sfc), FT(0), FT(0), Δz, d, (FT(2), FT(0)), (FT(0), FT(0)),
                nothing, cfg, SF.PointValueScheme(),
                SF.SolverOptions{FT}(maxiter = 40, tol = FT(1e-6), rtol = FT(1e-6),
                    forced_fixed_iters = false),
                nothing, nothing, nothing)
            out = sfc(290)
            @test out.converged && out.ζ > ζ_p
            @test out.L_eff ≈ (Δz - d) / ζ_p
            U_rec(L) = SF.compute_profile_value(param_set, L, FT(0.1), Δz - d, out.ustar,
                FT(0), UF.MomentumTransport())
            @test U_rec(out.L_eff) ≈ 2 rtol = 1e-3
            @test U_rec(out.L_MO) > FT(1.5) * 2
            # Below the cap (and in unstable conditions), L_eff is L_MO
            for T_sfc in (299.8, 305)
                out = sfc(T_sfc)
                @test out.ζ < ζ_p
                @test out.L_eff == out.L_MO
            end
        end

        @testset "Supercritical conditions: convergence and monotone flux ($FT)" begin
            for U in (1, 2, 4)
                shf_prev = zero(FT)
                for ΔT in (0.5, 1, 2, 4, 8, 16, 32)
                    out = solve(param_set, cfg_max, 300, 300 - ΔT, U, Δz, d)
                    @test out.converged
                    @test isfinite(out.shf) && out.shf < 0
                    @test out.shf < shf_prev  # more downward flux for larger ΔT
                    shf_prev = out.shf
                end
            end
            # Without the cap, the flux collapses in supercritical conditions
            a = solve(param_set, cfg_none, 300, 270, 1, Δz, d)
            b = solve(param_set, cfg_max, 300, 270, 1, Δz, d)
            @test abs(b.shf) > 5 * abs(a.shf)
        end

        @testset "Roots beyond ζ_max with a cap ($FT)" begin
            # Strong stratification at low wind over a deep layer: beyond the cap, Ri_b is
            # linear in ζ and the root lies beyond ζ_max = 100. It must be bracketed and
            # found (not saturated at ζ_max), with default and tight solver options.
            thermo_params = SFP.thermodynamics_params(param_set)
            for cfg in (cfg_const, cfg_max), (Δz_d, ΔT) in ((50, 10), (50, 30), (20, 30))
                inputs = SF.build_surface_flux_inputs(
                    FT(300), FT(0.005), FT(0), FT(0), FT(1.2), FT(300 - ΔT), FT(0.005),
                    FT(0), FT(Δz_d), FT(0), (FT(1), FT(0)), (FT(0), FT(0)), cfg,
                    nothing,
                    SF.FluxSpecs{FT}(), nothing, nothing,
                )
                inputs = SF.with_stability_cap(inputs, param_set, pv)
                residual = SF.ResidualFunction(param_set, inputs, pv, uf, thermo_params)
                for opts in (SF.SolverOptions{FT}(),
                    SF.SolverOptions{FT}(maxiter = 40, tol = FT(1e-6), rtol = FT(1e-6),
                        forced_fixed_iters = false))
                    out = @inferred SF.surface_fluxes(param_set, inputs, pv, opts)
                    @test out.converged
                    @test out.ζ > 100
                    @test abs(residual(out.ζ)) <= 1e-3 * abs(residual(zero(FT)))
                end
            end
        end

        @testset "Prescribed-flux modes honor the cap ($FT)" begin
            # Strongly stable prescribed fluxes (ζ ≈ 17 > ζ_cap): the diagnostic heat
            # conductance and L_eff are evaluated at the cap, as in the MOST solve.
            ustar, shf = FT(0.05), FT(-20)
            U = FT(2)
            for (cfg, ζ_cap) in ((cfg_const, FT(0.5)), (cfg_max, ζ_p)),
                specs in (SF.FluxSpecs{FT}(; shf, lhf = FT(0), ustar),
                    SF.FluxSpecs{FT}(; shf, lhf = FT(0), Cd = FT(ustar^2 / (U^2 + 1))))

                out = SF.surface_fluxes(param_set, FT(300), FT(0.005), FT(0), FT(0),
                    FT(1.2), FT(295), FT(0.005), FT(0), Δz, d, (U, FT(0)),
                    (FT(0), FT(0)),
                    nothing, cfg, pv, nothing, specs)
                ref = SF.surface_fluxes(param_set, FT(300), FT(0.005), FT(0), FT(0),
                    FT(1.2), FT(295), FT(0.005), FT(0), Δz, d, (U, FT(0)),
                    (FT(0), FT(0)),
                    nothing, cfg_none, pv, nothing, specs)
                @test out.ζ > ζ_cap
                @test out.ζ == ref.ζ && out.shf == ref.shf
                Ch = SF.heat_exchange_coefficient(
                    param_set,
                    ζ_cap,
                    FT(0.1),
                    FT(0.01),
                    Δz - d,
                )
                Ch_ref =
                    SF.heat_exchange_coefficient(
                        param_set,
                        ref.ζ,
                        FT(0.1),
                        FT(0.01),
                        Δz - d,
                    )
                @test out.g_h / ref.g_h ≈ Ch / Ch_ref
                @test out.g_h > ref.g_h
                @test out.L_eff ≈ (Δz - d) / ζ_cap
                @test ref.L_eff == ref.L_MO
            end
        end

        @testset "Convergence flag with default solver options ($FT)" begin
            # Beyond the cap, the residual is linear in ζ: regula falsi hits the root
            # and the far bracket endpoint never moves, so the bracket width alone
            # would report non-convergence of an accurate solve.
            opts = SF.SolverOptions{FT}()  # 7 fixed iterations, tol = rtol = 1e-2
            for cfg in (cfg_const, cfg_max), U in (1, 2, 4, 8), ΔT in (0.5, 2, 8, 32)
                ref = solve(param_set, cfg, 300, 300 - ΔT, U, Δz, d)
                out = SF.surface_fluxes(
                    param_set, FT(300), FT(0.005), FT(0), FT(0), FT(1.2), FT(300 - ΔT),
                    FT(0.005), FT(0), Δz, d, (FT(U), FT(0)), (FT(0), FT(0)), nothing,
                    cfg,
                    SF.PointValueScheme(), opts, nothing, nothing, nothing,
                )
                @test ref.converged
                @test out.converged
                @test abs(out.ζ - ref.ζ) <= max(opts.tol, opts.rtol * abs(ref.ζ))
            end
            # Linear residual with the root at an endpoint: step 0, bracket width 8
            ζ, conv = SF.monin_obukhov_finalize(
                FT(2), FT(10), FT(0), FT(8), FT(2), true, false, FT(100), opts,
            )
            @test ζ == 2 && conv
            # No sign change: never converged
            _, conv = SF.monin_obukhov_finalize(
                FT(2), FT(10), FT(1), FT(8), FT(2), false, false, FT(100), opts,
            )
            @test !conv
        end

        @testset "Callbacks see the capped conductance ($FT)" begin
            # Surface temperature relaxes toward a fixed skin temperature through a
            # conductance in series with the aerodynamic conductance (big-leaf canopy)
            T_leaf = FT(290)
            g_leaf = FT(0.05)
            function update_T_sfc(ζ, param_set, thermo_params, inputs, scheme, u_star,
                z0m, z0h)
                g_h = SF.heat_conductance(param_set, ζ, u_star, inputs, z0m, z0h, scheme)
                return (inputs.T_int + T_leaf * g_leaf / g_h) / (1 + g_leaf / g_h)
            end
            out = solve(param_set, cfg_max, 300, 290, 1, Δz, d; update_T_sfc)
            @test out.converged
            @test out.ζ > ζ_p
            Ch = SF.heat_exchange_coefficient(param_set, ζ_p, FT(0.1), FT(0.01), Δz - d)
            @test out.g_h ≈ Ch * FT(1)
            @test out.T_sfc ≈ (FT(300) + T_leaf * g_leaf / out.g_h) / (1 + g_leaf / out.g_h) rtol =
                1e-4
        end
    end

    @testset "Validation, constructors, and helpers with builder inputs" begin
        FT = Float64
        param_set = SFP.SurfaceFluxesParameters(FT, UF.GryanikParams)
        pv = SF.PointValueScheme()
        @test_throws ArgumentError SF.ConstantStabilityCap(0.0)
        @test_throws ArgumentError SF.ConstantStabilityCap(-1.0)
        @test SF.ConstantStabilityCap(0.5f0) isa SF.ConstantStabilityCap{Float32}
        # Positional constructor without L_eff (L_eff = L_MO)
        sfc =
            SF.SurfaceFluxConditions(1.0, 2.0, 3.0, 4.0, 5.0, 0.3, 0.1, 1e-3, 1e-2, 290.0,
                0.01, 50.0, true)
        @test sfc.L_eff == sfc.L_MO == 50.0 && sfc.converged
        @test sfc isa SF.SurfaceFluxConditions{Float64}
        # ζ_p stays within the search interval for extreme geometries
        for (Δz, z0m) in ((1.0, 0.999), (1e6, 1e-6), (1e-3, 1e-9))
            ζ_p = SF.max_heat_flux_stability(param_set, Δz, z0m)
            @test 1e-2 <= ζ_p <= 20
        end
        # Far probe: for a fixed surface state and a cap ≤ 10, the extrapolated root is
        # exact (Ri_b beyond the cap is a line through the origin, r0 = -Ri_b,state)
        k, Ri_state = 0.02, 5.0
        r0, r2 = -Ri_state, 10k - Ri_state
        @test SF.monin_obukhov_far_probe(0.5, 100.0, 10.0, r0, r2) ≈ 2 * Ri_state / k
        @test SF.monin_obukhov_far_probe(0.5, 100.0, 10.0, r0, 10k - 1.0) == 100.0  # root < ζ_max
        @test SF.monin_obukhov_far_probe(nothing, 100.0, 10.0, r0, r2) == 100.0
        # Helpers called with builder inputs (ζ_cap = nothing) agree with the solve
        rough = SF.ConstantRoughnessParams(0.1, 0.01)
        gust = SF.ConstantGustinessSpec(1.0)
        for cap in (SF.ConstantStabilityCap(0.5), SF.MaxHeatFluxStabilityCap(),
            SF.NoStabilityCap())
            cfg =
                SF.SurfaceFluxConfig(rough, gust, SF.MoistModel(), SF.NoRoughnessSubLayer(),
                    cap)
            args = (300.0, 0.005, 0.0, 0.0, 1.2, 290.0, 0.005, 0.0, 10.0, 0.0, (2.0, 0.0),
                (0.0, 0.0))
            out = SF.surface_fluxes(param_set, args..., nothing, cfg, pv)
            inputs =
                SF.build_surface_flux_inputs(args..., cfg, nothing, SF.FluxSpecs{FT}(),
                    nothing, nothing)
            @test inputs.ζ_cap === nothing
            @test SF.heat_conductance(param_set, out.ζ, out.ustar, inputs, 0.1, 0.01, pv) ≈
                  out.g_h
            @test SF.compute_ustar(param_set, out.ζ, 0.1, inputs, pv, 1.0) ≈ out.ustar rtol =
                1e-6
            ζ_c = SF.capped_stability(param_set, inputs, pv, out.ζ)
            @test ζ_c == SF.capped_stability(
                out.ζ,
                SF.stability_cap_value(cap, param_set,
                    inputs, pv),
            )
            @test SF.resolved_stability_cap(param_set, inputs, pv) ===
                  SF.resolved_stability_cap(
                param_set,
                SF.with_stability_cap(inputs, param_set, pv),
                pv,
            )
        end
        # neutral_momentum_roughness: the model's roughness length for roughness independent
        # of ustar (ConstantRoughnessParams, RaupachRoughnessParams) and the neutral solve
        # otherwise (COARE3RoughnessParams)
        @test !SF.depends_on_ustar(rough) && !SF.depends_on_ustar(gust)
        @test !SF.depends_on_ustar(SF.RaupachRoughnessParams{FT}())
        @test SF.depends_on_ustar(SF.COARE3RoughnessParams{FT}())
        @test SF.depends_on_ustar(SF.DeardorffGustinessSpec())
        args = (300.0, 0.005, 0.0, 0.0, 1.2, 290.0, 0.005, 0.0, 10.0, 0.0, (5.0, 0.0),
            (0.0, 0.0))
        for (rm, ri) in (
            (rough, nothing),
            (SF.RaupachRoughnessParams{FT}(), (; LAI = 3.0, h = 10.0)),
            (SF.COARE3RoughnessParams{FT}(), nothing),
        )
            cfg_rm = SF.SurfaceFluxConfig(
                rm, gust, SF.MoistModel(), SF.NoRoughnessSubLayer(),
                SF.MaxHeatFluxStabilityCap(),
            )
            inp = SF.build_surface_flux_inputs(
                args..., cfg_rm, ri, SF.FluxSpecs{FT}(), nothing, nothing,
            )
            z0m_n = @inferred SF.neutral_momentum_roughness(rm, param_set, inp, pv)
            inp_uncapped =
                (; inp..., stability_cap = SF.NoStabilityCap(), ζ_cap = nothing)
            _, z0m_ref, _ =
                SF.compute_ustar_and_roughness(param_set, 0.0, inp_uncapped, pv)
            @test z0m_n ≈ z0m_ref rtol = 1e-6
            @test SF.stability_cap_value(SF.MaxHeatFluxStabilityCap(), param_set, inp, pv) ≈
                  SF.max_heat_flux_stability(param_set, 10.0, z0m_ref, pv)
        end
        # ζ_max of a ConstantStabilityCap is differentiable (e.g., for calibration)
        opts = SF.SolverOptions{FT}(maxiter = 40, tol = 1e-10, rtol = 1e-10,
            forced_fixed_iters = false)
        function shf_of_cap(c)
            cfg =
                SF.SurfaceFluxConfig(rough, gust, SF.DryModel(), SF.NoRoughnessSubLayer(),
                    SF.ConstantStabilityCap(c))
            return SF.surface_fluxes(param_set, 300.0, 0.0, 0.0, 0.0, 1.2, 290.0, 0.0, 0.0,
                10.0, 0.0, (2.0, 0.0), (0.0, 0.0), nothing, cfg, pv, opts).shf
        end
        d_ad = ForwardDiff.derivative(shf_of_cap, 0.5)
        h = 1e-5
        d_fd = (shf_of_cap(0.5 + h) - shf_of_cap(0.5 - h)) / 2h
        @test isfinite(d_ad) && d_ad != 0
        @test d_ad ≈ d_fd rtol = 1e-5
    end

    @testset "Automatic differentiation with respect to Δz and z0m" begin
        # Derivatives with respect to Δz or z0m make the cap value ζ_cap a dual number,
        # while the solver's ζ probes are plain floats.
        FT = Float64
        param_set = SFP.SurfaceFluxesParameters(FT, UF.GryanikParams)
        opts = SF.SolverOptions{FT}(maxiter = 40, tol = 1e-10, rtol = 1e-10,
            forced_fixed_iters = false)
        function shf(Δz, z0m, T_sfc, cap)
            cfg = SF.SurfaceFluxConfig(SF.ConstantRoughnessParams(z0m, oftype(z0m, 0.01)),
                SF.ConstantGustinessSpec(1.0), SF.DryModel(), SF.NoRoughnessSubLayer(),
                cap)
            return SF.surface_fluxes(param_set, 300.0, 0.0, 0.0, 0.0, 1.2, T_sfc, 0.0, 0.0,
                Δz, 0.0, (2.0, 0.0), (0.0, 0.0), nothing, cfg, SF.PointValueScheme(), opts,
                nothing, nothing, nothing)
        end
        # One closure definition site, so that all gradients share the dual tag type
        shf_of(T_sfc, cap) = x -> shf(x[1], x[2], T_sfc, cap).shf
        x0 = [10.0, 0.1]
        # Weakly stable state below the cap
        T_sfc = 299.8
        @test shf(x0..., T_sfc, SF.MaxHeatFluxStabilityCap()).ζ <
              SF.max_heat_flux_stability(param_set, x0...)
        f = shf_of(T_sfc, SF.MaxHeatFluxStabilityCap())
        ad = ForwardDiff.gradient(f, x0)
        fd = [
            (f(x0 .+ h .* e) - f(x0 .- h .* e)) / 2h
            for (h, e) in ((1e-4, [1.0, 0.0]), (1e-6, [0.0, 1.0]))
        ]
        @test all(isfinite, ad)
        @test ad ≈ fd rtol = 1e-5

        # Capped state: the derivatives include the dependence of ζ_p on Δz and z0m.
        # Independent reference: ζ_p from the analytical point-value condition
        # F_m(ζ) = 3 [φ_m(ζ) - φ_m(ζ z0m/Δz)] (bisection), with the fluxes computed
        # at a constant cap equal to this ζ_p, differentiated by central differences.
        uf = SFP.uf_params(param_set)
        tm = UF.MomentumTransport()
        function ζ_p_exact(Δz, z0m)
            G(ζ) =
                UF.dimensionless_profile(uf, Δz, ζ, z0m, tm) -
                3 * (UF.phi(uf, ζ, tm) - UF.phi(uf, ζ * z0m / Δz, tm))
            lo, hi = 1e-2, 20.0
            for _ in 1:100
                mid = (lo + hi) / 2
                G(mid) > 0 ? (lo = mid) : (hi = mid)
            end
            return (lo + hi) / 2
        end
        for (i, h) in ((1, 1e-3), (2, 1e-5))
            e = i == 1 ? [1.0, 0.0] : [0.0, 1.0]
            dζ_ad =
                ForwardDiff.gradient(x -> SF.max_heat_flux_stability(param_set, x...), x0)[i]
            dζ_ref = (ζ_p_exact((x0 .+ h .* e)...) - ζ_p_exact((x0 .- h .* e)...)) / 2h
            @test dζ_ad ≈ dζ_ref rtol = 1e-2
        end
        T_sfc = 290.0
        @test shf(x0..., T_sfc, SF.MaxHeatFluxStabilityCap()).ζ >
              SF.max_heat_flux_stability(param_set, x0...)
        ad = ForwardDiff.gradient(shf_of(T_sfc, SF.MaxHeatFluxStabilityCap()), x0)
        f_ref(x) = shf(x..., T_sfc, SF.ConstantStabilityCap(ζ_p_exact(x...))).shf
        ref = [
            (f_ref(x0 .+ h .* e) - f_ref(x0 .- h .* e)) / 2h
            for (h, e) in ((1e-3, [1.0, 0.0]), (1e-5, [0.0, 1.0]))
        ]
        @test ad ≈ ref rtol = 1e-2
    end
end

end # module
