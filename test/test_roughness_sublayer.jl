# Roughness sublayer models: the RSL correction matches reference quadrature for point
# values and layer averages at all stabilities, is positive and vanishes above the RSL,
# keeps the exchange coefficients finite and positive over tall canopies, and enters
# u*, θ*, q*, the fluxes, and profile recovery consistently.

module TestRoughnessSubLayer

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import ClimaParams as CP

const MT = UF.MomentumTransport()
const HT = UF.HeatTransport()
const PV = UF.PointValueScheme()
const LA = UF.LayerAverageScheme()

μ_ref(::SF.LinearRSL, c, z, zR) = z < zR ? 1 - c * (1 - z / zR) : 1.0
μ_ref(::SF.ExponentialRSL, c, z, zR) = z < zR ? exp(-c * (1 - z / zR)) : 1.0
coef(m, ::UF.MomentumTransport) = m.c_m
coef(m, ::UF.HeatTransport) = m.c_h

# Midpoint rule in ln z for ∫_{a}^{b} g(z) dz / z
function ∫logz(g, a, b; n = 4000)
    a >= b && return 0.0
    la, lb = log(a), log(b)
    h = (lb - la) / n
    return sum(g(exp(la + (i - 0.5) * h)) for i in 1:n) * h
end

# Reference point-value correction (Float64)
function P_point_ref(uf, m, Δz, ζ, z0, tr)
    c = coef(m, tr)
    zR = m.z_RSL
    g(z) = UF.phi(uf, ζ * z / Δz, tr) * (1 - μ_ref(m, c, z, zR))
    return ∫logz(g, max(min(Δz, zR), z0), zR)
end

# Reference layer-averaged correction: (1/Δz) ∫_{z0}^{Δz} P_point(z; ζ z / Δz) dz
function P_layer_ref(uf, m, Δz, ζ, z0, tr; n = 400)
    h = (Δz - z0) / n
    return sum(
        P_point_ref(uf, m, z, ζ * z / Δz, z0, tr) for z in (z0 .+ ((1:n) .- 0.5) .* h)
    ) *
           h / Δz
end

# Same, with the order of integration swapped:
# (1/Δz) ∫_{z0}^{z_RSL} φ(z'/L) (1 - μ(z')) (min(z', Δz) - z0) dz'/z'
function P_layer_ref_swapped(uf, m, Δz, ζ, z0, tr; n = 40_000)
    c = coef(m, tr)
    zR = m.z_RSL
    g(z) = UF.phi(uf, ζ * z / Δz, tr) * (1 - μ_ref(m, c, z, zR)) * (min(z, Δz) - z0) / Δz
    zc = min(Δz, zR)
    return ∫logz(g, z0, zc; n) + ∫logz(g, zc, zR; n)
end

models(zR) = (
    SF.LinearRSL(c_m = 0.5, c_h = 0.4, z_RSL = zR),
    SF.ExponentialRSL(c_m = 0.7, c_h = 0.9, z_RSL = zR),
)

@testset "Roughness sublayer models" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    uf = SFP.uf_params(param_set)
    Pr_0 = SFP.Pr_0(param_set)

    @testset "NoRoughnessSubLayer — zero correction" begin
        rsl = SF.NoRoughnessSubLayer()
        for (Δz, ζ, z0, tr, sch) in ((10.0, 0.0, 0.1, MT, PV), (30.0, -2.0, 1.0, HT, LA))
            @test SF.rsl_profile_correction(uf, rsl, Δz, ζ, z0, tr, sch) == 0
            @test SF.rsl_corrected_profile(uf, rsl, Δz, ζ, z0, tr, sch) ==
                  UF.dimensionless_profile(uf, Δz, ζ, z0, tr, sch)
        end
    end

    @testset "Neutral closed form of the linear RSL (incl. Pr_0 for scalars)" begin
        z0, zR = 0.5, 30.0
        for Δz in (5.0, 20.0, 60.0), tr in (MT, HT)
            c = tr isa UF.MomentumTransport ? 0.5 : 0.4
            ϕ0 = tr isa UF.MomentumTransport ? 1.0 : Pr_0
            zc = min(Δz, zR)
            m_t = SF.LinearRSL(c_m = 0.5, c_h = 0.4, z_RSL = zR)
            P_t = ϕ0 * c * (log(zR / zc) - (zR - zc) / zR)
            @test SF.rsl_profile_correction(uf, m_t, Δz, 0.0, z0, tr) ≈ P_t rtol = 1e-4 atol =
                1e-12
        end
    end

    @testset "Stability dependence and layer averages vs reference quadrature" begin
        for (Δz, z0, zR) in ((12.0, 1.0, 30.0), (40.0, 2.0, 30.0), (8.0, 0.05, 20.0)),
            m in models(zR),
            tr in (MT, HT),
            ζ in (-10.0, -1.0, 0.0, 0.3, 2.0)

            P = SF.rsl_profile_correction(uf, m, Δz, ζ, z0, tr, PV)
            @test P ≈ P_point_ref(uf, m, Δz, ζ, z0, tr) rtol = 2e-3 atol = 1e-6
        end
        for (Δz, z0, zR) in ((12.0, 1.0, 30.0), (40.0, 2.0, 20.0)),
            m in models(zR),
            ζ in (-2.0, 0.0, 1.0)

            P = SF.rsl_profile_correction(uf, m, Δz, ζ, z0, MT, LA)
            @test P ≈ P_layer_ref(uf, m, Δz, ζ, z0, MT) rtol = 5e-3 atol = 1e-5
        end
        # Small roughness lengths (z_c / z0 up to 3·10⁴), where the weight (z' - z0)/z'
        # of the lower integral varies over many e-folds in ln z'
        for (Δz, z0, zR) in ((40.0, 1e-3, 30.0), (12.0, 1e-3, 30.0), (25.0, 0.01, 30.0)),
            m in models(zR),
            tr in (MT, HT),
            ζ in (-5.0, 0.0, 0.5)

            P = SF.rsl_profile_correction(uf, m, Δz, ζ, z0, tr, LA)
            @test P ≈ P_layer_ref_swapped(uf, m, Δz, ζ, z0, tr) rtol = 1e-2
        end
    end

    @testset "Sign and behavior above the RSL" begin
        z0, zR = 1.0, 20.0
        for m in models(zR), tr in (MT, HT)
            P(Δz, ζ = 0.0) = SF.rsl_profile_correction(uf, m, Δz, ζ, z0, tr)
            @test P(5.0) > 0
            @test P(5.0) > P(15.0) > 0          # decreases toward the RSL top
            @test P(20.0) == 0 && P(50.0) == 0   # MOST above the RSL
            @test SF.rsl_profile_correction(uf, m, 50.0, 0.0, z0, tr, LA) > 0  # layer includes the RSL
        end
        # The RSL reduces Cd within the RSL (for given apparent z0 and d)
        Cd(m) = SF.drag_coefficient(param_set, 0.0, z0, 10.0, PV, m)
        Cd0 = Cd(SF.NoRoughnessSubLayer())
        @test Cd(SF.ExponentialRSL(z_RSL = zR)) < Cd0
        @test Cd(SF.LinearRSL(z_RSL = zR)) < Cd0
        # Above the RSL, the corrected profile is exactly MOST
        @test SF.drag_coefficient(
            param_set,
            0.3,
            z0,
            40.0,
            PV,
            SF.ExponentialRSL(z_RSL = zR),
        ) ==
              SF.drag_coefficient(param_set, 0.3, z0, 40.0, PV)
    end

    @testset "Continuity and monotonicity in height" begin
        for m in models(20.0), tr in (MT, HT), sch in (PV, LA)
            zs = range(1.2, 60.0; length = 3000)
            F̂ = [SF.rsl_corrected_profile(uf, m, z, 0.0, 1.0, tr, sch) for z in zs]
            @test all(>(0), diff(F̂))
            @test maximum(abs, diff(F̂)) < 0.05
        end
    end

    @testset "Exact bounds and positivity for tall canopies" begin
        # z0 = 0.1 h, d = 0.67 h, z0h = z0/10, forcing 0.5–30 m above the canopy top
        for h in (5.0, 10.0, 20.0, 30.0, 50.0), dzabove in (0.5, 1.0, 2.0, 5.0, 30.0)
            z0m, d = 0.1h, 0.67h
            z0h = z0m / 10
            Δz_eff = h + dzabove - d
            for m in (
                    SF.LinearRSL(c_m = 0.6, c_h = 0.6, z_RSL = 2h - d),
                    SF.ExponentialRSL(c_m = 0.7, c_h = 0.7, z_RSL = h),
                    SF.ExponentialRSL(c_m = 0.7, c_h = 0.7, z_RSL = 2h - d),
                    SF.LinearRSL(c_m = 0.9, c_h = 0.9, z_RSL = 3h),
                ),
                sch in (PV, LA),
                ζ in (-100.0, -10.0, -1.0, 0.0, 1.0, 10.0)

                for (z0, tr) in ((z0m, MT), (z0h, HT))
                    F = UF.dimensionless_profile(uf, Δz_eff, ζ, z0, tr, sch)
                    F̂ = SF.rsl_corrected_profile(uf, m, Δz_eff, ζ, z0, tr, sch)
                    @test F̂ > 0
                    @test F̂ >= F
                end
                Cd = SF.drag_coefficient(param_set, ζ, z0m, Δz_eff, sch, m)
                Ch = SF.heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, sch, m)
                @test isfinite(Cd) && Cd > 0
                @test isfinite(Ch) && Ch > 0
            end
        end
    end

    @testset "Full solves over tall canopies are finite and bounded" begin
        gust = SF.ConstantGustinessSpec(FT(0.5))
        for h in (10.0, 30.0, 50.0), dzabove in (1.0, 2.0, 10.0),
            (ΔT, U) in ((2.0, 5.0), (2.0, 1.0), (0.0, 1.0), (8.0, 1.0)),
            sch in (PV, LA)

            z0m, d = 0.1h, 0.67h
            rough = SF.ConstantRoughnessParams{FT}(z0m = z0m, z0s = z0m / 10)
            solve(rsl) = SF.surface_fluxes(
                param_set, 290.0, 0.0, 0.0, 0.0, 1.2, 290.0 + ΔT, 0.0, 0.0, h + dzabove,
                d,
                (U, 0.0), (0.0, 0.0), nothing,
                SF.SurfaceFluxConfig(rough, gust, SF.DryModel(), rsl), sch,
            )
            ref = solve(SF.NoRoughnessSubLayer())
            for rsl in (
                SF.ExponentialRSL(c_m = 0.7, c_h = 0.7, z_RSL = h),
                SF.ExponentialRSL(c_m = 0.7, c_h = 0.7, z_RSL = 2h - d),
                SF.LinearRSL(c_m = 0.4, c_h = 0.4, z_RSL = 2h - d),
            )
                r = solve(rsl)
                @test isfinite(r.shf) && isfinite(r.ustar) && isfinite(r.Cd) &&
                      isfinite(r.g_h)
                @test r.ustar > 0 && r.Cd > 0 && r.g_h > 0
                if ΔT > 0  # unstable: no supercritical collapse
                    @test r.converged && ref.converged
                    # Within the RSL, the coefficients are reduced relative to MOST with
                    # the same apparent z0 and d, but remain bounded away from zero
                    @test 0.1 < r.Cd / ref.Cd <= 1 + 1e-6
                    @test 0.1 < r.g_h / ref.g_h <= 1 + 1e-6
                    @test r.ustar > 1e-4  # no collapse to the lower u* bracket limit
                end
            end
        end
    end

    @testset "Consistency of scales, fluxes, and profile recovery" begin
        rough = SF.ConstantRoughnessParams{FT}(z0m = FT(3), z0s = FT(0.3))
        gust = SF.ConstantGustinessSpec(FT(0))
        thermo_params = SFP.thermodynamics_params(param_set)
        for rsl in models(40.0), sch in (PV, LA), T_sfc in (288.0, 294.0)
            cfg = SF.SurfaceFluxConfig(rough, gust, SF.MoistModel(), rsl)
            Δz, d, U = 40.0, 20.0, 5.0
            args =
                (290.0, 0.008, 0.0, 0.0, 1.2, T_sfc, 0.01, 0.0, Δz, d, (U, 0.0), (0.0, 0.0))
            r = SF.surface_fluxes(param_set, args..., nothing, cfg, sch)
            inputs =
                SF.build_surface_flux_inputs(args..., cfg, nothing, SF.FluxSpecs{FT}(),
                    nothing, nothing)
            @test SF.compute_ustar(param_set, r.ζ, 3.0, inputs, sch, 0.0) ≈ r.ustar rtol =
                1e-3
            # θ* and q* from the similarity coefficients vs. from the fluxes
            ρ_sfc =
                SF.surface_density(param_set, 290.0, 1.2, T_sfc, Δz, 0.008, 0.0, 0.0, 0.01)
            cp_m = TD.cp_m(thermo_params, 0.008, 0.0, 0.0)
            θ_star = SF.compute_theta_star(param_set, r.ζ, 0.3, inputs, sch)
            q_star = SF.compute_q_star(param_set, r.ζ, 0.3, inputs, sch)
            @test q_star ≈ -r.evaporation / (ρ_sfc * r.ustar) rtol = 1e-3
            @test θ_star * r.ustar * ρ_sfc * cp_m ≈ -(r.shf) rtol = 0.05
            # Profile recovery returns the forcing wind speed at the forcing height
            U_rec = SF.compute_profile_value(
                param_set, r.L_MO, 3.0, Δz - d, r.ustar, 0.0, MT, sch, rsl,
            )
            @test U_rec ≈ U rtol = 1e-3
        end
    end

    @testset "Degenerate parameters and validation" begin
        for m in (SF.LinearRSL(z_RSL = 0.0), SF.ExponentialRSL(z_RSL = 0.0),
            SF.LinearRSL(z_RSL = 0.5), SF.ExponentialRSL(c_m = 0.0, c_h = 0.0))
            for sch in (PV, LA), Δz in (0.5, 10.0)
                @test SF.rsl_profile_correction(uf, m, Δz, -1.0, 1.0, MT, sch) == 0
            end
        end
        @test_throws ArgumentError SF.LinearRSL(c_m = 1.0)
        @test_throws ArgumentError SF.LinearRSL(c_h = -0.1)
        @test_throws ArgumentError SF.ExponentialRSL(c_m = -0.5)
        @test_throws ArgumentError SF.ExponentialRSL(z_RSL = -1.0)
    end

    @testset "Constructors and floating-point types" begin
        @test SF.ExponentialRSL() isa SF.ExponentialRSL{Float64}
        @test SF.LinearRSL(Float32) isa SF.LinearRSL{Float32}
        @test SF.LinearRSL(0.3f0, 0.3f0, 10.0f0) isa SF.LinearRSL{Float32}
        cfg = SF.SurfaceFluxConfig(SF.ConstantRoughnessParams{FT}(),
            SF.ConstantGustinessSpec(FT(1)), SF.DryModel(), SF.ExponentialRSL())
        @test cfg.rsl_model isa SF.ExponentialRSL

        ps32 = SFP.SurfaceFluxesParameters(Float32, UF.BusingerParams)
        uf32 = SFP.uf_params(ps32)
        for m in (SF.ExponentialRSL(Float32), SF.LinearRSL(Float32),
                SF.LinearRSL(Float32; c_m = 0.8)),
            sch in (PV, LA)

            # Float32 models and inputs give Float32 results
            @test SF.rsl_corrected_profile(uf32, m, 8.0f0, -0.5f0, 1.0f0, HT, sch) isa
                  Float32
            @test SF.drag_coefficient(ps32, 0.2f0, 1.0f0, 8.0f0, sch, m) isa Float32
            rough = SF.ConstantRoughnessParams{Float32}(z0m = 3.0f0, z0s = 0.3f0)
            cfg32 = SF.SurfaceFluxConfig(rough, SF.ConstantGustinessSpec(0.5f0),
                SF.DryModel(), m)
            r = SF.surface_fluxes(ps32, 290.0f0, 0.0f0, 0.0f0, 0.0f0, 1.2f0, 292.0f0, 0.0f0,
                0.0f0, 40.0f0, 20.0f0, (5.0f0, 0.0f0), (0.0f0, 0.0f0), nothing, cfg32,
                sch)
            @test r.shf isa Float32 && r.Cd isa Float32 && r.ustar isa Float32
            @test r.converged
        end
    end
end

end # module
