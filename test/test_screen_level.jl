# Screen-level reconstruction from a solve: neutral profiles are logarithmic, the
# reference level returns the interior state and the roughness length the surface
# state for every stability, including under a stability cap, where the capped
# `L_eff` and `ζ_eff` keep the profile consistent with the fluxes.

module TestScreenLevel

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaParams as CP
import ForwardDiff
using Main: @test_allocs_and_ts

const Δz = 20
const d = 2
const z0m = 0.1
const z0h = 0.01

# Allocation and type-stability checks need the arguments passed in, since locals of
# a testset body captured by the check's closure are boxed
screen_level_checked(ps, sc, inputs) =
    @test_allocs_and_ts SF.screen_level_values(ps, sc, inputs, 2.0f0, 10.0f0)

function solve(FT, param_set; T_sfc, T_int = 290, u = 3, q_int = 0.005, q_sfc = 0.008,
    cap = SF.NoStabilityCap(), gust = SF.ConstantGustinessSpec(FT(1)), opts = nothing,
    scheme = SF.PointValueScheme())
    config = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(z0m), FT(z0h)), gust, SF.MoistModel(),
        SF.NoRoughnessSubLayer(), cap,
    )
    inputs = SF.build_surface_flux_inputs(
        FT(T_int), FT(q_int), FT(0), FT(0), FT(1.2), FT(T_sfc), FT(q_sfc), FT(0),
        FT(Δz), FT(d), (FT(u), FT(0)), (FT(0), FT(0)), config, nothing,
        SF.FluxSpecs{FT}(), nothing, nothing,
    )
    sc = SF.surface_fluxes(param_set, inputs, scheme, opts)
    return sc, inputs
end

@testset "Screen-level values" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)
    κ = SFP.von_karman_const(param_set)
    g = SFP.grav(param_set)
    cp_d = TD.Parameters.cp_d(thermo_params)
    Δz_eff = FT(Δz - d)
    opts = SF.SolverOptions{FT}(
        maxiter = 50,
        tol = FT(1e-8),
        rtol = FT(1e-8),
        forced_fixed_iters = false,
    )

    @testset "Near-neutral: logarithmic profiles" begin
        # Equal dry static energies and humidities leave a small residual stability
        # from the virtual temperature, and the wind profile is logarithmic to that
        # order
        T_sfc = FT(290) + g * Δz_eff / cp_d
        sc, inputs = solve(FT, param_set; T_sfc, T_int = 290, q_int = 0.005, q_sfc = 0.005)
        @test abs(sc.ζ) < 1e-3
        s = SF.screen_level_values(param_set, sc, inputs, FT(2), FT(10))
        @test s.u ≈ sc.ustar / κ * log((z0m + 10) / z0m) rtol = 1e-3
        @test s.q ≈ FT(0.005)
        # The screen temperature follows the dry adiabat from the surface state at the
        # displacement height
        @test s.T ≈ T_sfc - g / cp_d * (z0h + 2) rtol = 1e-6
    end

    @testset "Reference level returns the interior state" begin
        for T_sfc in (FT(280), FT(290), FT(300)),
            cap in (SF.NoStabilityCap(), SF.MaxHeatFluxStabilityCap())

            sc, inputs = solve(FT, param_set; T_sfc, cap, opts)
            @test sc.ζ_eff == SF.capped_stability(
                sc.ζ,
                SF.stability_cap_value(cap, param_set, inputs, SF.PointValueScheme()),
            )
            @test sc.ζ_eff ≈ Δz_eff / sc.L_eff
            s = SF.screen_level_values(param_set, sc, inputs, Δz_eff - z0h, Δz_eff - z0m)
            @test s.T ≈ inputs.T_int rtol = 1e-8
            @test s.q ≈ inputs.q_tot_int rtol = 1e-8
            # At the reference level, the wind is the effective wind speed of the solve
            @test s.u ≈ sc.ustar / sqrt(sc.Cd) rtol = 1e-8
            # Beyond the reference level, the interior values are held
            above = SF.screen_level_values(param_set, sc, inputs, FT(100), FT(100))
            @test above.T ≈ inputs.T_int rtol = 1e-8
            @test above.q ≈ inputs.q_tot_int rtol = 1e-8
            @test above.u ≈ s.u rtol = 1e-8
            # At the roughness lengths, the apparent sinks at d + z0 above the surface,
            # the dry static energy and humidity are those of the surface state at d
            dse(T, z) = cp_d * T + g * z
            at_z0 = SF.screen_level_values(param_set, sc, inputs, FT(0), FT(0))
            @test dse(at_z0.T, z0h + d) ≈ dse(sc.T_sfc, d)
            @test at_z0.q ≈ sc.q_vap_sfc
            @test at_z0.u == 0
            # Between them, the screen values lie between the surface and interior
            # values of the dry static energy and humidity
            mid = SF.screen_level_values(param_set, sc, inputs, FT(2), FT(10))
            lo, hi = minmax(dse(sc.T_sfc, d), dse(inputs.T_int, Δz))
            @test lo - 1e-9 <= dse(mid.T, z0h + 2 + d) <= hi + 1e-9
            lo_q, hi_q = minmax(sc.q_vap_sfc, inputs.q_tot_int)
            @test lo_q <= mid.q <= hi_q
            @test 0 < mid.u < s.u
        end
    end

    @testset "Capped stable case: the uncapped length overestimates the differences" begin
        sc, inputs =
            solve(FT, param_set; T_sfc = 275, cap = SF.MaxHeatFluxStabilityCap(), opts)
        @test sc.ζ_eff < sc.ζ
        @test sc.L_eff > sc.L_MO
        uncapped = SF.SurfaceFluxConditions(
            sc.shf, sc.lhf, sc.evaporation, sc.ρτxz, sc.ρτyz, sc.ustar, sc.ζ, sc.Cd,
            sc.g_h,
            sc.T_sfc, sc.q_vap_sfc, sc.L_MO, sc.L_MO, sc.ζ, sc.converged,
        )
        s = SF.screen_level_values(param_set, sc, inputs, FT(2), FT(10))
        s_uncapped = SF.screen_level_values(param_set, uncapped, inputs, FT(2), FT(10))
        # With the uncapped length, the profile misses the reference-level state
        ref = SF.screen_level_values(
            param_set, uncapped, inputs, Δz_eff - z0h, Δz_eff - z0m,
        )
        @test !(ref.u ≈ sc.ustar / sqrt(sc.Cd))
        @test s.u != s_uncapped.u
    end

    @testset "Layer-average solve: point values at the screen and anemometer" begin
        # The interior state is a layer average, but the screen and anemometer values
        # are point values of the profile with the fluxes of the solve
        scheme = SF.LayerAverageScheme()
        sc, inputs = solve(FT, param_set; T_sfc = 295, scheme, opts)
        point_u(z) =
            sc.ustar / κ * SF.dimensionless_profile_value(
                param_set, sc.L_eff, z0m, z0m + z, Δz_eff, UF.MomentumTransport(),
                SF.PointValueScheme(), SF.NoRoughnessSubLayer(),
            )
        s = SF.screen_level_values(param_set, sc, inputs, FT(2), FT(10), scheme)
        @test s.u ≈ point_u(10) rtol = 1e-8
        # The layer-averaged wind of the solve is the point value low in the layer, near
        # Δz_eff / e for a logarithmic profile
        u_mean = sc.ustar / sqrt(sc.Cd)
        @test point_u(2) < u_mean < point_u(Δz_eff - z0m)
        # Beyond the reference level, the point values at the reference level are held
        above = SF.screen_level_values(param_set, sc, inputs, FT(100), FT(100), scheme)
        @test above.u ≈ point_u(Δz_eff - z0m) rtol = 1e-8
        lo_q, hi_q = minmax(sc.q_vap_sfc, inputs.q_tot_int)
        @test lo_q <= s.q <= hi_q
    end

    @testset "Legacy constructors derive ζ_eff" begin
        sc, _ = solve(FT, param_set; T_sfc = 275, cap = SF.MaxHeatFluxStabilityCap(), opts)
        legacy = SF.SurfaceFluxConditions(
            sc.shf, sc.lhf, sc.evaporation, sc.ρτxz, sc.ρτyz, sc.ustar, sc.ζ, sc.Cd,
            sc.g_h,
            sc.T_sfc, sc.q_vap_sfc, sc.L_MO, sc.L_eff, sc.converged,
        )
        @test legacy.ζ_eff ≈ sc.ζ_eff
        no_cap = SF.SurfaceFluxConditions(
            sc.shf, sc.lhf, sc.evaporation, sc.ρτxz, sc.ρτyz, sc.ustar, sc.ζ, sc.Cd,
            sc.g_h,
            sc.T_sfc, sc.q_vap_sfc, sc.L_MO, sc.converged,
        )
        @test no_cap.ζ_eff == sc.ζ
        @test no_cap.L_eff == sc.L_MO
        neutral = SF.SurfaceFluxConditions(
            FT(0), FT(0), FT(0), FT(0), FT(0), FT(0.3), FT(0), FT(1e-3), FT(0.01),
            FT(290),
            FT(0.005), FT(Inf), true,
        )
        @test neutral.ζ_eff == 0
        # Neutral stability with a finite effective length stays finite
        neutral_capped = SF.SurfaceFluxConditions(
            FT(0), FT(0), FT(0), FT(0), FT(0), FT(0.3), FT(0), FT(1e-3), FT(0.01),
            FT(290), FT(0.005), FT(Inf), FT(100), true,
        )
        @test neutral_capped.ζ_eff == 0
    end

    @testset "Interior vapor" begin
        _, inputs = solve(FT, param_set; T_sfc = 290)
        @test SF.interior_vapor_specific_humidity(inputs) == inputs.q_tot_int
        wet = (; inputs..., q_liq_int = FT(1e-4), q_ice_int = FT(2e-4))
        @test SF.interior_vapor_specific_humidity(wet) ≈ inputs.q_tot_int - FT(3e-4)
    end

    @testset "Float32, allocations, type stability, duals" begin
        param_set32 = SFP.SurfaceFluxesParameters(Float32, UF.BusingerParams)
        sc32, inputs32 = solve(Float32, param_set32; T_sfc = 295)
        s32 = SF.screen_level_values(param_set32, sc32, inputs32, 2.0f0, 10.0f0)
        @test s32.T isa Float32 && s32.q isa Float32 && s32.u isa Float32
        screen_level_checked(param_set32, sc32, inputs32)

        # Sensitivity of the screen temperature to the surface temperature propagates
        sc, inputs = solve(FT, param_set; T_sfc = 295)
        dT = ForwardDiff.derivative(FT(295)) do T_sfc
            scd = SF.SurfaceFluxConditions(
                sc.shf, sc.lhf, sc.evaporation, sc.ρτxz, sc.ρτyz, sc.ustar, sc.ζ,
                sc.Cd,
                sc.g_h, T_sfc, sc.q_vap_sfc, sc.L_MO, sc.L_eff, sc.ζ_eff,
                sc.converged,
            )
            SF.screen_level_values(param_set, scd, inputs, FT(2), FT(10)).T
        end
        @test 0 < dT < 1
    end
end

end # module
