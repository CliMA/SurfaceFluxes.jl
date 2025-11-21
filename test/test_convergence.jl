using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

import Thermodynamics as TD
using RootSolvers
const RS = RootSolvers

const SYNTH_T_SURFACE = (261.0, 282.0, 310.0)
const SYNTH_DELTA_T = (-8.0, -2.0, 3.5)
const SYNTH_RH_IN = (0, 0.85, 0.99)
const SYNTH_DELTA_QT = (-2e-3, 0.0, 3e-3)
const SYNTH_HEIGHTS = (5.0, 20.0, 80.0)
const SYNTH_SPEEDS = (1.0, 12.0)
const SYNTH_WIND_DIRS = ((1.0, 0.0), (-0.7, 0.4))
const SYNTH_Z0M = (1e-5, 8e-4, 3e-3)
const SYNTH_Z0B_FACTORS = (0.1, 1.0, 10.0)
const SYNTH_PRESSURES = (9.5e4, 1.0e5)

const TEMP_NEUTRAL_THRESHOLD = 0.3
const HUMIDITY_NEUTRAL_THRESHOLD = 1e-4

# Helper function to compute specific humidity from relative humidity
function qt_from_RH(::Type{FT}, RH, T, p) where {FT}
    # Create a minimal param_set to get thermo_params for Thermodynamics calculations
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    # Approximate saturation specific humidity (neglecting density effect of water vapor)
    qt_sat = TD.q_vap_saturation_from_pressure(thermo_params, FT(0), FT(p), FT(T), TD.PhaseEquil)

    # Multiply by relative humidity to get actual specific humidity
    return FT(RH * qt_sat)
end

function synthetic_cases(::Type{FT}) where {FT}
    cases = NamedTuple[]
    for T_sfc in SYNTH_T_SURFACE,
        ΔT in SYNTH_DELTA_T,
        RH_in in SYNTH_RH_IN,
        Δqt in SYNTH_DELTA_QT,
        z in SYNTH_HEIGHTS,
        speed in SYNTH_SPEEDS,
        wind_dir in SYNTH_WIND_DIRS,
        z0m in SYNTH_Z0M,
        z0b_factor in SYNTH_Z0B_FACTORS,
        p in SYNTH_PRESSURES

        T_in = FT(T_sfc + ΔT)
        qt_in = qt_from_RH(FT, RH_in, T_in, p)
        push!(
            cases,
            (
                T_sfc = FT(T_sfc),
                T_in = T_in,
                qt_in = qt_in,
                qt_sfc = FT(clamp(qt_in + Δqt, 0.0, 1.0)),
                z = FT(z),
                wind = (
                    FT(speed * wind_dir[1]),
                    FT(speed * wind_dir[2]),
                ),
                z0m = FT(z0m),
                z0b = FT(max(z0m * z0b_factor, 1e-6)),
                pressure = FT(p),
            ),
        )
    end
    return cases
end

function density_from_state(thermo_params, T, pressure, qt)
    q_pt = TD.PhasePartition_equil_given_p(
        thermo_params,
        T,
        pressure,
        qt,
        TD.PhaseEquil,
    )
    return TD.air_density(thermo_params, T, pressure, q_pt)
end

function build_surface_condition(param_set, case, roughness_model)
    FT = eltype(case.T_sfc)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = density_from_state(thermo_params, case.T_sfc, case.pressure, case.qt_sfc)
    ρ_in = density_from_state(thermo_params, case.T_in, case.pressure, case.qt_in)
    ts_sfc = TD.PhaseEquil_ρTq(thermo_params, ρ_sfc, case.T_sfc, case.qt_sfc)
    ts_in = TD.PhaseEquil_ρTq(thermo_params, ρ_in, case.T_in, case.qt_in)
    state_sfc = SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc)
    state_in = SF.StateValues(case.z, case.wind, ts_in)
    return SF.ValuesOnly(state_in, state_sfc, case.z0m, case.z0b; roughness_model)
end

function assert_flux_expectations(result, case, FT, param_set, sc)
    Δqt = case.qt_in - case.qt_sfc
    ΔDSEᵥ = SF.ΔDSEᵥ(param_set, sc)
    heat_tolerance = FT(SFP.cp_d(param_set) * TEMP_NEUTRAL_THRESHOLD)
    if abs(ΔDSEᵥ) > heat_tolerance
        expected_heat_sign = -sign(ΔDSEᵥ)
        if expected_heat_sign == 0
            @test isapprox(result.shf, FT(0); atol = FT(5e-2))
        else
            @test sign(result.shf) == expected_heat_sign
            @test sign(result.buoy_flux) == expected_heat_sign
        end
    else
        @test isapprox(result.shf, FT(0); atol = FT(5e-2))
    end
    if abs(Δqt) > FT(HUMIDITY_NEUTRAL_THRESHOLD)
        expected_evap_sign = sign(case.qt_sfc - case.qt_in)
        @test sign(result.evaporation) == expected_evap_sign
    end
    for (stress, wind_component) in zip((result.ρτxz, result.ρτyz), case.wind)
        if abs(wind_component) > FT(1e-3)
            @test sign(stress) == -sign(wind_component)
        else
            @test isapprox(stress, FT(0); atol = FT(1e-6))
        end
    end
    @test result.ustar >= FT(0)
    @test result.Cd > FT(0)
    @test result.Ch > FT(0)
end

@testset "SurfaceFluxes convergence matrix" begin
    schemes = (SF.PointValueScheme(), SF.LayerAverageScheme())
    roughness_models = (SF.ScalarRoughness(), SF.CharnockRoughness())
    for FT in (Float32, Float64)
        cases = synthetic_cases(FT)
        for uf_params in (UF.BusingerParams, UF.GryanikParams, UF.GrachevParams)
            param_set = SFP.SurfaceFluxesParameters(FT, uf_params)
            tol_neutral = FT(SF.Parameters.cp_d(param_set) / 10)
            scheme_set = uf_params === UF.GrachevParams ? (SF.PointValueScheme(),) : schemes
            for case in cases, roughness_model in roughness_models, scheme in scheme_set
                sc = build_surface_condition(param_set, case, roughness_model)
                result = SF.surface_conditions(
                    param_set,
                    sc,
                    scheme;
                    maxiter = 15,
                    tol_neutral = tol_neutral,
                    soltype = RS.CompactSolution(),
                )
                assert_flux_expectations(result, case, FT, param_set, sc)
            end
        end
    end
    @info "Convergence matrix exercised Businger/Gryanik/Grachev UFs, Float32/Float64, synthetic dry/moist gradients, Scalar/Charnock roughness"
end
