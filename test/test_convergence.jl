using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
import Thermodynamics as TD

import RootSolvers as RS

const SYNTH_T_SURFACE = (261.0, 282.0, 310.0)
const SYNTH_DELTA_T = (-8.0, -2.0, 3.5)
const SYNTH_RH_INT = (0, 0.85, 0.99)
const SYNTH_DELTA_QT = (-2e-3, 0.0, 3e-3)
const SYNTH_HEIGHTS = (5.0, 20.0, 80.0)
const SYNTH_SPEEDS = (1.0, 12.0)
const SYNTH_WIND_DIRS = ((1.0, 0.0), (-0.7, 0.4))
const SYNTH_Z0M = (1e-5, 8e-4, 3e-3)
const SYNTH_Z0h_FACTORS = (0.1, 1.0, 10.0)
const SYNTH_PRESSURES = (9.5e4, 1.0e5)

const TEMP_NEUTRAL_THRESHOLD = 0.3
const HUMIDITY_NEUTRAL_THRESHOLD = 1e-4

# Helper function to compute specific humidity from relative humidity
function qt_from_RH(::Type{FT}, RH, T, p) where {FT}
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    # Saturation vapor pressure over liquid water [Pa]
    e_sat = TD.saturation_vapor_pressure(thermo_params, FT(T), TD.Liquid())

    # Saturation specific humidity: q_sat = ε * e_sat / (p - (1 - ε) * e_sat)
    ε = FT(0.622)  # Rd / Rv
    q_sat = ε * e_sat / (p - (FT(1) - ε) * e_sat)

    return FT(RH * q_sat)
end

function synthetic_cases(::Type{FT}) where {FT}
    cases = NamedTuple[]
    for T_sfc in SYNTH_T_SURFACE,
        ΔT in SYNTH_DELTA_T,
        RH_int in SYNTH_RH_INT,
        Δqt in SYNTH_DELTA_QT,
        z in SYNTH_HEIGHTS,
        speed in SYNTH_SPEEDS,
        wind_dir in SYNTH_WIND_DIRS,
        z0m in SYNTH_Z0M,
        z0h_factor in SYNTH_Z0h_FACTORS,
        p in SYNTH_PRESSURES

        T_int = FT(T_sfc + ΔT)
        qt_int = qt_from_RH(FT, RH_int, T_int, p)
        push!(
            cases,
            (
                T_sfc = FT(T_sfc),
                T_int = T_int,
                qt_int = qt_int,
                qt_sfc = FT(clamp(qt_int + Δqt, 0.0, 1.0)),
                z = FT(z),
                wind = (
                    FT(speed * wind_dir[1]),
                    FT(speed * wind_dir[2]),
                ),
                z0m = FT(z0m),
                z0h = FT(max(z0m * z0h_factor, 1e-6)),
                pressure = FT(p),
            ),
        )
    end
    return cases
end

function density_from_state(thermo_params, T, pressure, qt)
    FT = eltype(thermo_params)
    R_m = TD.gas_constant_air(thermo_params, qt, FT(0), FT(0))
    return pressure / (R_m * T)
end

function build_surface_inputs(param_set, case, config)
    FT = eltype(case.T_sfc)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_int = density_from_state(thermo_params, case.T_int, case.pressure, case.qt_int)

    Φs = FT(0)  # Surface geopotential
    Δz = case.z

    return SF.build_surface_flux_inputs(
        case.T_int,
        case.qt_int,
        FT(0), # q_liq_int
        FT(0), # q_ice_int
        ρ_int,
        case.T_sfc,
        case.qt_sfc,
        Φs,
        Δz,
        FT(0), # d
        case.wind,
        (FT(0), FT(0)), # u_sfc
        config,
        nothing, # roughness_inputs
        SF.FluxSpecs{FT}(),
        nothing, # update_T_sfc
        nothing, # update_q_vap_sfc
    )
end

function compute_ΔDSEᵥ(param_set, T_sfc, q_sfc, T_int, q_int, Φ_sfc, Φ_int)
    thermo_params = SFP.thermodynamics_params(param_set)

    function local_ΔDSEᵥ(param_set, T_int, q_int, Φ_int, T_sfc, q_sfc, Φ_sfc)
        cp_d = SFP.cp_d(param_set)
        Tv_int = TD.virtual_temperature(thermo_params, T_int, q_int)
        Tv_sfc = TD.virtual_temperature(thermo_params, T_sfc, q_sfc)
        DSEv_int = cp_d * Tv_int + Φ_int
        DSEv_sfc = cp_d * Tv_sfc + Φ_sfc
        return DSEv_int - DSEv_sfc
    end

    return local_ΔDSEᵥ(param_set, T_int, q_int, Φ_int, T_sfc, q_sfc, Φ_sfc)
end

function assert_flux_expectations(result, case, FT, param_set, inputs)
    Δqt = case.qt_int - case.qt_sfc

    grav = SFP.grav(param_set)
    Φ_int = grav * inputs.Δz
    Φ_sfc = inputs.Φ_sfc # 0

    ΔDSEᵥ = compute_ΔDSEᵥ(
        param_set,
        case.T_sfc, case.qt_sfc,
        case.T_int, case.qt_int,
        Φ_sfc, Φ_int,
    )

    heat_tolerance = FT(SFP.cp_d(param_set) * TEMP_NEUTRAL_THRESHOLD)
    if abs(ΔDSEᵥ) > heat_tolerance
        expected_heat_sign = -sign(ΔDSEᵥ)
        if expected_heat_sign == 0
            @test isapprox(result.shf, FT(0); atol = FT(5e-2))
        else
            @test sign(result.shf) == expected_heat_sign
        end
    else
        @test isapprox(result.shf, FT(0); atol = FT(5e-2))
    end
    if abs(Δqt) > FT(HUMIDITY_NEUTRAL_THRESHOLD)
        expected_evap_sign = sign(case.qt_sfc - case.qt_int)
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

@testset "SurfaceFluxes Convergence Matrix" begin
    schemes = (SF.PointValueScheme(), SF.LayerAverageScheme())

    # Counter for convergence statistics (per UF type)
    converged_count = Dict{Type, Int}()
    total_count = Dict{Type, Int}()

    # Storage for Ri_b of failed cases (per UF type)
    failed_Rib = Dict{Type, Vector{Float64}}()

    for uf_params in (UF.BusingerParams, UF.GryanikParams, UF.GrachevParams)
        converged_count[uf_params] = 0
        total_count[uf_params] = 0
        failed_Rib[uf_params] = Float64[]
    end

    for FT in (Float32, Float64)
        # Define roughness configs as functions of (z0m, z0h)
        roughness_config_factories = (
            (z0m, z0h) -> SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(FT(z0m), FT(z0h)),
                SF.ConstantGustinessSpec(FT(1.0)),
            ),
            (z0m, z0h) -> SF.SurfaceFluxConfig(
                SF.COARE3RoughnessParams{FT}(),
                SF.ConstantGustinessSpec(FT(1.0)),
            ),
        )
        cases = synthetic_cases(FT)
        for uf_params in (UF.BusingerParams, UF.GryanikParams, UF.GrachevParams)
            param_set = SFP.SurfaceFluxesParameters(FT, uf_params)
            thermo_params = SFP.thermodynamics_params(param_set)
            scheme_set = uf_params === UF.GrachevParams ? (SF.PointValueScheme(),) : schemes
            for case in cases, config_factory in roughness_config_factories,
                scheme in scheme_set

                config = config_factory(case.z0m, case.z0h)
                inputs = build_surface_inputs(param_set, case, config)

                # Solver options for convergence testing
                solver_opts = SF.SolverOptions{FT}(
                    maxiter = 15,
                    tol = 1e-2,
                    forced_fixed_iters = false,
                )

                result = SF.surface_fluxes(
                    param_set,
                    inputs.T_int, inputs.q_tot_int, FT(0), FT(0), inputs.ρ_int,
                    inputs.T_sfc_guess, inputs.q_vap_sfc_guess,
                    inputs.Φ_sfc, inputs.Δz, inputs.d,
                    inputs.u_int, inputs.u_sfc,
                    nothing, # roughness_inputs
                    config,
                    scheme,
                    solver_opts,
                )

                total_count[uf_params] += 1
                if result.converged
                    converged_count[uf_params] += 1
                else
                    # Compute expected Ri_b for failed case
                    ΔU = SF.windspeed(inputs, FT(0))
                    ρ_sfc = density_from_state(
                        thermo_params,
                        inputs.T_sfc_guess,
                        case.pressure,
                        inputs.q_vap_sfc_guess,
                    )

                    Rib = SF.state_bulk_richardson_number(
                        param_set,
                        inputs,
                        inputs.T_sfc_guess,
                        ρ_sfc,
                        ΔU,
                        inputs.q_vap_sfc_guess,
                    )
                    push!(failed_Rib[uf_params], Float64(Rib))
                end

                assert_flux_expectations(result, case, FT, param_set, inputs)
            end
        end
    end

    @info "Convergence matrix exercised Businger/Gryanik/Grachev UFs, Float32/Float64, synthetic dry/moist gradients, Scalar/Charnock roughness"

    # Report statistics for each UF type
    total_all = sum(values(total_count))
    converged_all = sum(values(converged_count))
    @info "Overall Convergence Statistics" Total = total_all Converged = converged_all Failed =
        (total_all - converged_all) FailurePercentage =
        round((1 - converged_all / total_all) * 100; digits = 2)

    import Statistics
    for uf_params in (UF.BusingerParams, UF.GryanikParams, UF.GrachevParams)
        uf_name = string(nameof(uf_params))
        tc = total_count[uf_params]
        cc = converged_count[uf_params]
        fc = tc - cc
        fp = tc > 0 ? round((1 - cc / tc) * 100; digits = 2) : 0.0
        @info "Convergence: $uf_name" Total = tc Converged = cc Failed = fc FailurePercentage =
            fp

        rib_vec = failed_Rib[uf_params]
        if !isempty(rib_vec)
            @info "  Failed Ri_b ($uf_name)" Min = minimum(rib_vec) Max =
                maximum(rib_vec) Median = Statistics.median(rib_vec) Count = length(rib_vec)
        end
    end
end
