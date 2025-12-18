module TestRegressions

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP
import SurfaceFluxes.Parameters as SFP
using Test
include("test_utils.jl")

"""
    case_definitions(::Type{FT}) -> Vector

Numerical regression cases derived from challenging configurations.
Each case stores the state specification and the expected `SurfaceFluxConditions`
returned by `surface_fluxes`. The expected values were generated using 
SurfaceFluxes v0.13.1 and provide a regression target for changes.

State values are stored as direct primitives (T, q_tot, ρ) for simplicity.
"""
function case_definitions(::Type{FT}) where {FT}
    return [
        (
            name = "Strong wind shear, strongly stable #1",
            Δz = FT(15.000001),
            u_int = (FT(-19.07545), FT(16.88031)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0h = FT(1e-5),
            # Interior state: ρ=1.2595116, T=275.624, q_tot=0.0044478197
            T_int = FT(275.624),
            q_tot_int = FT(0.0044478197),
            ρ_int = FT(1.2595116),
            # Surface state: ρ=1.2544012, T=275.1768, q_tot=0.0044396375
            T_sfc = FT(275.1768),
            q_sfc = FT(0.0044396375),
            expected = (;  # Regression targets
                L_MO = FT(906302.9),
                shf = FT(-20.373146),
                lhf = FT(-0.69898504),
                buoy_flux = FT(-0.0005755669),
                ρτxz = FT(0.48220262),
                ρτyz = FT(-0.4267123),
                ustar = FT(0.71645653),
                Cd = FT(0.00079114654),
                Ch = FT(0.0010691151),
                evaporation = FT(-2.7950458e-7),
            ),
        ),
        (
            name = "Weak wind, weakly stratified #2",
            Δz = FT(15.000001),
            u_int = (FT(-0.168524), FT(-0.000566946)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0h = FT(1e-5),
            T_int = FT(276.95068),
            q_tot_int = FT(0.002202735),
            ρ_int = FT(1.2605726),
            T_sfc = FT(276.01752),
            q_sfc = FT(0.0047157165),
            expected = (;
                L_MO = FT(-5.2631125),
                shf = FT(-1.9265472),
                lhf = FT(11.302928),
                buoy_flux = FT(-3.274816e-5),
                ρτxz = FT(0.00021387689),
                ρτyz = FT(7.195215e-7),
                ustar = FT(0.031864032),
                Cd = FT(0.0010153166),
                Ch = FT(0.0014388718),
                evaporation = FT(4.519725e-6),
            ),
        ),
        (
            name = "Strong wind shear, strongly stable #3",
            Δz = FT(15.000001),
            u_int = (FT(-14.154735), FT(-5.1905923)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0h = FT(1e-5),
            T_int = FT(290.8733),
            q_tot_int = FT(0.012817842),
            ρ_int = FT(1.1730341),
            T_sfc = FT(290.97592),
            q_sfc = FT(0.012941063),
            expected = (;
                L_MO = FT(-22733.092),
                shf = FT(-0.75684166),
                lhf = FT(5.833838),
                buoy_flux = FT(-9.548945e-6),
                ρτxz = FT(0.19829303),
                ρτyz = FT(0.07271477),
                ustar = FT(0.42413536),
                Cd = FT(0.00079142884),
                Ch = FT(0.0010695357),
                evaporation = FT(2.3327887e-6),
            ),
        ),
        (
            name = "Strong wind shear, strongly stable #4",
            Δz = FT(15.000001),
            u_int = (FT(-13.526638), FT(-8.794365)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0h = FT(1e-5),
            T_int = FT(291.46088),
            q_tot_int = FT(0.013289474),
            ρ_int = FT(1.1698402),
            T_sfc = FT(291.5569),
            q_sfc = FT(0.013432522),
            expected = (;
                L_MO = FT(-22413.842),
                shf = FT(-0.9256648),
                lhf = FT(7.227552),
                buoy_flux = FT(-1.1460169e-5),
                ρτxz = FT(0.20222539),
                ρτyz = FT(0.13147715),
                ustar = FT(0.4538926),
                Cd = FT(0.00079143274),
                Ch = FT(0.0010695414),
                evaporation = FT(2.8900959e-6),
            ),
        ),
        (
            name = "High wind, near-neutral #5",
            Δz = FT(15.000001),
            u_int = (FT(-41.34482), FT(-23.609104)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0h = FT(1e-5),
            T_int = FT(276.25174),
            q_tot_int = FT(0.0048752176),
            ρ_int = FT(1.2182463),
            T_sfc = FT(276.38446),
            q_sfc = FT(0.004953378),
            expected = (;
                L_MO = FT(-351333.28),
                shf = FT(-0.8296894),
                lhf = FT(12.135787),
                buoy_flux = FT(-2.043232e-7),
                ρτxz = FT(1.8995668),
                ρτyz = FT(1.0847083),
                ustar = FT(1.3391849),
                Cd = FT(0.00079117215),
                Ch = FT(0.0010691539),
                evaporation = FT(4.852762e-6),
            ),
        ),
        (
            name = "Weak wind, stable #6",
            Δz = FT(15.000001),
            u_int = (FT(-0.75088084), FT(-0.09317328)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0h = FT(1e-5),
            T_int = FT(282.31366),
            q_tot_int = FT(0.0026684932),
            ρ_int = FT(1.2317619),
            T_sfc = FT(280.76575),
            q_sfc = FT(0.006637053),
            expected = (;
                L_MO = FT(-3.385903),
                shf = FT(-3.034892),
                lhf = FT(18.299839),
                buoy_flux = FT(-5.0476297e-5),
                ρτxz = FT(0.0009687537),
                ρτyz = FT(0.000120208104),
                ustar = FT(0.032587063),
                Cd = FT(0.0010619167),
                Ch = FT(0.0015176912),
                evaporation = FT(7.317594e-6),
            ),
        ),
    ]
end

const CASE_NUMERIC_FIELDS = (
    :L_MO,
    :shf,
    :lhf,
    # :buoy_flux, # Removed from output struct
    Symbol("ρτxz"),
    Symbol("ρτyz"),
    :ustar,
    :Cd,
    :Ch,
    :evaporation,
)

function inputs_from_case(case, ::Type{FT}, param_set) where {FT}
    # Use direct primitives from case definition
    return (;
        T_int = case.T_int,
        q_tot_int = case.q_tot_int,
        ρ_int = case.ρ_int,
        Ts = case.T_sfc,
        qs = case.q_sfc,
        Φs = FT(0),
        Δz = case.Δz,
        d = FT(0),
        u_int = case.u_int,
        u_sfc = case.u_sfc,
        roughness = SF.roughness_lengths(case.z0m, case.z0h),
        gustiness = SF.ConstantGustinessSpec(FT(1.0)),
        moisture_model = SF.MoistModel(),
    )
end

function assert_coefficient_reasonableness(result, ::Type{FT}) where {FT}
    @test result.ustar >= FT(0)
    for coeff in (result.Cd, result.Ch)
        @test coeff > FT(0)
        @test coeff < FT(0.01)
    end
end


@testset "Numerical regression cases" begin
    FT = Float32
    REGRESSION_RTOL = FT(0.01)
    REGRESSION_ATOL = FT(0.05)

    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    for case in case_definitions(FT)
        @testset "$(case.name)" begin
            inputs = inputs_from_case(case, FT, param_set)
            config = SF.SurfaceFluxConfig(inputs.roughness, inputs.gustiness, inputs.moisture_model)

            result = SF.surface_fluxes(
                param_set,
                inputs.T_int, inputs.q_tot_int, inputs.ρ_int,
                inputs.Ts, inputs.qs,
                inputs.Φs, inputs.Δz, inputs.d,
                inputs.u_int, inputs.u_sfc,
                nothing,
                config,
            )

            for field in CASE_NUMERIC_FIELDS
                expected_value = getfield(case.expected, field)
                actual_value = getproperty(result, field)
                # Infinite L_MO is permissible and would 
                # result in failing ATOL, RTOL checks.
                if isinf(actual_value)
                    @test expected_value > FT(1 / eps(FT))
                else
                    @test isapprox(
                        actual_value,
                        expected_value;
                        rtol = REGRESSION_RTOL,
                        atol = REGRESSION_ATOL,
                    )
                end
            end

            assert_coefficient_reasonableness(result, FT)
        end
    end
end

@testset "Surface condition container consistency" begin
    FT = Float32
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    base_case = case_definitions(FT)[1]
    inputs = inputs_from_case(base_case, FT, param_set)
    config = SF.SurfaceFluxConfig(inputs.roughness, inputs.gustiness, inputs.moisture_model)

    base_result = SF.surface_fluxes(
        param_set,
        inputs.T_int, inputs.q_tot_int, inputs.ρ_int,
        inputs.Ts, inputs.qs,
        inputs.Φs, inputs.Δz, inputs.d,
        inputs.u_int, inputs.u_sfc,
        nothing,
        config,
    )

    z0m, z0h = base_case.z0m, base_case.z0h

    @testset "Flux-prescribed" begin
        # Pass fluxes via FluxSpecs
        flux_specs = SF.FluxSpecs(FT; shf = base_result.shf, lhf = base_result.lhf)

        flux_result = SF.surface_fluxes(
            param_set,
            inputs.T_int, inputs.q_tot_int, inputs.ρ_int,
            inputs.Ts, inputs.qs,
            inputs.Φs, inputs.Δz, inputs.d,
            inputs.u_int, inputs.u_sfc,
            nothing,
            config,
            SF.PointValueScheme(),
            nothing, # solver_opts
            flux_specs,
        )
        @test isapprox(flux_result.L_MO, base_result.L_MO; rtol = FT(1e-3))
    end

    @testset "Flux+ustar prescribed" begin
        flux_specs = SF.FluxSpecs(FT; shf = base_result.shf, lhf = base_result.lhf, ustar = base_result.ustar)

        result_fluxustar = SF.surface_fluxes(
            param_set,
            inputs.T_int, inputs.q_tot_int, inputs.ρ_int,
            inputs.Ts, inputs.qs,
            inputs.Φs, inputs.Δz, inputs.d,
            inputs.u_int, inputs.u_sfc,
            nothing,
            config,
            SF.PointValueScheme(),
            nothing, # solver_opts
            flux_specs,
        )
        @test isapprox(
            result_fluxustar.ustar,
            base_result.ustar;
            rtol = sqrt(eps(FT)),
            atol = FT(1e-9),
        )
        @test isapprox(
            result_fluxustar.shf,
            base_result.shf;
            rtol = sqrt(eps(FT)),
            atol = FT(1e-9),
        )
        @test isapprox(
            result_fluxustar.lhf,
            base_result.lhf;
            rtol = sqrt(eps(FT)),
            atol = FT(1e-9),
        )
    end

    @testset "Coefficient-prescribed" begin
        flux_specs = SF.FluxSpecs(FT; Cd = base_result.Cd, Ch = base_result.Ch)

        coeff_result = SF.surface_fluxes(
            param_set,
            inputs.T_int, inputs.q_tot_int, inputs.ρ_int,
            inputs.Ts, inputs.qs,
            inputs.Φs, inputs.Δz, inputs.d,
            inputs.u_int, inputs.u_sfc,
            nothing,
            config,
            SF.PointValueScheme(),
            nothing, # solver_opts
            flux_specs,
        )
        @test isapprox(
            coeff_result.ustar,
            base_result.ustar;
            rtol = FT(1e-3),
        )
        @test isapprox(
            coeff_result.Cd,
            base_result.Cd;
            rtol = sqrt(eps(FT)),
            atol = FT(1e-9),
        )
        @test isapprox(
            coeff_result.Ch,
            base_result.Ch;
            rtol = sqrt(eps(FT)),
            atol = FT(1e-9),
        )
    end
end

end # module
