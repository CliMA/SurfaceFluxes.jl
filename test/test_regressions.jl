import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import ClimaParams as CP
import SurfaceFluxes.Parameters as SFP
include("test_utils.jl")

"""
    case_definitions(::Type{FT}) -> Vector

Numerical regression cases derived from challenging configurations.
Each case stores the state specification and the expected `SurfaceFluxConditions`
returned by `surface_fluxes`. The expected values were generated using 
SurfaceFluxes v0.13.1 and provide a regression target for changes.
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
            ts_int = TD.PhaseEquil{FT}(
                1.2595116f0,
                99902.82f0,
                12337.749f0,
                0.0044478197f0,
                275.624f0,
            ),
            ts_sfc = TD.PhaseEquil{FT}(
                1.2544012f0,
                99335.55f0,
                11996.086f0,
                0.0044396375f0,
                275.1768f0,
            ),
            expected = (;  # Regression targets
                L_MO = FT(1.0154548e6),
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
            ts_int = TD.PhaseEquil{FT}(
                1.2605726f0,
                100331.47f0,
                7956.4053f0,
                0.002202735f0,
                276.95068f0,
            ),
            ts_sfc = TD.PhaseEquil{FT}(
                1.2499729f0,
                99303.92f0,
                13258.002f0,
                0.0047157165f0,
                276.01752f0,
            ),
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
            ts_int = TD.PhaseEquil{FT}(
                1.1730341f0,
                98689.72f0,
                43302.703f0,
                0.012817842f0,
                290.8733f0,
            ),
            ts_sfc = TD.PhaseEquil{FT}(
                1.1740736f0,
                98819.375f0,
                43671.336f0,
                0.012941063f0,
                290.97592f0,
            ),
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
            ts_int = TD.PhaseEquil{FT}(
                1.1698402f0,
                98647.89f0,
                44855.285f0,
                0.013289474f0,
                291.46088f0,
            ),
            ts_sfc = TD.PhaseEquil{FT}(
                1.1708081f0,
                98770.55f0,
                45266.523f0,
                0.013432522f0,
                291.5569f0,
            ),
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
            ts_int = TD.PhaseEquil{FT}(
                1.2182463f0,
                96874.9f0,
                13805.914f0,
                0.0048752176f0,
                276.25174f0,
            ),
            ts_sfc = TD.PhaseEquil{FT}(
                1.2197124f0,
                97042.68f0,
                14087.365f0,
                0.004953378f0,
                276.38446f0,
            ),
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
            ts_int = TD.PhaseEquil{FT}(
                1.2317619f0,
                99965.086f0,
                12921.355f0,
                0.0026684932f0,
                282.31366f0,
            ),
            ts_sfc = TD.PhaseEquil{FT}(
                1.214932f0,
                98294.87f0,
                21252.703f0,
                0.006637053f0,
                280.76575f0,
            ),
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
    :buoy_flux,
    Symbol("ρτxz"),
    Symbol("ρτyz"),
    :ustar,
    :Cd,
    :Ch,
    :evaporation,
)

function build_values_only_case(case, ::Type{FT}) where {FT}
    state_sfc = SF.StateValues(FT(0), case.u_sfc, case.ts_sfc)
    state_int = SF.StateValues(case.Δz, case.u_int, case.ts_int)
    sc = SF.ValuesOnly(state_int, state_sfc, case.z0m, case.z0h)
    return sc, state_int, state_sfc
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
            sc, state_int, state_sfc = build_values_only_case(case, FT)
            result = surface_fluxes_wrapper(param_set, sc)

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
    base_sc, state_int, state_sfc = build_values_only_case(base_case, FT)
    base_result = surface_fluxes_wrapper(param_set, base_sc)

    z0m, z0h = base_case.z0m, base_case.z0h
    @testset "Flux-prescribed container" begin
        flux_sc = SF.Fluxes(
            state_int,
            state_sfc,
            base_result.shf,
            base_result.lhf,
            z0m,
            z0h,
        )
        flux_result = surface_fluxes_wrapper(param_set, flux_sc)
        @test isapprox(flux_result.L_MO, base_result.L_MO; rtol = FT(1e-3))
    end

    @testset "Flux+ustar container" begin
        fluxustar_sc = SF.FluxesAndFrictionVelocity(
            state_int,
            state_sfc,
            base_result.shf,
            base_result.lhf,
            base_result.ustar,
            z0m,
            z0h,
        )
        result_fluxustar = surface_fluxes_wrapper(param_set, fluxustar_sc)
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

    @testset "Coefficient-prescribed container" begin
        coeff_sc = SF.Coefficients(
            state_int,
            state_sfc,
            base_result.Cd,
            base_result.Ch,
            z0m,
            z0h,
        )
        coeff_result = surface_fluxes_wrapper(param_set, coeff_sc)
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
