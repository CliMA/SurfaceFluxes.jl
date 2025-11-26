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
            u_in = (FT(-19.07545), FT(16.88031)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0b = FT(1e-5),
            ts_in = TD.PhaseEquil{FT}(
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
                L_MO = FT(2156.5898),
                shf = FT(-20.321943),
                lhf = FT(-0.69523025),
                buoy_flux = FT(-5.7318626f-4),
                ρτxz = FT(0.47999802),
                ρτyz = FT(-0.42476144),
                ustar = FT(0.71481687),
                Cd = FT(7.875295f-4),
                Ch = FT(0.0010633721),
                evaporation = FT(-2.7800314f-7),
            ),
        ),
        (
            name = "Weak wind, weakly stratified #2",
            Δz = FT(15.000001),
            u_in = (FT(-0.168524), FT(-0.000566946)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0b = FT(1e-5),
            ts_in = TD.PhaseEquil{FT}(
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
                L_MO = FT(1.1920929f-7),
                shf = FT(-6.100155f-16),
                lhf = FT(3.5936135f-15),
                buoy_flux = FT(-1.0295639f-20),
                ρτxz = FT(9.6366225f-20),
                ρτyz = FT(3.2419386f-22),
                ustar = FT(1.1920929f-7),
                Cd = FT(4.5746986f-19),
                Ch = FT(4.5746986f-19),
                evaporation = FT(1.4369855f-21),
            ),
        ),
        (
            name = "Strong wind shear, strongly stable #3",
            Δz = FT(15.000001),
            u_in = (FT(-14.154735), FT(-5.1905923)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0b = FT(1e-5),
            ts_in = TD.PhaseEquil{FT}(
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
                L_MO = FT(22423.342),
                shf = FT(-0.6993717),
                lhf = FT(5.82858),
                buoy_flux = FT(-7.745661f-6),
                ρτxz = FT(0.19813664),
                ρτyz = FT(0.07265742),
                ustar = FT(0.42396805),
                Cd = FT(7.908046f-4),
                Ch = FT(0.0010685718),
                evaporation = FT(2.3306861f-6),
            ),
        ),
        (
            name = "Strong wind shear, strongly stable #4",
            Δz = FT(15.000001),
            u_in = (FT(-13.526638), FT(-8.794365)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0b = FT(1e-5),
            ts_in = TD.PhaseEquil{FT}(
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
                L_MO = FT(22546.48),
                shf = FT(-0.8574441),
                lhf = FT(7.22102),
                buoy_flux = FT(-9.300392f-6),
                ρτxz = FT(0.20206536),
                ρτyz = FT(0.1313731),
                ustar = FT(0.45371297),
                Cd = FT(7.908065f-4),
                Ch = FT(0.0010685748),
                evaporation = FT(2.887484f-6),
            ),
        ),
        (
            name = "High wind, near-neutral #5",
            Δz = FT(15.000001),
            u_in = (FT(-41.34482), FT(-23.609104)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0b = FT(1e-5),
            ts_in = TD.PhaseEquil{FT}(
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
                L_MO = FT(2.065302f7),
                shf = FT(-0.5786611),
                lhf = FT(8.980261),
                buoy_flux = FT(9.084079f-7),
                ρτxz = FT(1.899524),
                ρτyz = FT(1.0846839),
                ustar = FT(1.3391694),
                Cd = FT(7.9115434f-4),
                Ch = FT(7.9115434f-4),
                evaporation = FT(3.5909552f-6),
            ),
        ),
        (
            name = "Weak wind, stable #6",
            Δz = FT(15.000001),
            u_in = (FT(-0.75088084), FT(-0.09317328)),
            u_sfc = (FT(0), FT(0)),
            z0m = FT(1e-5),
            z0b = FT(1e-5),
            ts_in = TD.PhaseEquil{FT}(
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
                L_MO = FT(1.1920929f-7),
                shf = FT(-9.024182f-16),
                lhf = FT(5.5160268f-15),
                buoy_flux = FT(-1.4714911f-20),
                ρτxz = FT(4.1733562f-19),
                ρτyz = FT(5.178522f-20),
                ustar = FT(1.1920929f-7),
                Cd = FT(4.5746986f-19),
                Ch = FT(4.5746986f-19),
                evaporation = FT(2.205705f-21),
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
    state_in = SF.StateValues(case.Δz, case.u_in, case.ts_in)
    sc = SF.ValuesOnly(state_in, state_sfc, case.z0m, case.z0b)
    return sc, state_in, state_sfc
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
            sc, state_in, state_sfc = build_values_only_case(case, FT)
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
    base_sc, state_in, state_sfc = build_values_only_case(base_case, FT)
    base_result = surface_fluxes_wrapper(param_set, base_sc)

    z0m, z0b = base_case.z0m, base_case.z0b
    @testset "Flux-prescribed container" begin
        flux_sc = SF.Fluxes(
            state_in,
            state_sfc,
            base_result.shf,
            base_result.lhf,
            z0m,
            z0b,
        )
        flux_result = surface_fluxes_wrapper(param_set, flux_sc)
        @test isapprox(flux_result.L_MO, base_result.L_MO; rtol = FT(1e-3))
    end

    @testset "Flux+ustar container" begin
        fluxustar_sc = SF.FluxesAndFrictionVelocity(
            state_in,
            state_sfc,
            base_result.shf,
            base_result.lhf,
            base_result.ustar,
            z0m,
            z0b,
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
            state_in,
            state_sfc,
            base_result.Cd,
            base_result.Ch,
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
