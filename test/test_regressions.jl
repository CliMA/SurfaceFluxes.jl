module TestRegressions

using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
import Thermodynamics as TD

# Include the generated case definitions
include(joinpath(@__DIR__, "regression_cases.jl"))

const CASE_NUMERIC_FIELDS = (
    :shf,
    :lhf,
    :ustar,
    :Cd,
    :Ch,
    :evaporation,
)

function inputs_from_case(case, ::Type{FT}) where {FT}
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
        roughness = case.roughness_config, # Updated to use the struct directly
        gustiness = SF.ConstantGustinessSpec(FT(1.0)),
        moisture_model = SF.MoistModel(),
    )
end

function assert_coefficient_reasonableness(result, ::Type{FT}) where {FT}
    @test result.ustar >= FT(0)
    @test result.Cd >= FT(0)
    @test result.Ch >= FT(0)
end

@testset "Numerical regression cases" begin
    FT = Float32
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    for case in case_definitions(FT)
        @testset "$(case.name)" begin
            inputs = inputs_from_case(case, FT)
            config = SF.SurfaceFluxConfig(
                inputs.roughness,
                inputs.gustiness,
                inputs.moisture_model,
            )

            #options = SF.SolverOptions{FT}(maxiter = 10, tol = FT(1e-2))
            options = SF.SolverOptions{FT}() # use default options

            result = SF.surface_fluxes(
                param_set,
                inputs.T_int, inputs.q_tot_int, inputs.ρ_int,
                inputs.Ts, inputs.qs,
                inputs.Φs, inputs.Δz, inputs.d,
                inputs.u_int, inputs.u_sfc,
                nothing,
                config,
                SF.PointValueScheme(),
                options,
            )

            for field in CASE_NUMERIC_FIELDS
                expected_value = getfield(case.expected, field)
                actual_value = getproperty(result, field)

                if isinf(actual_value)
                    @test expected_value > FT(1 / eps(FT))
                else
                    # Use physically reasonable tolerances (10% rtol, 0.1 atol, which is mostly relevant for ustar)
                    is_approx = isapprox(
                        actual_value,
                        expected_value;
                        rtol = FT(0.1),
                        atol = FT(1e-1),
                    )
                    if !is_approx
                        println(
                            "Field: $field, Expected: $expected_value, Actual: $actual_value",
                        )
                    end
                    @test is_approx
                end
            end

            assert_coefficient_reasonableness(result, FT)
        end
    end
end

end # module
