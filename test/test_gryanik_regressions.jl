# Numerical regression tests for the Gryanik (2020) universal function parameterization.
#
# These cases mirror the structure of `test_regressions.jl` (which uses Businger) to
# ensure that the Gryanik parameterization also has pinned numerical values guarding
# against accidental changes.

module TestGryanikRegressions

using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

include(joinpath(@__DIR__, "gryanik_regression_cases.jl"))

const CASE_NUMERIC_FIELDS = (
    :shf,
    :lhf,
    :ustar,
    :Cd,
    :g_h,
    :evaporation,
)

@testset "Gryanik Numerical Regression Cases" begin
    FT = Float32
    param_set = SFP.SurfaceFluxesParameters(FT, UF.GryanikParams)

    for case in gryanik_case_definitions(FT)
        @testset "$(case.name)" begin
            config = SF.SurfaceFluxConfig(
                case.roughness_config,
                SF.ConstantGustinessSpec(FT(1.0)),
                SF.MoistModel(),
            )

            options = SF.SolverOptions{FT}()

            result = SF.surface_fluxes(
                param_set,
                case.T_int,
                case.q_tot_int,
                FT(0),
                FT(0),
                case.ρ_int,
                case.T_sfc,
                case.q_sfc,
                FT(0),
                case.Δz,
                FT(0),
                case.u_int,
                case.u_sfc,
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
                    rtol = FT(0.025)
                    atol = if field in (:shf, :lhf)
                        FT(0.5)
                    elseif field == :ustar
                        FT(0.05)
                    elseif field == :Cd
                        FT(5e-7)
                    elseif field == :g_h
                        FT(5e-6)
                    elseif field == :evaporation
                        FT(5e-10)
                    else
                        FT(0.05)
                    end

                    @test isapprox(actual_value, expected_value; rtol, atol)
                end
            end

            # Basic coefficient sanity
            @test result.ustar >= FT(0)
            @test result.Cd >= FT(0)
            @test result.g_h >= FT(0)
        end
    end
end

end # module
