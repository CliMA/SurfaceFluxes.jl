# Accessors for the floor of a gustiness model and for the model without it.

module TestGustinessFloor

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

@testset "Gustiness floor accessors" begin
    for FT in (Float32, Float64)
        c = SF.ConstantGustinessSpec(FT(2))
        f = SF.FlooredDeardorffGustinessSpec(FT(1))
        dd = SF.DeardorffGustinessSpec()
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
        @test SF.minimum_wind_speed(c, param_set) === FT(2)
        @test SF.minimum_wind_speed(f, param_set) === FT(1)
        @test SF.minimum_wind_speed(dd, param_set) === FT(0)
        # The floor takes the floating-point type of the parameter set
        @test SF.minimum_wind_speed(SF.ConstantGustinessSpec(2.0), param_set) === FT(2)
        @test SF.without_floor(c) === SF.ConstantGustinessSpec(FT(0))
        @test SF.without_floor(f) === SF.FlooredDeardorffGustinessSpec(FT(0))
        @test SF.without_floor(dd) === dd
        @test SF.minimum_wind_speed(SF.without_floor(f), param_set) === FT(0)

        # The floor is the gustiness in stable conditions, and the model without it
        # has none there
        B = FT(-0.01)
        @test SF.gustiness_value(f, param_set, B) == SF.minimum_wind_speed(f, param_set)
        @test SF.gustiness_value(SF.without_floor(f), param_set, B) == 0
        @test SF.gustiness_value(SF.without_floor(c), param_set, B) == 0
    end
end

end # module
