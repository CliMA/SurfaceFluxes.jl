# Tests for gustiness parameterizations
#
# Verifies:
# 1. ConstantGustinessSpec returns constant values
# 2. DeardorffGustinessSpec scales with buoyancy flux
# 3. Gustiness affects windspeed calculation correctly

module TestGustiness

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

@testset "Gustiness Parameterizations" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    @testset "ConstantGustinessSpec" begin
        gust_val = FT(2.0)
        spec = SF.ConstantGustinessSpec(gust_val)

        # Should return constant regardless of buoyancy flux
        @test SF.gustiness_value(spec, param_set, FT(0)) == gust_val
        @test SF.gustiness_value(spec, param_set, FT(0.1)) == gust_val
        @test SF.gustiness_value(spec, param_set, FT(-0.01)) == gust_val
    end

    @testset "DeardorffGustinessSpec" begin
        spec = SF.DeardorffGustinessSpec()
        β = SFP.gustiness_coeff(param_set)
        zi = SFP.gustiness_zi(param_set)

        # Stable: B <= 0 => gustiness = 0
        @test SF.gustiness_value(spec, param_set, FT(0)) == FT(0)
        @test SF.gustiness_value(spec, param_set, FT(-0.1)) == FT(0)

        # Unstable: B > 0 => gustiness = β * (B * zi)^(1/3)
        B_unstable = FT(0.01)
        expected = β * cbrt(B_unstable * zi)
        @test SF.gustiness_value(spec, param_set, B_unstable) ≈ expected

        # Larger B => larger gustiness
        B_large = FT(0.1)
        @test SF.gustiness_value(spec, param_set, B_large) >
              SF.gustiness_value(spec, param_set, B_unstable)
    end

    @testset "Gustiness in Windspeed" begin
        # Test that gustiness affects effective windspeed
        Δu = (FT(5.0), FT(0.0))  # 5 m/s in x direction
        mean_speed = hypot(Δu...)

        # With gustiness > mean speed, gustiness wins
        gustiness_high = FT(10.0)
        @test SF.windspeed(Δu, gustiness_high) == gustiness_high

        # With gustiness < mean speed, mean speed wins
        gustiness_low = FT(1.0)
        @test SF.windspeed(Δu, gustiness_low) == mean_speed

        # Zero wind uses gustiness as floor
        Δu_zero = (FT(0.0), FT(0.0))
        @test SF.windspeed(Δu_zero, gustiness_high) == gustiness_high
    end
end

end # module
