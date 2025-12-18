# Tests for Raupach (1994) canopy roughness parameterization

module TestRaupachRoughness

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

# Mock parameter set for roughness tests
struct MockRoughnessParamSet{FT} <: SFP.AbstractSurfaceFluxesParameters{FT} end
SFP.grav(::MockRoughnessParamSet{FT}) where {FT} = FT(9.81)
SFP.von_karman_const(::MockRoughnessParamSet{FT}) where {FT} = FT(0.4)
SFP.z0m_fixed(::MockRoughnessParamSet{FT}) where {FT} = FT(1e-4)  # Fallback for zero LAI

@testset "Raupach Roughness Parameterization" begin
    FT = Float64
    param_set = MockRoughnessParamSet{FT}()

    # Default Raupach parameters
    spec = SF.RaupachRoughnessParams{FT}()
    u_star = FT(0.3)

    @testset "Basic functionality" begin
        # Normal canopy: LAI = 3, h = 10m
        inputs = (LAI = FT(3.0), h = FT(10.0))

        z0m = SF.momentum_roughness(spec, u_star, param_set, inputs)
        z0s = SF.scalar_roughness(spec, u_star, param_set, inputs)

        @test z0m > 0
        @test z0s > 0
        @test z0s < z0m  # Scalar roughness < momentum roughness (Stanton number < 1)
        @test z0m < inputs.h  # Roughness length < canopy height
    end

    @testset "LAI dependence" begin
        h = FT(10.0)

        # Raupach formula is non-monotonic in LAI due to exponential term
        # At low LAI, z0m increases with LAI
        inputs_low = (LAI = FT(0.5), h = h)
        inputs_mid = (LAI = FT(2.0), h = h)

        z0m_low = SF.momentum_roughness(spec, u_star, param_set, inputs_low)
        z0m_mid = SF.momentum_roughness(spec, u_star, param_set, inputs_mid)

        # At low-to-mid LAI range, z0m increases
        @test z0m_mid > z0m_low

        # All values should be reasonable fractions of canopy height
        @test z0m_low < h
        @test z0m_mid < h
    end


    @testset "Canopy height dependence" begin
        LAI = FT(3.0)

        # Taller canopy => Higher z0m
        inputs_short = (LAI = LAI, h = FT(5.0))
        inputs_tall = (LAI = LAI, h = FT(20.0))

        z0m_short = SF.momentum_roughness(spec, u_star, param_set, inputs_short)
        z0m_tall = SF.momentum_roughness(spec, u_star, param_set, inputs_tall)

        @test z0m_tall > z0m_short
    end

    @testset "Zero/very low LAI fallback" begin
        # With zero LAI, should fall back to fixed roughness
        inputs_zero = (LAI = FT(0.0), h = FT(10.0))
        z0m_zero = SF.momentum_roughness(spec, u_star, param_set, inputs_zero)

        @test z0m_zero ≈ SFP.z0m_fixed(param_set)
    end

    @testset "Stanton number scaling" begin
        inputs = (LAI = FT(3.0), h = FT(10.0))

        z0m, z0s = SF.momentum_and_scalar_roughness(spec, u_star, param_set, inputs)

        # z0s = z0m * stanton_number
        @test z0s ≈ z0m * spec.stanton_number
    end
end

end # module
