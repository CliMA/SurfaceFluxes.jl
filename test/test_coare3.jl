# Tests for COARE 3.0 roughness parameterization (Fairall et al. 2003)
#
# Tests cover:
# 1. Charnock parameter piecewise behavior
# 2. Smooth and rough flow limits
# 3. Scalar roughness scaling

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
import SurfaceFluxes:
    COARE3RoughnessParams, momentum_roughness, scalar_roughness, charnock_parameter

@testset "COARE 3.0 Roughness" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    spec = COARE3RoughnessParams{FT}()

    @testset "Charnock Parameter Piecewise Behavior" begin
        α_low = spec.α_low
        α_high = spec.α_high
        u_low = spec.u_low
        u_high = spec.u_high

        # 1. Below low threshold: α = α_low
        @test charnock_parameter(FT(5), α_low, α_high, u_low, u_high) == α_low
        @test charnock_parameter(FT(0), α_low, α_high, u_low, u_high) == α_low
        @test charnock_parameter(u_low, α_low, α_high, u_low, u_high) == α_low

        # 2. Above high threshold: α = α_high
        @test charnock_parameter(FT(20), α_low, α_high, u_low, u_high) == α_high
        @test charnock_parameter(FT(100), α_low, α_high, u_low, u_high) == α_high
        @test charnock_parameter(u_high, α_low, α_high, u_low, u_high) == α_high

        # 3. Linear interpolation in between
        u_mid = (u_low + u_high) / 2
        α_mid = (α_low + α_high) / 2
        @test charnock_parameter(u_mid, α_low, α_high, u_low, u_high) ≈ α_mid

        # 4. Monotonicity
        u_grid = range(FT(0), FT(25), length = 50)
        α_vals = [charnock_parameter(u, α_low, α_high, u_low, u_high) for u in u_grid]
        @test issorted(α_vals)
    end

    @testset "Momentum Roughness Behavior" begin
        grav = SFP.grav(param_set)
        κ = SFP.von_karman_const(param_set)
        ν = spec.kinematic_visc

        # Low u_star: smooth flow dominates
        u_star_low = FT(0.05)
        z0m_low = momentum_roughness(spec, u_star_low, param_set, nothing)

        # Smooth flow term: 0.11 * ν / u_star
        z0_smooth = FT(0.11) * ν / u_star_low
        @test z0m_low > 0
        # For very low u_star, smooth component dominates
        @test z0m_low > z0_smooth * 0.9

        # High u_star: rough flow dominates
        u_star_high = FT(1.0)
        z0m_high = momentum_roughness(spec, u_star_high, param_set, nothing)

        # Rough flow term: α * u_star^2 / g
        # For high u_star, rough dominates
        @test z0m_high > z0m_low

        # Very low u_star: avoid division by zero
        u_star_tiny = eps(FT)
        z0m_tiny = momentum_roughness(spec, u_star_tiny, param_set, nothing)
        @test isfinite(z0m_tiny)
        @test z0m_tiny > 0
    end

    @testset "Scalar Roughness Scaling" begin
        u_star = FT(0.3)
        z0m = momentum_roughness(spec, u_star, param_set, nothing)
        z0s = scalar_roughness(spec, u_star, param_set, nothing)

        @test z0s > 0
        @test z0s < FT(1.1e-4)  # Upper bound from COARE formula

        # Re_star dependence
        ν = spec.kinematic_visc
        Re_star = z0m * u_star / ν

        # z0s = min(1.1e-4, 5.5e-5 * Re_star^-0.6)
        z0s_expected = min(FT(1.1e-4), FT(5.5e-5) * Re_star^FT(-0.6))
        @test z0s ≈ z0s_expected
    end

    @testset "Combined Momentum and Scalar Roughness" begin
        u_star = FT(0.5)
        z0m, z0s = SF.momentum_and_scalar_roughness(spec, u_star, param_set, nothing)

        @test z0m == momentum_roughness(spec, u_star, param_set, nothing)
        @test z0s == scalar_roughness(spec, u_star, param_set, nothing)
    end

    @testset "Wind Speed Dependence via 10m Neutral Profile" begin
        # Higher u_star => higher 10m wind => higher Charnock α => higher z0m
        u_star_range = FT.(0.1:0.2:1.5)
        z0m_vals = [momentum_roughness(spec, u, param_set, nothing) for u in u_star_range]

        # z0m should generally increase with u_star
        # (though not strictly monotonic due to smooth/rough transition)
        @test z0m_vals[end] > z0m_vals[1]
    end
end
