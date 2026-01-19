module TestExchangeCoefficients

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.Parameters.SurfaceFluxesParameters

# Use Float64 for tests
const FT = Float64

@testset "Exchange Coefficients" begin
    # Setup parameters
    param_set = SurfaceFluxesParameters(FT, UF.BusingerParams)

    # 1. Neutral Limit Check: Cd = (κ / ln(z/z0))^2
    @testset "Neutral Limit" begin
        z = FT(10)
        z0m = FT(0.01)
        z0h = FT(0.01)
        Δz = z
        ζ = FT(0) # Neutral

        κ = SFP.von_karman_const(param_set)
        Pr_0 = SFP.Pr_0(param_set)

        expected_Cd = (κ / log(z / z0m))^2
        # F_h includes Pr_0 for Businger
        expected_Ch = κ^2 / (log(z / z0m) * (Pr_0 * log(z / z0h)))

        Cd = SF.drag_coefficient(param_set, ζ, z0m, Δz)
        Ch = SF.heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz)

        @test Cd ≈ expected_Cd atol = 1e-6
        @test Ch ≈ expected_Ch atol = 1e-6

        # Verify positivity
        @test Cd > 0
        @test Ch > 0
    end

    # 2. Roughness Dependence
    @testset "Roughness Dependence" begin
        z = FT(10)
        Δz = z
        ζ = FT(0)

        z0m_small = FT(0.001)
        z0m_large = FT(0.1)

        Cd_small = SF.drag_coefficient(param_set, ζ, z0m_small, Δz)
        Cd_large = SF.drag_coefficient(param_set, ζ, z0m_large, Δz)

        # Larger roughness -> Larger drag coefficient
        @test Cd_large > Cd_small
    end

    # 3. Stability Dependence
    @testset "Stability Dependence" begin
        z = FT(10)
        z0m = FT(0.01)
        Δz = z

        # Unstable (ζ < 0) -> Enhanced mixing -> Higher Cd
        ζ_unstable = FT(-1.0)
        Cd_unstable = SF.drag_coefficient(param_set, ζ_unstable, z0m, Δz)

        # Stable (ζ > 0) -> Suppressed mixing -> Lower Cd
        ζ_stable = FT(1.0)
        Cd_stable = SF.drag_coefficient(param_set, ζ_stable, z0m, Δz)

        # Neutral reference
        Cd_neutral = SF.drag_coefficient(param_set, FT(0), z0m, Δz)

        @test Cd_unstable > Cd_neutral
        @test Cd_stable < Cd_neutral
    end

    # 4. Heat Conductance Consistency
    @testset "Heat Conductance" begin
        z = FT(10)
        z0m = FT(0.01)
        z0h = FT(0.01)
        Δz = z
        ζ = FT(-0.5)
        # Calculate ustar consistent with speed=10.0 and ζ=-0.5
        speed_target = FT(10.0)
        Cd_target = SF.drag_coefficient(param_set, ζ, z0m, Δz)
        ustar_consistent = speed_target * sqrt(Cd_target)

        # Create minimal inputs using builder
        inputs = SF.build_surface_flux_inputs(
            FT(300), FT(0), FT(0), FT(0), FT(1.2), # T_int, q...
            FT(300), FT(0), FT(0), FT(10), FT(0), # T_sfc, q_vap, Φ, Δz, d 
            (speed_target, FT(0)), (FT(0), FT(0)), # u_int, u_sfc
            SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(z0m, z0h),
                SF.ConstantGustinessSpec(FT(0.0)),
                SF.DryModel(),
            ),
            nothing,
            SF.FluxSpecs{FT}(ustar = ustar_consistent), # Consistently prescribed ustar
            nothing, nothing,
        )

        # heat_conductance(param_set, ζ, ustar, inputs, z0m, z0h, scheme)
        gh = SF.heat_conductance(param_set, ζ, ustar_consistent, inputs, z0m, z0h)

        # Direct calculation
        Ch = SF.heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz)

        # Back-calculate speed implied by conductance
        implied_speed = gh / Ch

        @test implied_speed ≈ speed_target atol = 1e-5
    end

    # 5. Drag Coefficient from Speed (Inputs/Speed signature)
    @testset "Drag Coefficient (Inputs/Speed)" begin
        # Test signature: drag_coefficient(inputs, speed)
        # Cd = (ustar / speed)^2
        ustar = FT(0.5)
        speed = FT(10.0)

        # Mock inputs with just ustar (Corrected args)
        inputs = SF.build_surface_flux_inputs(
            FT(300), FT(0), FT(0), FT(0), FT(1.2),
            FT(300), FT(0), FT(0), FT(10), FT(0),
            (FT(10), FT(0)), (FT(0), FT(0)),
            SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(FT(0.01), FT(0.01)),
                SF.ConstantGustinessSpec(FT(0.0)),
                SF.DryModel(),
            ),
            nothing,
            SF.FluxSpecs{FT}(ustar = ustar),
            nothing, nothing,
        )

        Cd = SF.drag_coefficient(inputs, speed)
        expected = (ustar / speed)^2
        @test Cd ≈ expected atol = 1e-8
    end

end

end # module
