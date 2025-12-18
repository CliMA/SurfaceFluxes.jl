# Tests for Obukhov length and stability parameter helper functions

module TestObukhovHelpers

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

@testset "Obukhov Helper Functions" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    κ = SFP.von_karman_const(param_set)

    @testset "obukhov_length" begin
        ustar = FT(0.3)

        # Unstable: B > 0 => L < 0
        B_unstable = FT(0.01)
        L_unstable = SF.obukhov_length(param_set, ustar, B_unstable)
        @test L_unstable < 0
        @test L_unstable ≈ -ustar^3 / (κ * B_unstable)

        # Stable: B < 0 => L > 0
        B_stable = FT(-0.005)
        L_stable = SF.obukhov_length(param_set, ustar, B_stable)
        @test L_stable > 0
        @test L_stable ≈ -ustar^3 / (κ * B_stable)

        # ustar = 0 => L = 0
        @test SF.obukhov_length(param_set, FT(0), B_unstable) == FT(0)

        # Near-neutral: very small B triggers non_zero guard
        # The result will be a very large L (not exactly 0)
        B_tiny = eps(FT) / 10
        L_tiny = SF.obukhov_length(param_set, ustar, B_tiny)
        @test isfinite(L_tiny)  # Should be finite, not Inf
    end

    @testset "obukhov_stability_parameter" begin
        ustar = FT(0.3)
        Δz = FT(10.0)

        # Unstable: ζ < 0
        B_unstable = FT(0.01)
        ζ_unstable = SF.obukhov_stability_parameter(param_set, Δz, ustar, B_unstable)
        @test ζ_unstable < 0

        # Should equal Δz / L_MO
        L_MO = SF.obukhov_length(param_set, ustar, B_unstable)
        @test ζ_unstable ≈ Δz / L_MO

        # Stable: ζ > 0
        B_stable = FT(-0.005)
        ζ_stable = SF.obukhov_stability_parameter(param_set, Δz, ustar, B_stable)
        @test ζ_stable > 0
    end

    @testset "Consistency with buoyancy_flux" begin
        # Round-trip: B -> L -> ζ -> B_calc
        ustar = FT(0.5)
        Δz = FT(15.0)
        B_true = FT(0.02)

        L = SF.obukhov_length(param_set, ustar, B_true)
        ζ = SF.obukhov_stability_parameter(param_set, Δz, ustar, B_true)

        # Create minimal inputs for buoyancy_flux(param_set, ζ, ustar, inputs)
        roughness = SF.ConstantRoughnessParams(FT(1e-4), FT(1e-4))
        gustiness = SF.ConstantGustinessSpec(FT(1))
        config = SF.SurfaceFluxConfig(roughness, gustiness)

        inputs = SF.build_surface_flux_inputs(
            param_set,
            FT(300), FT(0.01), FT(0), FT(0), FT(1.2),
            FT(300), FT(0.01), FT(0), Δz, FT(0),
            (FT(10), FT(0)), (FT(0), FT(0)),
            config,
            nothing, SF.FluxSpecs(FT), nothing, nothing,
        )

        B_calc = SF.buoyancy_flux(param_set, ζ, ustar, inputs)
        @test isapprox(B_calc, B_true; rtol = 1e-10)
    end
end

end # module
