# A reference level at or below the momentum roughness length cannot be connected to
# the surface by a Monin-Obukhov profile: the host-side check throws, and the solve,
# which runs in kernels and cannot throw, returns `NaN` with `converged = false`.

module TestReferenceHeight

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

function solve(FT, param_set, Δz, d, z0m; flux_specs = SF.FluxSpecs{FT}())
    config = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(z0m), FT(z0m / 10)),
        SF.ConstantGustinessSpec(FT(1)),
    )
    inputs = SF.build_surface_flux_inputs(
        FT(290), FT(0.005), FT(0), FT(0), FT(1.2), FT(295), FT(0.008), FT(0), FT(Δz),
        FT(d), (FT(3), FT(0)), (FT(0), FT(0)), config, nothing, flux_specs, nothing,
        nothing,
    )
    return SF.surface_fluxes(param_set, inputs), inputs
end

@testset "Reference height validity" begin
    for FT in (Float32, Float64)
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
        @test SF.reference_height_valid((; Δz = FT(10), d = FT(2)), FT(0.1))
        @test !SF.reference_height_valid((; Δz = FT(10), d = FT(9.95)), FT(0.1))
        @test !SF.reference_height_valid((; Δz = FT(10), d = FT(12)), FT(0.1))
        # A scalar roughness length above the momentum one also bounds the level
        @test SF.reference_height_valid((; Δz = FT(10), d = FT(9.95)), FT(0.01), FT(0.04))
        @test !SF.reference_height_valid((; Δz = FT(10), d = FT(9.95)), FT(0.01), FT(0.1))
        @test_throws ArgumentError SF.check_reference_height(
            FT(10),
            FT(9.95),
            FT(0.01),
            FT(0.1),
        )
        @test SF.check_reference_height(FT(10), FT(2), FT(0.1)) === nothing
        @test_throws ArgumentError SF.check_reference_height(FT(10), FT(9.95), FT(0.1))
        @test_throws ArgumentError SF.check_reference_height(FT(10), FT(12), FT(0.1))

        sc, inputs = solve(FT, param_set, 10, 2, 0.1)
        @test sc.converged
        @test isfinite(sc.shf)
        @test SF.reference_height_valid(inputs, FT(0.1))
        # The same point with the reference level below the roughness length
        for d in (9.95, 12)
            bad, _ = solve(FT, param_set, 10, d, 0.1)
            @test !bad.converged
            @test isnan(bad.shf) && isnan(bad.lhf) && isnan(bad.ustar) && isnan(bad.ζ_eff)
            @test isnan(bad.T_sfc) && isnan(bad.L_eff)
        end
        # Every solver mode flags it
        specs = SF.FluxSpecs{FT}(; Cd = FT(1e-3), Ch = FT(1e-3))
        bad, _ = solve(FT, param_set, 10, 12, 0.1; flux_specs = specs)
        @test !bad.converged && isnan(bad.shf)
        good, _ = solve(FT, param_set, 10, 2, 0.1; flux_specs = specs)
        @test good.converged && isfinite(good.shf)
    end
end

end # module
