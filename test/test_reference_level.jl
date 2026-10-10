# The surface state applies at the displacement height, so the fluxes depend on the
# geometry only through the effective height `Δz - d`, and a reference height measured
# from the apparent sink (`ReferenceAboveApparentSink`) is converted to one measured
# from the surface before the solve.

module TestReferenceLevel

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
using Main: @test_allocs_and_ts

const z0m = 0.1
const z0h = 0.01

function inputs_for(FT, Δz, d; T_sfc = 295,
    roughness = SF.ConstantRoughnessParams(FT(z0m), FT(z0h)),
    roughness_inputs = nothing, reference_level = SF.ReferenceAboveSurface())
    config = SF.SurfaceFluxConfig(
        roughness, SF.ConstantGustinessSpec(FT(1)), SF.MoistModel(),
        SF.NoRoughnessSubLayer(), SF.NoStabilityCap(), reference_level,
    )
    return SF.build_surface_flux_inputs(
        FT(290), FT(0.005), FT(0), FT(0), FT(1.2), FT(T_sfc), FT(0.008), FT(0), FT(Δz),
        FT(d), (FT(3), FT(0)), (FT(0), FT(0)), config, roughness_inputs,
        SF.FluxSpecs{FT}(), nothing, nothing,
    )
end

# Allocation and type-stability checks need the arguments passed in, since locals of
# a testset body captured by the check's closure are boxed
convert_checked(ps, inputs) = @test_allocs_and_ts SF.reference_above_surface(ps, inputs)

same_conditions(a, b; rtol, except = (:converged,)) = all(
    isapprox(getproperty(a, f), getproperty(b, f); rtol) for
    f in fieldnames(typeof(a)) if !(f in except)
)

@testset "Reference level" begin
    for FT in (Float32, Float64)
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
        g = SFP.grav(param_set)
        rtol = FT == Float32 ? 1e-5 : 1e-10

        @testset "Surface state at the displacement height ($FT)" begin
            inputs = inputs_for(FT, 20, 5)
            @test SF.surface_geopotential(param_set, inputs) == g * FT(5)
            # The deprecated one-argument form is the geopotential of the ground
            @test SF.surface_geopotential(inputs) == inputs.Φ_sfc
            @test SF.interior_geopotential(param_set, inputs) == g * FT(20)
            @test SF.surface_density(param_set, inputs, FT(295), FT(0.008)) ==
                  SF.surface_density(
                param_set, FT(290), FT(1.2), FT(295), FT(15), FT(0.005), FT(0), FT(0),
                FT(0.008),
            )
            # Raising the displacement height and the reference level together leaves
            # the exchange unchanged for every stability; the sensible heat flux gains
            # the potential energy of the vapor, which leaves the surface state higher
            for T_sfc in (288, 290, 300), δ in (3, 15)
                sc = SF.surface_fluxes(param_set, inputs_for(FT, 20, 2; T_sfc))
                lifted = SF.surface_fluxes(param_set, inputs_for(FT, 20 + δ, 2 + δ; T_sfc))
                @test sc.converged && lifted.converged
                @test same_conditions(sc, lifted; rtol, except = (:converged, :shf))
                atol = 200 * eps(FT) * max(abs(sc.shf), one(FT))
                @test lifted.shf - sc.shf ≈ g * δ * sc.evaporation rtol = 1e-3 atol = atol
            end
        end

        @testset "Reference height above the apparent sink ($FT)" begin
            above_sink =
                inputs_for(FT, 10, 15; reference_level = SF.ReferenceAboveApparentSink())
            above_surface = inputs_for(FT, 10 + 15 + z0m, 15)
            converted = convert_checked(param_set, above_sink)
            @test converted.Δz ≈ above_surface.Δz
            @test converted.reference_level === SF.ReferenceAboveSurface()
            @test SF.reference_above_surface(param_set, above_surface) === above_surface
            # Inputs without the field follow the default convention
            bare = Base.structdiff(above_surface, NamedTuple{(:reference_level,)})
            @test SF.reference_above_surface(param_set, bare) === bare
            @test same_conditions(
                SF.surface_fluxes(param_set, bare),
                SF.surface_fluxes(param_set, above_surface);
                rtol,
            )
            # A tall canopy below a low forcing height is valid under this convention
            sc = SF.surface_fluxes(param_set, above_sink)
            @test sc.converged
            @test same_conditions(sc, SF.surface_fluxes(param_set, above_surface); rtol)
            @test !SF.surface_fluxes(param_set, inputs_for(FT, 10, 15)).converged
            # Screen-level values apply the same conversion
            s = SF.screen_level_values(param_set, sc, above_sink, FT(2), FT(10))
            s_ref = SF.screen_level_values(param_set, sc, above_surface, FT(2), FT(10))
            @test s.T ≈ s_ref.T && s.q ≈ s_ref.q && s.u ≈ s_ref.u
            # Direct helper functions also convert inputs under ReferenceAboveApparentSink
            @test SF.effective_height(param_set, above_sink) ≈
                  SF.effective_height(above_surface)
            @test SF.interior_geopotential(param_set, above_sink) ≈
                  SF.interior_geopotential(param_set, above_surface)
            @test SF.surface_density(param_set, above_sink, FT(295), FT(0.008)) ≈
                  SF.surface_density(param_set, above_surface, FT(295), FT(0.008))
            sch = SF.PointValueScheme()
            @test SF.heat_conductance(
                param_set, sc.ζ, sc.ustar, above_sink, FT(z0m), FT(z0h), sch,
            ) ≈ SF.heat_conductance(
                param_set, sc.ζ, sc.ustar, above_surface, FT(z0m), FT(z0h), sch,
            )
            @test SF.compute_ustar(
                param_set, sc.ζ, FT(z0m), above_sink, sch, FT(1),
            ) ≈ SF.compute_ustar(
                param_set, sc.ζ, FT(z0m), above_surface, sch, FT(1),
            )
            @test SF.compute_theta_star(
                param_set, sc.ζ, FT(z0h), above_sink, sch, FT(295),
            ) ≈ SF.compute_theta_star(
                param_set, sc.ζ, FT(z0h), above_surface, sch, FT(295),
            )
            @test SF.compute_q_star(
                param_set, sc.ζ, FT(z0h), above_sink, sch, FT(0.008),
            ) ≈ SF.compute_q_star(
                param_set, sc.ζ, FT(z0h), above_surface, sch, FT(0.008),
            )
            @test SF.buoyancy_flux(param_set, sc.ζ, sc.ustar, above_sink) ≈
                  SF.buoyancy_flux(param_set, sc.ζ, sc.ustar, above_surface)
            # Canopy roughness from the Raupach model
            raupach = SF.RaupachRoughnessParams{FT}()
            canopy = (PAI = FT(3), h = FT(20))
            d = SF.displacement_height(raupach, canopy)
            z0m_canopy = SF.momentum_roughness(raupach, FT(0), param_set, canopy)
            sink = inputs_for(
                FT, 10, d; roughness = raupach, roughness_inputs = canopy,
                reference_level = SF.ReferenceAboveApparentSink(),
            )
            surface = inputs_for(
                FT, 10 + d + z0m_canopy, d; roughness = raupach,
                roughness_inputs = canopy,
            )
            @test SF.reference_above_surface(param_set, sink).Δz ≈ surface.Δz
            @test same_conditions(
                SF.surface_fluxes(param_set, sink), SF.surface_fluxes(param_set, surface);
                rtol,
            )
            # A roughness length that depends on the friction velocity is rejected
            coare = inputs_for(
                FT, 10, 0; roughness = SF.COARE3RoughnessParams{FT}(),
                reference_level = SF.ReferenceAboveApparentSink(),
            )
            @test_throws ArgumentError SF.reference_above_surface(param_set, coare)
        end

        @testset "Float64 roughness lengths keep the solve in $FT" begin
            rough64 = SF.ConstantRoughnessParams(z0m, z0h)
            sc = SF.surface_fluxes(param_set, inputs_for(FT, 20, 5; roughness = rough64))
            @test sc.shf isa FT && sc.ustar isa FT
            sink = inputs_for(
                FT, 10, 15; roughness = rough64,
                reference_level = SF.ReferenceAboveApparentSink(),
            )
            @test SF.reference_above_surface(param_set, sink).Δz isa FT
        end
    end
end

end # module
