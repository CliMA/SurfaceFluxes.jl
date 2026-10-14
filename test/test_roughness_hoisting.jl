# For roughness models independent of the friction velocity, `surface_fluxes` evaluates
# the roughness lengths once per solve and stores them in the inputs
# (`with_roughness_lengths`). The tests count the evaluations and check that the fluxes
# and their derivatives equal those of a solve that evaluates the roughness model
# wherever it needs the roughness lengths.

module TestRoughnessHoisting

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
import ForwardDiff
using Main: @test_allocs_and_ts

# Constant roughness lengths that count their evaluations; the type parameter `D` sets
# `depends_on_ustar` at compile time, as for the package's models, and so decides whether
# the solver evaluates them once
struct CountingRoughness{D, FT} <: SF.AbstractRoughnessParams
    z0m::FT
    z0s::FT
end
CountingRoughness(z0m::FT, z0s::FT, ustar_dependent::Bool) where {FT} =
    CountingRoughness{ustar_dependent, FT}(z0m, z0s)
const EVALUATIONS = Ref(0)
SF.depends_on_ustar(::CountingRoughness{D}) where {D} = D
function SF.momentum_and_scalar_roughness(spec::CountingRoughness, u★, ps, roughness_inputs)
    EVALUATIONS[] += 1
    return (spec.z0m, spec.z0s)
end
SF.momentum_roughness(spec::CountingRoughness, u★, ps, roughness_inputs) =
    first(SF.momentum_and_scalar_roughness(spec, u★, ps, roughness_inputs))

function inputs_for(
    FT,
    roughness;
    T_sfc = 288,
    Δz = 20,
    d = 5,
    roughness_inputs = (PAI = FT(3), h = FT(8)),
    gustiness = SF.ConstantGustinessSpec(FT(1)),
    stability_cap = SF.NoStabilityCap(),
    reference_level = SF.ReferenceAboveSurface(),
    flux_specs = SF.FluxSpecs{FT}(),
    update_T_sfc = nothing,
)
    config = SF.SurfaceFluxConfig(
        roughness, gustiness, SF.MoistModel(), SF.NoRoughnessSubLayer(), stability_cap,
        reference_level,
    )
    return SF.build_surface_flux_inputs(
        FT(290), FT(0.005), FT(0), FT(0), FT(1.2), FT(T_sfc), FT(0.008), FT(0), FT(Δz),
        FT(d), (FT(3), FT(0)), (FT(0), FT(0)), config, roughness_inputs, flux_specs,
        update_T_sfc, nothing,
    )
end

# The solve without the roughness lengths in the inputs, which evaluates the roughness
# model wherever it needs the roughness lengths
function unhoisted_surface_fluxes(param_set, inputs, scheme)
    inputs = SF.reference_above_surface(param_set, inputs)
    sc = SF.surface_fluxes_by_mode(param_set, inputs, scheme, nothing)
    z0m, z0h = SF.momentum_and_scalar_roughness(
        inputs.roughness_model, sc.ustar, param_set, inputs.roughness_inputs,
    )
    return SF.invalidate_unless(sc, SF.reference_height_valid(inputs, z0m, z0h))
end

same_conditions(a, b) =
    all(isequal(getproperty(a, f), getproperty(b, f)) for f in fieldnames(typeof(a)))

function evaluations(f)
    EVALUATIONS[] = 0
    f()
    return EVALUATIONS[]
end

# Allocation and type-stability checks need the arguments passed in, since locals of a
# testset body captured by the check's closure are boxed
hoist_checked(inputs, ps) = @test_allocs_and_ts SF.with_roughness_lengths(inputs, ps)
solve_checked(ps, inputs) = @test_allocs_and_ts SF.surface_fluxes(ps, inputs)

# Surface temperature relaxed toward the air temperature by the friction velocity
relax_T_sfc(ζ, ps, thermo_params, inputs, scheme, u★, z0m, z0h) =
    (inputs.T_sfc_guess + u★ * inputs.T_int) / (1 + u★)

@testset "Roughness hoisting" begin
    for FT in (Float32, Float64)
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
        raupach = SF.RaupachRoughnessParams{FT}()
        constant = SF.ConstantRoughnessParams{FT}(FT(0.1), FT(0.01))
        coare = SF.COARE3RoughnessParams{FT}()

        @testset "Roughness lengths in the inputs ($FT)" begin
            inputs = inputs_for(FT, raupach)
            @test inputs.z0m === nothing && inputs.z0h === nothing
            hoisted = hoist_checked(inputs, param_set)
            @test (hoisted.z0m, hoisted.z0h) === SF.momentum_and_scalar_roughness(
                raupach, zero(FT), param_set, inputs.roughness_inputs,
            )
            @test hoisted.z0m isa FT && hoisted.z0h isa FT
            @test (@inferred SF.roughness_lengths(param_set, hoisted, FT(0.3))) ===
                  (hoisted.z0m, hoisted.z0h)
            @test SF.roughness_lengths(param_set, inputs, FT(0.3)) ===
                  (hoisted.z0m, hoisted.z0h)
            # Set roughness lengths are kept
            @test SF.with_roughness_lengths(hoisted, param_set) === hoisted
            # Roughness lengths that depend on u★ are evaluated where they are needed
            coare_inputs = inputs_for(FT, coare; roughness_inputs = nothing)
            @test SF.with_roughness_lengths(coare_inputs, param_set) === coare_inputs
            @test SF.roughness_lengths(param_set, coare_inputs, FT(0.3)) ===
                  SF.momentum_and_scalar_roughness(coare, FT(0.3), param_set, nothing)
            # Inputs without the fields evaluate the roughness model
            bare = Base.structdiff(inputs, NamedTuple{(:z0m, :z0h)})
            @test SF.roughness_lengths(param_set, bare, FT(0.3)) ===
                  (hoisted.z0m, hoisted.z0h)
            @test SF.with_roughness_lengths(bare, param_set).z0m === hoisted.z0m
            # The reference height conversion reads the roughness length from the inputs
            sink =
                inputs_for(FT, raupach; reference_level = SF.ReferenceAboveApparentSink())
            @test SF.reference_above_surface(param_set, sink).Δz ===
                  SF.reference_above_surface(
                param_set,
                SF.with_roughness_lengths(sink, param_set),
            ).Δz
            sc = solve_checked(param_set, inputs)
            @test sc.converged
        end

        @testset "One evaluation per solve ($FT)" begin
            counting = CountingRoughness(FT(0.1), FT(0.01), false)
            dependent = CountingRoughness(FT(0.1), FT(0.01), true)
            ζ_cap = SF.MaxHeatFluxStabilityCap()
            sink = SF.ReferenceAboveApparentSink()
            deardorff = SF.DeardorffGustinessSpec()
            prescribed = SF.FluxSpecs{FT}(shf = FT(-20), lhf = FT(10), ustar = FT(0.3))
            heat_and_drag = SF.FluxSpecs{FT}(shf = FT(-20), lhf = FT(10), Cd = FT(0.01))
            for kwargs in (
                    (;),
                    (; T_sfc = 295),
                    (; stability_cap = ζ_cap),
                    (; reference_level = sink, stability_cap = ζ_cap),
                    (; gustiness = deardorff, T_sfc = 295),
                    (; update_T_sfc = relax_T_sfc),
                    (; flux_specs = prescribed, stability_cap = ζ_cap),
                    (; flux_specs = heat_and_drag),
                ),
                scheme in (SF.PointValueScheme(), SF.LayerAverageScheme())

                inputs = inputs_for(FT, counting; kwargs...)
                @test evaluations(() -> SF.surface_fluxes(param_set, inputs, scheme)) == 1
            end
            # Screen-level recovery under the apparent-sink convention, which needs z0m
            # for the reference height and both roughness lengths for the profiles
            inputs = inputs_for(FT, counting; reference_level = sink)
            sc = SF.surface_fluxes(param_set, inputs)
            @test evaluations(
                () -> SF.screen_level_values(param_set, sc, inputs, FT(2), FT(10)),
            ) == 1
            # A model that depends on u★ is evaluated in every residual evaluation
            inputs = inputs_for(FT, dependent)
            @test evaluations(() -> SF.surface_fluxes(param_set, inputs)) > 10
        end

        @testset "Fluxes unchanged ($FT)" begin
            for roughness in (raupach, constant),
                T_sfc in (280, 288, 295, 305),
                kwargs in (
                    (;),
                    (; stability_cap = SF.MaxHeatFluxStabilityCap()),
                    (; reference_level = SF.ReferenceAboveApparentSink()),
                    (; gustiness = SF.DeardorffGustinessSpec()),
                    (; gustiness = SF.FlooredDeardorffGustinessSpec(FT(0.5))),
                    (; update_T_sfc = relax_T_sfc),
                    (;
                        flux_specs = SF.FluxSpecs{FT}(
                            shf = FT(-20), lhf = FT(10), ustar = FT(0.3),
                        ),
                        stability_cap = SF.MaxHeatFluxStabilityCap(),
                    ),
                ),
                scheme in (SF.PointValueScheme(), SF.LayerAverageScheme())

                inputs = inputs_for(FT, roughness; T_sfc, kwargs...)
                @test same_conditions(
                    SF.surface_fluxes(param_set, inputs, scheme),
                    unhoisted_surface_fluxes(param_set, inputs, scheme),
                )
            end
            # Reference level below the stored roughness lengths
            z0m_canopy, _ = SF.momentum_and_scalar_roughness(
                raupach, zero(FT), param_set, (PAI = FT(3), h = FT(8)),
            )
            low = inputs_for(FT, raupach; Δz = 5 + z0m_canopy / 2, d = 5)
            sc = SF.surface_fluxes(param_set, low)
            @test !sc.converged && isnan(sc.shf)
            @test same_conditions(
                sc,
                unhoisted_surface_fluxes(param_set, low, SF.PointValueScheme()),
            )
        end

        @testset "Differentiation with respect to the canopy ($FT)" begin
            # The stored roughness lengths are dual numbers that hold the derivatives with
            # respect to the canopy height
            canopy(h) = inputs_for(FT, raupach; roughness_inputs = (PAI = FT(3), h))
            shf(h) = SF.surface_fluxes(param_set, canopy(h)).shf
            shf_unhoisted(h) =
                unhoisted_surface_fluxes(param_set, canopy(h), SF.PointValueScheme()).shf
            h = FT(8)
            δ = FT == Float32 ? FT(1e-2) : FT(1e-5)
            fd = (shf(h + δ) - shf(h - δ)) / (2δ)
            ad = ForwardDiff.derivative(shf, h)
            @test isfinite(ad) && ad != 0
            @test ad ≈ fd rtol = FT == Float32 ? 1e-2 : 1e-4
            @test ad == ForwardDiff.derivative(shf_unhoisted, h)
        end
    end
end

end # module
