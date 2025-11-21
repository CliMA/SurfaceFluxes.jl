# Quality assurance tests using Aqua.jl
#
# Aqua.jl provides automated quality assurance checks for Julia packages.
# This test file ensures that SurfaceFluxes.jl maintains good code quality
# by checking for common issues like unbound type parameters, method ambiguities,
# and other package health indicators.
#
# For more information about Aqua.jl, see:
# https://github.com/JuliaTesting/Aqua.jl

using Test
using SurfaceFluxes
using Aqua

@testset "Aqua tests (performance)" begin
    # Test for unbound type parameters.
    #
    # Unbound type parameters can cause performance issues by preventing type
    # inference and forcing runtime dispatch. This checks that we don't have
    # any unbound arguments that would trigger Julia issue #29393.
    # See: https://github.com/JuliaLang/julia/issues/29393
    unbound_args = Aqua.detect_unbound_args_recursively(SurfaceFluxes)
    @test length(unbound_args) == 0

    # Test for method ambiguities within SurfaceFluxes.
    #
    # Method ambiguities can cause runtime errors or performance degradation.
    # This test checks for ambiguities only within SurfaceFluxes (not across
    # dependencies) to ensure our API is unambiguous.
    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    all_ambiguities = Aqua.detect_ambiguities(SurfaceFluxes; recursive = true)

    # Filter to only include ambiguities where at least one method is from SurfaceFluxes.
    # This excludes ambiguities between dependencies (which are not our responsibility).
    pkg_dir_matches(pkgname::String, pkdir::Nothing) = false
    pkg_dir_matches(pkgname::String, pkdir::AbstractString) = occursin(pkgname, pkdir)
    surfacefluxes_ambiguities = filter(
        amb -> pkg_dir_matches("SurfaceFluxes", pkgdir(last(amb).module)),
        all_ambiguities,
    )

    # Require zero ambiguities within SurfaceFluxes.
    # If ambiguities are found, uncomment the debug block below to see details.
    @test length(surfacefluxes_ambiguities) == 0

    # Debug output (uncomment to see details of ambiguities):
    # if length(surfacefluxes_ambiguities) > 0
    #     @info "Method ambiguities found:" surfacefluxes_ambiguities
    #     for amb in surfacefluxes_ambiguities
    #         @show amb
    #     end
    # end
end

@testset "Aqua tests - remaining" begin
    # Run all other Aqua quality checks (project dependencies, version bounds, etc.).
    # We exclude ambiguities and unbound_args since those are tested above with
    # custom filtering/checks.
    Aqua.test_all(SurfaceFluxes; ambiguities = false, unbound_args = false)
end
