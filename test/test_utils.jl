using SurfaceFluxes
using JET

# macro based on test_all from ClimaCore.jl
# tests an expression for allocations and type stability,
# and then returns the evaluated expression
macro test_allocs_and_ts(expression)
    return quote
        local ex() = $(esc(expression))
        ex() # compile
        @test (@allocated ex()) <= 32 # allocations
        @test_opt ex() # type instabilities
        ex() # return evaluated expression
    end
end

# wrapper for surface_conditions that checks allocations and type stability
function surface_conditions_wrapper(sf_params, sc, scheme = SurfaceFluxes.PointValueScheme(); kwargs...)
    @test_allocs_and_ts SurfaceFluxes.surface_conditions(sf_params, sc; kwargs...)
end
