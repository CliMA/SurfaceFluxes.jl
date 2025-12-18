# Common testing utilities for SurfaceFluxes.jl

using SurfaceFluxes
using JET
using Test

# Test an expression for allocations and type stability, then return the result.
#
# This macro:
#   - Compiles the expression once
#   - Checks that allocations are ≤ 32 bytes (heap allocations should be minimal)
#   - Uses JET's `@test_opt` to verify type stability
#   - Returns the evaluated expression so it can be used in tests
#
# Example:
#   
#   result = @test_allocs_and_ts compute_something(x, y)
#   @test result ≈ expected
macro test_allocs_and_ts(expression)
    return quote
        local ex() = $(esc(expression))
        ex()  # compile
        @test (@allocated ex()) == 0  # allocations
        @test_opt ex()  # type instabilities
        ex()  # return evaluated expression
    end
end

# Wrapper for `SurfaceFluxes.surface_fluxes` that checks allocations and
# type stability.
#
# Arguments:
#   - `sf_params`: SurfaceFluxes parameters
#   - `sc`: State container (e.g., `SF.ValuesOnly`, `SF.Fluxes`, etc.)
#   - `scheme`: Computation scheme (default: `SurfaceFluxes.PointValueScheme()`)
#
# Returns:
#   The result of `surface_fluxes`, after verifying allocations and type stability.
function surface_fluxes_wrapper(args...)
    @test_allocs_and_ts SurfaceFluxes.surface_fluxes(args...)
end
