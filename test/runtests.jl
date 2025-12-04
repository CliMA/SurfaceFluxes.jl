# Central test runner for SurfaceFluxes.jl:
#   - sets up common test configuration (random seed, array type, float type),
#   - defines small utilities used by the tests,
#   - and groups the individual test files into named `@testset`s.

using Test

import Random
Random.seed!(1234)  # Ensure deterministic tests

import Thermodynamics as TD
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import ClimaParams as CP  # Load CreateParametersExt so `SurfaceFluxesParameters` accepts types

const ArrayType = Array
const FloatType = Float32  # Default scalar type used in several tests

include("test_utils.jl")  # common testing utilities (e.g. @test_allocs_and_ts)

@info "CPU Tests"
@info ArrayType

include("test_finite_difference_vs_finite_volume.jl")  # Finite-difference vs FV velocity-scale checks
include("test_floating_point_consistency.jl")         # Float32 vs Float64 and near-zero L_MO consistency
include("test_surface_fluxes_api.jl")                 # Primitive API regressions

@testset "Regression tests" begin
    # Regression tests with predefined (mostly stable) test cases.
    include("test_regressions.jl")
end

@testset "Test universal functions" begin
    # Unit tests for the asymptotic behavior and type stability of the universal 
    # functions used in SurfaceFluxes.
    include("test_universal_functions.jl")
end

@testset "Convergence tests" begin
    # Convergence tests for a broad variety of test cases.
    include("test_convergence.jl")
end

@testset "Quality assurance" begin
    # Meta-tests (Aqua, etc.) that check for common Julia package issues.
    include("aqua.jl")
end
