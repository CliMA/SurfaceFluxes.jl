# Central test runner for SurfaceFluxes.jl (Subset for UniversalFunctions)
using Test
import Random
Random.seed!(1234)

import Thermodynamics as TD
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import ClimaParams as CP

const ArrayType = Array
const FloatType = Float32 

include("test_utils.jl")

@testset "Test universal functions" begin
    include("test_universal_functions.jl")
end
