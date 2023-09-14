if !("." in LOAD_PATH) # for ease of local testing
    push!(LOAD_PATH, ".")
end

import Random
Random.seed!(1234)
using Test
import Thermodynamics
using SurfaceFluxes
using CLIMAParameters
const SF = SurfaceFluxes
const SFP = SF.Parameters
using StaticArrays
import KernelAbstractions: CPU

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
FloatType = Float32;
toml_dict = CLIMAParameters.create_toml_dict(FloatType; dict_type = "alias")
param_set = create_parameters(toml_dict, UF.BusingerType())
thermo_params = SFP.thermodynamics_params(param_set)
uft = SFP.universal_func_type(param_set)

const TD = Thermodynamics
device(::T) where {T <: Array} = CPU()

# FIXME: Refactor tests to work on GPUs as in `Thermodynamics.jl`
# if get(ARGS, 1, "Array") == "CuArray"
#    using CUDA
#    import CUDAKernels: CUDADevice
#    ArrayType = CUDA.CuArray
#    device(::T) where {T <: CuArray} = CUDADevice()
#    CUDA.allowscalar(false)
# else
ArrayType = Array
# end

@show ArrayType

# Parameter set generated in ClimaAtmos GCM run
const sf_params = SurfaceFluxes.Parameters.SurfaceFluxesParameters{
    FloatType,
    SurfaceFluxes.UniversalFunctions.BusingerParams{FloatType},
    Thermodynamics.Parameters.ThermodynamicsParameters{FloatType},
}(
    0.4f0,
    SurfaceFluxes.UniversalFunctions.BusingerParams{FloatType}(0.74f0, 4.7f0, 4.7f0, 2.5f0, 4.45f0),
    Thermodynamics.Parameters.ThermodynamicsParameters{FloatType}(
        273.16f0,
        100000.0f0,
        100000.0f0,
        1859.0f0,
        4181.0f0,
        2100.0f0,
        2.5008f6,
        2.8344f6,
        611.657f0,
        273.16f0,
        273.15f0,
        150.0f0,
        1000.0f0,
        298.15f0,
        6864.8f0,
        10513.6f0,
        0.2857143f0,
        8.31446f0,
        0.02897f0,
        0.01801528f0,
        290.0f0,
        220.0f0,
        9.80616f0,
        233.0f0,
        1.0f0,
    ),
)

@testset "Near-zero Obukhov length (Floating Point Consistency)" begin
    FloatTypes = (Float32, Float64)
    z_levels = [1, 5, 10, 20, 40, 80, 160, 320, 640] # [m] level of first interior grid point
    for (i, FloatType) in enumerate(FloatTypes)
        for (jj, z_int) in enumerate(z_levels)
            ts_int_test = Thermodynamics.PhaseEquil{FloatType}(1.1751807f0, 97086.64f0, 10541.609f0, 0.0f0, 287.85202f0)
            ts_sfc_test =
                Thermodynamics.PhaseEquil{FloatType}(1.2176297f0, 102852.51f0, 45087.812f0, 0.013232904f0, 291.96683f0)
            sc = SF.ValuesOnly{FloatType}(;
                state_in = SF.InteriorValues(FloatType(z_int), (FloatType(0), FloatType(0)), ts_int_test),
                state_sfc = SF.SurfaceValues(FloatType(0), (FloatType(0), FloatType(0)), ts_sfc_test),
                z0m = FloatType(1e-5),
                z0b = FloatType(1e-5),
            )

            sfc_output = SF.surface_conditions(sf_params, sc; maxiter = 20)
            @test abs(SF.obukhov_length(sfc_output)) > FloatType(0)
            @test sign(SF.non_zero(1.0)) == 1
            @test sign(SF.non_zero(-1.0)) == -1
            @test sign(SF.non_zero(-0.0)) == 1
            @test sign(SF.non_zero(0.0)) == 1
        end
    end
end

#@info "Identical thermodynamic states"
#@testset "Identical thermodynamic states (Floating Point Consistency)" begin
#    FloatTypes = (Float32, Float64)
#    z_levels = [1, 5, 10, 20, 40, 80, 160, ] # [m] level of first interior grid point
#    z0_m = [1e-5, 1e-4, 1e-3] # roughness length [momentum]
#    z0_b = [1e-5, 1e-4, 1e-3] # roughness length [heat] 
#    sol_mat = Array{Any, 4}(undef, 2, length(z_levels), length(z0_m), length(z0_b))
#    for (ii, FloatType) in enumerate(FloatTypes)
#        for (jj, z_int) in enumerate(z_levels)
#            for (kk, z0m) in enumerate(z0_m)
#                for (ll, z0b) in enumerate(z0_b)
#                    # Test case with identical interior and surface states
#                    ts_int_test =
#                        Thermodynamics.PhaseEquil{FloatType}(1.1751807f0, 97086.64f0, 10541.609f0, 0.0f0, 287.85202f0)
#                    ts_sfc_test = ts_int_test
#                    sc = SF.ValuesOnly{FloatType}(;
#                        state_in = SF.InteriorValues(FloatType(z_int), (FloatType(0), FloatType(0)), ts_int_test),
#                        state_sfc = SF.SurfaceValues(FloatType(0), (FloatType(0), FloatType(0)), ts_sfc_test),
#                        z0m = FloatType(z0m),
#                        z0b = FloatType(z0b),
#                    )
#                    sfc_output = SF.surface_conditions(sf_params, sc; maxiter = 20)
#                    sol_mat[ii, jj, kk, ll] = isinf(sfc_output.L_MO) ? FloatType(1e6) : sfc_output.L_MO
#                end
#            end
#        end
#    end
#    rdiff_sol = (sol_mat[1, :, :, :] .- sol_mat[2, :, :, :]) ./ sol_mat[2, :, :, :]
#    @test all(x -> x <= FloatType(0.005), abs.(rdiff_sol))
#end

@info "Container + Solver Method Comparison"
@testset "Exercise container structs, evaluate and compare Lₘₒ across all available solver methods" begin
    include("test_profiles.jl")
end
@info "Universal Functions"
@testset "Test universal functions" begin
    include("test_universal_functions.jl")
end
@info "Generated Thermodynamic States"
@testset "Test generated thermodynamic states" begin
    include("test_convergence.jl")
end

@testset "Quality assurance" begin
    include("aqua.jl")
end
