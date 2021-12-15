if !("." in LOAD_PATH) # for ease of local testing
    push!(LOAD_PATH, ".")
end

import Random
Random.seed!(1234)
using Test
import Thermodynamics
using SurfaceFluxes
const SF = SurfaceFluxes
const UF = SF.UniversalFunctions
using StaticArrays
import KernelAbstractions: CPU



import CLIMAParameters
const CP = CLIMAParameters
const CPP = CP.Planet
const APS = CP.AbstractEarthParameterSet
const CPSGS = CP.SubgridScale
const CPSGS = CP.SubgridScale
const TD = Thermodynamics

# TODO: add GPU tests back in, similar to Thermodynamics
# if get(ARGS, 1, "Array") == "CuArray"
#     using CUDA
#     import CUDAKernels: CUDADevice
#     ArrayType = CUDA.CuArray
#     CUDA.allowscalar(false)
#     device(::T) where {T <: CuArray} = CUDADevice()
# else
ArrayType = Array
device(::T) where {T <: Array} = CPU()
# end

@show ArrayType

using CLIMAParameters: AbstractEarthParameterSet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

@testset "SurfaceFluxes - Recovery Profiles" begin
    FT = Float32
    ρ_sfc = FT(1.15)
    ρ_in = FT(1.13)
    qt_sfc = FT(0.01)
    qt_in = FT(0.009)
    ## Discretisation altitude z
    z = ArrayType(FT[29.432779269303, 30.0497139076724, 31.6880000418153, 34.1873479240475])
    ## Virtual Potential Temperature at height z
    θ = ArrayType(FT[268.559120403867, 269.799228886728, 277.443023238556, 295.79192777341])
    ## Surface Pottemp
    θ_sfc = ArrayType(FT[273.42369841804, 272.551410044203, 278.638168565727, 298.133068766049])
    ## Roughness lengths
    z0 = ArrayType(FT[5.86144925739178e-05, 0.0001, 0.000641655193293549, 3.23383768877187e-05])
    ## Speed
    speed = ArrayType(FT[2.9693638452068, 2.43308757772094, 5.69418282305367, 9.5608693754561])
    ## Scale velocity and moisture
    u_star = ArrayType(FT[0.109462510724615, 0.0932942802513508, 0.223232887323184, 0.290918439028557])
    # No explicit buoyancy terms in ClimateMachine
    b_star = ArrayType([0.00690834676781433, 0.00428178089592372, 0.00121229800895103, 0.00262353784027441])

    κ = CPSGS.von_karman_const(param_set)
    for ii in 1:length(b_star)
        # Compute L_MO given u_star and b_star
        L_MO = u_star[ii]^2 / κ / b_star[ii]

        ts_sfc = TD.PhaseEquil_ρθq(param_set, ρ_sfc, θ_sfc[ii], qt_sfc)
        ts_in = TD.PhaseEquil_ρθq(param_set, ρ_in, θ[ii], qt_in)

        state_sfc = SF.SurfaceValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
        state_in = SF.InteriorValues(z[ii], SVector{2, FT}(speed[ii], 0), ts_in)

        # State containers
        kwargs = (; state_in, state_sfc, z0m = z0[ii], z0b = FT(0.001))
        sc = (
            SF.Fluxes{FT}(; kwargs..., shf = FT(0), lhf = FT(0)),
            SF.FluxesAndFrictionVelocity{FT}(; kwargs..., shf = FT(0), lhf = FT(0), ustar = u_star[ii]),
            SF.Coefficients{FT}(; kwargs..., Cd = FT(0.001), Ch = FT(0.001)),
            SF.ValuesOnly{FT}(; kwargs...),
        )
        uf = UF.Businger()
        for jj in 1:length(sc)
            u_scale_fd =
                SF.compute_physical_scale_coeff(param_set, sc[jj], L_MO, UF.MomentumTransport(), uf, SF.FDScheme())
            Δu_fd = u_star[ii] / u_scale_fd
            u_scale_fv =
                SF.compute_physical_scale_coeff(param_set, sc[jj], L_MO, UF.MomentumTransport(), uf, SF.FVScheme())
            Δu_fv = u_star[ii] / u_scale_fv
            @test (Δu_fd - Δu_fv) ./ Δu_fd * 100 <= FT(50)
        end
    end
end

include("test_profiles.jl")
include("test_universal_functions.jl")
