if !("." in LOAD_PATH) # for ease of local testing
    push!(LOAD_PATH, ".")
end

import Random
Random.seed!(1234)
using Test
import Thermodynamics
using SurfaceFluxes
const SF = SurfaceFluxes
const SFP = SF.Parameters
using StaticArrays
import KernelAbstractions: CPU

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
const FT = Float32;
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
param_set = create_parameters(toml_dict, UF.BusingerType())
thermo_params = SFP.thermodynamics_params(param_set)
uft = SFP.universal_func_type(param_set)

const TD = Thermodynamics

if get(ARGS, 1, "Array") == "CuArray"
    using CUDA
    import CUDAKernels: CUDADevice
    ArrayType = CUDA.CuArray
    CUDA.allowscalar(false)
    device(::T) where {T <: CuArray} = CUDADevice()
else
    ArrayType = Array
    device(::T) where {T <: Array} = CPU()
end

@show ArrayType

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

    κ = SFP.von_karman_const(param_set)
    for ii in 1:length(b_star)
        # Compute L_MO given u_star and b_star
        L_MO = u_star[ii]^2 / κ / b_star[ii]

        ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc[ii], qt_sfc)
        ts_in = TD.PhaseEquil_ρθq(thermo_params, ρ_in, θ[ii], qt_in)

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
        for jj in 1:length(sc)
            u_scale_fd =
                SF.compute_physical_scale_coeff(param_set, sc[jj], L_MO, UF.MomentumTransport(), uft, SF.FDScheme())
            Δu_fd = u_star[ii] / u_scale_fd
            u_scale_fv =
                SF.compute_physical_scale_coeff(param_set, sc[jj], L_MO, UF.MomentumTransport(), uft, SF.FVScheme())
            Δu_fv = u_star[ii] / u_scale_fv
            @test (Δu_fd - Δu_fv) ./ Δu_fd * 100 <= FT(50)
        end
    end
end

# Parameter set generated in ClimaAtmos GCM run
const sf_params = SurfaceFluxes.Parameters.SurfaceFluxesParameters{
    FT,
    SurfaceFluxes.UniversalFunctions.BusingerParams{FT},
    Thermodynamics.Parameters.ThermodynamicsParameters{FT},
}(
    0.4f0,
    SurfaceFluxes.UniversalFunctions.BusingerParams{FT}(0.74f0, 4.7f0, 4.7f0),
    Thermodynamics.Parameters.ThermodynamicsParameters{FT}(
        273.16f0,
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
    z_levels = [0.01, 0.1, 1, 5, 10, 20, 40, 80, 160, 320, 640] # [m] level of first interior grid point
    for (i, FT) in enumerate(FloatTypes)
        for (jj, z_int) in enumerate(z_levels)
            ts_int_test = Thermodynamics.PhaseEquil{FT}(1.1751807f0, 97086.64f0, 10541.609f0, 0.0f0, 287.85202f0)
            ts_sfc_test =
                Thermodynamics.PhaseEquil{FT}(1.2176297f0, 102852.51f0, 45087.812f0, 0.013232904f0, 291.96683f0)
            sc = SF.ValuesOnly{FT}(;
                state_in = SF.InteriorValues(FT(z_int), (FT(0), FT(0)), ts_int_test),
                state_sfc = SF.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test),
                z0m = FT(1e-5),
                z0b = FT(1e-5),
            )

            sfc_output = SF.surface_conditions(sf_params, sc; maxiter = 20)
            @test abs(SF.obukhov_length(sfc_output)) > FT(0)
            @test sign(SF.non_zero(1.0)) == 1
            @test sign(SF.non_zero(-1.0)) == -1
            @test sign(SF.non_zero(-0.0)) == 1
            @test sign(SF.non_zero(0.0)) == 1
        end
    end
end

@testset "Identical thermodynamic states (Floating Point Consistency)" begin
    FloatTypes = (Float32, Float64)
    z_levels = [0.1, 1, 5, 10, 20, 40, 80, 160, 320, 640] # [m] level of first interior grid point
    z0_m = [1e-5, 1e-3] # roughness length [momentum]
    z0_b = [1e-5, 1e-3] # roughness length [heat] 
    sol_mat = zeros(2, length(z_levels), length(z0_m), length(z0_b))
    for (ii, FT) in enumerate(FloatTypes)
        for (jj, z_int) in enumerate(z_levels)
            for (kk, z0m) in enumerate(z0_m)
                for (ll, z0b) in enumerate(z0_b)
                    # Test case with identical interior and surface states
                    ts_int_test =
                        Thermodynamics.PhaseEquil{FT}(1.1751807f0, 97086.64f0, 10541.609f0, 0.0f0, 287.85202f0)
                    ts_sfc_test = ts_int_test
                    sc = SF.ValuesOnly{FT}(;
                        state_in = SF.InteriorValues(FT(z_int), (FT(0), FT(0)), ts_int_test),
                        state_sfc = SF.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test),
                        z0m = FT(z0m),
                        z0b = FT(z0b),
                    )
                    sfc_output = SF.surface_conditions(sf_params, sc; maxiter = 20)
                    sol_mat[ii, jj, kk, ll] = isinf(sfc_output.L_MO) ? FT(1e6) : sfc_output.L_MO
                end
            end
        end
    end
    rdiff_sol = (sol_mat[1, :, :, :] .- sol_mat[2, :, :, :]) ./ sol_mat[1, :, :, :]
    @test all(x -> x <= sqrt(eps(FT)), rdiff_sol)
end

@testset "Test profiles" begin
    include("test_profiles.jl")
end
@testset "Test universal functions" begin
    include("test_universal_functions.jl")
end
@testset "Test generated thermodynamic states" begin
    include("test_convergence.jl")
end
