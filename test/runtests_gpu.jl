using Test
using CUDA

import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import ClimaParams as CP

if !CUDA.functional()
    @info "CUDA driver not available; skipping GPU tests"
else

    # GPU tests ensure that the GPU path produces the same physical outputs as the 
    # CPU reference implementation for the same states.

    const ArrayType = CUDA.CuArray

    @info "GPU Compatibility Tests"
    @info ArrayType

    # Small synthetic dataset for the broadcasted test. 
    # Temperature values are pre-computed from potential temperature using standard atmosphere.
    # (θ values for reference: 268.559, 269.799, 277.443, 295.792)
    # Using simple approximation: T ≈ θ for near-surface conditions
    const RAW_GPU_DATA = (
        z = (29.432779269303, 30.0497139076724, 31.6880000418153, 34.1873479240475),
        T_int = (268.559, 269.799, 277.443, 295.792),
        T_sfc = (273.424, 272.551, 278.638, 298.133),
        q_tot_int = (0.009, 0.009, 0.009, 0.009),
        q_sfc = (0.01, 0.01, 0.01, 0.01),
        ρ_int = (1.13, 1.13, 1.13, 1.13),
        z0 = (5.86144925739178e-05, 0.0001, 0.000641655193293549, 3.23383768877187e-05),
        speed = (2.9693638452068, 2.43308757772094, 5.69418282305367, 9.5608693754561),
        u_star = (0.109462510724615, 0.0932942802513508, 0.223232887323184, 0.290918439028557),
        b_star = (0.00690834676781433, 0.00428178089592372, 0.00121229800895103, 0.00262353784027441),
    )

    function problem_data(::Type{FT}) where {FT}
        return (; (name => FT.(collect(values)) for (name, values) in pairs(RAW_GPU_DATA))...)
    end

    # Compute a CPU reference solution for the supplied dataset so that GPU results
    # can be compared against identical physics.
    function cpu_reference_L_MO(FT, param_set, data)
        ρ_int = FT(1.13)
        z0h = FT(0.001)
        reference = Vector{FT}(undef, length(data.z))

        for ii in eachindex(data.z)
            T_int = data.T_int[ii]
            T_sfc_guess = data.T_sfc[ii]
            q_tot_int = data.q_tot_int[ii]
            q_vap_sfc_guess = data.q_sfc[ii]

            # Using data.z0[ii] for roughness config
            local_config =
                SF.SurfaceFluxConfig(SF.roughness_lengths(data.z0[ii], z0h), SF.ConstantGustinessSpec(FT(1.0)))

            reference[ii] =
                SF.surface_fluxes(
                    param_set,
                    T_int, q_tot_int, ρ_int,
                    T_sfc_guess, q_vap_sfc_guess,
                    FT(0), data.z[ii], zero(FT),
                    (data.speed[ii], FT(0)), (FT(0), FT(0)),
                    nothing,
                    local_config,
                ).L_MO
        end
        return reference
    end

    # Exercise the GPU broadcast path and compare with CPU.
    function run_gpu_broadcast_test(FT)
        param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
        data = problem_data(FT)
        reference = cpu_reference_L_MO(FT, param_set, data)

        z = ArrayType(data.z)
        T_int = ArrayType(data.T_int)
        T_sfc = ArrayType(data.T_sfc)
        q_tot_int = ArrayType(data.q_tot_int)
        q_sfc = ArrayType(data.q_sfc)
        z0 = ArrayType(data.z0)
        speed = ArrayType(data.speed)
        z0h = FT(0.001)
        ρ_int = FT(1.13)

        # We need to closure or map correctly.
        gpu_outputs = map(
            z,
            speed,
            T_int,
            T_sfc,
            q_tot_int,
            q_sfc,
            z0,
        ) do z_i, speed_i, T_int_i, T_sfc_i, q_tot_int_i, q_sfc_i, z0_i
            local_config = SF.SurfaceFluxConfig(SF.roughness_lengths(z0_i, z0h), SF.ConstantGustinessSpec(FT(1.0)))
            SF.surface_fluxes(
                param_set,
                T_int_i, q_tot_int_i, ρ_int,
                T_sfc_i, q_sfc_i,
                FT(0), z_i, zero(FT),
                (speed_i, FT(0)), (FT(0), FT(0)),
                nothing,
                local_config,
            ).L_MO
        end

        return Array(gpu_outputs), reference
    end

    @testset "GPU broadcast surface_fluxes" begin
        for FT in (Float32, Float64)
            gpu_vals, reference = run_gpu_broadcast_test(FT)
            @test all(isfinite, gpu_vals)
            @test isapprox(gpu_vals, reference; rtol = FT(1e-5), atol = FT(1e-6))
        end
    end

end  # if CUDA.functional()
