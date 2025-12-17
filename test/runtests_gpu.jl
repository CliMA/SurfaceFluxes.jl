using Test
using CUDA

import Thermodynamics as TD
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

    # Small synthetic dataset for the broadcasted test. Values are taken from the
    # CPU regression setup so the CPU reference solution is well defined.
    const RAW_GPU_DATA = (
        z = (29.432779269303, 30.0497139076724, 31.6880000418153, 34.1873479240475),
        theta = (268.559120403867, 269.799228886728, 277.443023238556, 295.79192777341),
        theta_sfc = (273.42369841804, 272.551410044203, 278.638168565727, 298.133068766049),
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
        thermo_params = SFP.thermodynamics_params(param_set)
        ρ_sfc = FT(1.15)
        ρ_int = FT(1.13)
        qt_sfc = FT(0.01)
        qt_int = FT(0.009)
        z0h = FT(0.001)
        reference = Vector{FT}(undef, length(data.z))

        config = SF.SurfaceFluxConfig(SF.roughness_lengths(FT(0.001), z0h), SF.ConstantGustinessSpec(FT(1.0)))

        for ii in eachindex(data.z)
            ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, data.theta_sfc[ii], qt_sfc)
            ts_int = TD.PhaseEquil_ρθq(thermo_params, ρ_int, data.theta[ii], qt_int)
            T_int = TD.air_temperature(thermo_params, ts_int)
            T_sfc_guess = TD.air_temperature(thermo_params, ts_sfc)
            q_tot_int = TD.total_specific_humidity(thermo_params, ts_int)
            q_vap_sfc_guess = TD.total_specific_humidity(thermo_params, ts_sfc)

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
        θ = ArrayType(data.theta)
        θ_sfc = ArrayType(data.theta_sfc)
        z0 = ArrayType(data.z0)
        speed = ArrayType(data.speed)
        z0h = FT(0.001)

        thermo_params = SFP.thermodynamics_params(param_set)
        ρ_sfc = FT(1.15)
        ρ_int = FT(1.13)
        qt_sfc = FT(0.01)
        qt_int = FT(0.009)

        ts_sfc = @. TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc, qt_sfc)
        ts_int = @. TD.PhaseEquil_ρθq(thermo_params, ρ_int, θ, qt_int)

        T_int = @. TD.air_temperature(thermo_params, ts_int)
        T_sfc_guess = @. TD.air_temperature(thermo_params, ts_sfc)
        q_tot_int = @. TD.total_specific_humidity(thermo_params, ts_int)
        q_vap_sfc_guess = @. TD.total_specific_humidity(thermo_params, ts_sfc)

        # We need to construct roughness models for broadcast
        config = SF.SurfaceFluxConfig(SF.roughness_lengths(FT(0.001), z0h), SF.ConstantGustinessSpec(FT(1.0)))

        # Currently surface_fluxes on GPU with varying roughness via config structs inside array inputs?
        # The main `surface_fluxes` doesn't broadcast well if inputs are arrays but function call is singular.
        # But we are mapping.

        # We need to closure or map correctly.
        gpu_outputs = map(
            z,
            speed,
            T_int,
            T_sfc_guess,
            q_tot_int,
            q_vap_sfc_guess,
            z0,
        ) do z_i, speed_i, T_int_i, T_sfc_guess_i, q_tot_int_i, q_vap_sfc_guess_i, z0_i
            local_config = SF.SurfaceFluxConfig(SF.roughness_lengths(z0_i, z0h), SF.ConstantGustinessSpec(FT(1.0)))
            SF.surface_fluxes(
                param_set,
                T_int_i, q_tot_int_i, ρ_int,
                T_sfc_guess_i, q_vap_sfc_guess_i,
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
