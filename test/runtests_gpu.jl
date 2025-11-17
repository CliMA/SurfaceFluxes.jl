using Test
using CUDA
using KernelAbstractions

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

    # Small synthetic dataset shared by both the kernel-launch test and the
    # broadcasted test. Values are taken from the CPU regression setup so the CPU
    # reference solution is well defined.
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
        ρ_in = FT(1.13)
        qt_sfc = FT(0.01)
        qt_in = FT(0.009)
        z0b = FT(0.001)
        reference = Vector{FT}(undef, length(data.z))
        for ii in eachindex(data.z)
            ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, data.theta_sfc[ii], qt_sfc)
            ts_in = TD.PhaseEquil_ρθq(thermo_params, ρ_in, data.theta[ii], qt_in)
            state_sfc = SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc)
            state_in = SF.StateValues(data.z[ii], (data.speed[ii], FT(0)), ts_in)
            sc = SF.ValuesOnly(state_in, state_sfc, data.z0[ii], z0b)
            reference[ii] = SF.surface_conditions(param_set, sc).L_MO
        end
        return reference
    end

    # Minimal wrapper kernel: each GPU thread reconstructs the thermodynamic states
    # and runs `surface_conditions`, storing only the Monin–Obukhov length.
    @kernel function surface_conditions_kernel!(
        param_set,
        output,
        z,
        θ,
        θ_sfc,
        z0,
        speed,
        κ,
        ρ_sfc,
        ρ_in,
        qt_sfc,
        qt_in,
    )
        ii = @index(Global)
        thermo_params = SFP.thermodynamics_params(param_set)
        ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc[ii], qt_sfc)
        ts_in = TD.PhaseEquil_ρθq(thermo_params, ρ_in, θ[ii], qt_in)
        state_sfc = SF.StateValues(zero(eltype(z0)), (zero(eltype(z0)), zero(eltype(z0))), ts_sfc)
        state_in = SF.StateValues(z[ii], (speed[ii], zero(eltype(z0))), ts_in)
        sc = SF.ValuesOnly(state_in, state_sfc, z0[ii], eltype(z0)(0.001))
        output[ii] = SF.surface_conditions(param_set, sc).L_MO
    end

    # Launch the KernelAbstractions kernel and compare its output to the CPU
    # reference.
    function run_gpu_kernel_test(FT)
        param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
        data = problem_data(FT)
        reference = cpu_reference_L_MO(FT, param_set, data)

        z = ArrayType(data.z)
        θ = ArrayType(data.theta)
        θ_sfc = ArrayType(data.theta_sfc)
        z0 = ArrayType(data.z0)
        speed = ArrayType(data.speed)

        output = ArrayType{FT}(undef, length(data.z))
        fill!(output, FT(-99999.99))

        backend = KernelAbstractions.get_backend(output)
        kernel! = surface_conditions_kernel!(backend)
        kernel!(
            param_set,
            output,
            z,
            θ,
            θ_sfc,
            z0,
            speed,
            SFP.von_karman_const(param_set),
            FT(1.15),
            FT(1.13),
            FT(0.01),
            FT(0.009);
            ndrange = (length(data.z),),
        )
        KernelAbstractions.synchronize(backend)
        return Array(output), reference
    end

    # Exercise the GPU broadcast path (no explicit kernel) and compare with CPU.
    function run_gpu_broadcast_test(FT)
        param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
        data = problem_data(FT)
        reference = cpu_reference_L_MO(FT, param_set, data)

        z = ArrayType(data.z)
        θ = ArrayType(data.theta)
        θ_sfc = ArrayType(data.theta_sfc)
        z0 = ArrayType(data.z0)
        speed = ArrayType(data.speed)
        z0b = FT(0.001)

        thermo_params = SFP.thermodynamics_params(param_set)
        ρ_sfc = FT(1.15)
        ρ_in = FT(1.13)
        qt_sfc = FT(0.01)
        qt_in = FT(0.009)

        ts_sfc = @. TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc, qt_sfc)
        ts_in = @. TD.PhaseEquil_ρθq(thermo_params, ρ_in, θ, qt_in)
        state_sfc = SF.StateValues.(Ref(FT(0)), Ref((FT(0), FT(0))), ts_sfc)
        state_in = map(z, speed, ts_in) do z_i, speed_i, ts_in_i
            SF.StateValues(z_i, (speed_i, zero(speed_i)), ts_in_i)
        end
        sc = SF.ValuesOnly.(state_in, state_sfc, z0, z0b)
        gpu_outputs = map(sfc -> SF.surface_conditions(param_set, sfc).L_MO, sc)
        return Array(gpu_outputs), reference
    end

    @testset "GPU kernel surface_conditions" begin
        for FT in (Float32, Float64)
            gpu_vals, reference = run_gpu_kernel_test(FT)
            @test all(isfinite, gpu_vals)
            @test isapprox(gpu_vals, reference; rtol = FT(1e-5), atol = FT(1e-6))
        end
    end

    @testset "GPU broadcast surface_conditions" begin
        for FT in (Float32, Float64)
            gpu_vals, reference = run_gpu_broadcast_test(FT)
            @test all(isfinite, gpu_vals)
            @test isapprox(gpu_vals, reference; rtol = FT(1e-5), atol = FT(1e-6))
        end
    end

end  # if CUDA.functional()
