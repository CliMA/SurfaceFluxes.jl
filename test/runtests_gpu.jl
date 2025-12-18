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
    const RAW_GPU_DATA = (
        z = (29.432779269303, 30.0497139076724, 31.6880000418153, 34.1873479240475),
        T_int = (268.559, 269.799, 277.443, 295.792),
        T_sfc = (273.424, 272.551, 278.638, 298.133),
        q_tot_int = (0.009, 0.009, 0.009, 0.009),
        q_sfc = (0.01, 0.01, 0.01, 0.01),
        ρ_int = (1.13, 1.13, 1.13, 1.13),
        z0 = (5.86144925739178e-05, 0.0001, 0.000641655193293549, 3.23383768877187e-05),
        speed = (2.9693638452068, 2.43308757772094, 5.69418282305367, 9.5608693754561),
    )

    function problem_data(::Type{FT}) where {FT}
        return (;
            (name => FT.(collect(values)) for (name, values) in pairs(RAW_GPU_DATA))...
        )
    end

    # Test broadcasting over homogeneous config arrays (same type, different values)
    @testset "GPU broadcast - Homogeneous Configs (Varying Roughness)" begin
        for FT in (Float32, Float64)
            param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
            data = problem_data(FT)
            n = length(data.z)

            # Build array of configs with varying z0 but same structure
            z0h = FT(0.001)
            cpu_configs = [
                SF.SurfaceFluxConfig(
                    SF.roughness_lengths(data.z0[i], z0h),
                    SF.ConstantGustinessSpec(FT(1.0)),
                ) for i in 1:n
            ]

            # CPU reference
            cpu_ρτxz = Vector{FT}(undef, n)
            cpu_ρτyz = Vector{FT}(undef, n)
            cpu_shf = Vector{FT}(undef, n)
            cpu_lhf = Vector{FT}(undef, n)

            for i in 1:n
                result = SF.surface_fluxes(
                    param_set,
                    data.T_int[i], data.q_tot_int[i], data.ρ_int[i],
                    data.T_sfc[i], data.q_sfc[i],
                    FT(0), data.z[i], FT(0),
                    (data.speed[i], FT(0)), (FT(0), FT(0)),
                    nothing,
                    cpu_configs[i],
                )
                cpu_ρτxz[i] = result.ρτxz
                cpu_ρτyz[i] = result.ρτyz
                cpu_shf[i] = result.shf
                cpu_lhf[i] = result.lhf
            end

            # GPU broadcast over config array
            gpu_configs = ArrayType(cpu_configs)
            T_int = ArrayType(data.T_int)
            T_sfc = ArrayType(data.T_sfc)
            q_tot_int = ArrayType(data.q_tot_int)
            q_sfc = ArrayType(data.q_sfc)
            ρ_int_array = ArrayType(fill(FT(1.13), n))
            z = ArrayType(data.z)
            speed = ArrayType(data.speed)
            Φ_sfc_array = ArrayType(fill(FT(0), n))
            d_array = ArrayType(fill(FT(0), n))
            u_sfc_array = ArrayType([(FT(0), FT(0)) for _ in 1:n])

            # Broadcast surface_fluxes over the config array
            # Broadcast surface_fluxes once to get all results
            u_int_cpu = [(data.speed[i], FT(0)) for i in 1:n]
            u_int_array = ArrayType(u_int_cpu)

            gpu_results =
                SF.surface_fluxes.(
                    Ref(param_set), T_int, q_tot_int, ρ_int_array, T_sfc, q_sfc,
                    Φ_sfc_array, z, d_array, u_int_array, u_sfc_array,
                    Ref(nothing), gpu_configs,
                    Ref(SF.PointValueScheme()),
                    Ref(SF.SolverOptions{FT}(FT(1e-2), 15, true)),
                    Ref(SF.FluxSpecs{FT}()),
                )

            # Extract individual fields
            gpu_ρτxz = map(x -> x.ρτxz, gpu_results)
            gpu_ρτyz = map(x -> x.ρτyz, gpu_results)
            gpu_shf = map(x -> x.shf, gpu_results)
            gpu_lhf = map(x -> x.lhf, gpu_results)

            # Verify
            @test all(isfinite, Array(gpu_ρτxz))
            @test all(isfinite, Array(gpu_ρτyz))
            @test all(isfinite, Array(gpu_shf))
            @test all(isfinite, Array(gpu_lhf))
            @test isapprox(Array(gpu_ρτxz), cpu_ρτxz; rtol = FT(1e-5))
            @test isapprox(Array(gpu_ρτyz), cpu_ρτyz; rtol = FT(1e-5))
            @test isapprox(Array(gpu_shf), cpu_shf; rtol = FT(1e-5))
            @test isapprox(Array(gpu_lhf), cpu_lhf; rtol = FT(1e-5))
        end
    end

    # Test broadcasting over configs with varying gustiness values
    @testset "GPU broadcast - Varying Gustiness Values" begin
        for FT in (Float32, Float64)
            param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
            data = problem_data(FT)
            n = length(data.z)

            # Varying gustiness values at different locations
            gustiness_values = FT.([0.5, 1.0, 1.5, 2.0])
            z0m = FT(0.0001)
            z0h = FT(0.001)

            cpu_configs = [
                SF.SurfaceFluxConfig(
                    SF.roughness_lengths(z0m, z0h),
                    SF.ConstantGustinessSpec(gustiness_values[i]),
                ) for i in 1:n
            ]

            # CPU reference
            cpu_shf = Vector{FT}(undef, n)
            cpu_lhf = Vector{FT}(undef, n)

            for i in 1:n
                result = SF.surface_fluxes(
                    param_set,
                    data.T_int[i], data.q_tot_int[i], data.ρ_int[i],
                    data.T_sfc[i], data.q_sfc[i],
                    FT(0), data.z[i], FT(0),
                    (data.speed[i], FT(0)), (FT(0), FT(0)),
                    nothing,
                    cpu_configs[i],
                )
                cpu_shf[i] = result.shf
                cpu_lhf[i] = result.lhf
            end

            # GPU broadcast
            gpu_configs = ArrayType(cpu_configs)
            T_int = ArrayType(data.T_int)
            T_sfc = ArrayType(data.T_sfc)
            q_tot_int = ArrayType(data.q_tot_int)
            q_sfc = ArrayType(data.q_sfc)
            ρ_int_array = ArrayType(fill(data.ρ_int[1], n))
            z = ArrayType(data.z)
            speed = ArrayType(data.speed)
            Φ_sfc_array = ArrayType(fill(FT(0), n))
            d_array = ArrayType(fill(FT(0), n))
            u_sfc_array = ArrayType([(FT(0), FT(0)) for _ in 1:n])

            u_int_cpu = [(data.speed[i], FT(0)) for i in 1:n]
            u_int_array = ArrayType(u_int_cpu)

            gpu_results =
                SF.surface_fluxes.(
                    Ref(param_set), T_int, q_tot_int, ρ_int_array, T_sfc, q_sfc,
                    Φ_sfc_array, z, d_array, u_int_array, u_sfc_array,
                    Ref(nothing), gpu_configs,
                )

            gpu_shf = map(x -> x.shf, gpu_results)
            gpu_lhf = map(x -> x.lhf, gpu_results)

            @test all(isfinite, Array(gpu_shf))
            @test all(isfinite, Array(gpu_lhf))
            @test isapprox(Array(gpu_shf), cpu_shf; rtol = FT(1e-5))
            @test isapprox(Array(gpu_lhf), cpu_lhf; rtol = FT(1e-5))
        end
    end

    # Test broadcasting over heterogeneous configs (different roughness types)
    @testset "GPU broadcast - Heterogeneous Configs (COARE3 + Raupach)" begin
        for FT in (Float32, Float64)
            param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
            data = problem_data(FT)
            n = length(data.z)

            # Mix COARE3 (ocean-like) and Raupach (land-like) roughness
            # Use COARE3 for odd indices, Raupach for even indices
            # Define concrete config types to avoid abstract container promotion
            # (which would make the struct non-isbits and incompatible with CuArray)
            ConfCOARE = SF.SurfaceFluxConfig{
                SF.COARE3RoughnessParams{FT},
                SF.ConstantGustinessSpec{FT},
                SF.MoistModel,
            }
            ConfRaupach = SF.SurfaceFluxConfig{
                SF.RaupachRoughnessParams{FT},
                SF.ConstantGustinessSpec{FT},
                SF.MoistModel,
            }
            ConfigType = Union{ConfCOARE, ConfRaupach}

            cpu_configs = Vector{ConfigType}(undef, n)
            for i in 1:n
                if isodd(i)
                    cpu_configs[i] = SF.SurfaceFluxConfig(
                        SF.COARE3RoughnessParams{FT}(),
                        SF.ConstantGustinessSpec(FT(1.0)),
                    )
                else
                    cpu_configs[i] = SF.SurfaceFluxConfig(
                        SF.RaupachRoughnessParams{FT}(),
                        SF.ConstantGustinessSpec(FT(1.0)),
                    )
                end
            end

            # Define roughness inputs (needed for Raupach)
            roughness_input = (LAI = FT(0.7), h = FT(8))

            # Displacement height: 0 for COARE, 2.5 for Raupach
            cpu_d = [isodd(i) ? FT(0) : FT(3.5) for i in 1:n]

            # CPU reference
            cpu_ρτxz = Vector{FT}(undef, n)
            cpu_shf = Vector{FT}(undef, n)
            cpu_lhf = Vector{FT}(undef, n)

            for i in 1:n
                result = SF.surface_fluxes(
                    param_set,
                    data.T_int[i], data.q_tot_int[i], data.ρ_int[i],
                    data.T_sfc[i], data.q_sfc[i],
                    FT(0), data.z[i], cpu_d[i],
                    (data.speed[i], FT(0)), (FT(0), FT(0)),
                    roughness_input,
                    cpu_configs[i],
                )
                cpu_ρτxz[i] = result.ρτxz
                cpu_shf[i] = result.shf
                cpu_lhf[i] = result.lhf
            end

            # GPU broadcast over heterogeneous config array
            gpu_configs = ArrayType(cpu_configs)
            T_int = ArrayType(data.T_int)
            T_sfc = ArrayType(data.T_sfc)
            q_tot_int = ArrayType(data.q_tot_int)
            q_sfc = ArrayType(data.q_sfc)
            ρ_int_array = ArrayType(fill(FT(1.13), n))
            z = ArrayType(data.z)
            speed = ArrayType(data.speed)
            Φ_sfc_array = ArrayType(fill(FT(0), n))
            d_array = ArrayType(fill(FT(0), n))
            # Construct tuple arrays on CPU first to avoid scalar indexing
            u_sfc_cpu = [(FT(0), FT(0)) for _ in 1:n]
            u_sfc_array = ArrayType(u_sfc_cpu)

            u_int_cpu = [(data.speed[i], FT(0)) for i in 1:n]
            u_int_array = ArrayType(u_int_cpu)
            d_array = ArrayType(cpu_d)
            gpu_roughness_inputs = Ref(roughness_input)

            gpu_results =
                SF.surface_fluxes.(
                    Ref(param_set), T_int, q_tot_int, ρ_int_array, T_sfc, q_sfc,
                    Φ_sfc_array, z, d_array, u_int_array, u_sfc_array,
                    gpu_roughness_inputs, gpu_configs,
                )

            gpu_ρτxz = map(x -> x.ρτxz, gpu_results)
            gpu_shf = map(x -> x.shf, gpu_results)
            gpu_lhf = map(x -> x.lhf, gpu_results)

            @test all(isfinite, Array(gpu_ρτxz))
            @test all(isfinite, Array(gpu_shf))
            @test all(isfinite, Array(gpu_lhf))
            @test isapprox(Array(gpu_ρτxz), cpu_ρτxz; rtol = FT(1e-5))
            @test isapprox(Array(gpu_shf), cpu_shf; rtol = FT(1e-5))
            @test isapprox(Array(gpu_lhf), cpu_lhf; rtol = FT(1e-5))
        end
    end

end  # if CUDA.functional()
