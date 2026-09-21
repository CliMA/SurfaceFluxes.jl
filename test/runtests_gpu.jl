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
                    SF.ConstantRoughnessParams(FT(data.z0[i]), FT(z0h)),
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
                    data.T_int[i], data.q_tot_int[i], FT(0), FT(0), data.ρ_int[i],
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
                    Ref(param_set), T_int, q_tot_int, Ref(FT(0)), Ref(FT(0)), ρ_int_array,
                    T_sfc, q_sfc,
                    Φ_sfc_array, z, d_array, u_int_array, u_sfc_array,
                    Ref(nothing), gpu_configs,
                    Ref(SF.PointValueScheme()),
                    Ref(SF.SolverOptions{FT}(tol = FT(1e-2), maxiter = 15)),
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
                    SF.ConstantRoughnessParams(z0m, z0h),
                    SF.ConstantGustinessSpec(gustiness_values[i]),
                ) for i in 1:n
            ]

            # CPU reference
            cpu_shf = Vector{FT}(undef, n)
            cpu_lhf = Vector{FT}(undef, n)

            for i in 1:n
                result = SF.surface_fluxes(
                    param_set,
                    data.T_int[i], data.q_tot_int[i], FT(0), FT(0), data.ρ_int[i],
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
                    Ref(param_set), T_int, q_tot_int, Ref(FT(0)), Ref(FT(0)), ρ_int_array,
                    T_sfc, q_sfc,
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
            # Fully concrete config types (including RSL) so the Union is an
            # isbitsunion that CuArray can allocate inline.
            ConfCOARE = SF.SurfaceFluxConfig{
                SF.COARE3RoughnessParams{FT},
                SF.ConstantGustinessSpec{FT},
                SF.MoistModel,
                SF.NoRoughnessSubLayer,
            }
            ConfRaupach = SF.SurfaceFluxConfig{
                SF.RaupachRoughnessParams{FT},
                SF.ConstantGustinessSpec{FT},
                SF.MoistModel,
                SF.NoRoughnessSubLayer,
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
                    data.T_int[i], data.q_tot_int[i], FT(0), FT(0), data.ρ_int[i],
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
                    Ref(param_set), T_int, q_tot_int, Ref(FT(0)), Ref(FT(0)), ρ_int_array,
                    T_sfc, q_sfc,
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

    # Surface-state callback broadcast on GPU.
    #
    # The `update_T_sfc` / `update_q_vap_sfc` callbacks (used by ClimaLand canopy
    # coupling) are passed as trailing positional arguments to `surface_fluxes`.
    # Because these test closures capture only isbits `FT` scalars, the closure
    # itself is isbits and can be broadcast over a `CuArray` via `Ref(callback)`.
    # These testsets mirror the CPU callback cases in
    # `test_supercritical_stability.jl`, verifying the GPU path reproduces the
    # CPU reference through the coupled canopy solve.

    @testset "GPU broadcast - Supercritical Canopy Callbacks" begin
        for FT in (Float32, Float64),
            uf_type in (BusingerParams, SF.UniversalFunctions.GryanikParams)

            param_set = SFP.SurfaceFluxesParameters(FT, uf_type)

            # Canopy coupling constants (point-independent; captured by closures)
            AI = FT(3.1223678081630233)
            leaf_Cd = FT(0.0726)
            g_stomata = FT(5.848035071201371e-6)

            # Canopy energy balance: T_sfc weighted between air and canopy states
            # by the ratio of canopy conductance to aerodynamic conductance.
            update_T_sfc =
                (ζ, ps, thermo_params, inputs, scheme, u_star, z0m, z0h) -> begin
                    Φ_sfc = SF.surface_geopotential(inputs)
                    Φ_int = SF.interior_geopotential(ps, inputs)
                    g_h = SF.heat_conductance(ps, ζ, u_star, inputs, z0m, z0h, scheme)
                    g_land = leaf_Cd * u_star * AI
                    cp_d = SFP.cp_d(ps)
                    r = g_land / g_h
                    return (
                        inputs.T_int + inputs.T_sfc_guess * r + (Φ_int - Φ_sfc) / cp_d
                    ) /
                           (1 + r)
                end

            # Canopy moisture balance with stomatal + leaf boundary-layer conductance.
            update_q_vap_sfc =
                (ζ, ps, thermo_params, inputs, scheme, T_sfc, u_star, z0m, z0h) -> begin
                    g_leaf = leaf_Cd * u_star * AI
                    g_land = g_stomata * g_leaf / (g_leaf + g_stomata)
                    g_h = SF.heat_conductance(ps, ζ, u_star, inputs, z0m, z0h, scheme)
                    q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
                    r = g_land / g_h
                    return (r * inputs.q_vap_sfc_guess + q_vap_int) / (1 + r)
                end

            # Base supercritical canopy state (strong inversion, weak wind),
            # replicated across points with the canopy temperature varied so the
            # broadcast exercises distinct per-point solves.
            T_int_c = FT(294.673095703125)
            q_tot_int_c = FT(0.009545918211858238)
            ρ_int_c = FT(1.157759649975361)
            q_canopy_c = FT(0.008910620278696149)
            Δz_c = FT(10)
            d_c = FT(0.0573)
            speed_c = FT(1.3673065900802612)

            T_canopy_base = FT(284.6341213016535)
            n = 4
            T_canopy_vals = [T_canopy_base + FT(δ) for δ in (-2, -1, 0, 1)]

            config_c = SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(FT(0.359), FT(0.0544)),
                SF.ConstantGustinessSpec(FT(1)),
            )
            cpu_configs = [config_c for _ in 1:n]

            # CPU reference (scalar solves)
            cpu_shf = Vector{FT}(undef, n)
            cpu_lhf = Vector{FT}(undef, n)
            cpu_ustar = Vector{FT}(undef, n)
            cpu_ζ = Vector{FT}(undef, n)
            cpu_T_sfc = Vector{FT}(undef, n)
            for i in 1:n
                result = SF.surface_fluxes(
                    param_set,
                    T_int_c, q_tot_int_c, FT(0), FT(0), ρ_int_c,
                    T_canopy_vals[i], q_canopy_c,
                    FT(0), Δz_c, d_c,
                    (speed_c, FT(0)), (FT(0), FT(0)),
                    nothing, config_c,
                    SF.PointValueScheme(), nothing, nothing,
                    update_T_sfc, update_q_vap_sfc,
                )
                cpu_shf[i] = result.shf
                cpu_lhf[i] = result.lhf
                cpu_ustar[i] = result.ustar
                cpu_ζ[i] = result.ζ
                cpu_T_sfc[i] = result.T_sfc
            end

            # GPU broadcast over the point array
            gpu_configs = ArrayType(cpu_configs)
            T_int_array = ArrayType(fill(T_int_c, n))
            q_tot_int_array = ArrayType(fill(q_tot_int_c, n))
            ρ_int_array = ArrayType(fill(ρ_int_c, n))
            T_canopy_array = ArrayType(T_canopy_vals)
            q_canopy_array = ArrayType(fill(q_canopy_c, n))
            u_int_array = ArrayType([(speed_c, FT(0)) for _ in 1:n])
            u_sfc_array = ArrayType([(FT(0), FT(0)) for _ in 1:n])

            gpu_results =
                SF.surface_fluxes.(
                    Ref(param_set), T_int_array, q_tot_int_array, Ref(FT(0)), Ref(FT(0)),
                    ρ_int_array,
                    T_canopy_array, q_canopy_array,
                    Ref(FT(0)), Ref(Δz_c), Ref(d_c),
                    u_int_array, u_sfc_array,
                    Ref(nothing), gpu_configs,
                    Ref(SF.PointValueScheme()), Ref(nothing), Ref(nothing),
                    Ref(update_T_sfc), Ref(update_q_vap_sfc),
                )

            gpu_shf = Array(map(x -> x.shf, gpu_results))
            gpu_lhf = Array(map(x -> x.lhf, gpu_results))
            gpu_ustar = Array(map(x -> x.ustar, gpu_results))
            gpu_ζ = Array(map(x -> x.ζ, gpu_results))
            gpu_T_sfc = Array(map(x -> x.T_sfc, gpu_results))

            @test all(isfinite, gpu_shf)
            @test all(isfinite, gpu_lhf)
            @test all(isfinite, gpu_ustar)
            @test all(isfinite, gpu_ζ)

            @test isapprox(gpu_shf, cpu_shf; rtol = FT(1e-3), atol = FT(1e-4))
            @test isapprox(gpu_lhf, cpu_lhf; rtol = FT(1e-3), atol = FT(1e-4))
            @test isapprox(gpu_ustar, cpu_ustar; rtol = FT(1e-3), atol = FT(1e-4))
            @test isapprox(gpu_ζ, cpu_ζ; rtol = FT(1e-3), atol = FT(1e-4))
            @test isapprox(gpu_T_sfc, cpu_T_sfc; rtol = FT(1e-3), atol = FT(1e-4))
        end
    end

    @testset "GPU broadcast - Opposite-branch Callback" begin
        for FT in (Float32, Float64)
            param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)

            grav = SFP.grav(param_set)
            cp_d = SFP.cp_d(param_set)
            T_int_o = FT(295)
            q_o = FT(0.005)
            ρ_o = FT(1.15)
            Δz_o = FT(10)
            T_neutral = T_int_o + grav * Δz_o / cp_d

            # Synthetic callback: weakly stable (supercritical, no stable root)
            # for ζ >= 0, strongly unstable for ζ < 0. The only root lives on the
            # unstable branch and must be found via the opposite-branch probe.
            δT = FT(0.1)
            update_T_sfc_flip =
                (ζ, ps, thermo_params, inputs, scheme, u_star, z0m, z0h) -> begin
                    return ζ >= 0 ? T_neutral - δT - 3 * ζ / (1 + ζ) :
                           T_neutral - δT + 5 * (-ζ) / (1 - ζ)
                end

            config_o = SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(FT(0.01), FT(0.001)),
                SF.ConstantGustinessSpec(FT(1)),
            )

            # Vary wind speed per point to exercise distinct broadcast solves.
            n = 4
            speed_vals = FT.([0.3, 0.5, 0.7, 0.9])
            T_sfc_guess = T_neutral - FT(0.1)
            cpu_configs = [config_o for _ in 1:n]

            cpu_shf = Vector{FT}(undef, n)
            cpu_ustar = Vector{FT}(undef, n)
            cpu_ζ = Vector{FT}(undef, n)
            for i in 1:n
                result = SF.surface_fluxes(
                    param_set,
                    T_int_o, q_o, FT(0), FT(0), ρ_o,
                    T_sfc_guess, q_o,
                    FT(0), Δz_o, FT(0),
                    (speed_vals[i], FT(0)), (FT(0), FT(0)),
                    nothing, config_o,
                    SF.PointValueScheme(), nothing, nothing,
                    update_T_sfc_flip, nothing,
                )
                cpu_shf[i] = result.shf
                cpu_ustar[i] = result.ustar
                cpu_ζ[i] = result.ζ
            end

            gpu_configs = ArrayType(cpu_configs)
            T_int_array = ArrayType(fill(T_int_o, n))
            q_array = ArrayType(fill(q_o, n))
            ρ_array = ArrayType(fill(ρ_o, n))
            T_sfc_array = ArrayType(fill(T_sfc_guess, n))
            u_int_array = ArrayType([(speed_vals[i], FT(0)) for i in 1:n])
            u_sfc_array = ArrayType([(FT(0), FT(0)) for _ in 1:n])

            gpu_results =
                SF.surface_fluxes.(
                    Ref(param_set), T_int_array, q_array, Ref(FT(0)), Ref(FT(0)), ρ_array,
                    T_sfc_array, q_array,
                    Ref(FT(0)), Ref(Δz_o), Ref(FT(0)),
                    u_int_array, u_sfc_array,
                    Ref(nothing), gpu_configs,
                    Ref(SF.PointValueScheme()), Ref(nothing), Ref(nothing),
                    Ref(update_T_sfc_flip), Ref(nothing),
                )

            gpu_shf = Array(map(x -> x.shf, gpu_results))
            gpu_ustar = Array(map(x -> x.ustar, gpu_results))
            gpu_ζ = Array(map(x -> x.ζ, gpu_results))

            @test all(isfinite, gpu_shf)
            @test all(isfinite, gpu_ustar)
            @test all(isfinite, gpu_ζ)

            # Root should sit on the unstable branch (warm surface -> upward SHF)
            @test all(gpu_shf .> 0)
            @test all(gpu_ustar .> 0)

            @test isapprox(gpu_shf, cpu_shf; rtol = FT(1e-3), atol = FT(1e-4))
            @test isapprox(gpu_ustar, cpu_ustar; rtol = FT(1e-3), atol = FT(1e-4))
            @test isapprox(gpu_ζ, cpu_ζ; rtol = FT(1e-3), atol = FT(1e-4))
        end
    end

    # Roughness-sublayer configs on GPU: canopy-like states with PG / HF RSL.
    @testset "GPU broadcast - Roughness Sublayer (PG + HF)" begin
        for FT in (Float32, Float64)
            param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)

            # Forest-like column: Δz = 40 m, d = 7 m, z_RSL = 20 m above d.
            n = 4
            T_int_vals = FT.([288.0, 290.0, 292.0, 295.0])
            T_sfc_vals = FT.([286.0, 289.0, 293.0, 298.0])
            speed_vals = FT.([2.0, 3.5, 5.0, 7.0])
            Δz = FT(40)
            d = FT(7)
            q = FT(0.008)
            ρ = FT(1.2)
            roughness = SF.ConstantRoughnessParams{FT}(z0m = FT(1.0), z0s = FT(0.1))
            gustiness = SF.ConstantGustinessSpec(FT(0.001))

            rsl_cases = (
                (
                    "PhysickGarrattRSL",
                    SF.PhysickGarrattRSL{FT}(
                        c_m = FT(0.4),
                        c_h = FT(0.4),
                        z_RSL = FT(20.0),
                    ),
                ),
                (
                    "HarmanFinniganRSL",
                    SF.HarmanFinniganRSL{FT}(
                        c1_m = FT(0.5),
                        c1_h = FT(0.5),
                        z_RSL = FT(20.0),
                    ),
                ),
            )

            for (rsl_name, rsl_model) in rsl_cases
                @testset "$rsl_name ($FT)" begin
                    config =
                        SF.SurfaceFluxConfig(roughness, gustiness, SF.DryModel(), rsl_model)
                    cpu_configs = [config for _ in 1:n]

                    cpu_shf = Vector{FT}(undef, n)
                    cpu_ustar = Vector{FT}(undef, n)
                    cpu_Cd = Vector{FT}(undef, n)
                    for i in 1:n
                        result = SF.surface_fluxes(
                            param_set,
                            T_int_vals[i],
                            q,
                            FT(0),
                            FT(0),
                            ρ,
                            T_sfc_vals[i],
                            q,
                            FT(0),
                            Δz,
                            d,
                            (speed_vals[i], FT(0)),
                            (FT(0), FT(0)),
                            nothing,
                            config,
                        )
                        cpu_shf[i] = result.shf
                        cpu_ustar[i] = result.ustar
                        cpu_Cd[i] = result.Cd
                    end

                    gpu_configs = ArrayType(cpu_configs)
                    T_int_array = ArrayType(T_int_vals)
                    T_sfc_array = ArrayType(T_sfc_vals)
                    q_array = ArrayType(fill(q, n))
                    ρ_array = ArrayType(fill(ρ, n))
                    u_int_array = ArrayType([(speed_vals[i], FT(0)) for i in 1:n])
                    u_sfc_array = ArrayType([(FT(0), FT(0)) for _ in 1:n])
                    Φ_array = ArrayType(fill(FT(0), n))
                    Δz_array = ArrayType(fill(Δz, n))
                    d_array = ArrayType(fill(d, n))

                    gpu_results =
                        SF.surface_fluxes.(
                            Ref(param_set),
                            T_int_array,
                            q_array,
                            Ref(FT(0)),
                            Ref(FT(0)),
                            ρ_array,
                            T_sfc_array,
                            q_array,
                            Φ_array,
                            Δz_array,
                            d_array,
                            u_int_array,
                            u_sfc_array,
                            Ref(nothing),
                            gpu_configs,
                            Ref(SF.PointValueScheme()),
                            Ref(SF.SolverOptions{FT}(tol = FT(1e-2), maxiter = 15)),
                            Ref(SF.FluxSpecs{FT}()),
                        )

                    gpu_shf = Array(map(x -> x.shf, gpu_results))
                    gpu_ustar = Array(map(x -> x.ustar, gpu_results))
                    gpu_Cd = Array(map(x -> x.Cd, gpu_results))

                    @test all(isfinite, gpu_shf)
                    @test all(isfinite, gpu_ustar)
                    @test all(isfinite, gpu_Cd)
                    @test all(gpu_ustar .> 0)
                    @test all(gpu_Cd .> 0)
                    @test isapprox(gpu_shf, cpu_shf; rtol = FT(1e-5), atol = FT(1e-5))
                    @test isapprox(gpu_ustar, cpu_ustar; rtol = FT(1e-5), atol = FT(1e-5))
                    @test isapprox(gpu_Cd, cpu_Cd; rtol = FT(1e-5), atol = FT(1e-5))
                end
            end
        end
    end

end  # if CUDA.functional()
