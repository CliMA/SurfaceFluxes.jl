using Test
using SurfaceFluxes
import DifferentiationInterface as DI
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
using ForwardDiff: ForwardDiff
using Enzyme: Enzyme

# `SolverOptions` and RootSolvers' solution structs mix FT fields with Int/Bool
# fields; Enzyme cannot statically classify all memory as float vs. non-float.
Enzyme.API.looseTypeAnalysis!(true)

function ad_backends()
    return (
        ("ForwardDiff", DI.AutoForwardDiff()),
        (
            "Enzyme (Forward)",
            DI.AutoEnzyme(
                mode = Enzyme.set_runtime_activity(
                    Enzyme.set_strong_zero(Enzyme.Forward),
                ),
            ),
        ),
        (
            "Enzyme (Reverse)",
            DI.AutoEnzyme(
                mode = Enzyme.set_runtime_activity(
                    Enzyme.set_strong_zero(Enzyme.Reverse),
                ),
            ),
        ),
    )
end

@testset verbose = true "AD Compatibility - Finite Difference Validation" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    # 1. Define physical state
    T_int = FT(300)
    q_int = FT(0.01)
    ρ_int = FT(1.1)
    q_sfc = FT(0.02)
    z = FT(20)
    u_int = (FT(10), FT(0))
    u_sfc = (FT(0), FT(0))
    dist = FT(0)

    # 2. Define wrapper function: Surface Flux (SHF) vs T_sfc
    #    Uses full MOST solving inside `surface_fluxes`.
    function compute_shf(T_sfc_val)
        result = SF.surface_fluxes(
            param_set,
            T_int,
            q_int,
            FT(0),
            FT(0),
            ρ_int,
            T_sfc_val,
            q_sfc,
            FT(0),
            z,
            dist,
            u_int,
            u_sfc,
        )
        return result.shf
    end

    # 3. Test across different stability regimes
    # T_int = 300 K.
    # T_sfc = 295 K => Stable (Ri > 0)
    # T_sfc = 300.1 K => Near Neutral / Slightly Unstable
    # T_sfc = 305 K => Unstable (Ri < 0)
    T_sfc_range = FT[295, 300.1, 305]

    for (backend_name, backend) in ad_backends()
        @testset "AD backend: $backend_name" begin
            for T_sfc_base in T_sfc_range
                # 4. AD Derivative
                dSHF_dT_ad = DI.derivative(compute_shf, backend, T_sfc_base)

                # 5. Finite Difference Approximation (Central Difference)
                ϵ = FT(1e-4)
                shf_plus = compute_shf(T_sfc_base + ϵ)
                shf_minus = compute_shf(T_sfc_base - ϵ)
                dSHF_dT_fd = (shf_plus - shf_minus) / (2ϵ)

                # 6. Comparison
                @info "Comparing AD and FD derivatives" T_sfc = T_sfc_base dSHF_dT_ad dSHF_dT_fd

                # Use relatively generous O(ϵ) tolerance
                @test isapprox(dSHF_dT_ad, dSHF_dT_fd, rtol = ϵ)
            end
        end
    end
end

@testset verbose = true "AD Compatibility - Roughness Sublayer" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    # Canopy-like state so the RSL correction is active (Δz_eff within / above z_RSL).
    T_int = FT(290)
    q_int = FT(0.008)
    ρ_int = FT(1.2)
    q_sfc = FT(0.01)
    Δz = FT(40)
    d = FT(7)
    u_int = (FT(5), FT(0))
    u_sfc = (FT(0), FT(0))
    roughness = SF.ConstantRoughnessParams{FT}(z0m = FT(1.0), z0s = FT(0.1))
    gustiness = SF.ConstantGustinessSpec(FT(0.001))
    T_sfc_range = FT[288, 290.1, 293]

    rsl_models = (
        (
            "PhysickGarrattRSL",
            SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(20.0)),
        ),
        (
            "HarmanFinniganRSL",
            SF.HarmanFinniganRSL{FT}(c1_m = FT(0.5), c1_h = FT(0.5), z_RSL = FT(20.0)),
        ),
    )

    for (rsl_name, rsl_model) in rsl_models
        config = SF.SurfaceFluxConfig(roughness, gustiness, SF.DryModel(), rsl_model)

        function compute_shf_rsl(T_sfc_val)
            result = SF.surface_fluxes(
                param_set,
                T_int,
                q_int,
                FT(0),
                FT(0),
                ρ_int,
                T_sfc_val,
                q_sfc,
                FT(0),
                Δz,
                d,
                u_int,
                u_sfc,
                nothing,
                config,
            )
            return result.shf
        end

        @testset "$rsl_name — surface_fluxes d(SHF)/d(T_sfc)" begin
            for (backend_name, backend) in ad_backends()
                @testset "AD backend: $backend_name" begin
                    for T_sfc_base in T_sfc_range
                        dSHF_dT_ad = DI.derivative(compute_shf_rsl, backend, T_sfc_base)
                        ϵ = FT(1e-4)
                        dSHF_dT_fd =
                            (
                                compute_shf_rsl(T_sfc_base + ϵ) -
                                compute_shf_rsl(T_sfc_base - ϵ)
                            ) / (2ϵ)
                        @test isfinite(dSHF_dT_ad)
                        @test isapprox(dSHF_dT_ad, dSHF_dT_fd, rtol = ϵ, atol = FT(1e-4))
                    end
                end
            end
        end
    end

    # Direct AD through the HF quadrature kernel (ForwardDiff): exercises
    # gauss_legendre4 + the log-mapped integrand under Dual arithmetic.
    @testset "HarmanFinniganRSL — rsl_profile_correction dP/dc1" begin
        Δz_eff = FT(40)
        z0m = FT(0.5)
        z_RSL = FT(30)
        function P_of_c1(c1)
            rsl = SF.HarmanFinniganRSL(c1_m = c1, c1_h = c1, z_RSL = typeof(c1)(z_RSL))
            return SF.rsl_profile_correction(rsl, Δz_eff, z0m, UF.MomentumTransport())
        end
        c1_base = FT(0.5)
        dP_ad = ForwardDiff.derivative(P_of_c1, c1_base)
        ϵ = FT(1e-6)
        dP_fd = (P_of_c1(c1_base + ϵ) - P_of_c1(c1_base - ϵ)) / (2ϵ)
        @test isfinite(dP_ad)
        @test dP_ad < 0  # stronger c1 → more negative P
        @test isapprox(dP_ad, dP_fd, rtol = FT(1e-4), atol = FT(1e-6))
    end

    @testset "PhysickGarrattRSL — rsl_profile_correction dP/dc" begin
        Δz_eff = FT(40)
        z0m = FT(0.5)
        z_RSL = FT(30)
        function P_of_c(c)
            rsl = SF.PhysickGarrattRSL(c_m = c, c_h = c, z_RSL = typeof(c)(z_RSL))
            return SF.rsl_profile_correction(rsl, Δz_eff, z0m, UF.MomentumTransport())
        end
        c_base = FT(0.4)
        dP_ad = ForwardDiff.derivative(P_of_c, c_base)
        ϵ = FT(1e-6)
        dP_fd = (P_of_c(c_base + ϵ) - P_of_c(c_base - ϵ)) / (2ϵ)
        @test isfinite(dP_ad)
        @test dP_ad < 0
        @test isapprox(dP_ad, dP_fd, rtol = FT(1e-4), atol = FT(1e-6))
    end
end
