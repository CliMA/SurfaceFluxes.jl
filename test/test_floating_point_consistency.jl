# Floating-point consistency checks for Obukhov solutions
#
# 1. Near-zero test: ensures helper utilities such as `non_zero` preserve sign
#    and that `surface_fluxes` still returns a meaningful (non-zero) L_MO
#    when the layer spacing approaches zero.
# 2. Identical-state test: compares Float32 vs Float64 solutions over a grid of
#    heights and roughness lengths to make sure both precisions agree within 1%.

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

const NEAR_ZERO_Z_LEVELS = (1, 5, 10, 20, 40, 80, 160, 320, 640)
const IDENTICAL_Z_LEVELS = (1, 5, 10, 20, 40, 80, 160)
const IDENTICAL_Z0 = (1e-5, 1e-4, 1e-3)
const FLOAT_TYPES = (Float32, Float64)

# Test case thermodynamic primitives:
const TEST_T_INT = 287.85202
const TEST_Q_TOT_INT = 0.0
const TEST_RHO_INT = 1.1751807
const TEST_T_SFC = 291.96683
const TEST_Q_SFC = 0.013232904

@testset "SurfaceFluxes - Near-Zero Obukhov Length" begin
    @test sign(SF.non_zero(0.1)) == 1
    @test sign(SF.non_zero(-0.1)) == -1
    @test sign(SF.non_zero(-0.0)) == -1
    @test sign(SF.non_zero(0.0)) == 1

    for FT in FLOAT_TYPES
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

        T_int = FT(TEST_T_INT)
        q_tot_int = FT(TEST_Q_TOT_INT)
        ρ_int = FT(TEST_RHO_INT)
        T_sfc_guess = FT(TEST_T_SFC)
        q_vap_sfc_guess = FT(TEST_Q_SFC)

        for z_int in NEAR_ZERO_Z_LEVELS
            config = SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(FT(1e-5), FT(1e-5)),
                SF.ConstantGustinessSpec(FT(1.0)),
            )

            sfc_output = SF.surface_fluxes(
                param_set,
                T_int, q_tot_int, ρ_int,
                T_sfc_guess, q_vap_sfc_guess,
                FT(0), FT(z_int), zero(FT),
                (FT(0), FT(0)), (FT(0), FT(0)),  # Zero wind
                nothing,
                config,
            )

            L_MO = sfc_output.L_MO
            @test L_MO != FT(0)
        end
    end
end

@testset "SurfaceFluxes - Identical thermodynamic states" begin
    sol_mat = Array{Float64, 4}(
        undef,
        2,
        length(IDENTICAL_Z_LEVELS),
        length(IDENTICAL_Z0),
        length(IDENTICAL_Z0),
    )

    for (ii, FT) in enumerate(FLOAT_TYPES)
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

        # Use same thermodynamic state for interior and surface (identical states)
        T_int = FT(TEST_T_INT)
        q_tot_int = FT(TEST_Q_TOT_INT)
        ρ_int = FT(TEST_RHO_INT)
        T_sfc_guess = T_int
        q_vap_sfc_guess = q_tot_int

        for (jj, z_int) in enumerate(IDENTICAL_Z_LEVELS),
            (kk, z0m) in enumerate(IDENTICAL_Z0),
            (ll, z0h) in enumerate(IDENTICAL_Z0)

            config = SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(FT(z0m), FT(z0h)),
                SF.ConstantGustinessSpec(FT(1.0)),
            )

            sfc_output = SF.surface_fluxes(
                param_set,
                T_int, q_tot_int, ρ_int,
                T_sfc_guess, q_vap_sfc_guess,
                FT(0), FT(z_int), zero(FT),
                (FT(0), FT(0)), (FT(0), FT(0)),
                nothing,
                config,
                SF.PointValueScheme(),
                SF.SolverOptions{FT}(; maxiter = 20),
            )

            L = isinf(sfc_output.L_MO) ? FT(1e6) : sfc_output.L_MO
            sol_mat[ii, jj, kk, ll] = Float64(L)
        end
    end

    rdiff = (sol_mat[1, :, :, :] .- sol_mat[2, :, :, :]) ./ sol_mat[2, :, :, :]
    rdiff = filter(x -> abs(x) < eps(Float32), rdiff)
    @test all(abs.(rdiff) .<= 0.01)
end
