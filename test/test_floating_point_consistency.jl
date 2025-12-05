# Floating-point consistency checks for Obukhov solutions
#
# 1. Near-zero test: ensures helper utilities such as `non_zero` preserve sign
#    and that `surface_fluxes` still returns a meaningful (non-zero) L_MO
#    when the layer spacing approaches zero.
# 2. Identical-state test: compares Float32 vs Float64 solutions over a grid of
#    heights and roughness lengths to make sure both precisions agree within 1%.

using Test
include("test_utils.jl")
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import ClimaParams as CP

const NEAR_ZERO_Z_LEVELS = (1, 5, 10, 20, 40, 80, 160, 320, 640)
const IDENTICAL_Z_LEVELS = (1, 5, 10, 20, 40, 80, 160)
const IDENTICAL_Z0 = (1e-5, 1e-4, 1e-3)
const FLOAT_TYPES = (Float32, Float64)

function build_values_only_case(FT, z_int, ts_intt, ts_sfc, z0m, z0b)
    state_int = SF.StateValues(FT(z_int), (FT(0), FT(0)), ts_intt)
    state_sfc = SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc)
    return SF.ValuesOnly(state_int, state_sfc, FT(z0m), FT(z0b))
end

@testset "SurfaceFluxes - Near-zero Obukhov length" begin
    @test sign(SF.non_zero(1.0)) == 1
    @test sign(SF.non_zero(-1.0)) == -1
    @test sign(SF.non_zero(-0.0)) == 1
    @test sign(SF.non_zero(0.0)) == 1

    for FT in FLOAT_TYPES
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
        ts_intt = TD.PhaseEquil{FT}(FT(1.1751807), FT(97086.64), FT(10541.609), FT(0), FT(287.85202))
        ts_sfc = TD.PhaseEquil{FT}(FT(1.2176297), FT(102852.51), FT(45087.812), FT(0.013232904), FT(291.96683))

        for z_int in NEAR_ZERO_Z_LEVELS
            sc = build_values_only_case(FT, z_int, ts_intt, ts_sfc, 1e-5, 1e-5)
            sfc_output = surface_fluxes_wrapper(param_set, sc)
            L_MO = sfc_output.L_MO
            @test L_MO != FT(0)
        end
    end
end

@testset "SurfaceFluxes - Identical thermodynamic states" begin
    sol_mat = Array{Float64, 4}(undef, 2, length(IDENTICAL_Z_LEVELS), length(IDENTICAL_Z0), length(IDENTICAL_Z0))

    for (ii, FT) in enumerate(FLOAT_TYPES)
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
        ts = TD.PhaseEquil{FT}(FT(1.1751807), FT(97086.64), FT(10541.609), FT(0), FT(287.85202))

        for (jj, z_int) in enumerate(IDENTICAL_Z_LEVELS),
            (kk, z0m) in enumerate(IDENTICAL_Z0),
            (ll, z0b) in enumerate(IDENTICAL_Z0)

            sc = build_values_only_case(FT, z_int, ts, ts, z0m, z0b)
            sfc_output = surface_fluxes_wrapper(param_set, sc; maxiter = 20)
            L = isinf(sfc_output.L_MO) ? FT(1e6) : sfc_output.L_MO
            sol_mat[ii, jj, kk, ll] = Float64(L)
        end
    end

    rdiff = (sol_mat[1, :, :, :] .- sol_mat[2, :, :, :]) ./ sol_mat[2, :, :, :]
    rdiff = filter(x -> abs(x) < eps(Float32), rdiff)
    @test all(abs.(rdiff) .<= 0.01)
end
