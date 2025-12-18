using Test
using SurfaceFluxes
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
using ForwardDiff

@testset "AD Compatibility - Finite Difference Validation" begin
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
            T_int, q_int, ρ_int,
            T_sfc_val, q_sfc,
            FT(0), z, dist,
            u_int, u_sfc,
        )
        return result.shf
    end

    # 3. Test across different stability regimes
    # T_int = 300 K.
    # T_sfc = 295 K => Stable (Ri > 0)
    # T_sfc = 300.1 K => Near Neutral / Slightly Unstable 
    # T_sfc = 305 K => Unstable (Ri < 0)
    T_sfc_range = FT[295, 300.1, 305]

    for T_sfc_base in T_sfc_range
        # 4. AD Derivative
        dSHF_dT_ad = ForwardDiff.derivative(compute_shf, T_sfc_base)

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
