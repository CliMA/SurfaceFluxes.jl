# `compute_profile_value` evaluates the dimensionless profile of the discretization
# scheme it is given: the closed forms of the Businger profiles in neutral and stable
# conditions (ψ = -a ζ, Ψ = -a ζ / 2) for point values and layer averages.

module TestProfileRecovery

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

@testset "Profile Recovery" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    uf = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)
    a_m, a_h, Pr_0 = UF.a_m(uf), UF.a_h(uf), UF.Pr_0(uf)

    z0 = FT(0.1)
    Δz_eff = FT(10)
    scale = FT(0.5) # u_star or theta_star
    val_sfc = FT(300)
    point = UF.PointValueScheme()
    layer = UF.LayerAverageScheme()
    momentum = UF.MomentumTransport()
    heat = UF.HeatTransport()
    profile(L, transport, args...) =
        SF.compute_profile_value(
            param_set,
            L,
            z0,
            Δz_eff,
            scale,
            val_sfc,
            transport,
            args...,
        )
    value(F) = val_sfc + scale / κ * F

    # Neutral: logarithmic point profile, and its layer average
    # (1/Δz) ∫_{z0}^{Δz} ln(z/z0) dz = ln(Δz/z0) - 1 + z0/Δz
    @test profile(FT(Inf), momentum, point) ≈ value(log(Δz_eff / z0))
    @test profile(FT(Inf), heat, layer) ≈
          value(Pr_0 * (log(Δz_eff / z0) - 1 + z0 / Δz_eff))

    # Stable (ζ = 1): the linear stable functions give closed forms
    L = FT(10)
    ζ = Δz_eff / L
    ζ0 = ζ * z0 / Δz_eff
    R = 1 - z0 / Δz_eff
    @test profile(L, momentum, point) ≈ value(log(Δz_eff / z0) + a_m * ζ * R)
    F_h_layer =
        Pr_0 * log(Δz_eff / z0) + a_h * ζ / 2 - (z0 / Δz_eff) * a_h * ζ0 / 2 +
        R * (-a_h * ζ0 - Pr_0)
    @test profile(L, heat, layer) ≈ value(F_h_layer)

    # The defaults are point values without a roughness-sublayer correction
    @test profile(L, momentum) == profile(L, momentum, point, SF.NoRoughnessSubLayer())
end

end # module
