using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions: BusingerParams

@testset "Profile Recovery" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    κ = SFP.von_karman_const(param_set)

    # Inputs
    L_MO = FT(10)
    z0 = FT(0.1)
    Δz = FT(10)
    scale = FT(0.5) # u_star or theta_star
    val_sfc = FT(300) # Surface value

    # Test Momentum with PointValueScheme
    scheme = UF.PointValueScheme()
    transport = UF.MomentumTransport()

    # Expected: val_sfc + (scale / k) * F_m
    uf_params = SFP.uf_params(param_set)
    ζ = Δz / L_MO
    F_m = UF.dimensionless_profile(uf_params, Δz, ζ, z0, transport, scheme)
    expected = val_sfc + (scale / κ) * F_m

    result = SF.compute_profile_value(
        param_set,
        L_MO,
        z0,
        Δz,
        scale,
        val_sfc,
        transport,
        scheme,
    )

    @test result ≈ expected

    # Test Heat with LayerAverageScheme
    scheme_fv = UF.LayerAverageScheme()
    transport_h = UF.HeatTransport()

    F_h = UF.dimensionless_profile(uf_params, Δz, ζ, z0, transport_h, scheme_fv)
    expected_h = val_sfc + (scale / κ) * F_h

    result_h = SF.compute_profile_value(
        param_set,
        L_MO,
        z0,
        Δz,
        scale,
        val_sfc,
        transport_h,
        scheme_fv,
    )

    @test result_h ≈ expected_h

    # Check that PointValueScheme is default
    result_default = SF.compute_profile_value(
        param_set,
        L_MO,
        z0,
        Δz,
        scale,
        val_sfc,
        transport,
    )
    @test result_default ≈ result
end
