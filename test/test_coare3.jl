
using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaParams as CP
import SurfaceFluxes: COARE3RoughnessParams, momentum_roughness, scalar_roughness, charnock_parameter

@testset "COARE 3.0 Roughness" begin
    FT = Float64

    # Helper to check if charnock_parameter behaves as expected
    α_low = FT(0.011)
    α_high = FT(0.018)
    u_low = FT(10)
    u_high = FT(18)

    # 1. Low wind speed (<= 10 m/s)
    @test charnock_parameter(FT(5), α_low, α_high, u_low, u_high) == FT(0.011)
    # 2. High wind speed (>= 18 m/s)
    @test charnock_parameter(FT(20), α_low, α_high, u_low, u_high) == FT(0.018)
    # 3. Interpolation region (14 m/s)
    @test charnock_parameter(FT(14), α_low, α_high, u_low, u_high) ≈ FT(0.0145)  # Mid wind (linear)

    # Mock parameters
    # Note: kinematics_viscosity_of_air might not be in the default set if we are running without full ClimaParams suite loaded properly in tests sometimes,
    # but we can try to construct our spec manually.

    spec = COARE3RoughnessParams{FT}(kinematic_visc = 1.5e-5)

    # Mock structs for simplified testing
    struct MockParamSet{FT} <: SFP.AbstractSurfaceFluxesParameters{FT} end

    SFP.grav(::MockParamSet{FT}) where {FT} = FT(9.81)
    SFP.von_karman_const(::MockParamSet{FT}) where {FT} = FT(0.4)
    # UF params required for profile recovery
    SFP.uf_params(::MockParamSet{FT}) where {FT} = UF.BusingerParams(FT)

    param_set = MockParamSet{FT}()

    # Mock Context - No longer needed for roughness calls

    u_star = FT(0.3)

    # Test momentum roughness
    # 10m wind will be estimated via profile recovery.
    # With u_star = 0.3, z0_proxy = 1e-4, L = 100
    # Neutral wind at 10m approx u_star/k * log(10/1e-4) = 0.3/0.4 * log(1e5) = 0.75 * 11.5 = 8.6 m/s
    # Charnock alpha should be low limit 0.011
    # Smooth part: 0.11 * 1.5e-5 / 0.3 = 5.5e-6
    # Rough part: 0.011 * 0.3^2 / 9.81 = 0.011 * 0.09 / 9.81 ≈ 1e-4
    # Total z0m ≈ 1e-4

    z0m = momentum_roughness(spec, u_star, param_set, nothing)
    @test z0m > 0
    @test z0m ≈ 1.06e-4 atol = 5e-5

    # Test scalar roughness
    # Re_star = z0m * u_star / nu ≈ 1e-4 * 0.3 / 1.5e-5 ≈ 2
    # z0s = min(1.1e-4, 5.5e-5 * Re_star^-0.6)
    # 5.5e-5 * 2^-0.6 ≈ 5.5e-5 * 0.66 ≈ 3.6e-5
    z0s = scalar_roughness(spec, u_star, param_set, nothing)
    @test z0s > 0
    @test z0s < 1.1e-4
    @test z0s ≈ 3.6e-5 atol = 1e-5

end
