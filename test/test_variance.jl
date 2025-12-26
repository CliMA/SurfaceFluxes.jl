module TestVariance

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import ClimaParams as CP

FT = Float32
param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

@testset "Variance Functions" begin
    # Construct standard inputs
    config = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(1e-3), FT(1e-3)),
        SF.ConstantGustinessSpec(FT(0.0)),
    )
    flux_specs = SF.FluxSpecs{FT}()

    inputs = SF.build_surface_flux_inputs(
        FT(300),    # T_int
        FT(0),      # q_tot_int
        FT(0),      # q_liq_int
        FT(0),      # q_ice_int
        FT(1.2),    # ρ_int
        FT(300),    # T_sfc_guess
        FT(0),      # q_vap_sfc_guess
        FT(0),      # Φ_sfc
        FT(10),     # Δz
        FT(0),      # d
        (FT(10), FT(0)), # u_int
        (FT(0), FT(0)),  # u_sfc
        config,
        nothing,    # roughness_inputs
        flux_specs,
        nothing,
        nothing,
    )
    ustar = FT(0.3)
    theta_star = FT(0.1)

    # 1. Unstable Conditions (L_MO < 0)
    L_MO = FT(-10.0) # ζ = 10 / -10 = -1

    # Calculate expected values
    # Calculate expected values for Tan et al. (2018)

    # Momentum: 3.75 + 0.2 * (w_*/u_*)^2 + (-ζ)^(2/3)
    # Calculate w_* first
    zi = SFP.gustiness_zi(param_set)
    κ = SFP.von_karman_const(param_set)
    # w_* = u_* * (zi / (-k * L))^(1/3)
    w_star = ustar * cbrt(zi / (-κ * L_MO))

    # Calculate effective height
    Δz_eff = inputs.Δz - inputs.d # d=0 here
    ζ = Δz_eff / L_MO
    # Expected phi^2
    expected_phi_u_sq = 3.75 + 0.2 * (w_star / ustar)^2 + cbrt(-ζ)^2

    var_u = SF.u_variance(param_set, Δz_eff, ustar, ζ)
    # var_u returns (u_* * phi)^2 = u_*^2 * phi^2
    @test isapprox(var_u, ustar^2 * expected_phi_u_sq; rtol = 1e-4)

    # Heat: Tan et al: 2 * (1 - 8.3ζ)^(-1/3)
    expected_phi_theta = 2.0 * cbrt(1 / (1 - 8.3 * ζ))

    var_theta = SF.scalar_variance(param_set, theta_star, ζ)
    sigma_theta = sqrt(var_theta)

    @test isapprox(sigma_theta, theta_star * expected_phi_theta; rtol = 1e-4)

    # 2. Stable Conditions (L_MO > 0)
    L_MO_stable = FT(10.0)

    # Momentum: Constant 3.75 (Tan et al. stable)
    ζ_stable = Δz_eff / L_MO_stable
    var_u_stable = SF.u_variance(param_set, Δz_eff, ustar, ζ_stable)
    @test isapprox(var_u_stable, (ustar * sqrt(3.75))^2; rtol = 1e-4)

    # Heat: Constant 2.0
    var_theta_stable = SF.scalar_variance(param_set, theta_star, ζ_stable)
    @test isapprox(sqrt(var_theta_stable), theta_star * 2.0; rtol = 1e-4)

    # 3. Neutral (L_MO -> Inf, falls into stable branch ζ >= 0)
    L_MO_neutral = FT(1e10)

    var_u_neutral = SF.u_variance(param_set, Δz_eff, ustar, Δz_eff / L_MO_neutral)
    @test isapprox(sqrt(var_u_neutral), ustar * sqrt(3.75); rtol = 1e-4)

    # 4. Theta Variance Wrapper (SHF input)
    shf = FT(100.0) # W/m2 
    rho = FT(1.2)
    thermo_params = SFP.thermodynamics_params(param_set)
    cp = TD.cp_m(thermo_params)
    theta_star_calc = -shf / (rho * cp * ustar)

    ζ_wrapper = Δz_eff / L_MO
    var_theta_wrapper = SF.theta_variance(param_set, inputs, shf, ustar, ζ_wrapper, rho)
    var_theta_explicit = SF.scalar_variance(param_set, theta_star_calc, ζ_wrapper)

    @test isapprox(var_theta_wrapper, var_theta_explicit; rtol = 1e-6)
    @test isapprox(var_theta_wrapper, var_theta_explicit; rtol = 1e-6)


end

@testset "Variance Functions - Tan et al. (2018) Unit Tests" begin
    # Direct unit tests for UF.phi to ensure correct implementation of Tan et al. (2018)
    uf_params = SFP.uf_params(param_set)

    # Momentum Variance
    # phi_m_var(p, ζ, u_star, w_star, MomentumVariance())
    # Eq: 3.75 + 0.2 (w_*/u_*)^2 + (-ζ)^(2/3) (Unstable)
    #     3.75 (Stable)

    u_star = FT(1.0)
    w_star = FT(2.0)
    mv = UF.MomentumVariance()

    # Stable case (ζ > 0)
    ζ_stable = FT(1.0)
    phi_stable = UF.phi(uf_params, ζ_stable, u_star, w_star, mv)
    @test phi_stable ≈ sqrt(FT(3.75))

    # Unstable case (ζ < 0)
    ζ_unstable = FT(-1.0)
    phi_unstable = UF.phi(uf_params, ζ_unstable, u_star, w_star, mv)

    # Expected: sqrt(3.75 + 0.2 * (2.0/1.0)^2 + (-(-1.0))^(2/3))
    expected_tke_norm = FT(3.75) + FT(0.2) * (w_star / u_star)^2 + cbrt(-ζ_unstable)^2
    @test phi_unstable ≈ sqrt(expected_tke_norm)

    # Heat Variance
    # phi_h_var(p, ζ, HeatVariance())
    # Eq: 2 * (1 - 8.3 * ζ)^(-1/3) (Unstable)
    #     2.0 (Stable)

    hv = UF.HeatVariance()

    # Stable case
    phi_h_stable = UF.phi(uf_params, ζ_stable, hv)
    @test phi_h_stable ≈ FT(2.0)

    # Unstable case
    phi_h_unstable = UF.phi(uf_params, ζ_unstable, hv)
    # Expected: 2 * (1 - 8.3 * (-1))^(-1/3)
    #         = 2 * (9.3)^(-1/3)
    expected_h = FT(2.0) * cbrt(1 / (1 - 8.3 * ζ_unstable))
    @test phi_h_unstable ≈ expected_h
end

end # module
