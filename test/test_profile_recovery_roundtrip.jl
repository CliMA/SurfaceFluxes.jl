# Profile recovery round-trip test.
#
# Verifies that solving for fluxes at height z, then recovering the profile value
# at z using compute_profile_value, returns the original input state.
#
# Round-trip:
#   1. Run surface_fluxes(T_int at z, T_sfc) -> get u_star, theta_star, L_MO
#   2. Recover U(z) = compute_profile_value(L_MO, z0m, Δz, u_star, 0, Momentum)
#   3. Verify U(z) ≈ |u_int - u_sfc|
#   4. Recover θ(z) = compute_profile_value(L_MO, z0h, Δz, θ_star, T_sfc, Heat)
#   5. Verify θ(z) ≈ T_int (approximately, modulo geopotential correction)

module TestProfileRecoveryRoundtrip

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaParams as CP

@testset "Profile Recovery Round-Trip" begin
    FT = Float64

    for uf_type in (UF.BusingerParams, UF.GryanikParams)
        @testset "$(nameof(uf_type))" begin
            param_set = SFP.SurfaceFluxesParameters(FT, uf_type)
            thermo_params = SFP.thermodynamics_params(param_set)
            κ = SFP.von_karman_const(param_set)

            z0m = FT(0.01)
            z0h = FT(0.001)
            Δz = FT(20)
            d = FT(0)

            config = SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(z0m, z0h),
                SF.ConstantGustinessSpec(FT(1)),
            )

            # Test across stability regimes
            test_cases = [
                (T_int = FT(295), T_sfc = FT(300), u = FT(10), label = "Unstable"),
                (T_int = FT(300), T_sfc = FT(300), u = FT(10), label = "Neutral"),
                (T_int = FT(305), T_sfc = FT(300), u = FT(10), label = "Stable"),
            ]

            for case in test_cases
                @testset "$(case.label)" begin
                    q_int = FT(0.01)
                    q_sfc = FT(0.012)
                    ρ_int = FT(1.1)
                    u_int = (case.u, FT(0))
                    u_sfc = (FT(0), FT(0))

                    opts = SF.SolverOptions{FT}(
                        maxiter = 50,
                        tol = FT(1e-10),
                        forced_fixed_iters = false,
                    )

                    result = SF.surface_fluxes(
                        param_set,
                        case.T_int, q_int, FT(0), FT(0), ρ_int,
                        case.T_sfc, q_sfc,
                        FT(0), Δz, d,
                        u_int, u_sfc,
                        nothing,
                        config,
                        SF.PointValueScheme(),
                        opts,
                    )

                    L_MO = result.L_MO
                    u_star = result.ustar
                    Δz_eff = Δz - d

                    # 1. Wind speed round-trip
                    # U(z) = (u_star / κ) * F_m(ζ, z0m)
                    # compute_profile_value returns: val_sfc + (scale/κ) * F
                    # For wind: val_sfc = 0 (surface wind = 0), scale = u_star
                    U_recovered = SF.compute_profile_value(
                        param_set,
                        L_MO,
                        z0m,
                        Δz_eff,
                        u_star,
                        FT(0), # val_sfc = |u_sfc| = 0
                        UF.MomentumTransport(),
                        UF.PointValueScheme(),
                    )

                    U_input = hypot(u_int[1] - u_sfc[1], u_int[2] - u_sfc[2])

                    # The effective wind speed used by the solver includes gustiness,
                    # so recovered U may differ if gustiness > |ΔU|.
                    # For these cases with u=10 m/s and gustiness=1, |ΔU| dominates.
                    @test isapprox(U_recovered, U_input; rtol = FT(0.02))

                    # 2. Temperature round-trip
                    # θ_star is defined via: SHF = -ρ_sfc * cp * u_star * θ_star
                    # (ignoring the VSE*E correction)
                    # compute_profile_value: T(z) = T_sfc + (θ_star / κ) * F_h
                    ρ_sfc = SF.surface_density(
                        param_set,
                        case.T_int,
                        ρ_int,
                        case.T_sfc,
                        Δz,
                        q_int,
                    )

                    # Compute theta_star from the result
                    inputs = SF.build_surface_flux_inputs(
                        case.T_int, q_int, FT(0), FT(0), ρ_int,
                        case.T_sfc, q_sfc,
                        FT(0), Δz, d,
                        u_int, u_sfc,
                        config,
                        nothing,
                        SF.FluxSpecs{FT}(),
                        nothing, nothing,
                    )

                    θ_star = SF.compute_theta_star(
                        param_set,
                        result.ζ,
                        z0h,
                        inputs,
                        SF.PointValueScheme(),
                        case.T_sfc,
                    )

                    # Recover DSE at z using similarity profile
                    # The profile function recovers θ(z) = θ_sfc + (θ*/κ) * F_h
                    # But the solver uses DSE = cp*T + gz, so we need to account for
                    # the geopotential difference.
                    grav = SFP.grav(param_set)
                    Φ_sfc = FT(0)
                    Φ_int = Φ_sfc + grav * Δz
                    cp = TD.cp_m(thermo_params, q_int, FT(0), FT(0))

                    DSE_sfc = TD.dry_static_energy(thermo_params, case.T_sfc, Φ_sfc)
                    DSE_int = TD.dry_static_energy(thermo_params, case.T_int, Φ_int)

                    # The similarity theory predicts: ΔDSE = (θ_star / κ) * F_h * cp
                    # (θ_star includes the 1/cp factor implicitly through the DSE formulation)
                    DSE_recovered = SF.compute_profile_value(
                        param_set,
                        L_MO,
                        z0h,
                        Δz_eff,
                        θ_star,
                        case.T_sfc,
                        UF.HeatTransport(),
                        UF.PointValueScheme(),
                    )

                    # DSE_recovered is θ(z) = T_sfc + (θ*/κ)*F_h
                    # The actual temperature is T_int, but θ_star uses DSE/cp scaling
                    # so T_int ≈ DSE_recovered + gz/cp (approximately)
                    T_recovered = DSE_recovered + grav * Δz / cp

                    # This should be close to T_int, within a few percent
                    # (tolerance accounts for the DSE vs T approximation)
                    @test isapprox(T_recovered, case.T_int; rtol = FT(0.02))
                end
            end
        end
    end
end

end # module
