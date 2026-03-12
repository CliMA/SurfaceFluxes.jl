# Energy budget closure test.
#
# Verifies internal consistency of the returned flux quantities:
#   SHF = -ρ_sfc * g_h * ΔDSE + VSE_sfc * E
#   LHF = LH_v0 * E
#   E   = -ρ_sfc * g_h * Δq_vap
#
# Since ρ_sfc is not returned in the result, we derive it from the evaporation
# equation and verify consistency with SHF. This tests the relationship:
#   SHF = E * (ΔDSE / Δq + VSE_sfc)    (when Δq ≠ 0)

module TestEnergyBudgetClosure

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP

@testset "Energy Budget Closure" begin
    FT = Float64

    for uf_type in (UF.BusingerParams, UF.GryanikParams)
        @testset "$(nameof(uf_type))" begin
            param_set = SFP.SurfaceFluxesParameters(FT, uf_type)
            thermo_params = SFP.thermodynamics_params(param_set)

            test_cases = [
                (T_int = FT(295), T_sfc = FT(305), q_int = FT(0.008), q_sfc = FT(0.020),
                    label = "Unstable, moist"),
                (T_int = FT(305), T_sfc = FT(295), q_int = FT(0.012), q_sfc = FT(0.008),
                    label = "Stable, dry surface"),
                (T_int = FT(300), T_sfc = FT(300), q_int = FT(0.010), q_sfc = FT(0.015),
                    label = "Neutral, moist surface"),
                (T_int = FT(290), T_sfc = FT(300), q_int = FT(0.005), q_sfc = FT(0.005),
                    label = "Unstable, dry"),
            ]

            for case in test_cases
                @testset "$(case.label)" begin
                    ρ_int = FT(1.1)
                    u_int = (FT(8), FT(0))
                    u_sfc = (FT(0), FT(0))
                    Δz = FT(15)
                    z0m = FT(1e-3)
                    z0h = FT(1e-3)

                    config = SF.SurfaceFluxConfig(
                        SF.ConstantRoughnessParams(z0m, z0h),
                        SF.ConstantGustinessSpec(FT(1)),
                    )

                    opts = SF.SolverOptions{FT}(
                        maxiter = 50,
                        tol = FT(1e-8),
                        forced_fixed_iters = false,
                    )

                    result = SF.surface_fluxes(
                        param_set,
                        case.T_int, case.q_int, FT(0), FT(0), ρ_int,
                        case.T_sfc, case.q_sfc,
                        FT(0), Δz, FT(0),
                        u_int, u_sfc,
                        nothing,
                        config,
                        SF.PointValueScheme(),
                        opts,
                    )

                    # Compute DSE and VSE from first principles
                    grav = SFP.grav(param_set)
                    Φ_sfc = FT(0)
                    Φ_int = Φ_sfc + grav * Δz

                    DSE_int = TD.dry_static_energy(thermo_params, case.T_int, Φ_int)
                    DSE_sfc = TD.dry_static_energy(thermo_params, case.T_sfc, Φ_sfc)
                    ΔDSE = DSE_int - DSE_sfc

                    VSE_sfc = TD.vapor_static_energy(thermo_params, case.T_sfc, Φ_sfc)

                    E = result.evaporation
                    g_h = result.g_h
                    LH_v0 = TP.LH_v0(thermo_params)

                    # Test 1: LHF = LH_v0 * E
                    lhf_reconstructed = LH_v0 * E
                    @test isapprox(result.lhf, lhf_reconstructed; rtol = FT(1e-6))

                    # Test 2: SHF/E relationship (eliminates ρ_sfc)
                    # From SHF = -ρ_sfc * g_h * ΔDSE + VSE_sfc * E
                    # and  E   = -ρ_sfc * g_h * Δq_vap
                    # we get SHF = E * (ΔDSE / Δq + VSE_sfc)  when Δq ≠ 0
                    q_vap_int = case.q_int  # No condensate
                    Δq = q_vap_int - case.q_sfc

                    if abs(Δq) > eps(FT)
                        shf_from_ratio = E * (ΔDSE / Δq + VSE_sfc)
                        @test isapprox(result.shf, shf_from_ratio; rtol = FT(1e-6))
                    end

                    # Test 3: Derive ρ_sfc from E, then verify SHF
                    if abs(Δq) > eps(FT) && abs(g_h) > eps(FT)
                        ρ_sfc_derived = -E / (g_h * Δq)
                        shf_reconstructed = -ρ_sfc_derived * g_h * ΔDSE + VSE_sfc * E
                        @test isapprox(result.shf, shf_reconstructed; rtol = FT(1e-6))
                    end

                    # Test 4: SHF decomposition (diffusive + vapor transport)
                    if abs(Δq) > eps(FT) && abs(g_h) > eps(FT)
                        ρ_sfc_derived = -E / (g_h * Δq)
                        shf_diffusive = -ρ_sfc_derived * g_h * ΔDSE
                        shf_transport = VSE_sfc * E
                        @test isapprox(
                            result.shf,
                            shf_diffusive + shf_transport;
                            rtol = FT(1e-10),
                        )
                    end
                end
            end
        end
    end
end

end # module
