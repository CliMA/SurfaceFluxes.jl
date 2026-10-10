using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

# Use ClimaParams for parameters
import ClimaParams as CP

@testset "Displacement Height Usage" begin
    FT = Float64
    # Create parameter set using Businger parameters
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    # Define common atmospheric state
    T_int = FT(300)
    p_int = FT(100000)
    q_tot = FT(0.01)
    u_int = (FT(10), FT(0))
    u_sfc = (FT(0), FT(0))

    # Surface state
    T_sfc = FT(301)
    q_sfc = FT(0.012)
    Φ_sfc = FT(0)
    grav = SFP.grav(param_set)

    # Common roughness/configs
    z0m = FT(0.1)
    z0h = FT(0.1)
    roughness = SF.ConstantRoughnessParams(z0m, z0h)
    gustiness = SF.ConstantGustinessSpec(FT(0.0))
    config = SF.SurfaceFluxConfig(roughness, gustiness)

    # Use prescribed fluxes to isolate MOST height dependency from thermodynamic feeback
    flux_specs = SF.FluxSpecs(
        shf = FT(50.0),
        lhf = FT(100.0),
        ustar = FT(0.5),
    )

    # Case 1: Reference case (d=0, Δz=8)
    Δz_1 = FT(8)
    d_1 = FT(0)
    ρ_int_1 = TD.air_density(thermo_params, T_int, p_int)

    inputs_1 = SF.build_surface_flux_inputs(
        T_int, q_tot, 0, 0, ρ_int_1,
        T_sfc, q_sfc, Φ_sfc, Δz_1, d_1,
        u_int, u_sfc,
        config, nothing, flux_specs, nothing, nothing,
    )

    sf_1 = SF.surface_fluxes(param_set, inputs_1)

    # Case 2: Displaced case (d=2, Δz=10) -> Effective height should be 8. The surface
    # state applies at the displacement height, so the hydrostatic extrapolation of the
    # surface density and the geopotential difference also span the effective height,
    # and the two cases are the same exchange problem.
    Δz_2 = FT(10)
    d_2 = FT(2)

    inputs_2 = SF.build_surface_flux_inputs(
        T_int, q_tot, 0, 0, ρ_int_1,
        T_sfc, q_sfc, Φ_sfc, Δz_2, d_2,
        u_int, u_sfc,
        config, nothing, flux_specs, nothing, nothing,
    )

    sf_2 = SF.surface_fluxes(param_set, inputs_2)

    # Verification
    # 1. Stability parameter ζ should be identical because ζ = (z-d)/L and L is identical (same fluxes)
    @test sf_1.ζ ≈ sf_2.ζ atol = 1e-10

    # 2. Drag/Exchange coefficients should be identical
    @test sf_1.Cd ≈ sf_2.Cd atol = 1e-10
    @test sf_1.g_h ≈ sf_2.g_h atol = 1e-10

    # 3. L_MO should be identical
    @test sf_1.L_MO ≈ sf_2.L_MO atol = 1e-10

    # Test full solver (no prescribed fluxes): d > 0 gives different results than
    # ignoring d (Case 2 would then behave like z=10), and the same results as Case 1,
    # apart from the potential energy g d E carried by the vapor, which leaves the
    # surface state d higher in Case 2.
    inputs_3 = SF.build_surface_flux_inputs(
        T_int, q_tot, 0, 0, ρ_int_1,
        T_sfc, q_sfc, Φ_sfc, Δz_2, d_1, # d=0, z=10
        u_int, u_sfc,
        config, nothing, SF.FluxSpecs{FT}(), nothing, nothing,
    )
    sf_3 = SF.surface_fluxes(param_set, inputs_3)

    inputs_2_solve = SF.build_surface_flux_inputs(
        T_int, q_tot, 0, 0, ρ_int_1,
        T_sfc, q_sfc, Φ_sfc, Δz_2, d_2, # d=2, z=10
        u_int, u_sfc,
        config, nothing, SF.FluxSpecs{FT}(), nothing, nothing,
    )
    sf_2_solve = SF.surface_fluxes(param_set, inputs_2_solve)

    @test sf_2_solve.ζ ≉ sf_3.ζ
    @test sf_2_solve.Cd ≉ sf_3.Cd

    inputs_1_solve = SF.build_surface_flux_inputs(
        T_int, q_tot, 0, 0, ρ_int_1,
        T_sfc, q_sfc, Φ_sfc, Δz_1, d_1, # d=0, z=8
        u_int, u_sfc,
        config, nothing, SF.FluxSpecs{FT}(), nothing, nothing,
    )
    sf_1_solve = SF.surface_fluxes(param_set, inputs_1_solve)

    @test sf_2_solve.ζ ≈ sf_1_solve.ζ rtol = 1e-10
    @test sf_2_solve.Cd ≈ sf_1_solve.Cd rtol = 1e-10
    @test sf_2_solve.lhf ≈ sf_1_solve.lhf rtol = 1e-10
    @test sf_2_solve.shf - sf_1_solve.shf ≈ grav * d_2 * sf_1_solve.evaporation rtol = 1e-6

    # =========================================================================
    # Test: Vertical profile equivalence
    # Verify that dimensionless_profile(Δz, d=0) == dimensionless_profile(Δz-d, d>0)
    # when using effective height properly.
    # =========================================================================
    @testset "Dimensionless Profile Equivalence" begin
        uf_params = SFP.uf_params(param_set)

        # Test parameters
        Δz_ref = FT(8)    # Reference effective height
        d_test = FT(2)    # Displacement height
        Δz_total = Δz_ref + d_test  # Total geometric height = 10
        ζ_test = FT(-0.1)  # Stability parameter (unstable)
        z0_test = FT(0.1)  # Roughness length

        # Case A: No displacement (d=0), Δz = 8
        # Effective height = 8 - 0 = 8
        F_m_no_disp = UF.dimensionless_profile(
            uf_params,
            Δz_ref,
            ζ_test,
            z0_test,
            UF.MomentumTransport(),
        )
        F_h_no_disp =
            UF.dimensionless_profile(uf_params, Δz_ref, ζ_test, z0_test, UF.HeatTransport())

        # Case B: With displacement (d=2), Δz = 10
        # Effective height = 10 - 2 = 8 (same as Case A)
        # When calling directly, we pass the effective height, not the inputs struct
        Δz_eff_with_disp = Δz_total - d_test  # = 8
        F_m_with_disp = UF.dimensionless_profile(
            uf_params,
            Δz_eff_with_disp,
            ζ_test,
            z0_test,
            UF.MomentumTransport(),
        )
        F_h_with_disp = UF.dimensionless_profile(
            uf_params,
            Δz_eff_with_disp,
            ζ_test,
            z0_test,
            UF.HeatTransport(),
        )

        # The profiles should be exactly equal
        @test F_m_no_disp ≈ F_m_with_disp atol = eps(FT)
        @test F_h_no_disp ≈ F_h_with_disp atol = eps(FT)

        # Additional verification: profiles at different effective heights should differ
        F_m_different = UF.dimensionless_profile(
            uf_params,
            Δz_total,
            ζ_test,
            z0_test,
            UF.MomentumTransport(),
        )
        @test F_m_no_disp ≉ F_m_different  # 8m vs 10m should give different profiles

        # Also test finite-volume scheme
        F_m_fv_no_disp = UF.dimensionless_profile(
            uf_params,
            Δz_ref,
            ζ_test,
            z0_test,
            UF.MomentumTransport(),
            UF.LayerAverageScheme(),
        )
        F_m_fv_with_disp = UF.dimensionless_profile(
            uf_params,
            Δz_eff_with_disp,
            ζ_test,
            z0_test,
            UF.MomentumTransport(),
            UF.LayerAverageScheme(),
        )
        @test F_m_fv_no_disp ≈ F_m_fv_with_disp atol = eps(FT)
    end

end
