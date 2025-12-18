using Test
using SurfaceFluxes
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaParams as CP

@testset "Physical Scales Consistency" begin
    FT = Float32
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    # Define range of conditions
    # Stable and Unstable
    # T_sfc, T_int pairs
    conditions = [
        # Unstable (T_sfc > T_int)
        (FT(300), FT(290), "Unstable"),
        # Stable (T_sfc < T_int)
        (FT(290), FT(300), "Stable"),
        # Neutral (T_sfc ≈ T_int) - might divide by zero in some implementations if not careful, but useful check
        (FT(300), FT(299.9), "Near Neutral"),
    ]

    # Fixed other parameters
    q_sfc = FT(0.02)
    q_int = FT(0.01) # q_sfc > q_int -> Evaporation -> Unstable for moisture
    u_int = (FT(10), FT(0))
    u_sfc = (FT(0), FT(0))
    z = FT(10)
    z0m = FT(0.1)  # Constant roughness
    z0h = FT(0.1)

    # Run test for each condition
    for (T_sfc, T_int, cond_name) in conditions
        @testset "$cond_name" begin
            # Create inputs
            ρ_int = FT(1.1)

            # Use API to compute fluxes
            config = SF.SurfaceFluxConfig(
                SF.ConstantRoughnessParams(z0m, z0h),
                SF.ConstantGustinessSpec(FT(0.0)), # No gustiness for simplicity
            )

            # Call surface_fluxes
            result = SF.surface_fluxes(
                param_set,
                T_int, q_int, ρ_int,
                T_sfc, q_sfc,
                FT(0), # rad level
                z,
                FT(0), # displacement
                u_int, u_sfc,
                nothing, # roughness_inputs
                config,
                SF.PointValueScheme(),
                SF.SolverOptions{FT}(),
                SF.FluxSpecs{FT}(),
            )

            # Extract basic scales
            ustar = result.ustar
            ζ = result.ζ

            inputs = SF.build_surface_flux_inputs(
                T_int, q_int, 0, 0, ρ_int, # q_liq, q_ice = 0
                T_sfc, q_sfc,
                FT(0), z, FT(0),
                u_int, u_sfc,
                config,
                nothing,
                SF.FluxSpecs{FT}(),
                nothing, nothing,
            )

            # Compute scales using new functions
            theta_star = SF.compute_theta_star(
                param_set,
                ζ,
                z0h,
                inputs,
                SF.PointValueScheme(),
                T_sfc,
            )

            q_star = SF.compute_q_star(
                param_set,
                ζ,
                z0h,
                inputs,
                SF.PointValueScheme(),
                q_sfc,
            )

            # --- Consistency Check ---

            # 1. Sensible Heat Flux
            # Compute correct ρ_sfc used by API
            ρ_sfc_calc = SF.surface_density(
                param_set, inputs.T_int, inputs.ρ_int, inputs.T_sfc_guess,
                q_int, FT(0), FT(0),
            )

            cp = TD.cp_m(thermo_params, q_int, FT(0), FT(0))

            # Diffusive SHF part
            shf_diffusive_scale = -ρ_sfc_calc * cp * ustar * theta_star

            # API SHF includes VSE * E term
            Φ_sfc = SF.surface_geopotential(inputs)
            VSE_sfc = TD.vapor_static_energy(thermo_params, T_sfc, Φ_sfc)
            shf_diffusive_api = result.shf - VSE_sfc * result.evaporation

            @test shf_diffusive_api ≈ shf_diffusive_scale rtol = 2e-2

            # 2. Evaporation / LHF
            evap_derived = -ρ_sfc_calc * ustar * q_star
            @test result.evaporation ≈ evap_derived rtol = 2e-2

            # 3. Geopotential Scale Consistency
            Φ_int = SF.interior_geopotential(param_set, inputs)
            Φ_sfc_check = SF.surface_geopotential(inputs)
            grav = SFP.grav(param_set)
            @test Φ_sfc_check == inputs.Φ_sfc
            @test Φ_int ≈ inputs.Φ_sfc + grav * inputs.Δz
        end
    end
end
