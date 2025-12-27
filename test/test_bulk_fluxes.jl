module TestBulkFluxes

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import ClimaParams as CP

# Helper to compute density from state using ideal gas law
function _density_from_state(thermo_params, T, pressure::FT, qt) where {FT}
    R_m = TD.gas_constant_air(thermo_params, qt, FT(0), FT(0))
    return pressure / (R_m * T)
end

@testset "Direct Bulk Flux Formulas" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    # Constants
    P_sfc = FT(101325) # Surface pressure [Pa]
    Φs = FT(0)         # Surface geopotential [m²/s²]
    z_int = FT(10)     # Reference height [m]
    d = FT(0)          # Displacement height [m]
    u_sfc = (FT(0), FT(0))

    # Dummy models for input construction
    roughness = SF.ConstantRoughnessParams(FT(1e-4), FT(1e-4))
    gustiness_spec = SF.ConstantGustinessSpec(FT(0))
    flux_specs = SF.FluxSpecs{FT}()

    # Helper to construct inputs for testing
    function make_inputs(T_int, q_tot_int, ρ_int, T_sfc, q_vap_sfc, u_int)
        config = SF.SurfaceFluxConfig(roughness, gustiness_spec)
        return SF.build_surface_flux_inputs(
            T_int,
            q_tot_int,
            FT(0), # q_liq_int
            FT(0), # q_ice_int
            ρ_int,
            T_sfc,
            q_vap_sfc,
            Φs,
            z_int, # Δz
            d,
            u_int,
            u_sfc,
            config,
            nothing, # roughness_inputs
            flux_specs,
            nothing, # update_T_sfc!
            nothing, # update_q_vap_sfc!
        )
    end

    # Common dummy coefficient values (must be positive)
    Cd = FT(0.0015)
    Ch = FT(0.0015)
    speed = FT(5) # Wind speed for conductance calc
    g_h = Ch * speed
    g_gustiness = FT(0) # No gustiness for these explicit checks

    @testset "Sensible Heat Flux (SHF) Signs" begin
        # 1. Unstable: T_sfc > T_int
        T_int = FT(290)
        T_sfc = FT(300)
        q_tot_int = FT(0.01)
        q_vap_sfc = FT(0.01)
        u_int = (FT(5), FT(0))

        ρ_int = _density_from_state(thermo_params, T_int, P_sfc, q_tot_int)
        ρ_sfc = _density_from_state(thermo_params, T_sfc, P_sfc, q_vap_sfc)

        inputs = make_inputs(T_int, q_tot_int, ρ_int, T_sfc, q_vap_sfc, u_int)
        E = FT(0) # Zero evaporation for isolation

        shf = SF.sensible_heat_flux(
            param_set, inputs, g_h, T_int, T_sfc, ρ_sfc, E,
        )
        @test shf > 0

        # 2. Stable: T_sfc < T_int
        T_sfc_stable = FT(280)
        shf_stable = SF.sensible_heat_flux(
            param_set, inputs, g_h, T_int, T_sfc_stable, ρ_sfc, E,
        )
        @test shf_stable < 0
    end

    @testset "Latent Heat Flux (LHF) Signs" begin
        # Fixed Temp
        T_int = FT(300)
        T_sfc = FT(300)
        u_int = (FT(5), FT(0))

        # 1. Evaporation: q_vap_sfc > q_tot_int
        q_vap_sfc = FT(0.02)
        q_tot_int = FT(0.01)

        ρ_int = _density_from_state(thermo_params, T_int, P_sfc, q_tot_int)
        ρ_sfc = _density_from_state(thermo_params, T_sfc, P_sfc, q_vap_sfc)

        inputs = make_inputs(T_int, q_tot_int, ρ_int, T_sfc, q_vap_sfc, u_int)

        E = SF.evaporation(
            param_set, inputs, g_h, q_tot_int, q_vap_sfc, ρ_sfc,
        )
        @test E > 0

        lhf = SF.latent_heat_flux(param_set, inputs, E)
        @test lhf > 0

        # 2. Condensation: q_sfc < q_tot_int
        q_sfc_cond = FT(0.005)

        E_cond = SF.evaporation(
            param_set, inputs, g_h, q_tot_int, q_sfc_cond, ρ_sfc,
        )
        @test E_cond < 0

        lhf_cond = SF.latent_heat_flux(param_set, inputs, E_cond)
        @test lhf_cond < 0
    end

    @testset "Buoyancy Flux Signs" begin
        # Verify virtual potential temperature flux sign matches stability
        T_int = FT(290)
        T_sfc = FT(300) # Unstable
        q_tot_int = FT(0.01)
        q_vap_sfc = FT(0.01)
        u_int = (FT(5), FT(0))

        ρ_int = _density_from_state(thermo_params, T_int, P_sfc, q_tot_int)
        ρ_sfc = _density_from_state(thermo_params, T_sfc, P_sfc, q_vap_sfc)

        inputs = make_inputs(T_int, q_tot_int, ρ_int, T_sfc, q_vap_sfc, u_int)

        # Positive fluxes
        shf = FT(100)
        lhf = FT(0) # Ignore moisture contribution for simple check

        buoy_flux = SF.buoyancy_flux(
            param_set, shf, lhf, T_sfc, ρ_sfc, q_vap_sfc,
        )
        @test buoy_flux > 0

        # Negative fluxes
        shf_neg = FT(-100)
        buoy_flux_neg = SF.buoyancy_flux(
            param_set, shf_neg, lhf, T_sfc, ρ_sfc, q_vap_sfc,
        )
        @test buoy_flux_neg < 0
    end

    @testset "Buoyancy Flux from Zeta" begin
        # Test consistency between B -> L -> ζ -> B
        B_true = FT(0.05)
        ustar = FT(0.3)
        κ = SFP.von_karman_const(param_set)
        L = SF.obukhov_length(param_set, ustar, B_true)
        ζ = z_int / L # Using z_int as Δz from global scope in test

        inputs = make_inputs(FT(300), FT(0.01), FT(1.2), FT(300), FT(0.01), (FT(10), FT(0)))

        B_calc = SF.buoyancy_flux(param_set, ζ, ustar, inputs)

        @test B_calc ≈ B_true

        # Test negative B (stable)
        B_stable = FT(-0.01)
        L_stable = SF.obukhov_length(param_set, ustar, B_stable)
        ζ_stable = z_int / L_stable
        B_calc_stable = SF.buoyancy_flux(param_set, ζ_stable, ustar, inputs)
        @test B_calc_stable ≈ B_stable
    end

    @testset "Momentum Flux Signs" begin
        T_int = FT(300)
        T_sfc = FT(300)
        q_tot_int = FT(0.01)
        q_vap_sfc = FT(0.01)

        # Wind +x
        u_int = (FT(10), FT(0))
        ρ_int = _density_from_state(thermo_params, T_int, P_sfc, q_tot_int)
        ρ_sfc = _density_from_state(thermo_params, T_sfc, P_sfc, q_vap_sfc)

        inputs = make_inputs(T_int, q_tot_int, ρ_int, T_sfc, q_vap_sfc, u_int)

        (ρτxz, ρτyz) = SF.momentum_fluxes(Cd, inputs, ρ_sfc, g_gustiness)

        # Surface stress opposes the flow difference (U_int - U_sfc)
        # Flow is +x. Stress should be -x.
        @test sign(ρτxz) == -sign(u_int[1] - u_sfc[1])
        @test ρτyz == 0

        # Wind +y
        u_int_y = (FT(0), FT(10))
        inputs_y = make_inputs(T_int, q_tot_int, ρ_int, T_sfc, q_vap_sfc, u_int_y)
        (ρτxz_y, ρτyz_y) = SF.momentum_fluxes(Cd, inputs_y, ρ_sfc, g_gustiness)

        @test ρτxz_y == 0
        @test sign(ρτyz_y) == -sign(u_int_y[2] - u_sfc[2])
    end

    @testset "Dry bulk flux defaults (omitted moisture args)" begin
        # Setup local vars
        T_sfc = FT(300)
        ρ_sfc = FT(1.2)
        shf = FT(50)
        lhf = FT(0)
        ΔU = FT(10)
        T_int = FT(290)
        ρ_int = FT(1.1)
        u_int = (FT(10), FT(0))

        # Verify that omitting moisture arguments defaults to 0 (dry)
        # 1. buoyancy_flux(param_set, thermo_params, shf, lhf, T_sfc, ρ_sfc, [q_vap=0, q_liq=0, q_ice=0])
        # We pass explicit 0s for comparison
        b_flux_explicit = SF.buoyancy_flux(
            param_set, shf, lhf, T_sfc, ρ_sfc, FT(0), FT(0), FT(0), SF.MoistModel(),
        )
        b_flux_implicit = SF.buoyancy_flux(
            param_set, shf, lhf, T_sfc, ρ_sfc,
        )
        @test b_flux_implicit ≈ b_flux_explicit

        # 2. state_bulk_richardson_number(..., T_sfc, ρ_sfc, ΔU, [q_vap_sfc=0])
        inputs_dry = make_inputs(T_int, FT(0), ρ_int, T_sfc, FT(0), u_int)

        Ri_explicit = SF.state_bulk_richardson_number(
            param_set, inputs_dry, T_sfc, ρ_sfc, ΔU, FT(0),
        )
        Ri_implicit = SF.state_bulk_richardson_number(
            param_set, inputs_dry, T_sfc, ρ_sfc, ΔU,
        )
        @test Ri_implicit ≈ Ri_explicit
    end
end

end # module
