using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

# Helper to compute density from state. 
# TODO Replace by Thermodynamics.jl function
function _density_from_state(thermo_params, T, pressure::FT, qt) where {FT}
    q_pt = TD.PhasePartition_equil_given_p(
        thermo_params,
        T,
        pressure,
        qt,
        TD.PhaseEquil,
    )
    return TD.air_density(thermo_params, T, pressure, q_pt)
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
    flux_specs = SF.FluxSpecs(FT)

    # Helper to construct inputs for testing
    function make_inputs(T_int, q_int, ρ_int, T_sfc, q_sfc, u_int)
        return SF.SurfaceFluxInputs(
            T_int,
            q_int,
            ρ_int,
            T_sfc,
            q_sfc,
            Φs,
            z_int,
            d,
            u_int,
            u_sfc,
            roughness,
            gustiness_spec,
            nothing, # roughness_inputs
            nothing, # update_Ts!
            nothing, # update_qs!
            flux_specs,
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
        q_int = FT(0.01)
        q_sfc = FT(0.01)
        u_int = (FT(5), FT(0))

        ρ_int = _density_from_state(thermo_params, T_int, P_sfc, q_int)
        ρ_sfc = _density_from_state(thermo_params, T_sfc, P_sfc, q_sfc)

        inputs = make_inputs(T_int, q_int, ρ_int, T_sfc, q_sfc, u_int)
        E = FT(0) # Zero evaporation for isolation

        shf = SF.sensible_heat_flux(
            param_set, thermo_params, inputs, g_h, T_int, T_sfc, ρ_sfc, E,
        )
        @test shf > 0

        # 2. Stable: T_sfc < T_int
        T_sfc_stable = FT(280)
        shf_stable = SF.sensible_heat_flux(
            param_set, thermo_params, inputs, g_h, T_int, T_sfc_stable, ρ_sfc, E,
        )
        @test shf_stable < 0
    end

    @testset "Latent Heat Flux (LHF) Signs" begin
        # Fixed Temp
        T_int = FT(300)
        T_sfc = FT(300)
        u_int = (FT(5), FT(0))

        # 1. Evaporation: q_sfc > q_int
        q_sfc = FT(0.02)
        q_int = FT(0.01)

        ρ_int = _density_from_state(thermo_params, T_int, P_sfc, q_int)
        ρ_sfc = _density_from_state(thermo_params, T_sfc, P_sfc, q_sfc)

        inputs = make_inputs(T_int, q_int, ρ_int, T_sfc, q_sfc, u_int)

        E = SF.evaporation(
            thermo_params, inputs, g_h, q_int, q_sfc, ρ_sfc,
        )
        @test E > 0

        lhf = SF.latent_heat_flux(thermo_params, inputs, E)
        @test lhf > 0

        # 2. Condensation: q_sfc < q_int
        q_sfc_cond = FT(0.005)
        # Note: using same inputs structure but passing different q values to formula 
        # is allowed by the function signature, but typically inputs.q_int matches q_int.
        # For strict correctness, let's remake inputs just in case function uses it.
        # Although `evaporation` function uses passed `q_vap_int`.

        E_cond = SF.evaporation(
            thermo_params, inputs, g_h, q_int, q_sfc_cond, ρ_sfc,
        )
        @test E_cond < 0

        lhf_cond = SF.latent_heat_flux(thermo_params, inputs, E_cond)
        @test lhf_cond < 0
    end

    @testset "Buoyancy Flux Signs" begin
        # Verify virtual potential temperature flux sign matches stability
        T_int = FT(290)
        T_sfc = FT(300) # Unstable
        q_int = FT(0.01)
        q_sfc = FT(0.01)
        u_int = (FT(5), FT(0))

        ρ_int = _density_from_state(thermo_params, T_int, P_sfc, q_int)
        ρ_sfc = _density_from_state(thermo_params, T_sfc, P_sfc, q_sfc)

        inputs = make_inputs(T_int, q_int, ρ_int, T_sfc, q_sfc, u_int)

        # Positive fluxes
        shf = FT(100)
        lhf = FT(0) # Ignore moisture contribution for simple check

        buoy_flux = SF.buoyancy_flux(
            param_set, thermo_params, shf, lhf, T_sfc, q_sfc, FT(0), FT(0), ρ_sfc,
        )
        @test buoy_flux > 0

        # Negative fluxes
        shf_neg = FT(-100)
        buoy_flux_neg = SF.buoyancy_flux(
            param_set, thermo_params, shf_neg, lhf, T_sfc, q_sfc, FT(0), FT(0), ρ_sfc,
        )
        @test buoy_flux_neg < 0
    end

    @testset "Momentum Flux Signs" begin
        T_int = FT(300)
        T_sfc = FT(300)
        q_int = FT(0.01)
        q_sfc = FT(0.01)

        # Wind +x
        u_int = (FT(10), FT(0))
        ρ_int = _density_from_state(thermo_params, T_int, P_sfc, q_int)
        ρ_sfc = _density_from_state(thermo_params, T_sfc, P_sfc, q_sfc)

        inputs = make_inputs(T_int, q_int, ρ_int, T_sfc, q_sfc, u_int)

        (ρτxz, ρτyz) = SF.momentum_fluxes(Cd, inputs, ρ_sfc, g_gustiness)

        # Surface stress opposes the flow difference (U_int - U_sfc)
        # Flow is +x. Stress should be -x.
        @test ρτxz < 0
        @test ρτyz == 0

        # Wind +y
        u_int_y = (FT(0), FT(10))
        inputs_y = make_inputs(T_int, q_int, ρ_int, T_sfc, q_sfc, u_int_y)
        (ρτxz_y, ρτyz_y) = SF.momentum_fluxes(Cd, inputs_y, ρ_sfc, g_gustiness)

        @test ρτxz_y == 0
        @test ρτyz_y < 0
    end
end
