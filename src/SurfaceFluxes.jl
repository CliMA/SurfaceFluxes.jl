"""
    SurfaceFluxes

Surface-layer flux calculations based on Monin-Obukhov similarity theory.
"""
module SurfaceFluxes

include("UniversalFunctions.jl")
include("Parameters.jl")

import Thermodynamics
const TD = Thermodynamics
const TP = Thermodynamics.Parameters

import RootSolvers
const RS = RootSolvers

import .UniversalFunctions
const UF = UniversalFunctions

import .Parameters
const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

const SolverScheme = UF.SolverScheme
const LayerAverageScheme = UF.LayerAverageScheme
const PointValueScheme = UF.PointValueScheme


include("types.jl")
include("roughness_lengths.jl")
include("input_builders.jl")
include("utilities.jl")
include("wind_and_gustiness.jl")
include("physical_scales.jl")
include("bulk_fluxes.jl")
include("exchange_coefficients.jl")
include("profile_recovery.jl")

@inline function default_solver_options(param_set::APS{FT}) where {FT}
    return SolverOptions(FT)
end

@inline function normalize_solver_options(param_set::APS{FT}, solver_opts) where {FT}
    if solver_opts === nothing
        return default_solver_options(param_set)
    elseif solver_opts isa Int
        return SolverOptions(FT; maxiter = solver_opts)
    elseif solver_opts isa SolverOptions{FT}
        return solver_opts
    elseif solver_opts isa SolverOptions
        return SolverOptions{FT}(
            solver_opts.tol,
            solver_opts.maxiter,
        )
    else
        throw(ArgumentError("Unsupported solver_opts specification: $(typeof(solver_opts))"))
    end
end

"""
    surface_fluxes(param_set, Tin, qin, ρin, T_sfc_guess, q_vap_sfc_guess, Φ_sfc, Δz, d,
                   u_int, u_sfc, config, scheme,
                   solver_opts, flux_specs, update_T_sfc, update_q_vap_sfc, q_liq_int=0, q_ice_int=0)

Compute near-surface fluxes using primitive inputs.
"""
function surface_fluxes(
    param_set::APS,
    T_int,
    q_tot_int,
    ρ_int,
    T_sfc_guess,
    q_vap_sfc_guess,
    Φ_sfc,
    Δz,
    d,
    u_int = nothing,
    u_sfc = nothing,
    roughness_inputs = nothing,
    config = nothing,
    scheme::SolverScheme = PointValueScheme(),
    solver_opts = nothing,
    flux_specs = nothing,
    update_T_sfc = nothing,
    update_q_vap_sfc = nothing,
    q_liq_int = 0,
    q_ice_int = 0,
)
    FT = eltype(param_set)

    # 1. Build inputs
    config_val = config === nothing ? default_surface_flux_config(FT) : config
    flux_specs_val = flux_specs === nothing ? FluxSpecs(FT) : flux_specs

    inputs = build_surface_flux_inputs(
        param_set,
        T_int,
        q_tot_int,
        q_liq_int,
        q_ice_int,
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        Φ_sfc,
        Δz,
        d,
        u_int,
        u_sfc,
        config_val,
        roughness_inputs,
        flux_specs_val,
        update_T_sfc,
        update_q_vap_sfc,
    )

    # Dispatching based on availability:
    # Case A: Coefficients known
    if inputs.Cd !== nothing && inputs.Ch !== nothing
        return compute_fluxes_given_coefficients(param_set, inputs, scheme)
    end

    # Case B: Fully Prescribed Fluxes (shf, lhf, ustar known)
    if inputs.shf !== nothing && inputs.lhf !== nothing && inputs.ustar !== nothing
        return compute_fluxes_from_prescribed(param_set, inputs, scheme)
    end

    # Case C: Standard MOST solve
    solver_opts_val = normalize_solver_options(param_set, solver_opts)
    return solve_monin_obukhov(param_set, inputs, scheme, solver_opts_val)
end

function default_surface_flux_config(::Type{FT}) where {FT}
    return SurfaceFluxConfig(
        ConstantRoughnessParams(FT(1e-3), FT(1e-3)),
        ConstantGustinessSpec(FT(1)),
    )
end

# ------------------------------------------------------------------------------
# Solver Logic
# ------------------------------------------------------------------------------

"""
    compute_fluxes_given_coefficients(param_set, inputs, scheme)

Computes fluxes when Cd and Ch are already known.
"""
function compute_fluxes_given_coefficients(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    scheme,
)

    FT = eltype(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)

    # If update functions exist, we assume they are not needed for prescribed coefficients.
    T_sfc = inputs.T_sfc_guess
    q_vap_sfc = inputs.q_vap_sfc_guess
    ρ_sfc = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc,
        inputs.q_tot_int,
    )

    # Coefficients
    Cd_in = inputs.Cd
    Ch_in = inputs.Ch
    if Cd_in === nothing || Ch_in === nothing
        # This function should only be called if both Cd and Ch are provided
        error(
            "Both Cd and Ch must be provided for compute_fluxes_given_coefficients",
        )
    end
    Cd = Cd_in
    Ch = Ch_in

    # We need gustiness to get ΔU.
    # For prescribed coefficients, we assume minimal or zero buoyancy flux for gustiness if not provided.
    buoyancy_flux_val = FT(0)

    # Fluxes
    (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc, q_vap_sfc, ρ_sfc, buoyancy_flux_val,
    )
    buoy_flux = buoyancy_flux(
        param_set,
        thermo_params,
        shf,
        lhf,
        T_sfc,
        q_vap_sfc,
        inputs.q_liq_int,
        inputs.q_ice_int,
        ρ_sfc,
        inputs.moisture_model,
    )

    # Helper to compute ustar from Cd
    # u*^2 = Cd * U^2
    ΔU = windspeed(inputs, param_set, buoyancy_flux_val)
    ustar = sqrt(Cd) * ΔU

    # Derived L_MO and stability parameter
    L_MO = obukhov_length(param_set, ustar, buoy_flux)
    ζ = obukhov_stability_parameter(param_set, inputs.Δz, ustar, buoy_flux)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        ustar, ζ, Cd, Ch,
        L_MO,
    )
end

"""
    compute_fluxes_from_prescribed(param_set, inputs, scheme)

Computes diagnostics when ustar, shf, lhf are all prescribed.
"""
function compute_fluxes_from_prescribed(param_set::APS, inputs::SurfaceFluxInputs, scheme)
    FT = eltype(param_set)
    # We assume state is known (or initial guess is sufficient).

    T_sfc = inputs.T_sfc_guess
    q_vap_sfc = inputs.q_vap_sfc_guess
    ρ_sfc = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc,
        inputs.q_tot_int,
    )

    thermo_params = SFP.thermodynamics_params(param_set)
    model = inputs.moisture_model

    # We need E for LHF calculation if LHF is prescribed but we need consistent E.
    # However, flux functions handle prescribed returning if set.
    g_h_dummy = FT(0)

    E = evaporation(thermo_params, inputs, g_h_dummy, inputs.q_tot_int, q_vap_sfc, ρ_sfc, model)
    lhf = latent_heat_flux(thermo_params, inputs, E, model)
    shf = sensible_heat_flux(param_set, thermo_params, inputs, g_h_dummy, inputs.T_int, T_sfc, ρ_sfc, E)

    # Compute buoyancy flux from prescribed fluxes
    buoy_flux = buoyancy_flux(
        param_set,
        thermo_params,
        shf,
        lhf,
        T_sfc,
        q_vap_sfc,
        inputs.q_liq_int,
        inputs.q_ice_int,
        ρ_sfc,
        model,
    )

    ustar = inputs.ustar

    # Compute L_MO and stability parameter
    L_MO = obukhov_length(param_set, ustar, buoy_flux)
    ζ = obukhov_stability_parameter(param_set, inputs.Δz, ustar, buoy_flux)

    # Compute Coefficients
    ΔU = windspeed(inputs, param_set, buoy_flux)
    Cd = (ustar / ΔU)^2

    # If ustar is prescribed, we compute roughness from ustar.
    z0m, z0s = momentum_and_scalar_roughness(inputs.roughness_model, ustar, param_set, inputs.roughness_inputs)
    z0h = z0s

    # Compute Ch
    Ch = heat_exchange_coefficient(param_set, ζ, z0m, z0h, inputs.Δz, scheme)

    (shf_out, lhf_out, E_out, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc, q_vap_sfc, ρ_sfc, buoy_flux,
    )

    return SurfaceFluxConditions(
        shf_out, lhf_out, E_out,
        ρτxz, ρτyz,
        ustar, ζ, Cd, Ch,
        L_MO,
    )
end

"""
    compute_flux_components(param_set, inputs, Ch, Cd, Ts, qs, ρ_sfc, b_flux)

Computes the individual flux components (sensible heat, latent heat, buoyancy, momentum)
given the exchange coefficients and surface state.
"""
function compute_flux_components(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    Ch,
    Cd,
    Ts,
    qs,
    ρ_sfc,
    b_flux,
)
    FT = eltype(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)

    g_h = Ch * windspeed(inputs, param_set, b_flux)
    g_q = g_h

    model = inputs.moisture_model
    E = evaporation(thermo_params, inputs, g_h, inputs.q_tot_int, qs, ρ_sfc, model)
    lhf = latent_heat_flux(thermo_params, inputs, E, model)
    shf = sensible_heat_flux(param_set, thermo_params, inputs, g_h, inputs.T_int, Ts, ρ_sfc, E)

    # Momentum fluxes
    gustiness = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    (ρτxz, ρτyz) = momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)

    return (shf, lhf, E, ρτxz, ρτyz)
end

struct ResidualFunction{PS, I, UF, TP, SCH} <: Function
    param_set::PS
    inputs::I
    scheme::SCH
    uf_params::UF
    thermo_params::TP
end

function (rf::ResidualFunction)(ζ)
    FT = eltype(rf.param_set)

    # Unpack parameters that do not change over iterations
    param_set = rf.param_set
    inputs = rf.inputs
    scheme = rf.scheme
    uf_params = rf.uf_params
    thermo_params = rf.thermo_params

    # 1. Compute u_star and roughness lengths, iteratively if they are mutually dependent
    u_star, z0m, z0s = compute_ustar_and_roughness(
        param_set,
        ζ,
        inputs,
        scheme,
    )
    z0h = z0s # Assuming scalar roughness = heat roughness

    # 2. Update T_sfc and q_vap_sfc via callbacks or use inputs
    if inputs.update_T_sfc !== nothing
        val = inputs.update_T_sfc(ζ, param_set, thermo_params, inputs)
        T_sfc_new = val isa FT ? val : inputs.T_sfc_guess
    else
        T_sfc_new = inputs.T_sfc_guess
    end

    if inputs.update_q_vap_sfc !== nothing
        val = inputs.update_q_vap_sfc(ζ, param_set, thermo_params, inputs)
        q_vap_sfc_new = val isa FT ? val : inputs.q_vap_sfc_guess
    else
        q_vap_sfc_new = inputs.q_vap_sfc_guess
    end

    # 3. Update density
    ρ_sfc = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc_new,
        inputs.q_tot_int,
    )

    # 4. Compute gustiness and ΔU
    # Use the buoyancy flux derived from the current ζ and ustar to calculate gustiness
    current_ΔU = windspeed(param_set, ζ, u_star, inputs)

    # 5. Calculate state bulk Richardson number
    Rib_state = state_bulk_richardson_number(
        param_set,
        thermo_params,
        inputs,
        T_sfc_new,
        q_vap_sfc_new,
        ρ_sfc,
        current_ΔU,
    )

    # 6. Evaluate residual
    Rib_theory = UF.bulk_richardson_number(uf_params, inputs.Δz, ζ, z0m, z0h, scheme)

    return Rib_theory - Rib_state
end

"""
    solve_monin_obukhov(param_set, inputs, scheme, options)

Solves the Monin-Obukhov Similarity Theory (MOST) equations for the surface fluxes.
Iterates to find the stability parameter `ζ` that satisfies the
surface layer profiles and surface balance equations.
"""
function solve_monin_obukhov(param_set::APS, inputs::SurfaceFluxInputs, scheme, options::SolverOptions)
    FT = eltype(param_set)

    T_sfc_val = inputs.T_sfc_guess
    q_vap_sfc_val = inputs.q_vap_sfc_guess
    ρ_sfc_val = inputs.ρ_int

    # Auxiliary variables constant during solve
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)
    grav = SFP.grav(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)

    root_function = ResidualFunction(
        param_set,
        inputs,
        scheme,
        uf_params,
        thermo_params,
    )

    # Solve
    Ri_b_0 = root_function(FT(0.001))

    x0, x1 = if Ri_b_0 < 0 # Stable
        (FT(0.01), FT(0.05))
    else # Unstable
        (FT(-0.01), FT(-0.05))
    end

    sol = RS.find_zero(
        root_function,
        RS.SecantMethod(x0, x1),
        RS.CompactSolution(),
        RS.SolutionTolerance(options.tol),
        options.maxiter,
    )

    # If failed, try other side?
    if !sol.converged
        x0_alt, x1_alt = (-x0, -x1)
        sol_alt = RS.find_zero(
            root_function,
            RS.SecantMethod(x0_alt, x1_alt),
            RS.CompactSolution(),
            RS.SolutionTolerance(options.tol),
            options.maxiter,
        )
        if sol_alt.converged
            sol = sol_alt
            ζ_final = sol.root
        else
            # Both failed.
            if Ri_b_0 < 0
                # Ri_b_state > 0 (Stable).
                # Non-convergence in stable regime often implies Ri_b > Ri_critical.
                # Collapse to laminar flow (large ζ).
                ζ_final = FT(100)
            else
                # Unstable non-convergence?
                # Use best guess from initial attempt
                ζ_final = sol.root
            end
        end
    else
        ζ_final = sol.root
    end

    # Finalize state 
    # 1. Compute u_star and roughness
    # Consistent with ζ_final
    u_star_curr, z0m, z0s = compute_ustar_and_roughness(
        param_set,
        ζ_final,
        inputs,
        scheme,
    )
    z0h = z0s

    if inputs.update_T_sfc !== nothing
        val = inputs.update_T_sfc(ζ_final, param_set, thermo_params, inputs)
        if val isa FT
            T_sfc_val = val
        end
    end
    if inputs.update_q_vap_sfc !== nothing
        val = inputs.update_q_vap_sfc(ζ_final, param_set, thermo_params, inputs)
        if val isa FT
            q_vap_sfc_val = val
        end
    end

    # Update ρ_sfc based on final state
    ρ_sfc_val = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc_val,
        inputs.q_tot_int,
    )

    # Consistent gustiness/fluxes
    b_flux_for_gustiness = buoyancy_flux(param_set, ζ_final, u_star_curr, inputs)

    # Use input coefficients if available, otherwise use MOST-derived ones
    Cd = inputs.Cd !== nothing ? inputs.Cd : drag_coefficient(param_set, ζ_final, z0m, inputs.Δz, scheme)
    Ch = inputs.Ch !== nothing ? inputs.Ch : heat_exchange_coefficient(param_set, ζ_final, z0m, z0h, inputs.Δz, scheme)

    (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc_val, q_vap_sfc_val, ρ_sfc_val, b_flux_for_gustiness,
    )

    L_MO = obukhov_length(param_set, u_star_curr, b_flux_for_gustiness)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        u_star_curr, ζ_final, Cd, Ch,
        L_MO,
    )
end

end # module
