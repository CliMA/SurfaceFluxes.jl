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
    eval_callback(callback, default, ζ, param_set, thermo_params, inputs, FT)

Evaluate a surface state callback, returning `default` if callback is nothing or returns a non-Real value.
"""
@inline function eval_callback(callback, default, ζ, param_set, thermo_params, inputs, FT)
    if callback !== nothing
        val = callback(ζ, param_set, thermo_params, inputs)
        return val isa Real ? FT(val) : default
    else
        return default
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
    thermo_params = SFP.thermodynamics_params(param_set)

    # Surface state from guesses (callbacks not used for prescribed coefficients)
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
        error(
            "Both Cd and Ch must be provided for compute_fluxes_given_coefficients",
        )
    end
    Cd = Cd_in
    Ch = Ch_in

    # First pass: compute fluxes with zero buoyancy flux for gustiness
    FT = eltype(param_set)
    b_flux_init = FT(0)
    (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc, q_vap_sfc, ρ_sfc, b_flux_init,
    )

    # Now compute actual buoyancy flux from the fluxes
    b_flux = buoyancy_flux(
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

    # Second pass: recompute with correct buoyancy flux for gustiness
    (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc, q_vap_sfc, ρ_sfc, b_flux,
    )

    # Compute ustar from Cd: u*^2 = Cd * ΔU^2
    ΔU = windspeed(inputs, param_set, b_flux)
    ΔU_safe = max(ΔU, eps(FT))
    ustar = sqrt(Cd) * ΔU_safe

    # Derived L_MO and stability parameter
    L_MO = obukhov_length(param_set, ustar, b_flux)
    ζ = obukhov_stability_parameter(param_set, inputs.Δz, ustar, b_flux)

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
    thermo_params = SFP.thermodynamics_params(param_set)
    model = inputs.moisture_model

    T_sfc = inputs.T_sfc_guess
    q_vap_sfc = inputs.q_vap_sfc_guess
    ρ_sfc = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc,
        inputs.q_tot_int,
    )

    # Use prescribed flux values directly
    shf = inputs.shf
    lhf = inputs.lhf
    ustar = inputs.ustar

    # Compute E consistent with prescribed LHF
    LH_v0 = TP.LH_v0(thermo_params)
    E = lhf / LH_v0

    # Compute buoyancy flux from prescribed fluxes
    b_flux = buoyancy_flux(
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

    # Compute L_MO and stability parameter
    L_MO = obukhov_length(param_set, ustar, b_flux)
    ζ = obukhov_stability_parameter(param_set, inputs.Δz, ustar, b_flux)

    # Compute Coefficients with division-by-zero guard
    ΔU = windspeed(inputs, param_set, b_flux)
    ΔU_safe = max(ΔU, eps(FT))
    Cd = (ustar / ΔU_safe)^2

    # Compute roughness from ustar
    z0m, z0s = momentum_and_scalar_roughness(
        inputs.roughness_model,
        ustar,
        param_set,
        inputs.roughness_inputs,
    )
    z0h = z0s

    # Compute Ch
    Ch = heat_exchange_coefficient(param_set, ζ, z0m, z0h, inputs.Δz, scheme)

    # Compute momentum fluxes using Cd
    gustiness = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    (ρτxz, ρτyz) = momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)

    return SurfaceFluxConditions(
        shf, lhf, E,
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
@inline function compute_flux_components(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    Ch,
    Cd,
    Ts,
    qs,
    ρ_sfc,
    b_flux,
)
    thermo_params = SFP.thermodynamics_params(param_set)

    g_h = Ch * windspeed(inputs, param_set, b_flux)

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
    T_sfc_new = eval_callback(inputs.update_T_sfc, inputs.T_sfc_guess, ζ, param_set, thermo_params, inputs, FT)
    q_vap_sfc_new =
        eval_callback(inputs.update_q_vap_sfc, inputs.q_vap_sfc_guess, ζ, param_set, thermo_params, inputs, FT)

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

    # Auxiliary variables constant during solve
    uf_params = SFP.uf_params(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)

    root_function = ResidualFunction(
        param_set,
        inputs,
        scheme,
        uf_params,
        thermo_params,
    )

    # Determine the stability regime by evaluating the residual at neutral limit (ζ = 0)
    Ri_b_0 = root_function(FT(0))

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

    T_sfc_val = eval_callback(inputs.update_T_sfc, inputs.T_sfc_guess, ζ_final, param_set, thermo_params, inputs, FT)
    q_vap_sfc_val =
        eval_callback(inputs.update_q_vap_sfc, inputs.q_vap_sfc_guess, ζ_final, param_set, thermo_params, inputs, FT)

    # Update ρ_sfc based on final state
    ρ_sfc_val = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc_val,
        inputs.q_tot_int,
    )

    # Consistent gustiness/fluxes
    b_flux = buoyancy_flux(param_set, ζ_final, u_star_curr, inputs)

    # Use input coefficients if available, otherwise use MOST-derived ones
    Cd = inputs.Cd !== nothing ? inputs.Cd : drag_coefficient(param_set, ζ_final, z0m, inputs.Δz, scheme)
    Ch = inputs.Ch !== nothing ? inputs.Ch : heat_exchange_coefficient(param_set, ζ_final, z0m, z0h, inputs.Δz, scheme)

    (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc_val, q_vap_sfc_val, ρ_sfc_val, b_flux,
    )

    L_MO = obukhov_length(param_set, u_star_curr, b_flux)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        u_star_curr, ζ_final, Cd, Ch,
        L_MO,
    )
end

end # module
