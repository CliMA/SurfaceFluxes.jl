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

import RootSolvers as RS

import .UniversalFunctions
const UF = UniversalFunctions

import .Parameters
const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

const SolverScheme = UF.SolverScheme
const LayerAverageScheme = UF.LayerAverageScheme
const PointValueScheme = UF.PointValueScheme

# Top-level API
export surface_fluxes

# From bulk_fluxes.jl
export sensible_heat_flux,
    evaporation,
    latent_heat_flux,
    buoyancy_flux,
    momentum_fluxes,
    state_bulk_richardson_number

# From exchange_coefficients.jl
export drag_coefficient, heat_exchange_coefficient, heat_conductance

# From physical_scales.jl
export compute_physical_scale_coeff,
    compute_ustar,
    compute_theta_star,
    compute_q_star,
    u_variance,
    scalar_variance,
    theta_variance,
    obukhov_length,
    obukhov_stability_parameter

# From types.jl
export SurfaceFluxConditions,
    SurfaceFluxConfig,
    FluxSpecs,
    SurfaceFluxInputs,
    SolverOptions,
    ConstantRoughnessParams,
    COARE3RoughnessParams,
    RaupachRoughnessParams,
    ConstantGustinessSpec,
    DeardorffGustinessSpec,
    MoistModel,
    DryModel

# From utilities.jl
export surface_density

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
    return SolverOptions{FT}()
end

@inline normalize_solver_options(param_set::APS{FT}, ::Nothing) where {FT} =
    default_solver_options(param_set)
@inline normalize_solver_options(param_set::APS{FT}, solver_opts::Int) where {FT} =
    SolverOptions{FT}(FT(1e-2), solver_opts, true)
@inline normalize_solver_options(
    param_set::APS{FT},
    solver_opts::SolverOptions{FT},
) where {FT} = solver_opts
@inline normalize_solver_options(
    param_set::APS{FT},
    solver_opts::SolverOptions,
) where {FT} =
    SolverOptions{FT}(solver_opts.tol, solver_opts.maxiter, solver_opts.forced_fixed_iters)
@inline normalize_solver_options(param_set::APS{FT}, solver_opts) where {FT} =
    default_solver_options(param_set)  # Fallback: use defaults for unsupported types

"""
    eval_callback(callback, default, ζ, param_set, thermo_params, inputs)

Evaluate a surface state callback, returning `default` if callback is nothing or returns a non-Real value.
"""
@inline function eval_callback(callback, default, args...)
    if callback !== nothing
        val = callback(args...)
        return val isa Real ? val : default
    else
        return default
    end
end

"""
    surface_fluxes(param_set, T_int, q_tot_int, q_liq_int, q_ice_int, ρ_int, T_sfc_guess, q_vap_sfc_guess, Φ_sfc, Δz, d, u_int, u_sfc, roughness_inputs=nothing, config=default, scheme=PointValueScheme, solver_opts=nothing, flux_specs=nothing)

Core entry point for calculating surface fluxes using Monin-Obukhov Similarity Theory (MOST).

# Functionality
Calculates sensible heat flux, latent heat flux, momentum flux (stress), and friction velocity.

Can operate in four modes depending on inputs:
1. **Prescribed Coefficients**: If `Cd` and `Ch` are provided in `flux_specs`, fluxes are computed directly.
2. **Fully Prescribed Fluxes**: If `shf`, `lhf`, and `ustar` are provided, they are validated and the fluxes are returned.
3. **Prescribed Heat and Drag**: If `shf`, `lhf`, and `Cd` are provided, `ustar` is derived from `Cd` and wind speed.
4. **Iterative Solver**: Otherwise, iterates to find the Obukhov stability parameter `ζ`. Optional functions can be provided to calculate the skin temperature and skin humidity during the iteration.

# Arguments
- `param_set`: SurfaceFluxes parameters (containing thermodynamics and universal function params).
- `T_int`: Interior (air) temperature [K] at height `z`.
- `q_tot_int`: Interior total specific humidity [kg/kg].
- `q_liq_int`, `q_ice_int`: Interior liquid/ice specific humidity [kg/kg].
- `ρ_int`: Interior air density [kg/m^3].
- `T_sfc_guess`: Initial guess for surface temperature [K], updated via callback if provided.
- `q_vap_sfc_guess`: Initial guess for surface vapor specific humidity [kg/kg], updated via callback if provided.
- `Φ_sfc`: Surface geopotential [m^2/s^2].
- `Δz`: Geometric height difference between the surface and the interior level [m], used for geopotential.
- `d`: Displacement height [m]. Aerodynamic calculations (MOST) use effective height `Δz - d`.
- `u_int`: Tuple of interior wind components `(u, v)` [m/s].
- `u_sfc`: Tuple of surface wind components `(u, v)` [m/s]. (Usually `(0, 0)`).
- `roughness_inputs`: Optional container of parameters (e.g., LAI, canopy height) that are passed
  directly to the specific roughness model (e.g., `RaupachRoughnessParams`).
- `config`: [`SurfaceFluxConfig`](@ref) struct containing:
    - `roughness`: Model for roughness lengths (e.g., `ConstantRoughnessParams`, `COARE3RoughnessSpec`).
      Note: This package currently assumes the roughness length for heat (`z0h`) is equal to the 
      roughness length for scalars (`z0s`).
    - `gustiness`: Model for gustiness (e.g., `ConstantGustinessSpec`).
    - `moisture_model`: `DryModel` or `WetModel`.
- `scheme`: Discretization scheme (`PointValueScheme` or `LayerAverageScheme`).
- `solver_opts`: Options for the root solver (`maxiter`, `tol`, `forced_fixed_iters`).
- `flux_specs`: Optional `FluxSpecs` to prescribe specific constraints (e.g., `ustar`, `shf`, `Cd`).
- `update_T_sfc`: Optional callback `f(T_sfc)` to update surface temperature during iteration.
- `update_q_vap_sfc`: Optional callback `f(q_vap)` to update surface humidity during iteration.

# Returns
A [`SurfaceFluxConditions`](@ref) struct containing:
- `shf`: Sensible Heat Flux [W/m^2].
- `lhf`: Latent Heat Flux [W/m^2].
- `evaporation`: Evaporation rate [kg/m^2/s].
- `ustar`: Friction velocity [m/s].
- `ρτxz`, `ρτyz`: Momentum flux components (stress) [N/m^2].
- `ζ`: Stability parameter (`(z-d)/L`).
- `Cd`, `Ch`: Drag and heat exchange coefficients.
- `T_sfc`, `q_vap_sfc`: Surface temperature [K] and vapor specific humidity [kg/kg].
- `L_MO`: Monin-Obukhov length [m].
- `converged`: Convergence status.
"""
function surface_fluxes(
    param_set::APS,
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
    u_int = (0, 0),
    u_sfc = (0, 0),
    roughness_inputs = nothing,
    config = default_surface_flux_config(eltype(param_set)),
    scheme::SolverScheme = PointValueScheme(),
    solver_opts = nothing,
    flux_specs = nothing,
    update_T_sfc = nothing,
    update_q_vap_sfc = nothing,
)
    FT = eltype(param_set)

    # Build inputs: Use explicit positional constructor for 
    # GPU compatibility (avoids kwargs)
    flux_specs_val = flux_specs === nothing ? FluxSpecs{FT}() : flux_specs

    inputs = build_surface_flux_inputs(
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
        config,
        roughness_inputs,
        flux_specs_val,
        update_T_sfc,
        update_q_vap_sfc,
    )

    return surface_fluxes(param_set, inputs, scheme, solver_opts)
end

"""
    surface_fluxes(param_set, inputs, scheme=PointValueScheme(), solver_opts=nothing)

Dispatch to the appropriate solver mode based on the availability of inputs (coefficients, fluxes, or state).
"""
function surface_fluxes(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    scheme::SolverScheme = PointValueScheme(),
    solver_opts::Union{SolverOptions, Nothing} = nothing,
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

    # Case C: Prescribed Heat Fluxes and Drag Coefficient
    if inputs.shf !== nothing && inputs.lhf !== nothing && inputs.Cd !== nothing
        return compute_fluxes_with_prescribed_heat_and_drag(param_set, inputs, scheme)
    end

    # Case D: Standard MOST solve
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
    # Use type of T_int to allow for ForwardDiff.Dual
    T_sfc_type = typeof(inputs.T_int)
    q_vap_sfc_type = typeof(inputs.q_tot_int)

    T_sfc::T_sfc_type =
        inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc::q_vap_sfc_type =
        inputs.q_vap_sfc_guess === nothing ? inputs.q_tot_int :
        inputs.q_vap_sfc_guess
    ρ_sfc = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc,
        inputs.Δz,
        inputs.q_tot_int,
        inputs.q_liq_int,
        inputs.q_ice_int,
        q_vap_sfc,
    )

    # Coefficients (caller must ensure both are provided)
    FT = eltype(param_set)
    Cd = FT(inputs.Cd)
    Ch = FT(inputs.Ch)

    # First pass: compute fluxes with zero buoyancy flux for gustiness
    b_flux_init = FT(0)
    (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc, q_vap_sfc, ρ_sfc, b_flux_init,
    )

    # Now compute actual buoyancy flux from the fluxes. (This and the following 
    # second pass computation could be skipped if gustiness does not depend on 
    # buoyancy flux.)
    b_flux = buoyancy_flux(
        param_set,
        shf,
        lhf,
        T_sfc,
        ρ_sfc,
        q_vap_sfc,
        inputs.q_liq_int,
        inputs.q_ice_int,
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
    Δz_eff = effective_height(inputs)
    ζ = obukhov_stability_parameter(param_set, Δz_eff, ustar, b_flux)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        ustar, ζ, Cd, Ch,
        T_sfc, q_vap_sfc,
        L_MO,
        true,
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
    # Use type of T_int to allow for ForwardDiff.Dual
    T_sfc_type = typeof(inputs.T_int)
    q_vap_sfc_type = typeof(inputs.q_tot_int)

    T_sfc::T_sfc_type =
        inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc::q_vap_sfc_type =
        inputs.q_vap_sfc_guess === nothing ? inputs.q_tot_int :
        inputs.q_vap_sfc_guess
    ρ_sfc = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc,
        inputs.Δz,
        inputs.q_tot_int,
        inputs.q_liq_int,
        inputs.q_ice_int,
        q_vap_sfc,
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
        shf,
        lhf,
        T_sfc,
        ρ_sfc,
        q_vap_sfc,
        inputs.q_liq_int,
        inputs.q_ice_int,
        model,
    )

    # Compute L_MO and stability parameter
    L_MO = obukhov_length(param_set, ustar, b_flux)
    Δz_eff = effective_height(inputs)
    ζ = obukhov_stability_parameter(param_set, Δz_eff, ustar, b_flux)

    # Compute Coefficients with division-by-zero guard
    ΔU = windspeed(inputs, param_set, b_flux)
    ΔU_safe = max(ΔU, eps(FT))
    Cd = (ustar / ΔU_safe)^2

    # Compute roughness from ustar
    # Note: We assume z0h = z0s (scalar roughness) for now
    z0m, z0h = momentum_and_scalar_roughness(
        inputs.roughness_model,
        ustar,
        param_set,
        inputs.roughness_inputs,
    )

    # Compute Ch
    Ch = heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, scheme)

    # Compute momentum fluxes using Cd
    gustiness = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    (ρτxz, ρτyz) = momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        ustar, ζ, Cd, Ch,
        T_sfc, q_vap_sfc,
        L_MO,
        true,
    )
end

"""
    compute_fluxes_with_prescribed_heat_and_drag(param_set, inputs, scheme)

Computes diagnostics when `shf`, `lhf`, and `Cd` are prescribed.
Returns a [`SurfaceFluxConditions`](@ref) struct.
"""
function compute_fluxes_with_prescribed_heat_and_drag(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    scheme,
)
    FT = eltype(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)
    model = inputs.moisture_model

    T_sfc = inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc =
        inputs.q_vap_sfc_guess === nothing ? inputs.q_tot_int : inputs.q_vap_sfc_guess
    ρ_sfc = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc,
        inputs.Δz,
        inputs.q_tot_int,
        inputs.q_liq_int,
        inputs.q_ice_int,
        q_vap_sfc,
    )

    # Use prescribed values
    shf = inputs.shf
    lhf = inputs.lhf
    Cd = inputs.Cd

    # Compute E consistent with prescribed LHF
    LH_v0 = TP.LH_v0(thermo_params)
    E = lhf / LH_v0

    # Compute buoyancy flux from prescribed fluxes
    b_flux = buoyancy_flux(
        param_set,
        shf,
        lhf,
        T_sfc,
        ρ_sfc,
        q_vap_sfc,
        inputs.q_liq_int,
        inputs.q_ice_int,
        model,
    )

    # Compute wind speed (accounting for gustiness via b_flux)
    ΔU = windspeed(inputs, param_set, b_flux)
    ΔU_safe = max(ΔU, eps(FT))

    # Compute ustar from Cd and wind speed
    # τ = ρ Cd ΔU^2 = ρ ustar^2  => ustar = sqrt(Cd) * ΔU
    ustar = sqrt(Cd) * ΔU_safe

    # Compute L_MO and stability parameter
    L_MO = obukhov_length(param_set, ustar, b_flux)
    Δz_eff = effective_height(inputs)
    ζ = obukhov_stability_parameter(param_set, Δz_eff, ustar, b_flux)

    # Compute roughness from ustar
    # Note: We assume z0h = z0s (scalar roughness) for now
    z0m, z0h = momentum_and_scalar_roughness(
        inputs.roughness_model,
        ustar,
        param_set,
        inputs.roughness_inputs,
    )

    # Compute Ch
    Ch = heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, scheme)

    # Compute momentum fluxes using Cd
    gustiness = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    (ρτxz, ρτyz) = momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        ustar, ζ, Cd, Ch,
        T_sfc, q_vap_sfc,
        L_MO,
        true,
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
    g_h = Ch * windspeed(inputs, param_set, b_flux)

    model = inputs.moisture_model
    E = evaporation(param_set, inputs, g_h, inputs.q_tot_int, qs, ρ_sfc, model)
    lhf = latent_heat_flux(param_set, inputs, E, model)
    shf = sensible_heat_flux(param_set, inputs, g_h, inputs.T_int, Ts, ρ_sfc, E)

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
    u_star, z0m, z0h = compute_ustar_and_roughness(
        param_set,
        ζ,
        inputs,
        scheme,
    )

    # Ensure type stability for default values (strip Union{Nothing, FT})
    # If guess is nothing, use interior values as safe dummy defaults
    # Use type of T_int to allow for Dual numbers during AD
    T_sfc_type = typeof(inputs.T_int)
    q_vap_sfc_type = typeof(inputs.q_tot_int)

    T_sfc_guess_safe::T_sfc_type =
        inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc_guess_safe::q_vap_sfc_type =
        inputs.q_vap_sfc_guess === nothing ? inputs.q_tot_int :
        inputs.q_vap_sfc_guess

    # 2. Update T_sfc and q_vap_sfc via callbacks or use inputs
    T_sfc_new = eval_callback(
        inputs.update_T_sfc,
        T_sfc_guess_safe,
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star,
        z0m,
        z0h,
    )
    q_vap_sfc_new = eval_callback(
        inputs.update_q_vap_sfc,
        q_vap_sfc_guess_safe,
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        T_sfc_new,
        u_star,
        z0m,
        z0h,
    )

    # 3. Update density
    ρ_sfc = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc_new,
        inputs.Δz,
        inputs.q_tot_int,
        inputs.q_liq_int,
        inputs.q_ice_int,
        q_vap_sfc_new,
    )

    # 4. Compute gustiness and ΔU
    # Use the buoyancy flux derived from the current ζ and ustar to calculate gustiness
    current_ΔU = windspeed(param_set, ζ, u_star, inputs)

    # 5. Calculate state bulk Richardson number
    Rib_state = state_bulk_richardson_number(
        param_set,
        inputs,
        T_sfc_new,
        ρ_sfc,
        current_ΔU,
        q_vap_sfc_new,
    )

    # 6. Evaluate residual
    Δz_eff = effective_height(inputs)
    Rib_theory = UF.bulk_richardson_number(uf_params, Δz_eff, ζ, z0m, z0h, scheme)

    return Rib_theory - Rib_state
end

"""
    solve_monin_obukhov(param_set, inputs, scheme, options)

Solves the Monin-Obukhov Similarity Theory (MOST) equations for the surface fluxes.
Iterates to find the stability parameter `ζ` that satisfies the
surface layer profiles and surface balance equations. Convergence is controlled
by `options.maxiter` and `options.tol`. If `options.forced_fixed_iters` is true, 
ignores tolerance and iterates for exactly `maxiter`.
"""
function solve_monin_obukhov(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    scheme,
    options::SolverOptions,
)
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

    # Define physical limits for the stability parameter ζ
    ζ_min = FT(-100)
    ζ_max = FT(100)

    sol = RS.find_zero(
        root_function,
        RS.BrentsMethod(ζ_min, ζ_max),
        RS.CompactSolution(),
        RS.SolutionTolerance(options.forced_fixed_iters ? FT(0) : options.tol),
        options.maxiter,
    )
    ζ_final = sol.root
    converged = sol.converged

    # Finalize state 
    # 1. Compute u_star and roughness
    # Consistent with ζ_final
    u_star_curr, z0m, z0h = compute_ustar_and_roughness(
        param_set,
        ζ_final,
        inputs,
        scheme,
    )

    # Ensure type stability for default values (strip Union{Nothing, FT})
    # Use type of T_int to allow for Dual numbers during AD
    T_sfc_type = typeof(inputs.T_int)
    q_vap_sfc_type = typeof(inputs.q_tot_int)

    T_sfc_guess_safe::T_sfc_type =
        inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc_guess_safe::q_vap_sfc_type =
        inputs.q_vap_sfc_guess === nothing ? inputs.q_tot_int :
        inputs.q_vap_sfc_guess

    T_sfc_val = eval_callback(
        inputs.update_T_sfc,
        T_sfc_guess_safe,
        ζ_final,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star_curr,
        z0m,
        z0h,
    )
    q_vap_sfc_val = eval_callback(
        inputs.update_q_vap_sfc,
        q_vap_sfc_guess_safe,
        ζ_final,
        param_set,
        thermo_params,
        inputs,
        scheme,
        T_sfc_val,
        u_star_curr,
        z0m,
        z0h,
    )

    # Update ρ_sfc based on final state
    ρ_sfc_val = surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc_val,
        inputs.Δz,
        inputs.q_tot_int,
        inputs.q_liq_int,
        inputs.q_ice_int,
        q_vap_sfc_val,
    )

    # Consistent gustiness/fluxes
    b_flux = buoyancy_flux(param_set, ζ_final, u_star_curr, inputs)

    # Use input coefficients if available, otherwise use MOST-derived ones
    Δz_eff = effective_height(inputs)
    Cd =
        inputs.Cd !== nothing ? inputs.Cd :
        drag_coefficient(param_set, ζ_final, z0m, Δz_eff, scheme)
    Ch =
        inputs.Ch !== nothing ? inputs.Ch :
        heat_exchange_coefficient(param_set, ζ_final, z0m, z0h, Δz_eff, scheme)

    (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc_val, q_vap_sfc_val, ρ_sfc_val, b_flux,
    )

    L_MO = obukhov_length(param_set, u_star_curr, b_flux)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        u_star_curr, ζ_final, Cd, Ch,
        T_sfc_val, q_vap_sfc_val,
        L_MO,
        converged,
    )
end

end # module
