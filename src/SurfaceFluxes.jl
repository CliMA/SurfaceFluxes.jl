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
    surface_tke,
    scalar_variance,
    theta_variance,
    obukhov_length,
    obukhov_stability_parameter

# From types.jl
export SurfaceFluxConditions,
    SurfaceFluxConfig,
    FluxSpecs,
    SolverOptions,
    ConstantRoughnessParams,
    COARE3RoughnessParams,
    RaupachRoughnessParams,
    ConstantGustinessSpec,
    DeardorffGustinessSpec,
    MoistModel,
    DryModel

# From roughness_sublayer.jl
export NoRoughnessSubLayer, PhysickGarrattRSL, HarmanFinniganRSL, rsl_profile_correction

# From utilities.jl
export surface_density

# From profile_recovery.jl
export compute_profile_value

# From UniversalFunctions.jl (solver schemes)
export PointValueScheme, LayerAverageScheme

include("types.jl")
include("roughness_sublayer.jl")
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
    SolverOptions{FT}(tol = FT(1e-2), maxiter = solver_opts)
@inline normalize_solver_options(
    param_set::APS{FT},
    solver_opts::SolverOptions{FT},
) where {FT} = solver_opts
@inline normalize_solver_options(
    param_set::APS{FT},
    solver_opts::SolverOptions,
) where {FT} =
    SolverOptions{FT}(
        tol = solver_opts.tol,
        rtol = solver_opts.rtol,
        maxiter = solver_opts.maxiter,
        forced_fixed_iters = solver_opts.forced_fixed_iters,
    )
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
    surface_fluxes(
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
    - `roughness`: Model for roughness lengths (e.g., `ConstantRoughnessParams`, `COARE3RoughnessParams`).
      Note: This package currently assumes the roughness length for heat (`z0h`) is equal to the
      roughness length for scalars (`z0s`).
    - `gustiness`: Model for gustiness (e.g., `ConstantGustinessSpec`).
    - `moisture_model`: `DryModel` or `MoistModel`.
- `scheme`: Discretization scheme (`PointValueScheme` or `LayerAverageScheme`).
- `solver_opts`: Options for the root solver (`maxiter`, `tol`, `rtol`, `forced_fixed_iters`).
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
- `ζ`: Stability parameter `(z-d)/L` [-].
- `Cd`: Drag coefficient [-].
- `g_h`: Heat conductance [m/s].
- `T_sfc`, `q_vap_sfc`: Final iterated surface temperature [K] and vapor specific humidity [kg/kg].
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
    inputs,
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
    # Generic fallback used when the caller omits `config`. The roughness lengths come from
    # the `ConstantRoughnessParams` keyword defaults (z0m = 2e-4 m, z0s = 2e-5 m), which is
    # the single source of truth for the default roughness. Real applications pass an explicit
    # `config` or load roughness lengths from ClimaParams.
    return SurfaceFluxConfig(
        ConstantRoughnessParams{FT}(),
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
    inputs,
    scheme,
)
    thermo_params = SFP.thermodynamics_params(param_set)

    # Surface state from guesses (callbacks not used for prescribed coefficients)
    # Don't use type annotations here to allow for Dual numbers during AD
    T_sfc =
        inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc =
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

    # Compute g_h
    g_h = Ch * ΔU_safe

    # Derived L_MO and stability parameter
    L_MO = obukhov_length(param_set, ustar, b_flux)
    Δz_eff = effective_height(inputs)
    ζ = obukhov_stability_parameter(param_set, Δz_eff, ustar, b_flux)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        ustar, ζ, Cd, g_h,
        T_sfc, q_vap_sfc,
        L_MO,
        true,
    )
end

"""
    compute_fluxes_from_prescribed(param_set, inputs, scheme)

Computes diagnostics when ustar, shf, lhf are all prescribed.
"""
function compute_fluxes_from_prescribed(param_set::APS, inputs, scheme)
    FT = eltype(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)
    model = inputs.moisture_model
    T_sfc =
        inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc =
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

    # Compute g_h
    g_h = heat_conductance(param_set, ζ, ustar, inputs, z0m, z0h, scheme)

    # Compute momentum fluxes using Cd
    gustiness = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    (ρτxz, ρτyz) = momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        ustar, ζ, Cd, g_h,
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
    inputs,
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

    # Compute g_h
    g_h = heat_conductance(param_set, ζ, ustar, inputs, z0m, z0h, scheme)

    # Compute momentum fluxes using Cd
    gustiness = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    (ρτxz, ρτyz) = momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        ustar, ζ, Cd, g_h,
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
    inputs,
    Ch,
    Cd,
    Ts,
    qs,
    ρ_sfc,
    b_flux,
)
    g_h = Ch * windspeed(inputs, param_set, b_flux)

    model = inputs.moisture_model
    q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
    E = evaporation(param_set, inputs, g_h, q_vap_int, qs, ρ_sfc, model)
    lhf = latent_heat_flux(param_set, inputs, E, model)
    shf = sensible_heat_flux(param_set, inputs, g_h, inputs.T_int, Ts, ρ_sfc, E)

    # Momentum fluxes
    gustiness = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    (ρτxz, ρτyz) = momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)

    return (shf, lhf, E, ρτxz, ρτyz)
end

"""
    bulk_richardson_number_rsl(uf_params, rsl_model, Δz_eff, ζ, z0m, z0h, scheme)

RSL-corrected bulk Richardson number used inside the stability solver residual.

Computes `ζ · F̂_h / F̂_m²` where `F̂ = F + P` includes the roughness sublayer
correction from [`rsl_profile_correction`](@ref). Reduces to the standard
[`UF.bulk_richardson_number`](@ref) when `rsl_model` is [`NoRoughnessSubLayer`](@ref).
"""
@inline function bulk_richardson_number_rsl(uf_params, rsl_model, Δz_eff, ζ, z0m, z0h, scheme)
    F_m = UF.dimensionless_profile(uf_params, Δz_eff, ζ, z0m, UF.MomentumTransport(), scheme)
    F_h = UF.dimensionless_profile(uf_params, Δz_eff, ζ, z0h, UF.HeatTransport(), scheme)
    P_m = rsl_profile_correction(rsl_model, Δz_eff, z0m, UF.MomentumTransport())
    P_h = rsl_profile_correction(rsl_model, Δz_eff, z0h, UF.HeatTransport())
    return ζ * (F_h + P_h) / (F_m + P_m)^2
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
    T_sfc_guess_safe =
        inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc_guess_safe =
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

    # 6. Evaluate residual (RSL-corrected theoretical Ri_b)
    Δz_eff = effective_height(inputs)
    Rib_theory = bulk_richardson_number_rsl(uf_params, inputs.rsl_model, Δz_eff, ζ, z0m, z0h, scheme)

    return Rib_theory - Rib_state
end

"""
    solve_stability_param(root_function, ζ_max, options)

GPU-friendly solver for the stability parameter ζ. It uses a fixed number
of residual evaluations (`options.maxiter` + 5) and no data-dependent control
flow when `options.forced_fixed_iters` is true: every point performs identical
work, so warps execute uniformly.

# Stage 1: branch detection and bracketing (5 evaluations)
The residual at neutral stability is `f(0) = -Ri_b_state` (the theoretical
bulk Richardson number `Ri_b(0)` vanishes), so its sign indicates the stability
branch on which a root is expected: stable (`ζ > 0`) for `f(0) <= 0`, unstable
(`ζ < 0`) for `f(0) > 0`.

For a fixed surface state, this indication is exact: `Ri_b(ζ)` is monotonic with
`sign(Ri_b) = sign(ζ)`, so the root lies on the branch matching the sign
of the state Richardson number. With surface-state callbacks
(`update_T_sfc`/`update_q_vap_sfc`), however, `Ri_b_state` itself varies with
ζ, and in transitional (near-neutral) conditions, its sign can differ between
the neutral evaluation and evaluations on the branch; therefore, the root can lie 
on the opposite branch from what `f(0)` indicates. To cover this, the residual is
probed at `ζ = ±1` on *both* branches, and at `|ζ| = 10, ζ_max` on the
indicated branch. The bracketing interval is selected from the sign changes in
priority order, innermost first: `(0, ±1)` on the indicated branch, `(0, ∓1)` 
on the opposite branch, then `(±1, ±10)` and `(±10, ±ζ_max)` on the indicated 
branch. (For a fixed surface state, the opposite-branch interval can never bracket, 
so this probe changes nothing in callback-free solves; an opposite-branch root 
at `|ζ| > 1` would require the callbacks to swing the state Richardson number 
by O(1) along the branch and is not searched for.)

If no probed interval brackets a sign change, no root exists within the
physical range (e.g., supercritical `Ri_b`, where the state is more stable
than the universal functions can support), and the solve saturates at the
indicated branch limit; this preserves the expected stability regime with
(near-)minimal fluxes.

# Stage 2: fixed-count refinement (`options.maxiter` evaluations)
The refinement is delegated to RootSolvers' `RegulaFalsiMethod` (safeguarded
regula falsi, Illinois variant), passing the stage-1 bracket with its
already-evaluated endpoint residuals so that no residual evaluation is
repeated. The bracket is preserved at every step, so iterates are bounded by
construction. In the saturated (no-root) case, the refinement runs on the
outermost same-sign interval (which RootSolvers rejects without iterating),
and its result is discarded by the final `ifelse`.

When `options.forced_fixed_iters` is true, `RootSolvers.NoTolerance` runs
exactly `maxiter` iterations with no data-dependent early exit. Otherwise,
`RootSolvers.RelativeOrAbsoluteSolutionTolerance(rtol, tol)` allows an early
exit once the step between iterates satisfies the tolerances. In both modes,
the returned root is sharpened by a final interpolation of the last bracket
(`RootSolvers.TwoPointSolution`) at no extra residual cost.

Returns `(ζ, converged)`. `converged` is `true` when a sign change was found
and either the solver's early-exit criterion fired (tolerance-checked mode)
or the final bracket width satisfies the tolerances. The saturated case
reports `false`. This flag is meaningful regardless of the `forced_fixed_iters` setting.
"""
function solve_stability_param(
    root_function::F,
    ζ_max::FT,
    options::SolverOptions,
) where {F, FT}
    r0 = root_function(zero(FT))

    # Indicated branch: stable (ζ > 0) iff Ri_b_state = -f(0) >= 0.
    # Exact for a fixed surface state; a heuristic when callbacks make
    # Ri_b_state vary with ζ (see the docstring).
    # `sgn` inherits the numeric type of the residual (e.g., Dual) so that
    # all iterates promote consistently under automatic differentiation.
    sgn = ifelse(r0 <= zero(r0), one(r0), -one(r0))

    # Near-neutral probes on both branches, log-spaced probes outward on the
    # indicated branch
    p1 = sgn
    m1 = -sgn
    p2 = FT(10) * sgn
    p3 = ζ_max * sgn
    r1 = root_function(p1)
    rm1 = root_function(m1)
    r2 = root_function(p2)
    r3 = root_function(p3)

    # Sign-change interval, selected innermost first; the opposite-branch
    # near-neutral interval (`cm`) covers transitional states whose callbacks
    # put the root on the branch opposite to the `f(0)` indication. If no
    # interval brackets (`bracketed == false`), the refinement below runs on
    # `(p2, p3)` and its result is discarded.
    c1 = r0 * r1 <= 0
    cm = r0 * rm1 <= 0
    c2 = r1 * r2 <= 0
    c3 = r2 * r3 <= 0
    bracketed = c1 | cm | c2 | c3

    a = ifelse(c1 | cm, zero(r0), ifelse(c2, p1, p2))
    fa = ifelse(c1 | cm, r0, ifelse(c2, r1, r2))
    b = ifelse(c1, p1, ifelse(cm, m1, ifelse(c2, p2, p3)))
    fb = ifelse(c1, r1, ifelse(cm, rm1, ifelse(c2, r2, r3)))

    # Stage 2: fixed-count safeguarded refinement, delegated to RootSolvers'
    # regula falsi (Illinois variant with a bisection fallback). The bracket
    # endpoints are passed pre-evaluated, so no residual evaluation is repeated
    # and the total count stays `5 + maxiter`. `TwoPointSolution` returns the
    # final bracket state, from which the convergence flag and a sharpened root
    # are computed.
    sol = if options.forced_fixed_iters
        # No data-dependent early exit: exactly `maxiter` iterations per point
        RS.find_zero(
            root_function,
            RS.RegulaFalsiMethod,
            a,
            b,
            fa,
            fb,
            RS.TwoPointSolution(),
            RS.NoTolerance(),
            options.maxiter,
        )
    else
        # Tolerance-checked (CPU) mode: early exit on the step between iterates
        RS.find_zero(
            root_function,
            RS.RegulaFalsiMethod,
            a,
            b,
            fa,
            fb,
            RS.TwoPointSolution(),
            RS.RelativeOrAbsoluteSolutionTolerance(options.rtol, options.tol),
            options.maxiter,
        )
    end

    # Final regula falsi interpolant of the last bracket: improves on the last
    # evaluated point at no extra residual cost (discarded, like the rest of
    # the refinement, when the interval does not bracket a root). The bracket
    # residuals may be Illinois-damped (halved), which preserves their signs
    # and hence the interpolant's validity.
    x = sol.root
    x_last = (sol.x0 * sol.y1 - sol.x1 * sol.y0) / (sol.y1 - sol.y0)
    lo = min(sol.x0, sol.x1)
    hi = max(sol.x0, sol.x1)
    use_last = isfinite(x_last) & (lo <= x_last) & (x_last <= hi)
    x = ifelse(use_last, x_last, x)

    # Converged when a sign change was found and either the solver's early-exit
    # criterion was met (tolerance-checked mode) or the final bracket width satisfies 
    # the tolerances.
    width = abs(sol.x1 - sol.x0)
    converged =
        bracketed &&
        (sol.converged || width < options.tol || width < options.rtol * abs(x))
    ζ = ifelse(bracketed, x, p3)
    return ζ, converged
end

"""
    solve_monin_obukhov(param_set, inputs, scheme, options)

Solves the Monin-Obukhov Similarity Theory (MOST) equations for the surface fluxes.
Iterates to find the stability parameter `ζ` that satisfies the
surface layer profiles and surface balance equations. Convergence is controlled
by `options.maxiter` and `options.tol`. If `options.forced_fixed_iters` is true,
ignores tolerance and iterates for exactly `maxiter`.

The ζ-iteration is performed by the internal `solve_stability_param`, a
fixed-evaluation-count and branchless bracketed solve suitable for GPU execution:
it detects the stability branch from the residual at neutral stability, brackets
the root with log-spaced probes within the physical range `|ζ| <= ζ_max = 100`, and refines
it with a fixed number of safeguarded regula falsi iterations. When no
root exists in the physical range (supercritical `Ri_b`), `ζ` saturates at the
limit of the appropriate stability branch (± ζ_max) and `converged = false` is reported.
"""
function solve_monin_obukhov(
    param_set::APS,
    inputs,
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

    # Physical limit for |ζ|. For supercritical Ri_b (e.g., very stable
    # stratification), no solution exists within this limit (e.g., for
    # Businger profiles, whose Ri_b(ζ) saturates at a critical value), and the
    # solve saturates at the limit of the appropriate stability branch.
    ζ_max = FT(100)

    ζ_final, converged = solve_stability_param(root_function, ζ_max, options)

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
    T_sfc_guess_safe =
        inputs.T_sfc_guess === nothing ? inputs.T_int : inputs.T_sfc_guess
    q_vap_sfc_guess_safe =
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

    # Use input coefficients if available, otherwise use MOST-derived ones (with RSL)
    Δz_eff = effective_height(inputs)
    ΔU = windspeed(inputs, param_set, b_flux)
    ΔU_safe = max(ΔU, eps(FT))
    Cd =
        inputs.Cd !== nothing ? inputs.Cd :
        inputs.ustar !== nothing ? (inputs.ustar / ΔU_safe)^2 :
        drag_coefficient(param_set, ζ_final, z0m, Δz_eff, scheme, inputs.rsl_model)
    Ch =
        inputs.Ch !== nothing ? inputs.Ch :
        heat_exchange_coefficient(param_set, ζ_final, z0m, z0h, Δz_eff, scheme, inputs.rsl_model)

    (shf, lhf, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, T_sfc_val, q_vap_sfc_val, ρ_sfc_val, b_flux,
    )

    g_h = Ch * ΔU_safe
    L_MO = obukhov_length(param_set, u_star_curr, b_flux)

    return SurfaceFluxConditions(
        shf, lhf, E,
        ρτxz, ρτyz,
        u_star_curr, ζ_final, Cd, g_h,
        T_sfc_val, q_vap_sfc_val,
        L_MO,
        converged,
    )
end

# ------------------------------------------------------------------------------
# Deprecations
# ------------------------------------------------------------------------------

# `u_variance` was renamed to `surface_tke` (it returns the surface-layer TKE, not the
# streamwise variance σ_u²). `Base.@deprecate` forwards `u_variance` to `surface_tke` and
# re-exports it, so existing `using SurfaceFluxes; u_variance(...)` keeps working (with a
# deprecation warning under `--depwarn=yes`).
Base.@deprecate u_variance surface_tke

end # module
