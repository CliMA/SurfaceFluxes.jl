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

@inline float_type(::APS{FT}) where {FT} = FT

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
    surface_fluxes(param_set, Tin, qin, ρin, Ts_guess, qs_guess, Φs, Δz, d,
                   u_int, u_sfc, config, scheme,
                   solver_opts, flux_specs, update_Ts!, update_qs!)

Compute near-surface fluxes using primitive inputs.
"""
function surface_fluxes(
    param_set::APS{FT},
    Tin::FT,
    qin::FT,
    ρin::FT,
    Ts_guess::FT,
    qs_guess::FT,
    Φs::FT,
    Δz::FT,
    d::FT,
    u_int = nothing,
    u_sfc = nothing,
    roughness_inputs = nothing,
    config = nothing,
    scheme::SolverScheme = PointValueScheme(),
    solver_opts = nothing,
    flux_specs = nothing,
    update_Ts! = nothing,
    update_qs! = nothing,
) where {FT}
    return surface_fluxes(
        param_set,
        Tin,
        qin,
        zero(FT), # ql_in
        zero(FT), # qi_in
        ρin,
        Ts_guess,
        qs_guess,
        Φs,
        Δz,
        d,
        u_int,
        u_sfc,
        roughness_inputs,
        config,
        scheme,
        solver_opts,
        flux_specs,
        update_Ts!,
        update_qs!,
    )
end

function surface_fluxes(
    param_set::APS{FT},
    Tin::FT,
    qin::FT,
    ql_in::FT,
    qi_in::FT,
    ρin::FT,
    Ts_guess::FT,
    qs_guess::FT,
    Φs::FT,
    Δz::FT,
    d::FT,
    u_int = nothing,
    u_sfc = nothing,
    roughness_inputs = nothing,
    config = nothing,
    scheme::SolverScheme = PointValueScheme(),
    solver_opts = nothing,
    flux_specs = nothing,
    update_Ts! = nothing,
    update_qs! = nothing,
) where {FT}
    # 1. Build inputs
    config_val = config === nothing ? default_surface_flux_config(FT) : config
    flux_specs_val = flux_specs === nothing ? FluxSpecs(FT) : flux_specs
    
    inputs = build_surface_flux_inputs(
        param_set,
        Tin,
        qin,
        ql_in,
        qi_in,
        ρin,
        Ts_guess,
        qs_guess,
        Φs,
        Δz,
        d,
        u_int,
        u_sfc,
        config_val,
        roughness_inputs,
        flux_specs_val,
        update_Ts!,
        update_qs!,
    )

    # 2. Check for fully prescribed coefficients (bypass MOST)
    if inputs.Cd !== nothing && inputs.Ch !== nothing
        # Compute generic bulk fluxes
        # Note: We need to ensure state consistency (Ts, qs) if update functions are present?
        # But usually prescribed Cd, Ch implies simple calculation.
        # We will assume fixed state unless we want to iterate just for Ts?
        # User said "bypass MOST calculations".
        # We calculate U, Δq, ΔT and flux = C * U * Δ
        return compute_bulk_fluxes_with_coefficients(param_set, inputs, scheme)
    end

    # 3. Solve Monin-Obukhov
    solver_opts_val = normalize_solver_options(param_set, solver_opts)
    return solve_monin_obukhov(param_set, inputs, scheme, solver_opts_val)
end


# Generic wrapper for sc-style inputs (legacy support / tests)
function surface_fluxes(
    param_set::APS,
    sc,
    scheme::SolverScheme = PointValueScheme();
    config = nothing,
    kwargs...,
)
    # Extract fields from sc (Fluxes, FluxesAndFrictionVelocity, etc)
    # This assumes sc has fields: state_int, state_sfc, z0m, z0h, etc.
    # We map them to the full arguments.
    
    FT = float_type(param_set)
    ts_int = sc.state_int.ts
    ts_sfc = sc.state_sfc.ts
    thermo_params = SFP.thermodynamics_params(param_set)

    Tin = TD.air_temperature(thermo_params, ts_int)
    qin = TD.total_specific_humidity(thermo_params, ts_int)
    ρin = TD.air_density(thermo_params, ts_int)

    Ts_guess = TD.air_temperature(thermo_params, ts_sfc)
    qs_guess = TD.total_specific_humidity(thermo_params, ts_sfc)

    grav = SFP.grav(param_set)
    Φs = grav * sc.state_sfc.z
    Δz = sc.state_int.z - sc.state_sfc.z
    d = zero(Δz)

    u_int = sc.state_int.u
    u_sfc = sc.state_sfc.u

    config_val = if config !== nothing
        config
    else
        roughness = roughness_lengths(sc.z0m, sc.z0h)
        SurfaceFluxConfig(roughness, ConstantGustinessSpec(FT(1.0)))
    end

    # Handle prescribed fluxes/ustar if present in sc
    shf = hasproperty(sc, :shf) ? sc.shf : nothing
    lhf = hasproperty(sc, :lhf) ? sc.lhf : nothing
    ustar = hasproperty(sc, :ustar) ? sc.ustar : nothing
    Cd = hasproperty(sc, :Cd) ? sc.Cd : nothing
    Ch = hasproperty(sc, :Ch) ? sc.Ch : nothing

    flux_spec_args = FluxSpecs(FT; shf=shf, lhf=lhf, ustar=ustar, Cd=Cd, Ch=Ch)

    return surface_fluxes(
        param_set,
        Tin, qin, zero(FT), zero(FT), ρin, Ts_guess, qs_guess, Φs, Δz, d,
        u_int, u_sfc, nothing, config_val, scheme,
        SolverOptions(FT; kwargs...),
        flux_spec_args,
        nothing, nothing
    )
end

function default_surface_flux_config(::Type{FT}) where {FT}
    return SurfaceFluxConfig(
        ConstantRoughnessParams(FT(1e-3), FT(1e-3)),
        ConstantGustinessSpec(FT(1))
    )
end

# ------------------------------------------------------------------------------
# Solver Logic
# ------------------------------------------------------------------------------

function solve_ustar_and_roughness(
    param_set::APS{FT},
    inputs::SurfaceFluxInputs,
    scheme,
    L_MO::FT,
    ustar_guess::FT,
    ΔU::FT,
) where {FT}
    u_star_curr = ustar_guess
    local z0m::FT, z0s::FT, z0h::FT, Cd::FT, Ch::FT
    
    for _ in 1:10
        z0m, z0s = momentum_and_scalar_roughness(inputs.roughness_model, u_star_curr, param_set, inputs.roughness_inputs)
        z0h = z0s
        Cd = drag_coefficient(param_set, L_MO, z0m, inputs.Δz, scheme)
        u_star_curr = sqrt(Cd) * ΔU
    end
    
    # Final synchronization
    z0m, z0s = momentum_and_scalar_roughness(inputs.roughness_model, u_star_curr, param_set, inputs.roughness_inputs)
    z0h = z0s
    Cd = drag_coefficient(param_set, L_MO, z0m, inputs.Δz, scheme)
    Ch = heat_exchange_coefficient(param_set, L_MO, z0m, z0h, inputs.Δz, scheme)
    
    return u_star_curr, z0m, z0h, Cd, Ch
end

function compute_flux_components(
    param_set::APS{FT},
    inputs::SurfaceFluxInputs,
    Ch::FT,
    Cd::FT,
    Ts::FT,
    qs::FT,
    ρ_sfc::FT,
    b_flux::FT,
) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    
    g_h = heat_conductance(inputs, Ch, param_set, b_flux)
    g_q = g_h
    
    E = evaporation(thermo_params, inputs, g_q, inputs.qin, qs, ρ_sfc)
    lhf = latent_heat_flux(thermo_params, inputs, E)
    shf = sensible_heat_flux(param_set, thermo_params, inputs, g_h, inputs.Tin, Ts, ρ_sfc, E)
    
    phase_sfc = TD.PhasePartition(qs)
    buoy_flux = buoyancy_flux(param_set, thermo_params, shf, lhf, Ts, phase_sfc.tot, phase_sfc.liq, phase_sfc.ice, ρ_sfc)
    
    # Momentum fluxes
    u_int_x, u_int_y = inputs.u_int
    u_sfc_x, u_sfc_y = inputs.u_sfc
    Δu_x = u_int_x - u_sfc_x
    Δu_y = u_int_y - u_sfc_y
    ΔU = windspeed(inputs, param_set, b_flux)
    ρτxz = -ρ_sfc * Cd * ΔU * Δu_x
    ρτyz = -ρ_sfc * Cd * ΔU * Δu_y
    
    return (shf, lhf, buoy_flux, E, ρτxz, ρτyz)
end

function update_state_with_callbacks!(
    inputs::SurfaceFluxInputs{FT},
    iter_state::SurfaceFluxIterationState{FT},
) where {FT}
    ctx = callable_context(inputs, iter_state)
    
    if inputs.update_Ts! !== nothing
        val = inputs.update_Ts!(iter_state, ctx)
        if val isa FT
            iter_state.Ts = val
        end
    end
    if inputs.update_qs! !== nothing
        val = inputs.update_qs!(iter_state, ctx)
        if val isa FT
            iter_state.qs = val
        end
    end
end

"""
    compute_bulk_fluxes_with_coefficients(param_set, inputs, scheme)

Computes fluxes when Cd and Ch are already known.
"""
function compute_bulk_fluxes_with_coefficients(param_set::APS{FT}, inputs::SurfaceFluxInputs, scheme) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    
    # We must use current state (Ts, qs). If update functions exist, 
    # we might technically need to update them? 
    # But for prescribed coefficients, we usually just compute straight.
    
    Ts = inputs.Ts_guess
    qs = inputs.qs_guess
    ρ_sfc = surface_density(
        param_set,
        inputs.Tin,
        inputs.ρin,
        Ts,
        inputs.qin
    )
    
    # Coefficients
    Cd_in = inputs.Cd
    Ch_in = inputs.Ch
    if Cd_in === nothing || Ch_in === nothing
        # This function should only be called if Cd and Ch are provided
        error("Cd and Ch must be provided for compute_bulk_fluxes_with_coefficients")
    end
    Cd::FT = Cd_in
    Ch::FT = Ch_in
    
    # Context
    iter_state = SurfaceFluxIterationState{FT}()
    iter_state.Ts = Ts
    iter_state.qs = qs
    iter_state.ρ_sfc = ρ_sfc
    iter_state.Cd = Cd
    iter_state.Ch = Ch
    
    
    # We need gustiness to get ΔU.
    # For compute_bulk_fluxes_with_coefficients, we don't have L_MO loop unless we iterate.
    # If we are using Deardorff, we need buoyancy flux.
    # We can try to use a "lagged" approach or just 0 if not stateful?
    # Or maybe we can't easily support Deardorff here without explicit iteration?
    # For now, we use 0.0 or the value coupled with Cd?
    # But Cd is prescribed.
    # We'll use 0.0 for buoyancy flux in this context as an approximation or initial state
    # unless we want to iterate just for gustiness.
    # Given "compute_bulk_fluxes_with_coefficients" implies "simple", 0.0 is safest default.
    buoyancy_flux_val = FT(0) # Or from inputs if we had history?
    
    # Fluxes
    (shf, lhf, buoy_flux, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, Ts, qs, ρ_sfc, buoyancy_flux_val
    )
    
    # Helper to compute ustar from Cd
    # u*^2 = Cd * U^2
    ΔU = windspeed(inputs, param_set, buoyancy_flux_val)
    ustar = sqrt(Cd) * ΔU
    
    # Derived L_MO (if ustar > 0)
    # L = -u*^3 / (k B_flux)
    L_MO = zero(FT)
    if non_zero(buoy_flux) != 0 && ustar > 0
        κ = SFP.von_karman_const(param_set)
        L_MO = -ustar^3 / (κ * buoy_flux)
    end
    
    return SurfaceFluxConditions(
        L_MO, shf, lhf, buoy_flux,
        ρτxz, ρτyz,
        ustar, Cd, Ch, E
    )
end


struct ResidualFunction{FT, PS, I, UF, TP, SCH, S} <: Function
    param_set::PS
    inputs::I
    scheme::SCH
    uf_params::UF
    thermo_params::TP
    iter_state::S
    κ::FT
    grav::FT
    Ts::FT
    qs::FT
    ρ_sfc::FT
    buoyancy_flux::FT
end

function (rf::ResidualFunction{FT})(ζ) where {FT}
    # Unpack
    inputs = rf.inputs
    param_set = rf.param_set
    uf_params = rf.uf_params
    thermo_params = rf.thermo_params
    scheme = rf.scheme
    κ = rf.κ
    grav = rf.grav

    # Determine current state
    if rf.iter_state !== nothing
        # Stateful mode: State may have evolved via callbacks
        current_buoyancy_flux = rf.iter_state.buoyancy_flux
    else
        # Pure mode: we need to compute gustiness from available info?
        current_buoyancy_flux = FT(0)  # Initial guess or zero for pure mode start
    end
    
    current_ΔU = windspeed(inputs, param_set, current_buoyancy_flux)
    
    # Neutral guess for u*
    u_star_curr = maximize(current_ΔU * κ / log(inputs.Δz / FT(1e-3)), FT(1e-6))
    
    local z0m::FT, z0s::FT, z0h::FT, Cd::FT, Ch::FT
    
    # Solve for ustar and roughness
    L_MO_loop = inputs.Δz / ζ
    u_star_curr, z0m, z0h, Cd, Ch = solve_ustar_and_roughness(
        param_set, inputs, scheme, L_MO_loop, u_star_curr, current_ΔU
    )
    
    # Update state with converged values if stateful
    if rf.iter_state !== nothing
        rf.iter_state.ustar = u_star_curr
        rf.iter_state.Cd = Cd
        rf.iter_state.Ch = Ch
        
        # Approximate buoyancy flux from u* and L for next iteration's gustiness
        # L = -u*^3 / (κ B)  => B = -u*^3 / (κ L)
        # Check singularity
        if abs(L_MO_loop) > eps(FT)
             B_approx = -u_star_curr^3 / (κ * L_MO_loop)
             rf.iter_state.buoyancy_flux = B_approx
        end

        update_state_with_callbacks!(inputs, rf.iter_state)

        Ts_curr = rf.iter_state.Ts
        qs_curr = rf.iter_state.qs
    else
        # Pure mode: no callbacks invoked
        Ts_curr = rf.Ts
        qs_curr = rf.qs
    end

    Rib_theory = UF.bulk_richardson_number(uf_params, inputs.Δz, ζ, z0m, z0h, scheme)

    # 4. Compute Actual Bulk Richardson Number
    ρ_sfc_curr = surface_density(
        param_set,
        inputs.Tin,
        inputs.ρin,
        Ts_curr,
        inputs.qin
    )

    if rf.iter_state !== nothing
        rf.iter_state.ρ_sfc = ρ_sfc_curr
    end

    ts_int = TD.PhaseEquil_ρTq(thermo_params, inputs.ρin, inputs.Tin, inputs.qin)
    phase_sfc = TD.PhasePartition(qs_curr)
    θv_sfc = TD.virtual_pottemp(thermo_params, Ts_curr, ρ_sfc_curr, phase_sfc)
    θv_int = TD.virtual_pottemp(thermo_params, ts_int)

    Δθv = θv_int - θv_sfc
    θv_ref = θv_int

    Rib_state = (grav * inputs.Δz * Δθv) / (θv_ref * current_ΔU^2)

    return Rib_theory - Rib_state
end

maximize(a, b) = a > b ? a : b

function solve_monin_obukhov(param_set::APS{FT}, inputs::SurfaceFluxInputs{FT}, scheme, options::SolverOptions) where {FT}

    # Prepare iteration state
    # We use local variables for the main state to avoid allocation of mutable struct in pure mode.
    Ts_val = inputs.Ts_guess
    qs_val = inputs.qs_guess
    ρ_sfc_val = inputs.ρin

    ustar_in = inputs.ustar
    ustar_val = if ustar_in isa FT
        ustar_in
    else
        FT(0.1)
    end

    # Decide if we need stateful iteration
    use_stateful = (inputs.update_Ts! !== nothing) || (inputs.update_qs! !== nothing)

    # If stateful, we must use the mutable struct
    iter_state = if use_stateful
        SurfaceFluxIterationState{FT}(;
            Ts = Ts_val,
            qs = qs_val,
            ρ_sfc = ρ_sfc_val,
            ustar = ustar_val,
            # Others default
        )
    else
        nothing
    end

    # Default gustiness (initial)
    # Default gustiness (initial)
    # We create a temporary context for initialization
    iter_state_init = SurfaceFluxIterationState{FT}(;
        Ts = Ts_val,
        qs = qs_val,
        ρ_sfc = ρ_sfc_val,
        ustar = ustar_val,
        Cd = FT(0),  # Initial guess
        Ch = FT(0),  # Initial guess
    )
    ctx_init = callable_context(inputs, iter_state_init)

    # Initial gustiness with B = 0 implicitly handled by B=0 in iter_state

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
        iter_state,
        κ,
        grav,
        Ts_val,
        qs_val,
        ρ_sfc_val,
        FT(0) # Initial buoyancy flux
    )
    
    # Solve
    Ri_b_0 = try 
        root_function(FT(0.001))
    catch
        FT(0)
    end
    
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
    # Use locals
    
    # If stateful, we use the last updated B from the solver loop
    # If pure, we use 0.0 (initial).
    
    current_buoyancy_flux = use_stateful ? iter_state.buoyancy_flux : FT(0)
    current_ΔU = windspeed(inputs, param_set, current_buoyancy_flux)
    u_star_curr = maximize(current_ΔU * κ / log(inputs.Δz / FT(1e-3)), FT(1e-6))
    
    L_MO = inputs.Δz / ζ_final
    u_star_curr, z0m, z0h, Cd, Ch = solve_ustar_and_roughness(
        param_set, inputs, scheme, L_MO, u_star_curr, current_ΔU
    )
    
    # Update B one last time derived from L_MO? Or compute real B?
    # We compute real B below.
    
    if use_stateful
        iter_state.ustar = u_star_curr
        iter_state.Cd = Cd
        iter_state.Ch = Ch
        
        update_state_with_callbacks!(inputs, iter_state)
        
        # Note: we don't update gustiness here again, as we are done iterating.
        
        # Read back
        Ts_val = iter_state.Ts
        qs_val = iter_state.qs
        current_buoyancy_flux = iter_state.buoyancy_flux
    end
    
    # Final calculations using locals
    u_star = u_star_curr
    
    # surface_density needed for compute_flux_components
    # Wait, compute_flux_components calculates shf, lhf...
    # But it takes ρ_sfc as input.
    # Re-calculate ρ_sfc here (as logic might have updated Ts)
    
    ρ_sfc = surface_density(
        param_set,
        inputs.Tin,
        inputs.ρin,
        Ts_val,
        inputs.qin
    )
    
    (shf, lhf, buoy_flux, E, ρτxz, ρτyz) = compute_flux_components(
        param_set, inputs, Ch, Cd, Ts_val, qs_val, ρ_sfc, current_buoyancy_flux
    )
    
    return SurfaceFluxConditions(
        L_MO, shf, lhf, buoy_flux,
        ρτxz, ρτyz,
        u_star, Cd, Ch, E
    )
end

@inline function callable_context(
    inputs::SurfaceFluxInputs,
    iter_state::SurfaceFluxIterationState,
)
    return CallableContext(
        inputs.Tin,
        inputs.qin,
        inputs.ρin,
        iter_state.Ts,
        iter_state.qs,
        inputs.Φs,
        inputs.Δz,
        inputs.d,
        inputs.u_int,
        inputs.u_sfc,
        iter_state.buoyancy_flux,
        iter_state.ustar,
        iter_state.Cd,
        iter_state.Ch,
        iter_state.ρ_sfc,
    )
end

# Backward compatibility for sc-based calls relying on old propertynames
obukhov_similarity_solution(sfc::SurfaceFluxConditions) = sfc.L_MO

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # module
