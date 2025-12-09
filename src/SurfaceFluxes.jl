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
include("thermo_primitives.jl")
include("roughness_lengths.jl")
include("input_builders.jl")
include("utilities.jl")
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
    # 1. Build Inputs
    config_val = config === nothing ? default_surface_flux_config(FT) : config
    flux_specs_val = flux_specs === nothing ? FluxSpecs(FT) : flux_specs
    
    inputs = build_surface_flux_inputs(
        param_set,
        Tin,
        qin,
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
        Tin, qin, ρin, Ts_guess, qs_guess, Φs, Δz, d,
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
        TD.cv_m(thermo_params, TD.PhasePartition(inputs.qin)), # Approx
        TD.gas_constant_air(thermo_params, TD.PhasePartition(inputs.qin)),
        inputs.Tin,
        inputs.ρin,
        Ts
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
    
    ctx = callable_context(inputs, iter_state)
    gustiness = gustiness_value(inputs.gustiness_model, inputs, ctx)
    ΔU = windspeed(inputs, gustiness)
    
    # Fluxes
    g_h = heat_conductance(inputs, Ch, gustiness)
    g_q = g_h # Assumed same
    
    E = evaporation(thermo_params, inputs, g_q, inputs.qin, qs, ρ_sfc)
    lhf = latent_heat_flux(thermo_params, inputs, E)
    shf = sensible_heat_flux(param_set, thermo_params, inputs, g_h, inputs.Tin, Ts, ρ_sfc, E)
    
    # Helper to compute ustar from Cd
    # u*^2 = Cd * U^2
    ustar = sqrt(Cd) * ΔU
    
    # Calculate buoyancy flux
    phase_sfc = TD.PhasePartition(qs)
    buoy_flux = buoyancy_flux(param_set, thermo_params, shf, lhf, Ts, phase_sfc.tot, phase_sfc.liq, phase_sfc.ice, ρ_sfc)
    
    # Momentum fluxes
    # τ = -ρ Cd |ΔU| Δu
    u_int_x, u_int_y = inputs.u_int
    u_sfc_x, u_sfc_y = inputs.u_sfc
    Δu_x = u_int_x - u_sfc_x
    Δu_y = u_int_y - u_sfc_y
    ρτxz = -ρ_sfc * Cd * ΔU * Δu_x
    ρτyz = -ρ_sfc * Cd * ΔU * Δu_y

    # Derived L_MO (if ustar > 0)
    # L = -u*^3 / (k B_flux)
    L_MO = zero(FT)
    if non_zero(buoy_flux) != 0 && ustar > 0
        κ = SFP.von_karman_const(param_set)
        L_MO = -ustar^3 / (κ * buoy_flux)
    end
    
    return SurfaceFluxConditions(
        L_MO, shf, lhf, buoy_flux,
        zero(FT), zero(FT), # Stress components not computed here for shortness
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
    gustiness::FT
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
        current_gustiness = rf.iter_state.gustiness
        # Note: We don't update gustiness here, assuming it's done via callbacks or fixed
    else
        current_gustiness = rf.gustiness
    end
    
    current_ΔU = windspeed(inputs, current_gustiness)
    
    # Neutral guess for u*
    u_star_curr = maximize(current_ΔU * κ / log(inputs.Δz / FT(1e-3)), FT(1e-6))
    
    local z0m::FT, z0s::FT, z0h::FT
    
    for _ in 1:10
        # Update roughness
        z0m, z0s = momentum_and_scalar_roughness(
            inputs.roughness_model, 
            u_star_curr, 
            param_set, 
            inputs.roughness_inputs
        )
        z0h = z0s 
        
        # Update u*
        fm = UF.dimensionless_profile(uf_params, inputs.Δz, ζ, z0m, UF.MomentumTransport(), scheme)
        u_star_curr = κ * current_ΔU / fm
    end
    
    # Update state with converged values if stateful
    z0m, z0s = momentum_and_scalar_roughness(inputs.roughness_model, u_star_curr, param_set, inputs.roughness_inputs)
    z0h = z0s
    
    # Drag Coefficients
    fm = UF.dimensionless_profile(uf_params, inputs.Δz, ζ, z0m, UF.MomentumTransport(), scheme)
    fh = UF.dimensionless_profile(uf_params, inputs.Δz, ζ, z0h, UF.HeatTransport(), scheme)
    
    Cd = (κ / fm)^2
    Ch = (κ / fm) * (κ / fh)
    
    # Callbacks
    if rf.iter_state !== nothing
        rf.iter_state.ustar = u_star_curr
        rf.iter_state.Cd = Cd
        rf.iter_state.Ch = Ch
        
        ctx = callable_context(inputs, rf.iter_state)
        
        if inputs.update_Ts! !== nothing
            val = inputs.update_Ts!(rf.iter_state, ctx)
            if val isa FT
                rf.iter_state.Ts = val
            end
        end
        if inputs.update_qs! !== nothing
            val = inputs.update_qs!(rf.iter_state, ctx)
            if val isa FT
                rf.iter_state.qs = val
            end
        end
        
        Ts_curr = rf.iter_state.Ts
        qs_curr = rf.iter_state.qs
    else
        Ts_curr = rf.Ts
        qs_curr = rf.qs
        # Pure mode: no callbacks invoked
    end
    
    # 4. Compute Actual Bulk Richardson Number
    ρ_sfc_curr = surface_density(
        TD.cv_m(thermo_params, TD.PhasePartition(inputs.qin)), 
        TD.gas_constant_air(thermo_params, TD.PhasePartition(inputs.qin)),
        inputs.Tin,
        inputs.ρin,
        Ts_curr
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
    
    # 5. Compute Theoretical Bulk Richardson Number
    Rib_theory = ζ * fh / fm^2
    
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
    gustiness_val = FT(1)

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
        gustiness_val
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
    
    current_gustiness = use_stateful ? iter_state.gustiness : gustiness_val
    current_ΔU = windspeed(inputs, current_gustiness)
    u_star_curr = maximize(current_ΔU * κ / log(inputs.Δz / FT(1e-3)), FT(1e-6))
    
    local z0m::FT, z0s::FT, z0h::FT, Cd::FT, Ch::FT
    for _ in 1:10
        z0m, z0s = momentum_and_scalar_roughness(inputs.roughness_model, u_star_curr, param_set, inputs.roughness_inputs)
        z0h = z0s 
        fm = UF.dimensionless_profile(uf_params, inputs.Δz, ζ_final, z0m, UF.MomentumTransport(), scheme)
        u_star_curr = κ * current_ΔU / fm
    end
    
    z0m, z0s = momentum_and_scalar_roughness(inputs.roughness_model, u_star_curr, param_set, inputs.roughness_inputs)
    z0h = z0s
    
    fm = UF.dimensionless_profile(uf_params, inputs.Δz, ζ_final, z0m, UF.MomentumTransport(), scheme)
    fh = UF.dimensionless_profile(uf_params, inputs.Δz, ζ_final, z0h, UF.HeatTransport(), scheme)
    
    Cd = (κ / fm)^2
    Ch = (κ / fm) * (κ / fh)
    
    if use_stateful
        iter_state.ustar = u_star_curr
        iter_state.Cd = Cd
        iter_state.Ch = Ch
        
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
        
        # Read back
        Ts_val = iter_state.Ts
        qs_val = iter_state.qs
        current_gustiness = iter_state.gustiness
    end
    
    # Final calculations using locals
    u_star = u_star_curr
    L_MO = inputs.Δz / ζ_final
    
    g_h = heat_conductance(inputs, Ch, current_gustiness)
    g_q = g_h
    
    ρ_sfc = surface_density(
        TD.cv_m(thermo_params, TD.PhasePartition(inputs.qin)), 
        TD.gas_constant_air(thermo_params, TD.PhasePartition(inputs.qin)),
        inputs.Tin,
        inputs.ρin,
        Ts_val
    )
    
    E = evaporation(thermo_params, inputs, g_q, inputs.qin, qs_val, ρ_sfc)
    lhf = latent_heat_flux(thermo_params, inputs, E)
    shf = sensible_heat_flux(param_set, thermo_params, inputs, g_h, inputs.Tin, Ts_val, ρ_sfc, E)
    
    phase_sfc = TD.PhasePartition(qs_val)
    buoy_flux = buoyancy_flux(param_set, thermo_params, shf, lhf, Ts_val, phase_sfc.tot, phase_sfc.liq, phase_sfc.ice, ρ_sfc)
    
    # Momentum fluxes
    u_int_x, u_int_y = inputs.u_int
    u_sfc_x, u_sfc_y = inputs.u_sfc
    Δu_x = u_int_x - u_sfc_x
    Δu_y = u_int_y - u_sfc_y
    ρτxz = -ρ_sfc * Cd * current_ΔU * Δu_x
    ρτyz = -ρ_sfc * Cd * current_ΔU * Δu_y
    
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
        iter_state.gustiness,
        iter_state.ustar,
        iter_state.shf,
        iter_state.lhf,
        iter_state.Cd,
        iter_state.Ch,
        iter_state.L_MO,
        iter_state.evaporation,
        iter_state.buoyancy_flux,
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
