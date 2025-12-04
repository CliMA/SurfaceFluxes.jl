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

import .UniversalFunctions
const UF = UniversalFunctions

import .Parameters
const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

include("types.jl")
include("thermo_primitives.jl")
include("input_builders.jl")
include("utilities.jl")
include("physical_scale_coefficient_methods.jl")
include("roughness_lengths.jl")
include("thermodynamic_fluxes.jl")
include("momentum_exchange_coefficient_methods.jl")
include("heat_exchange_coefficient_methods.jl")
include("friction_velocity_methods.jl")
include("conductances.jl")
include("profile_recovery.jl")

@inline float_type(::APS{FT}) where {FT} = FT

@inline function default_solver_options(param_set::APS{FT}) where {FT}
    return SolverOptions(
        FT;
        tol = sqrt(eps(FT)),
        tol_neutral = SFP.cp_d(param_set) / 100,
        maxiter = 30,
    )
end

@inline function default_flux_specs(::APS{FT}) where {FT}
    return FluxSpecs{FT}()
end

@inline function normalize_solver_options(param_set::APS{FT}, solver_opts) where {FT}
    if solver_opts === nothing
        return default_solver_options(param_set)
    elseif solver_opts isa Int
        return SolverOptions(param_set; maxiter = solver_opts)
    elseif solver_opts isa SolverOptions{FT}
        return solver_opts
    elseif solver_opts isa SolverOptions
        return SolverOptions(
            FT;
            tol = solver_opts.tol,
            tol_neutral = solver_opts.tol_neutral,
            maxiter = solver_opts.maxiter,
        )
    else
        throw(ArgumentError("Unsupported solver_opts specification: $(typeof(solver_opts))"))
    end
end

@inline function normalize_flux_specs(param_set::APS{FT}, specs) where {FT}
    if specs === nothing
        return default_flux_specs(param_set)
    elseif specs isa FluxSpecs{FT}
        return specs
    elseif specs isa FluxSpecs
        return FluxSpecs(
            FT;
            shf = specs.shf,
            lhf = specs.lhf,
            ustar = specs.ustar,
            Cd = specs.Cd,
            Ch = specs.Ch,
        )
    elseif specs isa NamedTuple
        return FluxSpecs(param_set; specs...)
    else
        throw(ArgumentError("Unsupported flux specification: $(typeof(specs))"))
    end
end

@inline FluxSpecs(param_set::APS{FT}; kwargs...) where {FT} = FluxSpecs(FT; kwargs...)
@inline SolverOptions(param_set::APS{FT}; kwargs...) where {FT} = SolverOptions(FT; kwargs...)

"""
    charnock_momentum(; α = 0.011, scalar = 1e-4)

Convenience helper returning a roughness specification that applies the
Charnock relation for the momentum roughness length while keeping a scalar
roughness length fixed.
"""
@inline function charnock_momentum(; α = 0.011, scalar = 1e-4)
    return CharnockRoughnessSpec(α, scalar)
end

"""
    surface_fluxes(param_set, Tin, qin, ρin, Ts_guess, qs_guess, Φs, Δz, d,
                   u_int, u_sfc, config, scheme,
                   solver_opts, flux_specs, update_Ts!, update_qs!)

Compute near-surface fluxes using primitive inputs. `Ts_guess` and `qs_guess`
seed the nonlinear MOST iteration, while the optional `update_Ts!` and
`update_qs!` hooks may adjust these values after every iteration. Roughness
and gustiness parameterizations are provided via `SurfaceFluxConfig`, ensuring
the solver only relies on module-defined models.

# Required arguments
- `param_set`: Surface flux parameter set
- `Tin`: Interior temperature [K]
- `qin`: Interior specific humidity [kg/kg]
- `ρin`: Interior air density [kg/m³]
- `Ts_guess`: Initial surface temperature guess [K]
- `qs_guess`: Initial surface specific humidity guess [kg/kg]
- `Φs`: Surface geopotential [m²/s²]
- `Δz`: Height difference between interior and surface reference levels [m]
- `d`: Displacement height [m]

All subsequent arguments are positional to remain GPU-friendly. Defaults are
provided for winds, the solver scheme, iteration tolerances, prescribed flux
specifications, and the surface parameterization config. Use `flux_spec` and
`gustiness_constant` for clarity when configuring optional specifications.
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
    config::SurfaceFluxConfig = SurfaceFluxConfig(),
    scheme::SolverScheme = PointValueScheme(),
    solver_opts = nothing,
    flux_specs = nothing,
    update_Ts! = nothing,
    update_qs! = nothing,
) where {FT}
    u_int_val = u_int === nothing ? (zero(FT), zero(FT)) : u_int
    u_sfc_val = u_sfc === nothing ? (zero(FT), zero(FT)) : u_sfc
    config_val = config === nothing ? SurfaceFluxConfig() : config
    solver_opts_val = normalize_solver_options(param_set, solver_opts)
    flux_specs_val = normalize_flux_specs(param_set, flux_specs)

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
        u_int_val,
        u_sfc_val,
        config_val,
        flux_specs_val,
        update_Ts!,
        update_qs!,
    )
    thermo_params = SFP.thermodynamics_params(param_set)
    solution = solve_surface_layer(
        param_set,
        thermo_params,
        inputs,
        scheme,
        solver_opts_val,
    )

    scales = solution.scales
    ρτxz, ρτyz = momentum_fluxes(
        param_set,
        solution.Cd,
        inputs,
        solution.ρ_sfc,
        solution.gustiness,
    )
    buoy_flux = solution.buoyancy_flux

    return SurfaceFluxConditions(
        scales.L_star,
        solution.shf,
        solution.lhf,
        buoy_flux,
        ρτxz,
        ρτyz,
        scales.u_star,
        solution.Cd,
        solution.Ch,
        solution.evaporation,
    )
end

function solve_surface_layer(
    param_set::APS,
    thermo_params,
    inputs::SurfaceFluxInputs,
    scheme::SolverScheme,
    solver_opts::SolverOptions,
)
    tol = solver_opts.tol
    tol_neutral = solver_opts.tol_neutral
    maxiter = solver_opts.maxiter
    T_int = inputs.Tin
    q_in = inputs.qin
    ρ_int = inputs.ρin
    Φ_int = interior_geopotential(param_set, inputs)
    Φ_sfc = surface_geopotential(inputs)
    phase_int = TD.PhasePartition(q_in)
    cp_m_int = TD.cp_m(thermo_params, phase_int)
    cv_m_int = TD.cv_m(thermo_params, phase_int)
    R_m_int = TD.gas_constant_air(thermo_params, phase_int)
    DSEᵥ_int_val = DSEᵥ(param_set, T_int, phase_int, Φ_int)
    uf_params = SFP.uf_params(param_set)
    grav = SFP.grav(param_set)
    FT = float_type(param_set)

    iter_state = SurfaceFluxIterationState{FT}()
    iter_state.Ts = inputs.Ts_guess
    iter_state.qs = inputs.qs_guess
    iter_state.ρ_sfc = ρ_int
    inputs.ustar !== nothing && (iter_state.ustar = inputs.ustar)
    inputs.shf !== nothing && (iter_state.shf = inputs.shf)
    inputs.lhf !== nothing && (iter_state.lhf = inputs.lhf)
    inputs.Cd !== nothing && (iter_state.Cd = inputs.Cd)
    inputs.Ch !== nothing && (iter_state.Ch = inputs.Ch)
    ctx0 = callable_context(inputs, iter_state, T_int, q_in, ρ_int)
    iter_state.gustiness = gustiness_value(inputs.gustiness_model, inputs, ctx0)

    prev_state = nothing
    last = nothing

    for _ in 1:maxiter
        ctx = callable_context(inputs, iter_state, T_int, q_in, ρ_int)
        did_update = false
        if inputs.update_Ts! !== nothing
            new_Ts = inputs.update_Ts!(iter_state, ctx)
            if new_Ts !== nothing
                iter_state.Ts = convert(FT, new_Ts)
            end
            did_update = true
        end
        if inputs.update_qs! !== nothing
            new_qs = inputs.update_qs!(iter_state, ctx)
            if new_qs !== nothing
                iter_state.qs = convert(FT, new_qs)
            end
            did_update = true
        end
        if did_update
            ctx = callable_context(inputs, iter_state, T_int, q_in, ρ_int)
        end

        ρ_sfc = dry_adiabatic_surface_density(cv_m_int, R_m_int, T_int, ρ_int, iter_state.Ts)
        iter_state.ρ_sfc = ρ_sfc
        phase_sfc = TD.PhasePartition(iter_state.qs)

        ctx = callable_context(inputs, iter_state, T_int, q_in, ρ_int)
        gustiness = gustiness_value(inputs.gustiness_model, inputs, ctx)
        iter_state.gustiness = gustiness

        ΔDSE_val = ΔDSEᵥ(param_set, T_int, phase_int, Φ_int, iter_state.Ts, phase_sfc, Φ_sfc)
        Δθᵥ_val = Δθᵥ(param_set, T_int, ρ_int, phase_int, iter_state.Ts, ρ_sfc, phase_sfc)
        Δq_val = Δqt(q_in, iter_state.qs)
        ΔU = windspeed(inputs, gustiness)

        if prev_state === nothing
            u_star₀ = iter_state.ustar
            ell_u₀ = compute_z0(u_star₀, param_set, inputs, UF.MomentumTransport(), ctx)
            ell_theta₀ = compute_z0(u_star₀, param_set, inputs, UF.HeatTransport(), ctx)
            ell_q₀ = ell_theta₀
            δ = iszero(ΔDSE_val) ? one(FT) : sign(ΔDSE_val)
            L_init = δ ≥ 0 ? FT(10) : FT(-10)
            prev_state = SimilarityScales(
                u_star₀,
                FT(δ),
                FT(δ),
                L_init,
                FT(δ),
                ell_u₀,
                ell_theta₀,
                ell_q₀,
            )
        end

        current = iterate_interface_fluxes(
            param_set,
            inputs,
            gustiness,
            prev_state,
            scheme,
            tol_neutral,
            uf_params,
            ctx,
            ΔU,
            ΔDSE_val,
            Δθᵥ_val,
            Δq_val,
            DSEᵥ_int_val,
            grav,
        )

        Cd_val =
            inputs.Cd === nothing ?
            momentum_exchange_coefficient(
                param_set,
                current.L_star,
                current.u_star,
                current.ell_u,
                inputs,
                scheme,
                tol_neutral,
                gustiness,
                ΔDSE_val,
            ) : inputs.Cd

        Ch_val =
            inputs.Ch === nothing ?
            heat_exchange_coefficient(
                param_set,
                current.L_star,
                current.u_star,
                current.ell_u,
                current.ell_theta,
                inputs,
                scheme,
                tol_neutral,
                gustiness,
                ΔDSE_val,
            ) : inputs.Ch

        g_h = heat_conductance(inputs, Ch_val, gustiness)
        g_q = heat_conductance(inputs, Ch_val, gustiness)
        E = evaporation(
            param_set,
            thermo_params,
            inputs,
            g_q,
            q_in,
            iter_state.qs,
            ρ_sfc,
        )
        lhf = latent_heat_flux(param_set, thermo_params, inputs, E)
        shf = sensible_heat_flux(
            param_set,
            thermo_params,
            inputs,
            g_h,
            T_int,
            iter_state.Ts,
            ρ_sfc,
            E,
        )
        buoy_flux = buoyancy_flux(
            param_set,
            thermo_params,
            shf,
            lhf,
            iter_state.Ts,
            phase_sfc.tot,
            phase_sfc.liq,
            phase_sfc.ice,
            ρ_sfc,
        )

        iter_state.ustar = current.u_star
        iter_state.L_MO = current.L_star
        iter_state.Cd = Cd_val
        iter_state.Ch = Ch_val
        iter_state.shf = shf
        iter_state.lhf = lhf
        iter_state.evaporation = E
        iter_state.buoyancy_flux = buoy_flux

        last = SolverSnapshot(
            current,
            ρ_sfc,
            gustiness,
            Cd_val,
            Ch_val,
            shf,
            lhf,
            E,
            buoyancy_flux,
        )

        has_converged(current, prev_state, tol) && break
        prev_state = current
    end

    return last
end

@inline function callable_context(
    inputs::SurfaceFluxInputs,
    iter_state::SurfaceFluxIterationState,
    Tin,
    qin,
    ρin,
)
    return CallableContext(
        Tin,
        qin,
        ρin,
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

@inline function compute_Fₘₕ(inputs::SurfaceFluxInputs, ufₛ, ζ, 𝓁, transport)
    return log(inputs.Δz / 𝓁) -
           UF.psi(ufₛ, ζ, transport) +
           UF.psi(ufₛ, 𝓁 * ζ / inputs.Δz, transport)
end

function iterate_interface_fluxes(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    gustiness,
    approximate_state,
    scheme::SolverScheme,
    tol_neutral,
    uf_params,
    ctx,
    ΔU,
    ΔDSE_val,
    Δθᵥ_val,
    Δq_val,
    DSEᵥ_int_val,
    grav,
)
    FT = typeof(approximate_state.u_star)
    u_star = approximate_state.u_star
    ell_u = compute_z0(u_star, param_set, inputs, UF.MomentumTransport(), ctx)
    ell_theta = compute_z0(u_star, param_set, inputs, UF.HeatTransport(), ctx)
    ell_q = compute_z0(u_star, param_set, inputs, UF.HeatTransport(), ctx)
    κ = SFP.von_karman_const(param_set)
    dsev_star = approximate_state.dsev_star
    b_star = dsev_star * grav / DSEᵥ_int_val
    if abs(b_star) <= eps(FT)
        sgn = iszero(b_star) ? one(FT) : sign(b_star)
        L_star = sgn * FT(Inf)
    else
        L_star = u_star^2 / (κ * b_star)
    end
    ζ = inputs.Δz / L_star
    χu = κ / compute_Fₘₕ(inputs, uf_params, ζ, ell_u, UF.MomentumTransport())
    χDSE = κ / compute_Fₘₕ(inputs, uf_params, ζ, ell_theta, UF.HeatTransport())
    χq = κ / compute_Fₘₕ(inputs, uf_params, ζ, ell_q, UF.HeatTransport())
    χθ = κ / compute_Fₘₕ(inputs, uf_params, ζ, ell_theta, UF.HeatTransport())
    u_star = χu * ΔU
    dsev_star = χDSE * ΔDSE_val
    q_star = iszero(Δq_val) ? zero(FT) : χq * Δq_val
    theta_v_star = χθ * Δθᵥ_val
    return SimilarityScales(u_star, dsev_star, q_star, L_star, theta_v_star, ell_u, ell_theta, ell_q)
end

@inline function has_converged(current, previous, tol)
    if previous === nothing
        return false
    end
    return abs(current.L_star - previous.L_star) <= tol &&
           abs(current.u_star - previous.u_star) <= tol &&
           abs(current.q_star - previous.q_star) <= tol &&
           abs(current.dsev_star - previous.dsev_star) <= tol
end

function momentum_fluxes(
    param_set::APS,
    Cd,
    inputs::SurfaceFluxInputs,
    ρ_sfc,
    gustiness,
)
    Δu = Δu_components(inputs)
    ΔU = windspeed(inputs, gustiness)
    ρτxz = -ρ_sfc * Cd * Δu[1] * ΔU
    ρτyz = -ρ_sfc * Cd * Δu[2] * ΔU
    return (ρτxz, ρτyz)
end

obukhov_similarity_solution(sfc::SurfaceFluxConditions) = sfc.L_MO

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # SurfaceFluxes module
