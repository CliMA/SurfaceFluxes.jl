"""
    SurfaceFluxes

Surface-layer flux calculations based on Monin-Obukhov similarity theory.
"""
module SurfaceFluxes

include("UniversalFunctions.jl")
include("Parameters.jl")

import Thermodynamics
const TD = Thermodynamics

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
include("evaporation_methods.jl")
include("latent_heat_methods.jl")
include("sensible_heat_methods.jl")
include("momentum_exchange_coefficient_methods.jl")
include("heat_exchange_coefficient_methods.jl")
include("friction_velocity_methods.jl")
include("buoyancy_flux_methods.jl")
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
    charnock_momentum(param_set; α = 0.011)

Convenience helper returning a type-stable callable that evaluates the Charnock
relation for the momentum roughness length.
"""
@inline function charnock_momentum(
    param_set::APS{FT};
    α = FT(0.011),
) where {FT}
    return SurfaceCallable(CharnockMomentum{FT}(convert(FT, α), SFP.grav(param_set)))
end

"""
    surface_fluxes(param_set, Tin, qin, ρin, Ts, qs, Φs, Δz, d,
                   u_in, u_sfc, roughness, gustiness, scheme,
                   solver_opts, flux_specs)

Compute near-surface fluxes using primitive inputs. `Ts` and `qs` may be
scalars or callables that are evaluated every iteration with the current flux
context. Roughness lengths and gustiness may likewise be scalars or the
module-defined functional forms exposed via `SurfaceFluxInputs`.

# Required arguments
- `param_set`: Surface flux parameter set
- `Tin`: Interior temperature [K]
- `qin`: Interior specific humidity [kg/kg]
- `ρin`: Interior air density [kg/m³]
- `Ts`: Surface temperature (scalar or callable)
- `qs`: Surface specific humidity (scalar or callable)
- `Φs`: Surface geopotential [m²/s²]
- `Δz`: Height difference between interior and surface reference levels [m]
- `d`: Displacement height [m]

All subsequent arguments are positional to remain GPU-friendly. Defaults are
provided for winds, roughness lengths, gustiness, the solver scheme,
iteration tolerances, and prescribed flux specifications. Use
`roughness_lengths`, `flux_spec`, and `gustiness_constant` for clarity when
configuring these optional specifications.
"""
function surface_fluxes(
    param_set::APS{FT},
    Tin::FT,
    qin::FT,
    ρin::FT,
    Ts,
    qs,
    Φs::FT,
    Δz::FT,
    d::FT,
    u_in = nothing,
    u_sfc = nothing,
    roughness = nothing,
    gustiness = nothing,
    scheme::SolverScheme = PointValueScheme(),
    solver_opts = nothing,
    flux_specs = nothing,
) where {FT}
    u_in_val = u_in === nothing ? (zero(FT), zero(FT)) : u_in
    u_sfc_val = u_sfc === nothing ? (zero(FT), zero(FT)) : u_sfc
    roughness_val = roughness === nothing ? default_roughness_lengths(param_set) : roughness
    gustiness_val = gustiness === nothing ? FT(1) : gustiness
    solver_opts_val = normalize_solver_options(param_set, solver_opts)
    flux_specs_val = normalize_flux_specs(param_set, flux_specs)

    inputs = build_surface_flux_inputs(
        param_set,
        Tin,
        qin,
        ρin,
        Ts,
        qs,
        Φs,
        Δz,
        d,
        u_in_val,
        u_sfc_val,
        roughness_val,
        gustiness_val,
        flux_specs_val,
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
    T_in = inputs.Tin
    q_in = inputs.qin
    ρ_in = inputs.ρin
    Φ_in = interior_geopotential(param_set, inputs)
    Φ_sfc = surface_geopotential(inputs)
    phase_in = TD.PhasePartition(q_in)
    cp_m_in = TD.cp_m(thermo_params, phase_in)
    cv_m_in = TD.cv_m(thermo_params, phase_in)
    R_m_in = TD.gas_constant_air(thermo_params, phase_in)
    DSEᵥ_in_val = DSEᵥ(param_set, T_in, phase_in, Φ_in)
    uf_params = SFP.uf_params(param_set)
    grav = SFP.grav(param_set)
    FT = float_type(param_set)

    iter_state = SurfaceFluxIterationState{FT}()
    iter_state.Ts = T_in
    iter_state.qs = q_in
    iter_state.ρ_sfc = ρ_in
    inputs.ustar !== nothing && (iter_state.ustar = inputs.ustar)
    inputs.shf !== nothing && (iter_state.shf = inputs.shf)
    inputs.lhf !== nothing && (iter_state.lhf = inputs.lhf)
    inputs.Cd !== nothing && (iter_state.Cd = inputs.Cd)
    inputs.Ch !== nothing && (iter_state.Ch = inputs.Ch)
    ctx0 = callable_context(inputs, iter_state, T_in, q_in, ρ_in)
    iter_state.gustiness = gustiness_value(inputs.gustiness_model, inputs, ctx0)

    prev_state = nothing
    last = nothing

    for _ in 1:maxiter
        ctx = callable_context(inputs, iter_state, T_in, q_in, ρ_in)
        iter_state.Ts = resolve_quantity(inputs.Ts_spec, ctx)
        iter_state.qs = resolve_quantity(inputs.qs_spec, ctx)

        ρ_sfc = dry_adiabatic_surface_density(cv_m_in, R_m_in, T_in, ρ_in, iter_state.Ts)
        iter_state.ρ_sfc = ρ_sfc
        phase_sfc = TD.PhasePartition(iter_state.qs)
        cp_m_sfc = TD.cp_m(thermo_params, phase_sfc)

        ctx = callable_context(inputs, iter_state, T_in, q_in, ρ_in)
        gustiness = gustiness_value(inputs.gustiness_model, inputs, ctx)
        iter_state.gustiness = gustiness

        ΔDSE_val = ΔDSEᵥ(param_set, T_in, phase_in, Φ_in, iter_state.Ts, phase_sfc, Φ_sfc)
        Δθᵥ_val = Δθᵥ(param_set, T_in, ρ_in, phase_in, iter_state.Ts, ρ_sfc, phase_sfc)
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
            DSEᵥ_in_val,
            grav,
        )

        Cd_val = inputs.Cd === nothing ? momentum_exchange_coefficient(
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

        Ch_val = inputs.Ch === nothing ? heat_exchange_coefficient(
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

        E = evaporation(
            param_set,
            inputs,
            Ch_val,
            q_in,
            iter_state.qs,
            ρ_sfc,
            gustiness,
        )
        lhf = latent_heat_flux(param_set, inputs, E)
        shf = sensible_heat_flux(
            param_set,
            thermo_params,
            inputs,
            Ch_val,
            cp_m_in,
            cp_m_sfc,
            T_in,
            iter_state.Ts,
            ρ_sfc,
            gustiness,
            E,
        )
        buoy_flux = compute_buoyancy_flux(
            param_set,
            thermo_params,
            shf,
            lhf,
            cp_m_in,
            T_in,
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
        inputs.u_in,
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
    DSEᵥ_in_val,
    grav,
) 
    FT = typeof(approximate_state.u_star)
    u_star = approximate_state.u_star
    ell_u = compute_z0(u_star, param_set, inputs, UF.MomentumTransport(), ctx)
    ell_theta = compute_z0(u_star, param_set, inputs, UF.HeatTransport(), ctx)
    ell_q = compute_z0(u_star, param_set, inputs, UF.HeatTransport(), ctx)
    κ = SFP.von_karman_const(param_set)
    dsev_star = approximate_state.dsev_star
    b_star = dsev_star * grav / DSEᵥ_in_val
    if abs(ΔDSE_val) <= tol_neutral
        sgn = iszero(ΔDSE_val) ? one(FT) : sign(ΔDSE_val)
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
