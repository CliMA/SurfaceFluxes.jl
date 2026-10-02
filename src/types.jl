# Surface flux configuration types: roughness, gustiness, moisture, and the
# containers (`SurfaceFluxConfig`, `FluxSpecs`, `SolverOptions`) that bundle the
# user-facing specifications consumed by `surface_fluxes`.

abstract type AbstractRoughnessParams end
abstract type AbstractGustinessSpec end
abstract type AbstractRoughnessSubLayerModel end
abstract type AbstractStabilityCap end



"""
    ConstantGustinessSpec{TG <: Real}

A gustiness model where the gustiness velocity is a constant value.
# Fields
- `value`: The constant gustiness velocity [m/s].
"""
struct ConstantGustinessSpec{TG <: Real} <: AbstractGustinessSpec
    value::TG
end

"""
    DeardorffGustinessSpec

A gustiness model based on Deardorff (1970) scaling with convective velocity scale ``w_*``.
The gustiness velocity is computed as:
``U_{gust} = \\beta w_*``
where ``w_* = (B z_i)^{1/3}``.
"""
struct DeardorffGustinessSpec <: AbstractGustinessSpec end

"""
    FlooredDeardorffGustinessSpec{FT <: Real}

The larger of a constant minimum wind speed `u_min` [m/s] and the convective
gustiness ``β w_*`` of [`DeardorffGustinessSpec`](@ref), with the convective velocity
scale ``w_* = (B z_i)^{1/3}`` of Deardorff (1970) for a positive surface buoyancy flux
``B``, the boundary layer depth ``z_i`` (`gustiness_zi`) and coefficient ``β``
(`gustiness_coeff`). The convective part vanishes in stable conditions, where the floor
applies; `FlooredDeardorffGustinessSpec(zero(FT))` is the pure convective gustiness.

Within the stability solve, the convective part is evaluated in closed form from the
surface and atmospheric state (see [`free_convection_wind_speed`](@ref)): at a stability
parameter ``ζ``, the bulk relations ``u_* = κ U / F_m`` and ``θ_{v*} = κ Δθ_v / F_h``
make the buoyancy flux ``B = (g/θ_v) u_* θ_{v*}`` linear in the effective wind speed
``U``, so that ``U = β w_*`` has the solution

```math
U^2 = β^3 κ^2 \\frac{g}{θ_v} z_i \\frac{Δθ_v}{F_m(ζ) F_h(ζ)}.
```

The effective wind speed is the fixed point of ``U = \\max(|Δu|, u_{min}, β w_*(U))``,
and the gustiness is independent of ``u_*`` (see [`depends_on_ustar`](@ref)), so the
friction velocity follows from ``ζ`` in closed form and the free-convection limit is well
posed at every ``ζ``. [`DeardorffGustinessSpec`](@ref) evaluates the gustiness from the
buoyancy flux implied by ``ζ`` and the current ``u_*``; at fixed ``ζ`` that gustiness is
proportional to ``u_*``, and a consistent ``u_*`` exists only up to the free-convection
limit.

# Fields
- `u_min`: Minimum wind speed [m/s].

# Examples
```julia
gustiness = FlooredDeardorffGustinessSpec(0.5)  # floor of 0.5 m/s
config = SurfaceFluxConfig(ConstantRoughnessParams(0.01, 0.001), gustiness)
```

# References
- Deardorff, J. W. (1970). Convective velocity and temperature scales for the unstable
  planetary boundary layer and for Rayleigh convection. J. Atmos. Sci., 27, 1211-1213.
- Beljaars, A. C. M. (1995). The parametrization of surface fluxes in large-scale models
  under free convection. Q. J. R. Meteorol. Soc., 121, 255-270.
"""
struct FlooredDeardorffGustinessSpec{FT <: Real} <: AbstractGustinessSpec
    u_min::FT
end

Base.broadcastable(p::AbstractRoughnessParams) = tuple(p)
Base.broadcastable(p::AbstractGustinessSpec) = tuple(p)

abstract type AbstractMoistureModel end

"""
    MoistModel

Indicates that moisture effects (latent heat, virtual temperature) should be included in the flux calculations.
"""
struct MoistModel <: AbstractMoistureModel end

"""
    DryModel

Indicates that moisture effects should be ignored (sensible heat and momentum only).
"""
struct DryModel <: AbstractMoistureModel end

Base.broadcastable(m::AbstractMoistureModel) = tuple(m)

"""
    SurfaceFluxConfig

Configuration for surface flux calculation components.

# Fields
- `roughness`: The roughness length parameterization to use (e.g., [`ConstantRoughnessParams`](@ref)).
- `gustiness`: The gustiness parameterization to use (e.g., [`ConstantGustinessSpec`](@ref)).
- `moisture_model`: The moisture model (e.g., [`MoistModel`](@ref) or [`DryModel`](@ref)).
- `rsl_model`: Roughness sublayer correction model (e.g., [`ExponentialRSL`](@ref)).
  Defaults to [`NoRoughnessSubLayer`](@ref) (standard MOST, no RSL correction).
- `stability_cap`: Cap on the stability parameter in stable conditions (e.g.,
  [`MaxHeatFluxStabilityCap`](@ref)). Defaults to [`NoStabilityCap`](@ref) (standard MOST).
"""
struct SurfaceFluxConfig{
    R <: AbstractRoughnessParams,
    G <: AbstractGustinessSpec,
    M <: AbstractMoistureModel,
    RSL <: AbstractRoughnessSubLayerModel,
    SC <: AbstractStabilityCap,
}
    roughness::R
    gustiness::G
    moisture_model::M
    rsl_model::RSL
    stability_cap::SC
end

function SurfaceFluxConfig(roughness, gustiness)
    return SurfaceFluxConfig(roughness, gustiness, MoistModel())
end

function SurfaceFluxConfig(roughness, gustiness, moisture_model)
    return SurfaceFluxConfig(roughness, gustiness, moisture_model, NoRoughnessSubLayer())
end

function SurfaceFluxConfig(roughness, gustiness, moisture_model, rsl_model)
    return SurfaceFluxConfig(
        roughness,
        gustiness,
        moisture_model,
        rsl_model,
        NoStabilityCap(),
    )
end



const FluxOption{FT} = Union{Nothing, FT}

"""
    FluxSpecs{FT}

Container for prescribed surface flux boundary conditions.

# Fields
- `shf`: Sensible Heat Flux [W/m^2].
- `lhf`: Latent Heat Flux [W/m^2].
- `ustar`: Friction velocity [m/s].
- `Cd`: Momentum exchange coefficient.
- `Ch`: Heat exchange coefficient.
"""
Base.@kwdef struct FluxSpecs{
    FT,
    A <: FluxOption{FT},
    B <: FluxOption{FT},
    C <: FluxOption{FT},
    D <: FluxOption{FT},
    E <: FluxOption{FT},
}
    shf::A = nothing
    lhf::B = nothing
    ustar::C = nothing
    Cd::D = nothing
    Ch::E = nothing
end

function FluxSpecs{FT}(;
    shf::A = nothing,
    lhf::B = nothing,
    ustar::C = nothing,
    Cd::D = nothing,
    Ch::E = nothing,
) where {FT, A, B, C, D, E}
    return FluxSpecs{FT, A, B, C, D, E}(shf, lhf, ustar, Cd, Ch)
end

"""
    SolverOptions{FT}

Options for the Monin-Obukhov similarity theory solver.

# Fields
- `tol`: Absolute tolerance on the stability parameter: the `converged` flag requires the
  final bracket width, or the step from the last iterate to the final regula falsi
  interpolant, to satisfy it, and in tolerance-checked mode it also bounds the step
  between iterates for the early exit.
- `rtol`: Relative tolerance on the stability parameter, used analogously to `tol`.
- `maxiter`: Number of bracket-refinement iterations. The ζ-solve performs
  `5 + maxiter` residual evaluations in total (branch detection + bracketing
  probes + refinement); see the internal `solve_stability_param`.
- `forced_fixed_iters`: If true (default), disables the early tolerance exit and runs
  exactly `maxiter` refinement iterations (via `RootSolvers.NoTolerance`), so every point
  performs identical work (uniform control flow on GPUs). The `converged` flag is still
  evaluated from the final bracket (width and final step) and the tolerances.
"""
Base.@kwdef struct SolverOptions{FT}
    tol::FT = FT(1e-2)
    rtol::FT = FT(1e-2)
    maxiter::Int = 7
    forced_fixed_iters::Bool = true
end


"""
    SurfaceFluxConditions{FT}

Surface flux conditions returned by [`surface_fluxes`](@ref).

All floating-point fields share the type `FT`, obtained by promoting the inputs.
Momentum-flux components are the kinematic stress times density, i.e. `ρτ = ρ u_* u_*`,
with units `[kg/(m·s²)] = [N/m²]`.

# Fields
- `shf`: Sensible heat flux [W/m²].
- `lhf`: Latent heat flux [W/m²].
- `evaporation`: Evaporation rate [kg/(m²·s)].
- `ρτxz`: Momentum flux, eastward component [kg/(m·s²)].
- `ρτyz`: Momentum flux, northward component [kg/(m·s²)].
- `ustar`: Friction velocity [m/s].
- `ζ`: Monin-Obukhov stability parameter `(z - d)/L` [-].
- `Cd`: Momentum exchange (drag) coefficient [-].
- `g_h`: Heat conductance `Ch * U_eff` [m/s].
- `T_sfc`: Surface temperature [K].
- `q_vap_sfc`: Surface air vapor specific humidity [kg/kg].
- `L_MO`: Monin-Obukhov length [m].
- `L_eff`: Effective Obukhov length for profile recovery, `Δz_eff / min(ζ, ζ_cap)` [m].
  It equals `L_MO` unless a stability cap (see [`MaxHeatFluxStabilityCap`](@ref)) is
  active, in which case the exchange coefficients and similarity scales were evaluated at
  the capped stability parameter. Pass `L_eff` (not `L_MO`) to
  [`compute_profile_value`](@ref) to recover profiles consistent with the fluxes.
- `converged`: Solver convergence status.

The positional constructor accepts the fields in this order, with or without `L_eff`
(without it, `L_eff = L_MO`).
"""
struct SurfaceFluxConditions{FT <: Real}
    shf::FT
    lhf::FT
    evaporation::FT
    ρτxz::FT
    ρτyz::FT
    ustar::FT
    ζ::FT
    Cd::FT
    g_h::FT
    T_sfc::FT
    q_vap_sfc::FT
    L_MO::FT
    L_eff::FT
    converged::Bool
end

SurfaceFluxConditions(
    shf,
    lhf,
    E,
    ρτxz,
    ρτyz,
    ustar,
    ζ,
    Cd,
    g_h,
    T_sfc,
    q_vap_sfc,
    L_MO,
    L_eff,
    converged,
) =
    let vars =
            promote(
                shf,
                lhf,
                E,
                ρτxz,
                ρτyz,
                ustar,
                ζ,
                Cd,
                g_h,
                T_sfc,
                q_vap_sfc,
                L_MO,
                L_eff,
            )
        SurfaceFluxConditions{eltype(vars)}(vars..., converged)
    end

# Without an effective Obukhov length (no stability cap): L_eff = L_MO
SurfaceFluxConditions(
    shf,
    lhf,
    E,
    ρτxz,
    ρτyz,
    ustar,
    ζ,
    Cd,
    g_h,
    T_sfc,
    q_vap_sfc,
    L_MO,
    converged::Bool,
) = SurfaceFluxConditions(
    shf, lhf, E, ρτxz, ρτyz, ustar, ζ, Cd, g_h, T_sfc, q_vap_sfc, L_MO, L_MO, converged,
)

function Base.show(io::IO, sfc::SurfaceFluxConditions)
    println(io, "----------------------- SurfaceFluxConditions")
    println(io, "Sensible Heat Flux                  = ", sfc.shf)
    println(io, "Latent Heat Flux                    = ", sfc.lhf)
    println(io, "Evaporation rate                    = ", sfc.evaporation)
    println(io, "Momentum Flux (x)                   = ", sfc.ρτxz)
    println(io, "Momentum Flux (y)                   = ", sfc.ρτyz)
    println(io, "Friction velocity u⋆                = ", sfc.ustar)
    println(io, "Obukhov stability ζ                 = ", sfc.ζ)
    println(io, "C_drag                              = ", sfc.Cd)
    println(io, "Heat conductance                    = ", sfc.g_h)
    println(io, "Surface temperature                 = ", sfc.T_sfc)
    println(io, "Surface air vapor specific humidity = ", sfc.q_vap_sfc)
    println(io, "Monin-Obukhov length                = ", sfc.L_MO)
    println(io, "Effective Obukhov length            = ", sfc.L_eff)
    println(io, "Converged                           = ", sfc.converged)
    println(io, "-----------------------")
end
