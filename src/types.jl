# Surface flux configuration types: roughness, gustiness, moisture, and the
# containers (`SurfaceFluxConfig`, `FluxSpecs`, `SolverOptions`) that bundle the
# user-facing specifications consumed by `surface_fluxes`.

abstract type AbstractRoughnessParams end
abstract type AbstractGustinessSpec end



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
"""
struct SurfaceFluxConfig{
    R <: AbstractRoughnessParams,
    G <: AbstractGustinessSpec,
    M <: AbstractMoistureModel,
}
    roughness::R
    gustiness::G
    moisture_model::M
end

function SurfaceFluxConfig(roughness, gustiness)
    return SurfaceFluxConfig(roughness, gustiness, MoistModel())
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
  final bracket width to satisfy it, and in tolerance-checked mode it also bounds the step
  between iterates for the early exit.
- `rtol`: Relative tolerance on the stability parameter, used analogously to `tol`.
- `maxiter`: Number of bracket-refinement iterations. The ζ-solve performs
  `5 + maxiter` residual evaluations in total (branch detection + bracketing
  probes + refinement); see the internal `solve_stability_param`.
- `forced_fixed_iters`: If true (default), disables the early tolerance exit and runs
  exactly `maxiter` refinement iterations (via `RootSolvers.NoTolerance`), so every point
  performs identical work (uniform control flow on GPUs). The `converged` flag is still
  evaluated from the final bracket width and the tolerances.
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
- `converged`: Solver convergence status.
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
    converged,
) =
    let vars = promote(shf, lhf, E, ρτxz, ρτyz, ustar, ζ, Cd, g_h, T_sfc, q_vap_sfc, L_MO)
        SurfaceFluxConditions{eltype(vars)}(vars..., converged)
    end

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
    println(io, "Converged                           = ", sfc.converged)
    println(io, "-----------------------")
end
