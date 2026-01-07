
"""
    Surface flux configuration specs
"""

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
``U_{gust} = \\beta w_*^{1/3}``
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
- `tol`: Absolute tolerance on the change in the stability parameter for determining convergence.
- `maxiter`: Maximum number of iterations.
- `forced_fixed_iters`: If true, disables the tolerance check and forces the solver to run for exactly `maxiter` iterations (or until machine precision is reached/bypassed). Default is `true`.
"""
Base.@kwdef struct SolverOptions{FT}
    tol::FT = FT(1e-2)
    maxiter::Int = 10
    forced_fixed_iters::Bool = true
end

"""
    SurfaceFluxInputs

Immutable container describing the atmospheric and surface state using primitive
quantities plus module-defined parameterizations. Instances of this type are
passed to the functional surface flux solver.

# Fields
- `T_int`, `q_tot_int`, `ρ_int`: Interior air temperature [K], specific humidity [kg/kg], and density [kg/m^3].
- `q_liq_int`, `q_ice_int`: Interior liquid and ice specific humidity [kg/kg].
- `T_sfc_guess`, `q_vap_sfc_guess`: Initial guesses for surface temperature [K] and vapor specific humidity [kg/kg]. Optional, can be `nothing` for default fallback.
- `Φ_sfc`: Surface geopotential [m^2/s^2].
- `Δz`: Height difference between interior and surface reference levels [m].
- `d`: Displacement height [m].
- `u_int`, `u_sfc`: Horizontal wind components (u, v) at interior and surface levels [m/s].
- `roughness_model`: Roughness parameterization (e.g., [`ConstantRoughnessParams`](@ref)).
- `gustiness_model`: Gustiness parameterization (e.g., [`ConstantGustinessSpec`](@ref)).
- `moisture_model`: Moisture model (e.g., [`MoistModel`](@ref) or [`DryModel`](@ref)).
- `roughness_inputs`: Optional inputs for roughness models.
- `update_T_sfc`, `update_q_vap_sfc`: Optional callbacks to update surface state during iteration.
- `shf`, `lhf`, `ustar`, `Cd`, `Ch`: Optional prescribed flux/scale quantities supplied via [`FluxSpecs`](@ref).
"""
struct SurfaceFluxInputs{
    FT,
    RM <: AbstractRoughnessParams,
    GM <: AbstractGustinessSpec,
    MM <: AbstractMoistureModel,
    R,
    F1,
    F2,
    SHF,
    LHF,
    UST,
    CD,
    CH,
}
    T_int::FT
    q_tot_int::FT
    q_liq_int::FT
    q_ice_int::FT
    ρ_int::FT
    T_sfc_guess::Union{FT, Nothing}
    q_vap_sfc_guess::Union{FT, Nothing}
    Φ_sfc::FT
    Δz::FT
    d::FT
    u_int::Tuple{FT, FT}
    u_sfc::Tuple{FT, FT}
    roughness_model::RM
    gustiness_model::GM
    moisture_model::MM
    roughness_inputs::R
    update_T_sfc::Union{F1, Nothing}
    update_q_vap_sfc::Union{F2, Nothing}
    shf::SHF
    lhf::LHF
    ustar::UST
    Cd::CD
    Ch::CH
end

"""
    SurfaceFluxConditions

Surface flux conditions, returned from `surface_fluxes`.

- `shf::FT`: Sensible heat flux [W/m²]
- `lhf::FT`: Latent heat flux [W/m²]
- `evaporation::FT`: Evaporation rate [kg/(m²·s)]
- `ρτxz::FT`: Momentum flux, eastward component [kg/(m·s²)]
- `ρτyz::FT`: Momentum flux, northward component [kg/(m·s²)]
- `ustar::FT`: Friction velocity [m/s]
- `ζ::FT`: Monin-Obukhov stability parameter (z/L)
- `Cd::FT`: Momentum exchange coefficient
- `Ch::FT`: Heat exchange coefficient
- `T_sfc::FT`: Surface temperature [K]
- `q_vap_sfc::FT`: Surface air vapor specific humidity [kg/kg]
- `L_MO::FT`: Monin-Obukhov lengthscale [m]
- `converged::Bool`: Solver convergence status
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
    Ch::FT
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
    Ch,
    T_sfc,
    q_vap_sfc,
    L_MO,
    converged,
) =
    let vars = promote(shf, lhf, E, ρτxz, ρτyz, ustar, ζ, Cd, Ch, T_sfc, q_vap_sfc, L_MO)
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
    println(io, "C_heat                              = ", sfc.Ch)
    println(io, "Surface temperature                 = ", sfc.T_sfc)
    println(io, "Surface air vapor specific humidity = ", sfc.q_vap_sfc)
    println(io, "Monin-Obukhov length                = ", sfc.L_MO)
    println(io, "Converged                           = ", sfc.converged)
    println(io, "-----------------------")
end
