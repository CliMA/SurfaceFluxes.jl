
"""
    Surface flux configuration specs
"""

abstract type AbstractRoughnessParams end
abstract type AbstractGustinessSpec end



struct ConstantGustinessSpec{TG <: Real} <: AbstractGustinessSpec
    value::TG
end

struct DeardorffGustinessSpec <: AbstractGustinessSpec end

Base.broadcastable(p::AbstractRoughnessParams) = tuple(p)
Base.broadcastable(p::AbstractGustinessSpec) = tuple(p)

abstract type AbstractMoistureModel end
struct MoistModel <: AbstractMoistureModel end
struct DryModel <: AbstractMoistureModel end

Base.broadcastable(m::AbstractMoistureModel) = tuple(m)

"""
    SurfaceFluxConfig

Configuration for surface flux calculation components.

# Fields
- `roughness`: The roughness length parameterization to use.
- `gustiness`: The gustiness parameterization to use.
- `moisture_model`: The moisture model (e.g. `MoistModel` or `DryModel`).
"""
struct SurfaceFluxConfig{R <: AbstractRoughnessParams, G <: AbstractGustinessSpec, M <: AbstractMoistureModel}
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
- `shf`: Sensible Heat Flux [W/m^2]
- `lhf`: Latent Heat Flux [W/m^2]
- `ustar`: Friction velocity [m/s]
- `Cd`: Momentum exchange coefficient
- `Ch`: Heat exchange coefficient
"""
Base.@kwdef struct FluxSpecs{FT}
    shf::FluxOption{FT} = nothing
    lhf::FluxOption{FT} = nothing
    ustar::FluxOption{FT} = nothing
    Cd::FluxOption{FT} = nothing
    Ch::FluxOption{FT} = nothing
end

FluxSpecs(::Type{FT}; kwargs...) where {FT} = FluxSpecs{FT}(; kwargs...)

"""
    SolverOptions{FT}

Options for the Monin-Obukhov similarity theory solver.

# Fields
- `tol`: Absolute tolerance for determining convergence.
- `maxiter`: Maximum number of iterations.
"""
Base.@kwdef struct SolverOptions{FT}
    tol::FT = sqrt(eps(FT))
    maxiter::Int = 100
end

SolverOptions(::Type{FT}; kwargs...) where {FT} = SolverOptions{FT}(; kwargs...)

"""
    SurfaceFluxInputs

Immutable container describing the atmospheric and surface state using primitive
quantities plus module-defined parameterizations. Instances of this type are
passed to the functional surface flux solver.

- `T_int`, `q_tot_int`, `ρ_int`: Interior air temperature [K], specific humidity [kg/kg], and density [kg/m³]
- `Ts_guess`, `qs_guess`: Scalar initial guesses for surface temperature and humidity
- `Φs`: Surface geopotential [m²/s²]
- `Δz`: Height difference between interior and surface reference levels [m]
- `d`: Displacement height [m]
- `u_int`, `u_sfc`: Horizontal wind components (u, v) at interior and surface levels [m/s]
- `roughness_model`: Module-defined roughness parameterization
- `gustiness_model`: Module-defined gustiness parameterization
- `update_T_sfc`, `update_q_vap_sfc`: Optional hooks invoked each solver iteration with signature `(ζ, param_set, thermo_params, inputs)`
- `shf`, `lhf`, `ustar`, `Cd`, `Ch`: Optional prescribed flux/scale quantities supplied via `FluxSpecs`
"""
struct SurfaceFluxInputs{
    FT,
    RM <: AbstractRoughnessParams,
    GM <: AbstractGustinessSpec,
    MM <: AbstractMoistureModel,
    R,
    F1,
    F2,
}
    T_int::FT
    q_tot_int::FT
    q_liq_int::FT
    q_ice_int::FT
    ρ_int::FT
    T_sfc_guess::FT
    q_vap_sfc_guess::FT
    Φ_sfc::FT
    Δz::FT
    d::FT
    u_int::Tuple{FT, FT}
    u_sfc::Tuple{FT, FT}
    roughness_model::RM
    gustiness_model::GM
    moisture_model::MM
    roughness_inputs::Union{R, Nothing}
    update_T_sfc::Union{F1, Nothing}
    update_q_vap_sfc::Union{F2, Nothing}
    shf::Union{Nothing, FT}
    lhf::Union{Nothing, FT}
    ustar::Union{Nothing, FT}
    Cd::Union{Nothing, FT}
    Ch::Union{Nothing, FT}
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
- `L_MO::FT`: Monin-Obukhov lengthscale [m]
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
    L_MO::FT
end

SurfaceFluxConditions(shf, lhf, E, ρτxz, ρτyz, ustar, ζ, Cd, Ch, L_MO) =
    let vars = promote(shf, lhf, E, ρτxz, ρτyz, ustar, ζ, Cd, Ch, L_MO)
        SurfaceFluxConditions{eltype(vars)}(vars...)
    end
function Base.show(io::IO, sfc::SurfaceFluxConditions)
    println(io, "----------------------- SurfaceFluxConditions")
    println(io, "Sensible Heat Flux     = ", sfc.shf)
    println(io, "Latent Heat Flux       = ", sfc.lhf)
    println(io, "Evaporation rate       = ", sfc.evaporation)
    println(io, "Momentum Flux (x)      = ", sfc.ρτxz)
    println(io, "Momentum Flux (y)      = ", sfc.ρτyz)
    println(io, "Friction velocity u⋆   = ", sfc.ustar)
    println(io, "Obukhov stability ζ    = ", sfc.ζ)
    println(io, "C_drag                 = ", sfc.Cd)
    println(io, "C_heat                 = ", sfc.Ch)
    println(io, "Monin-Obukhov length   = ", sfc.L_MO)
    println(io, "-----------------------")
end
