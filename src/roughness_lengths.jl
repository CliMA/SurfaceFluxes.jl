# Roughness-length evaluation
#
# These helpers resolve the active momentum and scalar roughness specifications.
# They are intentionally lightweight so they can be inlined inside GPU kernels.
#

"""
    RoughnessLengths{FT}

Struct to hold momentum and scalar roughness lengths.
"""
struct RoughnessLengths{FT}
    momentum::FT
    scalar::FT
end

"""
    ConstantRoughnessSpec{FT} <: AbstractRoughnessSpec

Roughness lengths fixed to constant values.

# Fields
- `z0m`: Momentum roughness length [m]
- `z0s`: Scalar roughness length [m]
"""
Base.@kwdef struct ConstantRoughnessSpec{FT} <: AbstractRoughnessSpec
    z0m::FT = 1e-3
    z0s::FT = 1e-4
end

"""
    CharnockRoughnessSpec{FT} <: AbstractRoughnessSpec

Charnock relationship for momentum roughness; scalar roughness fixed.

# Fields
- `α`: Charnock coefficient
- `z0s`: Scalar roughness length [m]
"""
Base.@kwdef struct CharnockRoughnessSpec{FT} <: AbstractRoughnessSpec
    α::FT = 0.011
    z0s::FT = 1e-4
end

"""
    RaupachRoughnessSpec <: AbstractRoughnessSpec

Raupach (1994) canopy roughness model.
"""
Base.@kwdef struct RaupachRoughnessSpec{FT} <: AbstractRoughnessSpec
    C_R::FT = 0.3
    C_S::FT = 0.003
    c_d1::FT = 7.5
end

SurfaceFluxConfig(; roughness = ConstantRoughnessSpec(), gustiness = ConstantGustinessSpec(1.0)) =
    SurfaceFluxConfig(roughness, gustiness)

"""
    roughness_lengths(z0m; scalar)

Helper to construct `ConstantRoughnessSpec`.
"""
function roughness_lengths(z0m; scalar)
    return ConstantRoughnessSpec(z0m = z0m, z0s = scalar)
end

"""
    charnock_momentum(; α, scalar)

Helper to construct `CharnockRoughnessSpec`.
"""
function charnock_momentum(; α, scalar)
    return CharnockRoughnessSpec(α = α, z0s = scalar)
end

# Accessors
@inline function momentum_roughness(spec::ConstantRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    return spec.z0m
end

@inline function scalar_roughness(spec::ConstantRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    return spec.z0s
end

@inline function momentum_roughness(spec::CharnockRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    grav = SFP.grav(sfc_param_set)
    return spec.α * u★^2 / grav
end

@inline function scalar_roughness(spec::CharnockRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    return spec.z0s
end

"""
    momentum_roughness(spec::RaupachRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)

Calculate momentum roughness length using the Raupach (1994) canopy roughness model.
Requires `roughness_inputs` to contain `LAI` (Leaf Area Index) and `h` (canopy height).
"""
@inline function momentum_roughness(spec::RaupachRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    FT = float_type(sfc_param_set)
    k = SFP.von_karman_const(sfc_param_set)
    
    LAI = roughness_inputs.LAI
    h = roughness_inputs.h
    
    # Roughness density lambda (frontal area index) approximated as LAI / 2
    λ = LAI / FT(2)
    
    C_R = spec.C_R
    C_S = spec.C_S
    c_d1 = spec.c_d1
    
    # Avoid division by zero if LAI is extremely small
    if λ <= eps(FT)
        return SFP.z0m_fixed(sfc_param_set)
    end
    
    # d/h (Eq 15 in Raupach 1994)
    sqrt_c_lambda = sqrt(c_d1 * λ)
    d_over_h = FT(1) - (FT(1) - exp(-sqrt_c_lambda)) / sqrt_c_lambda
    
    # u_star / U(h) (Eq 8 in Raupach 1994)
    ustar_over_Uh = sqrt(C_S + C_R * λ)
    Uh_over_ustar = FT(1) / ustar_over_Uh
    
    # Psi_h (Roughness sublayer influence function), approximated as 0.193
    Ψ_h = FT(0.193)
    
    z0m = h * (FT(1) - d_over_h) * exp(-k * Uh_over_ustar - Ψ_h)
    
    return z0m
end

@inline function scalar_roughness(::RaupachRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    return SFP.z0s_fixed(sfc_param_set) # Fallback to fixed scalar roughness
end

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.MomentumTransport,
    ctx,
)
    return momentum_roughness(inputs.roughness_model, u★, sfc_param_set, ctx, inputs.roughness_inputs)
end

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.HeatTransport,
    ctx,
)
    return scalar_roughness(inputs.roughness_model, u★, sfc_param_set, ctx, inputs.roughness_inputs)
end
