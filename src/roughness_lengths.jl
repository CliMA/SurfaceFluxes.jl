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
    roughness_lengths(momentum; scalar=momentum)

Convenience constructor that pairs the momentum and scalar roughness
specifications.
"""
roughness_lengths(momentum; scalar = momentum) = FixedRoughnessSpec(momentum, scalar)

@inline function default_roughness_lengths(::APS{FT}) where {FT}
    return RoughnessLengths{FT}(FT(1e-3), FT(1e-4))
end

struct FixedRoughnessSpec{TM <: Real, TS <: Real} <: AbstractRoughnessSpec
    momentum::TM
    scalar::TS
end

FixedRoughnessSpec(momentum::Real; scalar = momentum) = FixedRoughnessSpec(momentum, scalar)

struct CharnockRoughnessSpec{TA <: Real, TS <: Real} <: AbstractRoughnessSpec
    α::TA
    scalar::TS
end


SurfaceFluxConfig(; roughness = DefaultRoughnessSpec(), gustiness = ConstantGustinessSpec(1.0)) =
    SurfaceFluxConfig(roughness, gustiness)

struct FixedRoughnessModel{FT} <: AbstractRoughnessModel{FT}
    momentum::FT
    scalar::FT
end

struct CharnockRoughnessModel{FT} <: AbstractRoughnessModel{FT}
    α::FT
    scalar::FT
    grav::FT
end

@inline function instantiate_roughness(param_set::APS{FT}, ::DefaultRoughnessSpec) where {FT}
    defaults = default_roughness_lengths(param_set)
    return FixedRoughnessModel{FT}(defaults.momentum, defaults.scalar)
end

@inline function instantiate_roughness(param_set::APS{FT}, spec::FixedRoughnessSpec) where {FT}
    return FixedRoughnessModel{FT}(convert(FT, spec.momentum), convert(FT, spec.scalar))
end

@inline function instantiate_roughness(param_set::APS{FT}, spec::CharnockRoughnessSpec) where {FT}
    return CharnockRoughnessModel{FT}(convert(FT, spec.α), convert(FT, spec.scalar), SFP.grav(param_set))
end

"""
    charnock_momentum(; α = 0.011, scalar = 1e-4)

Convenience helper returning a roughness specification that applies the
Charnock relation for the momentum roughness length while keeping a scalar
roughness length fixed.
"""
@inline function charnock_momentum(; α = 0.011, scalar = 1e-4)
    return CharnockRoughnessSpec(α, scalar)
end
@inline momentum_roughness(model::FixedRoughnessModel{FT}, u★, sfc_param_set, ctx, roughness_inputs) where {FT} =
    model.momentum

@inline scalar_roughness(model::FixedRoughnessModel{FT}, u★, sfc_param_set, ctx, roughness_inputs) where {FT} =
    model.scalar

@inline function momentum_roughness(model::CharnockRoughnessModel{FT}, u★, sfc_param_set, ctx, roughness_inputs) where {FT}
    return model.α * u★^2 / model.grav
end

@inline scalar_roughness(model::CharnockRoughnessModel{FT}, u★, sfc_param_set, ctx, roughness_inputs) where {FT} =
    model.scalar

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
