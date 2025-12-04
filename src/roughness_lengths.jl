# Roughness-length evaluation
#
# These helpers resolve the active momentum and scalar roughness specifications.
# They are intentionally lightweight so they can be inlined inside GPU kernels.
#

@inline momentum_roughness(model::FixedRoughnessModel{FT}, u★, sfc_param_set, ctx) where {FT} =
    model.momentum

@inline scalar_roughness(model::FixedRoughnessModel{FT}, u★, sfc_param_set, ctx) where {FT} =
    model.scalar

@inline function momentum_roughness(model::CharnockRoughnessModel{FT}, u★, sfc_param_set, ctx) where {FT}
    return model.α * u★^2 / model.grav
end

@inline scalar_roughness(model::CharnockRoughnessModel{FT}, u★, sfc_param_set, ctx) where {FT} =
    model.scalar

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.MomentumTransport,
    ctx,
)
    return momentum_roughness(inputs.roughness_model, u★, sfc_param_set, ctx)
end

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.HeatTransport,
    ctx,
)
    return scalar_roughness(inputs.roughness_model, u★, sfc_param_set, ctx)
end
