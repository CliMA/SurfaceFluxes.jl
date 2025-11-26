# Roughness-length evaluation
#
# These helpers resolve the active momentum and scalar roughness specifications,
# regardless of whether they were provided as scalars or module-defined
# callables. They are intentionally lightweight so they can be inlined inside
# GPU kernels.

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.MomentumTransport,
    ctx,
)
    return resolve_quantity(inputs.z0m_spec, ctx)
end

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.HeatTransport,
    ctx,
)
    return resolve_quantity(inputs.z0b_spec, ctx)
end
