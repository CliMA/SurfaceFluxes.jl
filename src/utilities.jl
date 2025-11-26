function non_zero(v::FT) where {FT}
    sign_of_v = v == 0 ? 1 : sign(v)
    return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
end

@inline resolve_quantity(value::Number, ctx) = value
@inline resolve_quantity(::Nothing, ctx) = nothing
@inline resolve_quantity(spec::SurfaceScalar, ctx) = spec.value
@inline function resolve_quantity(spec::SurfaceCallable, ctx)
    return spec.fn(ctx)
end
@inline function resolve_quantity(callable, ctx)
    return callable(ctx)
end

@inline gustiness_value(model::ConstantGustiness, inputs, ctx) = model.value
@inline function gustiness_value(model::FunctionalGustiness, inputs, ctx)
    return model.fn(inputs, ctx)
end

@inline z_in(inputs::SurfaceFluxInputs) = inputs.d + inputs.Δz
@inline z_sfc(inputs::SurfaceFluxInputs) = inputs.d
@inline Δz(inputs::SurfaceFluxInputs) = inputs.Δz

@inline function Δu_components(inputs::SurfaceFluxInputs)
    return (
        inputs.u_in[1] - inputs.u_sfc[1],
        inputs.u_in[2] - inputs.u_sfc[2],
    )
end

@inline function windspeed(Δu::NTuple{2, FT}, gustiness::FT) where {FT}
    return max(hypot(Δu[1], Δu[2]), gustiness)
end

@inline function windspeed(inputs::SurfaceFluxInputs, gustiness::FT) where {FT}
    return windspeed(Δu_components(inputs), gustiness)
end

@inline function interior_geopotential(param_set::APS, inputs::SurfaceFluxInputs)
    return inputs.Φs + SFP.grav(param_set) * inputs.Δz
end

@inline surface_geopotential(inputs::SurfaceFluxInputs) = inputs.Φs
