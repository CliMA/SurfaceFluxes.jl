function non_zero(v::FT) where {FT}
    sign_of_v = v == 0 ? 1 : sign(v)
    return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
end

@inline gustiness_value(spec::ConstantGustinessSpec, inputs, ctx) = spec.value

@inline z_in(inputs::SurfaceFluxInputs) = inputs.d + inputs.Δz
@inline z_sfc(inputs::SurfaceFluxInputs) = inputs.d
@inline Δz(inputs::SurfaceFluxInputs) = inputs.Δz
@inline Δz(inputs::Fluxes) = inputs.state_int.z - inputs.state_sfc.z
@inline Δz(inputs::FluxesAndFrictionVelocity) = inputs.state_int.z - inputs.state_sfc.z
@inline Δz(inputs::ValuesOnly) = inputs.state_int.z - inputs.state_sfc.z

@inline function Δu_components(inputs::SurfaceFluxInputs)
    return (
        inputs.u_int[1] - inputs.u_sfc[1],
        inputs.u_int[2] - inputs.u_sfc[2],
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
