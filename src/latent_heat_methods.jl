function latent_heat_flux(param_set::APS, inputs::SurfaceFluxInputs, E)
    if inputs.lhf !== nothing
        return inputs.lhf
    end
    LH_v0 = SFP.LH_v0(param_set)
    return LH_v0 * E
end
