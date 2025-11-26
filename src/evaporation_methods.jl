function evaporation(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    Ch,
    q_in,
    q_sfc,
    ρ_sfc,
    gustiness::FT,
) where {FT}
    if inputs.lhf !== nothing
        LH_v0 = SFP.LH_v0(param_set)
        return inputs.lhf / LH_v0
    end
    Δqt_val = Δqt(q_in, q_sfc)
    return -ρ_sfc * Ch * windspeed(inputs, gustiness) * Δqt_val
end
