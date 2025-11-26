
function recover_profile(
    param_set::APS,
    L_MO,
    𝓁,
    Z,
    X_star,
    X_sfc,
    transport,
    scheme::Union{LayerAverageScheme, PointValueScheme},
)
    uf = SFP.uf_params(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    num1 = log(Z / 𝓁)
    num2 = -UF.psi(uf, Z / L_MO, transport)
    num3 = UF.psi(uf, 𝓁 / L_MO, transport)
    Σnum = num1 + num2 + num3
    return Σnum * X_star / 𝜅 + X_sfc
end
