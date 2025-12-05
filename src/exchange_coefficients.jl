@inline heat_conductance(inputs::SurfaceFluxInputs, Ch, gustiness) =
    Ch * windspeed(inputs, gustiness)

"""
    drag_coefficient(param_set, L_MO, z0m, Δz)

Compute and return `Cd`, the drag coefficient, for the current similarity state.
"""
function drag_coefficient(
    param_set::APS,
    L_MO,
    z0m,
    Δz,
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)
    ζ = Δz / L_MO

    F_m = UF.dimensionless_profile(uf_params, Δz, ζ, z0m, UF.MomentumTransport())
    Cd = (κ / F_m)^2
    return Cd
end

"""
    heat_exchange_coefficient(param_set, L_MO, z0m, z0b, Δz)

Compute and return `Ch`, the heat exchange coefficient, for the current
similarity state. Neutral and non-neutral regimes follow the log-law and MOST
formulations, respectively.
"""
function heat_exchange_coefficient(
    param_set::APS,
    L_MO,
    z0m,
    z0b,
    Δz,
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)
    ζ = Δz / L_MO

    F_m = UF.dimensionless_profile(uf_params, Δz, ζ, z0m, UF.MomentumTransport())
    F_h = UF.dimensionless_profile(uf_params, Δz, ζ, z0b, UF.HeatTransport())

    Ch = κ^2 / (F_m * F_h)
    return Ch
end
