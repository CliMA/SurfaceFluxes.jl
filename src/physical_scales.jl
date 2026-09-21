"""
    compute_physical_scale_coeff(
        param_set::APS,
        Δz_eff,
        ζ,
        z0,
        transport,
        scheme::SolverScheme,
        rsl_model = NoRoughnessSubLayer(),
    )

Compute the coefficient relating a bulk difference to its similarity scale.

Returns `ϕ` such that `scale = Δvalue * ϕ`; for example, `u★ = ΔU * ϕ_m`. It is given by

```math
ϕ = \\frac{κ}{\\hat{F}(Δz_{eff}, ζ, z_0)}
```

where `κ` is the von Kármán constant and `F̂ = F + P` is the RSL-corrected dimensionless
profile (`P ≤ 0` from [`rsl_profile_correction`](@ref), zero when `rsl_model` is
[`NoRoughnessSubLayer`](@ref)).

# Arguments
- `param_set`: Parameter set.
- `Δz_eff`: Effective aerodynamic height `Δz - d` [m].
- `ζ`: Monin-Obukhov stability parameter [-].
- `z0`: Roughness length for the transported variable [m].
- `transport`: Transport type (`MomentumTransport` or `HeatTransport`).
- `scheme`: Discretization scheme ([`PointValueScheme`](@ref) or [`LayerAverageScheme`](@ref)).
- `rsl_model`: Optional roughness sublayer model (default: [`NoRoughnessSubLayer`](@ref)).
"""
function compute_physical_scale_coeff(
    param_set::APS,
    Δz_eff,
    ζ,
    z0,
    transport,
    scheme::SolverScheme,
    rsl_model = NoRoughnessSubLayer(),
)
    κ = SFP.von_karman_const(param_set)
    uf = SFP.uf_params(param_set)

    profile = UF.dimensionless_profile(uf, Δz_eff, ζ, z0, transport, scheme)
    P = rsl_profile_correction(rsl_model, Δz_eff, z0, transport)
    return κ / (profile + P)
end

"""
    compute_ustar(param_set, ζ, z0, inputs, scheme, gustiness)

Return the friction velocity implied by the Monin-Obukhov solution.

If a friction velocity is prescribed via `inputs.ustar` (in the inputs container),
it is returned directly; otherwise it is recomputed from the similarity coefficients.

# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `z0`: Momentum roughness length [m].
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `scheme`: Discretization scheme.
- `gustiness`: Gustiness velocity scale [m/s].
"""
function compute_ustar(
    param_set::APS,
    ζ,
    z0,
    inputs,
    scheme::SolverScheme,
    gustiness,
)
    # Per-input ustar check
    if inputs.ustar !== nothing
        return inputs.ustar
    end

    ΔU = windspeed(inputs, gustiness)
    Δz_eff = effective_height(inputs)
    ϕ = compute_physical_scale_coeff(
        param_set,
        Δz_eff,
        ζ,
        z0,
        UF.MomentumTransport(),
        scheme,
        inputs.rsl_model,
    )
    return ΔU * ϕ
end

"""
    compute_theta_star(param_set, ζ, z0h, inputs, scheme, T_sfc)

Return the potential temperature scale implied by the Monin-Obukhov solution, 
where `z0h` is the roughness length for heat.


# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `z0h`: Thermal roughness length [m].
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `scheme`: Discretization scheme.
- `T_sfc`: Surface temperature [K]. Optional; defaults to `inputs.T_sfc_guess`, falling back
  to the interior temperature `inputs.T_int` when the guess is `nothing`.
"""
function compute_theta_star(
    param_set::APS,
    ζ,
    z0h,
    inputs,
    scheme::SolverScheme,
    T_sfc = something(inputs.T_sfc_guess, inputs.T_int),
)
    thermo_params = SFP.thermodynamics_params(param_set)
    Φ_int = interior_geopotential(param_set, inputs)
    Φ_sfc = surface_geopotential(inputs)

    DSE_int = TD.dry_static_energy(thermo_params, inputs.T_int, Φ_int)
    DSE_sfc = TD.dry_static_energy(thermo_params, T_sfc, Φ_sfc)
    ΔDSE = DSE_int - DSE_sfc

    c_p = TD.cp_m(thermo_params, inputs.q_tot_int, inputs.q_liq_int, inputs.q_ice_int)
    Δθ = ΔDSE / c_p

    Δz_eff = effective_height(inputs)
    ϕ = compute_physical_scale_coeff(
        param_set,
        Δz_eff,
        ζ,
        z0h,
        UF.HeatTransport(),
        scheme,
    )
    return Δθ * ϕ
end

"""
    compute_q_star(param_set, ζ, z0h, inputs, scheme, q_vap_sfc)

Return the specific humidity scale implied by the current Monin-Obukhov solution, 
where `z0h` is the roughness length for scalars (assumed equal to heat).


# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `z0h`: Thermal/scalar roughness length [m].
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `scheme`: Discretization scheme.
- `q_vap_sfc`: Surface vapor specific humidity [kg/kg]. Optional; defaults to
  `inputs.q_vap_sfc_guess`, falling back to the interior total specific humidity
  `inputs.q_tot_int` when the guess is `nothing`.
"""
function compute_q_star(
    param_set::APS,
    ζ,
    z0h,
    inputs,
    scheme::SolverScheme,
    q_vap_sfc = something(inputs.q_vap_sfc_guess, inputs.q_tot_int),
)
    # Δq = q_vap_int - q_vap_sfc
    q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
    Δq = q_vap_int - q_vap_sfc

    # Scalars use HeatTransport coefficients in MOST
    Δz_eff = effective_height(inputs)
    ϕ = compute_physical_scale_coeff(
        param_set,
        Δz_eff,
        ζ,
        z0h,
        UF.HeatTransport(),
        scheme,
    )
    return Δq * ϕ
end

"""
    surface_tke(param_set, Δz_eff, ustar, ζ)

Compute the surface-layer turbulent kinetic energy (TKE) `(u_* ϕ)^2` following Tan et al. (2018).

Returns `(u_* ϕ)^2` [m²/s²], where `ϕ = sqrt(TKE)/u_*` is the TKE-based velocity
similarity function. In unstable conditions this is
`TKE = 3.75 u_*^2 + 0.2 w_*^2 + u_*^2 (-ζ)^{2/3}`, and in stable/neutral conditions
it reduces to `3.75 u_*^2`. The convective (Deardorff) velocity scale `w_*` is computed
from the mixed-layer height `zi` (a fixed parameter) and the buoyancy flux implied by `ζ`.

!!! note
    This returns the full TKE, not the streamwise velocity variance `σ_u^2`. The streamwise
    similarity function `ϕ_σu = σ_u / u_*` (Panofsky et al. 1977) is available via
    `phi(uf, ζ, MomentumVariance())`, from which `σ_u^2 = (u_* ϕ_σu)^2`.

!!! warning "Range of validity"
    This closure is **independent of the flux-profile parameterization** in `param_set`
    (`Businger`/`Gryanik`/`Grachev` give the same result): Grachev et al. (2007) and Gryanik
    et al. (2020) define no variance functions. It is a convective surface-layer / surface-BC
    form; on the stable side it returns the constant `3.75 u_*^2`, which is not a validated
    stable-boundary-layer result. Use with care in stably stratified conditions.

# Arguments
- `param_set`: Parameter set.
- `Δz_eff`: Effective aerodynamic height `Δz - d` [m].
- `ustar`: Friction velocity [m/s].
- `ζ`: Monin-Obukhov stability parameter [-].
"""
function surface_tke(param_set::APS, Δz_eff, ustar, ζ)
    uf = SFP.uf_params(param_set)
    zi = SFP.gustiness_zi(param_set) # Mixed-layer height taken to be fixed

    κ = SFP.von_karman_const(param_set)

    FT = eltype(ustar)

    # Calculate convective velocity scale w_* for unstable conditions.
    # w_* = (B * zi)^(1/3)
    # B = -u_*^3 / (κ * L_MO) = -u_*^3 * ζ / (κ * Δz_eff)
    # => w_* = u_* * ( (zi * ζ) / (-κ * Δz_eff) )^(1/3)
    #
    # We use ifelse to avoid branching and potential warp divergence on GPUs.
    term = -(zi * ζ) / (κ * Δz_eff)
    w_star = ifelse(ζ < 0, ustar * cbrt(term), FT(0))

    ϕ = UF.phi(uf, ζ, ustar, w_star, UF.MomentumVariance())
    return (ustar * ϕ)^2
end

"""
    scalar_variance(param_set, scale, ζ)

Compute the scalar variance `σ_s^2 = (scale * ϕ_σs)^2`, using the temperature-variance
similarity `ϕ_σs = ϕ_σθ` (Wyngaard et al. 1971; Tan et al. 2018).

!!! warning "Range of validity"
    As for [`surface_tke`](@ref), this closure is **independent of the flux-profile
    parameterization** (Grachev/Gryanik define no variance functions) and returns the constant
    `2.0` on the stable side. The constant has some support in the very stable (z-less) limit
    but is not calibrated to stable-boundary-layer data.

# Arguments
- `param_set`: Parameter set.
- `scale`: Similarity scale of the scalar (e.g., `theta_star`, `q_star`).
- `ζ`: Monin-Obukhov stability parameter.
"""
function scalar_variance(param_set::APS, scale, ζ)
    uf = SFP.uf_params(param_set)

    transport = UF.HeatVariance()
    ϕ = UF.phi(uf, ζ, transport)
    return (scale * ϕ)^2
end

"""
    theta_variance(param_set, inputs, shf, ustar, ζ, rho_sfc)

Computes potential temperature variance from sensible heat flux `shf`.
Calculates `θ_* = -shf / (ρ * c_p * u_*)` and calls `scalar_variance`.

# Arguments
- `param_set`: Parameter set.
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `shf`: Sensible heat flux [W/m^2].
- `ustar`: Friction velocity [m/s].
- `ζ`: Monin-Obukhov stability parameter.
- `rho_sfc`: Surface density [kg/m^3].
"""
function theta_variance(param_set::APS, inputs, shf, ustar, ζ, rho_sfc)
    thermo_params = SFP.thermodynamics_params(param_set)

    c_p = TD.cp_m(thermo_params, inputs.q_tot_int, inputs.q_liq_int, inputs.q_ice_int)

    # Calculate scale θ_*
    # Flux convention: positive upward.
    theta_star = -shf / (rho_sfc * c_p * ustar)

    return scalar_variance(param_set, theta_star, ζ)
end

"""
    obukhov_length(param_set, ustar, buoy_flux)

Computes the Monin-Obukhov length [m].

Returns zero if `ustar` is zero.

# Arguments
- `param_set`: Parameter set.
- `ustar`: Friction velocity [m/s].
- `buoy_flux`: Surface buoyancy flux [m^2/s^3].
"""
@inline function obukhov_length(param_set::APS, ustar, buoy_flux)
    FT = eltype(param_set)
    κ = SFP.von_karman_const(param_set)
    L_MO_raw = -ustar^3 / (κ * non_zero(buoy_flux))
    return ifelse(ustar > 0, L_MO_raw, zero(FT))
end

"""
    obukhov_stability_parameter(param_set, Δz_eff, ustar, buoy_flux)

Computes the Monin-Obukhov stability parameter `ζ = Δz_eff / L_MO`, where
`Δz_eff` is the effective aerodynamic height (\$z-d\$).

Returns zero if `ustar` and hence `L_MO` are zero.

# Arguments
- `param_set`: Parameter set.
- `Δz_eff`: Effective aerodynamic height [m].
- `ustar`: Friction velocity [m/s].
- `buoy_flux`: Surface buoyancy flux [m^2/s^3].
"""
@inline function obukhov_stability_parameter(param_set::APS, Δz_eff, ustar, buoy_flux)
    FT = eltype(param_set)
    L_MO = obukhov_length(param_set, ustar, buoy_flux)
    return ifelse(L_MO != 0, Δz_eff / L_MO, zero(FT))
end
