"""
    compute_physical_scale_coeff(
        param_set::APS,
        Δz_eff,
        ζ,
        z0,
        transport,
        scheme::SolverScheme,
    )

Computes the coefficient for the physical scale of a variable.
Returns `ϕ` such that `scale = Δvalue * ϕ`. For example, `u★ = ΔU * ϕ_m`.

# Arguments
- `Δz_eff`: Effective aerodynamic height `Δz - d` [m]
This is computed as:
```math
ϕ = \\frac{κ}{F(Δz_eff, ζ, z_0)}
```
where `F` is the dimensionless profile function from `UniversalFunctions`.
"""
function compute_physical_scale_coeff(
    param_set::APS,
    Δz_eff,
    ζ,
    z0,
    transport,
    scheme::SolverScheme,
)
    κ = SFP.von_karman_const(param_set)
    uf = SFP.uf_params(param_set)

    profile = UF.dimensionless_profile(
        uf,
        Δz_eff,
        ζ,
        z0,
        transport,
        scheme,
    )
    return κ / profile
end

"""
    compute_ustar(param_set, ζ, z0, inputs, scheme, gustiness)

Return the friction velocity implied by the Monin-Obukhov solution.

If a friction velocity is prescribed via `inputs.ustar` (in [`SurfaceFluxInputs`](@ref)),
it is returned directly; otherwise it is recomputed from the similarity coefficients.

# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `z0`: Momentum roughness length [m].
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
- `scheme`: Discretization scheme.
- `gustiness`: Gustiness velocity scale [m/s].
"""
function compute_ustar(
    param_set::APS,
    ζ,
    z0,
    inputs::SurfaceFluxInputs,
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
    )
    return ΔU * ϕ
end

"""
    compute_theta_star(param_set, ζ, z0h, inputs, scheme, T_sfc)

Return the potential temperature scale implied by the Monin-Obukhov solution, 
where `z0h` is the roughness length for heat.

See [`SurfaceFluxInputs`](@ref).

# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `z0h`: Thermal roughness length [m].
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
- `scheme`: Discretization scheme.
- `T_sfc`: Surface temperature [K]. Optional, defaults to `inputs.T_sfc_guess`.
"""
function compute_theta_star(
    param_set::APS,
    ζ,
    z0h,
    inputs::SurfaceFluxInputs,
    scheme::SolverScheme,
    T_sfc = inputs.T_sfc_guess,
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

See [`SurfaceFluxInputs`](@ref).

# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `z0h`: Thermal/scalar roughness length [m].
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
- `scheme`: Discretization scheme.
- `q_vap_sfc`: Surface vapor specific humidity [kg/kg]. Optional, defaults to `inputs.q_vap_sfc_guess`.
"""
function compute_q_star(
    param_set::APS,
    ζ,
    z0h,
    inputs::SurfaceFluxInputs,
    scheme::SolverScheme,
    q_vap_sfc = inputs.q_vap_sfc_guess,
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
    u_variance(param_set, Δz_eff, ustar, ζ)

Compute the velocity variance `σ_u^2 = (u_star * ϕ_σu)^2`.

For unstable conditions, the convective (Deardorff) velocity scale `w_*` 
is calculated using the mixed-layer height `zi` from parameters, and 
passed to the universal function. `Δz_eff` is the effective aerodynamic height.

# Arguments
- `param_set`: Parameter set.
- `Δz_eff`: Effective aerodynamic height [m].
- `ustar`: Friction velocity [m/s].
- `ζ`: Monin-Obukhov stability parameter.
"""
function u_variance(param_set::APS, Δz_eff, ustar, ζ)
    uf = SFP.uf_params(param_set)
    zi = SFP.gustiness_zi(param_set) # Mixed-layer height taken to be fixed

    κ = SFP.von_karman_const(param_set)

    # Check for unstable conditions using ζ directly
    FT = eltype(ustar)
    w_star = FT(0)

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

Compute the scalar variance `σ_s^2 = (c_s * ϕ_σs)^2`.

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
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
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
- `buoy_flux`: Surace buoyancy flux [m^2/s^3].
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
- `buoy_flux`: Surace buoyancy flux [m^2/s^3].
"""
@inline function obukhov_stability_parameter(param_set::APS, Δz_eff, ustar, buoy_flux)
    FT = eltype(param_set)
    L_MO = obukhov_length(param_set, ustar, buoy_flux)
    return ifelse(L_MO != 0, Δz_eff / L_MO, zero(FT))
end
