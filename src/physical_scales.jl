"""
    compute_physical_scale_coeff(
        param_set::APS,
        Δz_layer,
        ζ,
        z0,
        transport,
        scheme::SolverScheme,
    )

Computes the coefficient for the physical scale of a variable.
Returns `ϕ` such that `scale = Δvalue * ϕ`.
For example, `u★ = ΔU * ϕ_m`.

# Arguments
- `Δz_layer`: Layer height [m]
This is computed as:
```math
ϕ = \\frac{κ}{F(z, ζ, z_0)}
```
where `F` is the dimensionless profile function from `UniversalFunctions`.
"""
function compute_physical_scale_coeff(
    param_set::APS,
    Δz_layer,
    ζ,
    z0,
    transport,
    scheme::SolverScheme,
)
    κ = SFP.von_karman_const(param_set)
    uf = SFP.uf_params(param_set)

    profile = UF.dimensionless_profile(
        uf,
        Δz_layer,
        ζ,
        z0,
        transport,
        scheme,
    )
    return κ / profile
end

"""
    compute_ustar(param_set, ζ, z0, inputs, scheme, gustiness)

Return the friction velocity implied by the current Monin-Obukhov solution.
If a friction velocity is prescribed via `inputs.ustar`, it is returned
directly; otherwise it is recomputed from the similarity coefficients.
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
    ϕ = compute_physical_scale_coeff(
        param_set,
        inputs.Δz,
        ζ,
        z0,
        UF.MomentumTransport(),
        scheme,
    )
    return ΔU * ϕ
end

"""
    u_variance(param_set, Δz_layer, ustar, ζ)

Compute the velocity variance `σ_u^2 = u_star^2 * ϕ_σu^2`.
For unstable conditions, the convective (Deardorff)velocity scale `w_*` 
is calculated using the mixed-layer height `zi` from parameters, and 
passed to the universal function.
"""
function u_variance(param_set::APS, Δz_layer, ustar, ζ)
    uf = SFP.uf_params(param_set)
    zi = SFP.gustiness_zi(param_set) # Mixed-layer height taken to be fixed

    κ = SFP.von_karman_const(param_set)

    # Check for unstable conditions using ζ directly
    FT = eltype(ustar)
    w_star = FT(0)

    # Calculate L_MO back from ζ to compute w_*? 
    # Or refactor w_* calculation to use ζ?
    # w_* = (B * zi)^(1/3)
    # B = -u_*^3 / (κ * L) = -u_*^3 * ζ / (κ * Δz)

    if ζ < 0
        L_MO = Δz_layer / ζ
        # w_* = (B * zi)^(1/3) where B = -u_*^3 / (κ * L_MO)
        term = zi / (-κ * L_MO)
        w_star = ustar * cbrt(term)
    end

    ϕ = UF.phi(uf, ζ, ustar, w_star, UF.MomentumVariance())
    return (ustar * ϕ)^2
end

"""
    scalar_variance(param_set, scale, ζ)

Compute the scalar variance `σ_c^2 = c_*^2 * ϕ_σc^2`.
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
Returns zero if `buoy_flux` is too small (using `non_zero`) or `ustar` is zero.
"""
function obukhov_length(param_set::APS, ustar, buoy_flux)
    FT = eltype(param_set)
    if non_zero(buoy_flux) != 0 && ustar > 0
        κ = SFP.von_karman_const(param_set)
        return -ustar^3 / (κ * buoy_flux)
    else
        return zero(FT)
    end
end

"""
    obukhov_stability_parameter(param_set, Δz, ustar, buoy_flux)

Computes the Monin-Obukhov stability parameter `ζ = Δz / L_MO`.
Returns zero if `buoy_flux` is too small (using `non_zero`) or `ustar` is zero.
"""
function obukhov_stability_parameter(param_set::APS, Δz, ustar, buoy_flux)
    FT = eltype(param_set)
    L_MO = obukhov_length(param_set, ustar, buoy_flux)
    if L_MO != 0
        return Δz / L_MO
    else
        return zero(FT)
    end
end
