"""
    compute_physical_scale_coeff(
        param_set::APS,
        inputs,
        L_MO,
        z0,
        transport,
        scheme::SolverScheme,
    )

Computes the coefficient for the physical scale of a variable.
Returns `ϕ` such that `scale = Δvalue * ϕ`.
For example, `u★ = ΔU * ϕ_m`.

This is computed as:
```math
ϕ = \\frac{κ}{F(z, L, z_0)}
```
where `F` is the dimensionless profile function from `UniversalFunctions`.
"""
function compute_physical_scale_coeff(
    param_set::APS,
    inputs,
    L_MO,
    z0,
    transport,
    scheme::SolverScheme,
)
    κ = SFP.von_karman_const(param_set)
    uf = SFP.uf_params(param_set)
    Δz_layer = Δz(inputs)
    ζ = Δz_layer / L_MO
    
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
    compute_ustar(param_set, L_MO, z0, inputs, scheme, gustiness)

Return the friction velocity implied by the current Monin-Obukhov solution.
If a friction velocity is prescribed via `SurfaceFluxInputs`, it is returned
directly; otherwise it is recomputed from the similarity coefficients.
"""
function compute_ustar(
    param_set::APS,
    L_MO,
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
        inputs,
        L_MO,
        z0,
        UF.MomentumTransport(),
        scheme,
    )
    return ΔU * ϕ
end

"""
    u_variance(param_set, inputs, ustar, L_MO)

Compute the velocity variance `σ_u^2 = u_star^2 * ϕ_σu^2`.
For unstable conditions, the convective velocity scale `w_*` is calculated
using the mixed-layer height `zi` from parameters, and passed to the
universal function.
"""
function u_variance(param_set::APS, inputs, ustar, L_MO)
    uf = SFP.uf_params(param_set)
    zi = SFP.gustiness_zi(param_set) # Mixed-layer height taken to be fixed

    κ = SFP.von_karman_const(param_set)
    Δz_layer = Δz(inputs)
    ζ = Δz_layer / L_MO 

    FT = eltype(ustar)
    w_star = FT(0)
    if L_MO < 0
        # Calculate convective velocity scale w_*
        # w_* = (B * zi)^(1/3) where B = -u_*^3 / (κ * L)
        term = zi / (-κ * L_MO)
        w_star = ustar * cbrt(term)
    end

    ϕ = UF.phi(uf, ζ, ustar, w_star, UF.MomentumVariance())
    return (ustar * ϕ)^2
end

"""
    scalar_variance(param_set, inputs, scale, L_MO)

Compute the scalar variance `σ_c^2 = c_*^2 * ϕ_σc^2`.
"""
function scalar_variance(param_set::APS, inputs, scale, L_MO)
    uf = SFP.uf_params(param_set)
    Δz_layer = Δz(inputs)
    ζ = Δz_layer / L_MO
    
    transport = UF.HeatVariance()
    ϕ = UF.phi(uf, ζ, transport)
    return (scale * ϕ)^2
end

"""
    theta_variance(param_set, inputs, shf, ustar, L_MO)

Convenience function to compute potential temperature variance from sensible heat flux `shf`.
Calculates `θ_* = -shf / (ρ * c_p * u_*)` and calls `scalar_variance`.
"""
function theta_variance(param_set::APS, inputs, shf, ustar, L_MO)
    thermo_params = SFP.thermodynamics_params(param_set)
    
    #TODO: Use current iteration state for density calculation, 
    # or get it from current iteration state
    FT = eltype(inputs.Tin) 
    ρ = surface_density(
        param_set,
        inputs.Tin,
        inputs.ρin,
        inputs.Ts_guess,
        inputs.qin,
        inputs.ql_in,
        inputs.qi_in,
    )

    c_p = TD.cp_m(thermo_params, inputs.qin, inputs.ql_in, inputs.qi_in)
    
    # Calculate scale θ_*
    # Flux convention: positive upward.
    theta_star = -shf / (ρ * c_p * ustar)
    
    return scalar_variance(param_set, inputs, theta_star, L_MO)
end
