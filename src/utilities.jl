"""
    non_zero(v)

Ensure that `v` is not zero, returning `eps(v)` (preserving sign) if `v` is too small.
"""
@inline function non_zero(v)
    FT = typeof(v)
    threshold = eps(FT)
    return ifelse(abs(v) < threshold, copysign(FT(threshold), v), v)
end

"""
    interior_geopotential(param_set, inputs)

Compute the geopotential at the interior (atmospheric) reference level.

# Arguments
- `param_set`: Parameter set containing gravitational constant.
- `inputs`: `SurfaceFluxInputs` struct with `Φ_sfc` and `Δz`.

Returns `Φ_sfc + g * Δz` [m²/s²].
"""
@inline function interior_geopotential(param_set::APS, inputs::SurfaceFluxInputs)
    return inputs.Φ_sfc + SFP.grav(param_set) * inputs.Δz
end

"""
    surface_geopotential(inputs)

Return the surface geopotential from the inputs.

# Arguments
- `inputs`: `SurfaceFluxInputs` struct.

Returns `inputs.Φ_sfc` [m²/s²].
"""
@inline surface_geopotential(inputs::SurfaceFluxInputs) = inputs.Φ_sfc

"""
    surface_density(param_set, T_int, ρ_int, T_sfc, Δz, q_tot_int=0, q_liq_int=0, q_ice_int=0, q_vap_sfc=nothing)

Estimates the surface air density assuming hydrostatic balance between the interior and surface.
It effectively extrapolates the interior pressure to the surface using the hydrostatic 
equation with an average virtual potential temperature, and then computes the surface 
density using the ideal gas law.

# Arguments
- `param_set`: AbstractSurfaceFluxesParameters.
- `T_int`: Interior temperature [K].
- `ρ_int`: Interior density [kg/m^3].
- `T_sfc`: Surface temperature [K].
- `Δz`: Height difference [m].
- `q_tot_int`: Interior total specific humidity.
- `q_liq_int`: Interior liquid specific humidity.
- `q_ice_int`: Interior ice specific humidity.
- `q_vap_sfc`: Surface vapor specific humidity (optional, defaults to `q_vap_int`).

Returns `ρ_sfc` [kg/m^3].
"""
@inline function surface_density(
    param_set::APS,
    T_int,
    ρ_int,
    T_sfc,
    Δz,
    q_tot_int = 0,
    q_liq_int = 0,
    q_ice_int = 0,
    q_vap_sfc = nothing,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)

    # Humidities for hydrostatic extrapolation
    q_liq_sfc = q_liq_int
    q_ice_sfc = q_ice_int
    q_vap_int = q_tot_int - q_liq_int - q_ice_int
    q_vap_sfc = q_vap_sfc === nothing ? q_vap_int : q_vap_sfc
    q_tot_sfc = q_vap_sfc + q_liq_sfc + q_ice_sfc

    # Gas constants 
    R_m_int = TD.gas_constant_air(thermo_params, q_tot_int, q_liq_int, q_ice_int)
    R_m_sfc = TD.gas_constant_air(thermo_params, q_tot_sfc, q_liq_sfc, q_ice_sfc)

    # Take average of R_m * T (correspondong to average virtual temperature)
    R_m_T_avg = (R_m_int * T_int + R_m_sfc * T_sfc) / 2

    # Using hydrostatic balance: p_sfc = p_int * exp(g * Δz / (R_m_T_avg)) together with ideal 
    # gas law ρ = p / (R_m * T), we get:
    ρ_sfc = ρ_int * (R_m_int * T_int) / (R_m_sfc * T_sfc) * exp(grav * Δz / R_m_T_avg)

    # Surface density
    return ρ_sfc
end

"""
    effective_height(inputs)

Compute the effective aerodynamic height `z_eff = Δz - d`.

# Arguments
- `inputs`: `SurfaceFluxInputs` struct with `Δz` and `d`.

Returns `Δz - d` [m].
"""
@inline function effective_height(inputs::SurfaceFluxInputs)
    FT = typeof(inputs.Δz)
    return max(inputs.Δz - inputs.d, eps(FT))
end
