"""
    shf = sensible_heat_flux(param_set, inputs, g_h, T_int, T_sfc, ρ_sfc, E)

Computes the sensible heat flux at the surface.

The sensible heat flux is given by

    SHF = -ρ_sfc * g_h * ΔDSE + VSE_sfc * E

where `ΔDSE = DSE_int - DSE_sfc` is the difference in dry static energy between
the interior and surface, `g_h` is the heat/moisture conductance, `VSE_sfc` is the dry 
static energy of water vapor at the surface temperature, and `E` is the evaporation rate. The 
second term, `VSE_sfc * E`, accounts for the vapor static energy 
`VSE_sfc = cp_v * (T_sfc - T_0) + Φ_sfc` (i.e., dry enthalpy `cp_v * (T_sfc - T_0)`, 
or sensible heat, plus potential energy Φ_sfc) carried by evaporating water.

If `inputs.shf` is provided (not `nothing`), the function returns that value directly,
allowing for prescribed sensible heat flux conditions. See [`SurfaceFluxInputs`](@ref).

# Arguments
- `param_set`: Parameter set.
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
- `g_h`: Heat/moisture conductance [m/s].
- `T_int`: Interior temperature [K].
- `T_sfc`: Surface temperature [K].
- `ρ_sfc`: Surface air density [kg/m^3].
- `E`: Evaporation rate [kg/m^2/s]. Optional, default 0.
"""
@inline function sensible_heat_flux(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    g_h,
    T_int,
    T_sfc,
    ρ_sfc,
    E = eltype(ρ_sfc)(0),
)
    shf_in = inputs.shf
    if shf_in !== nothing
        return shf_in
    end
    thermo_params = SFP.thermodynamics_params(param_set)
    Φ_sfc = surface_geopotential(inputs)
    Φ_int = interior_geopotential(param_set, inputs)
    DSE_int = TD.dry_static_energy(thermo_params, T_int, Φ_int)
    DSE_sfc = TD.dry_static_energy(thermo_params, T_sfc, Φ_sfc)
    ΔDSE = DSE_int - DSE_sfc
    VSE_sfc = TD.vapor_static_energy(thermo_params, T_sfc, Φ_sfc)
    return -ρ_sfc * g_h * ΔDSE + VSE_sfc * E
end

"""
    sensible_heat_flux(param_set, ζ, ustar, inputs, z0m, z0h, T_sfc, q_vap_sfc, ρ_sfc, scheme)

Computes the sensible heat flux given the Monin-Obukhov stability parameter `ζ`,
friction velocity `ustar`, roughness lengths, and surface state.
Useful for computing fluxes from variables available inside the solver loop.

# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `ustar`: Friction velocity [m/s].
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
- `z0m`: Momentum roughness length [m].
- `z0h`: Thermal roughness length [m].
- `T_sfc`: Surface temperature [K].
- `q_vap_sfc`: Surface vapor specific humidity [kg/kg].
- `ρ_sfc`: Surface air density [kg/m^3].
- `scheme`: Discretization scheme.
"""
@inline function sensible_heat_flux(
    param_set::APS,
    ζ,
    ustar,
    inputs::SurfaceFluxInputs,
    z0m,
    z0h,
    T_sfc,
    q_vap_sfc,
    ρ_sfc,
    scheme,
)

    # Compute conductance
    g_h = heat_conductance(param_set, ζ, ustar, inputs, z0m, z0h, scheme)
    q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int

    # Compute evaporation (needed for SHF correction)
    E = evaporation(
        param_set,
        inputs,
        g_h,
        q_vap_int,
        q_vap_sfc,
        ρ_sfc,
        inputs.moisture_model,
    )

    return sensible_heat_flux(
        param_set,
        inputs,
        g_h,
        inputs.T_int,
        T_sfc,
        ρ_sfc,
        E,
    )
end

"""
    E = evaporation(param_set, inputs, g_h, q_vap_int, q_vap_sfc, ρ_sfc, model)

Computes the evaporation rate at the surface.

The evaporation rate is given by

    E = -ρ_sfc * g_h * Δq_vap

where `Δq_vap = q_vap_int - q_vap_sfc` is the difference in vapor specific humidity between
the interior and surface, `g_h` is the heat/moisture conductance (heat and moisture 
conductances must be equal for energetic consistency), and `ρ_sfc` is the
surface air density. Here `q_vap_int` and `q_vap_sfc` are the vapor specific humidities
(not total specific humidity) at the interior and surface, respectively.

If `inputs.lhf` is provided (not `nothing`), the function returns the evaporation
rate computed from the prescribed latent heat flux: `E = LHF / LH_v0`, where
`LH_v0` is the latent heat of vaporization at the reference temperature.

Arguments:
- `param_set`: Parameter set.
- `inputs`: [`SurfaceFluxInputs`](@ref) struct.
- `g_h`: Heat conductance [m/s].
- `q_vap_int`: Interior vapor specific humidity [kg/kg].
- `q_vap_sfc`: Surface vapor specific humidity [kg/kg].
- `ρ_sfc`: Surface density [kg/m^3].
- `model`: Moisture model ([`MoistModel`](@ref) or [`DryModel`](@ref)).
"""
@inline function evaporation(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    g_h,
    q_vap_int,
    q_vap_sfc,
    ρ_sfc,
    model::AbstractMoistureModel = MoistModel(),
)
    lhf_in = inputs.lhf
    if lhf_in !== nothing
        thermo_params = SFP.thermodynamics_params(param_set)
        LH_v0 = TP.LH_v0(thermo_params)
        return lhf_in / LH_v0
    end
    Δq_vap = q_vap_int - q_vap_sfc
    return -ρ_sfc * g_h * Δq_vap
end

@inline function evaporation(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    g_h,
    q_vap_int,
    q_vap_sfc,
    ρ_sfc,
    model::DryModel,
)
    FT = eltype(param_set)
    return FT(0)
end

"""

    evaporation(param_set, ζ, ustar, inputs, z0m, z0h, q_vap_sfc, ρ_sfc, scheme)

Computes the evaporation rate given the Monin-Obukhov stability parameter `ζ`,
friction velocity `ustar`, roughness lengths, and surface state.

# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `ustar`: Friction velocity [m/s].
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
- `z0m`: Momentum roughness length [m].
- `z0h`: Thermal roughness length [m].
- `q_vap_sfc`: Surface vapor specific humidity [kg/kg].
- `ρ_sfc`: Surface air density [kg/m^3].
- `scheme`: Discretization scheme.
"""
@inline function evaporation(
    param_set::APS,
    ζ,
    ustar,
    inputs::SurfaceFluxInputs,
    z0m,
    z0h,
    q_vap_sfc,
    ρ_sfc,
    scheme,
)
    # Compute conductance
    g_h = heat_conductance(param_set, ζ, ustar, inputs, z0m, z0h, scheme)
    q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int

    return evaporation(
        param_set,
        inputs,
        g_h,
        q_vap_int,
        q_vap_sfc,
        ρ_sfc,
        inputs.moisture_model,
    )
end

"""
    latent_heat_flux(param_set, inputs, E, model)

Computes the latent heat flux at the surface.

The latent heat flux is given by

    LHF = LH_v0 * E

where `LH_v0` is the latent heat of vaporization at the reference temperature
and `E` is the evaporation rate.

If `inputs.lhf` is provided (not `nothing`), the function returns that value directly,
allowing for prescribed latent heat flux conditions.

Arguments:
- `param_set`: Parameter set.
- `inputs`: [`SurfaceFluxInputs`](@ref) struct.
- `E`: Evaporation rate [kg/m^2/s].
- `model`: Moisture model ([`MoistModel`](@ref) or [`DryModel`](@ref)).
"""
@inline function latent_heat_flux(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    E,
    model::AbstractMoistureModel = MoistModel(),
)
    lhf_in = inputs.lhf
    if lhf_in !== nothing
        return lhf_in
    end
    thermo_params = SFP.thermodynamics_params(param_set)
    LH_v0 = TP.LH_v0(thermo_params)
    return LH_v0 * E
end

@inline function latent_heat_flux(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    E,
    model::DryModel,
)
    FT = eltype(param_set)
    return FT(0)
end

"""

    latent_heat_flux(param_set, ζ, ustar, inputs, z0m, z0h, q_vap_sfc, ρ_sfc, scheme)

Computes the latent heat flux given the Monin-Obukhov stability parameter `ζ`,
friction velocity `ustar`, roughness lengths, and surface state.
Calculates conductance and evaporation internally.

# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `ustar`: Friction velocity [m/s].
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
- `z0m`: Momentum roughness length [m].
- `z0h`: Thermal roughness length [m].
- `q_vap_sfc`: Surface vapor specific humidity [kg/kg].
- `ρ_sfc`: Surface air density [kg/m^3].
- `scheme`: Discretization scheme.
"""
@inline function latent_heat_flux(
    param_set::APS,
    ζ,
    ustar,
    inputs::SurfaceFluxInputs,
    z0m,
    z0h,
    q_vap_sfc,
    ρ_sfc,
    scheme,
)
    E = evaporation(
        param_set,
        ζ,
        ustar,
        inputs,
        z0m,
        z0h,
        q_vap_sfc,
        ρ_sfc,
        scheme,
    )

    return latent_heat_flux(param_set, inputs, E, inputs.moisture_model)
end

"""
    buoyancy_flux(param_set, shf, lhf, T_sfc, ρ_sfc, q_vap_sfc, q_liq_sfc, q_ice_sfc, model)

Computes the buoyancy flux at the surface, accounting for the presence of liquid and ice condensate.

The buoyancy flux `B` is defined as the vertical flux of virtual potential temperature `θ_v`.
It is approximated by linearizing the density perturbations with respect to temperature and total water content:

    B ≈ (g / ρ_sfc) * ( SHF / (cp_m * T_sfc) + (ε_vd - 1) * LHF / LH_v0 )

Where:
 - `cp_m` is the specific heat of moist air, calculated using `q_tot_sfc`, `q_liq_sfc`, and `q_ice_sfc`.
 - `ε_vd` is the ratio of gas constants for water vapor and dry air.
 - The term `(ε_vd - 1)` (approximately 0.61) represents the density effect of water vapor relative to dry air 
    (the virtual temperature correction factor), ensuring the buoyancy flux accounts for the fact that moist 
    air is lighter than dry air.

Arguments:
 - `ρ_sfc`: Surface air density [kg/m³].
 - `q_vap_sfc`: Specific humidity of water vapor at the surface (default: 0).
 - `q_liq_sfc`: Specific humidity of liquid water at the surface (default: 0).
 - `q_ice_sfc`: Specific humidity of ice at the surface (default: 0).
 - `model`: Moisture model ([`MoistModel`](@ref) or [`DryModel`](@ref)).
"""
@inline function buoyancy_flux(
    param_set::APS,
    shf,
    lhf,
    T_sfc,
    ρ_sfc,
    q_vap_sfc = 0,
    q_liq_sfc = 0,
    q_ice_sfc = 0,
    model::AbstractMoistureModel = MoistModel(),
)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    ε_vd = TD.Parameters.Rv_over_Rd(thermo_params)

    # Calculate specific heat of moist air at the surface including condensate
    q_tot_sfc = q_vap_sfc + q_liq_sfc + q_ice_sfc
    cp_m_sfc = TD.cp_m(thermo_params, q_tot_sfc, q_liq_sfc, q_ice_sfc)

    LH_v0 = TP.LH_v0(thermo_params)

    # Term 1: Sensible heat flux contribution
    term_shf = shf / (cp_m_sfc * T_sfc)

    # Term 2: Latent heat flux contribution
    # This term captures the buoyancy effect of water vapor concentration changes.
    term_lhf = (ε_vd - 1) * lhf / LH_v0

    return (grav / ρ_sfc) * (term_shf + term_lhf)
end

@inline function buoyancy_flux(
    param_set::APS,
    shf,
    lhf,
    T_sfc,
    ρ_sfc,
    q_vap_sfc,
    q_liq_sfc,
    q_ice_sfc,
    model::DryModel,
)
    FT = eltype(param_set)
    return buoyancy_flux(
        param_set,
        shf,
        FT(0),
        T_sfc,
        ρ_sfc,
        FT(0),
        FT(0),
        FT(0),
        MoistModel(),
    )
end

"""
    buoyancy_flux(param_set, ζ, ustar, inputs)

Computes the buoyancy flux given the Monin-Obukhov stability parameter `ζ`,
friction velocity `ustar`, and geometric inputs via [`SurfaceFluxInputs`](@ref).

The relationship is derived from the definition of the Obukhov length:

    L = -u_*^3 / (κ * B)
    ζ = Δz / L
    => B = -(u_*^3 * ζ) / (κ * Δz)

# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `ustar`: Friction velocity [m/s].
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
"""
@inline function buoyancy_flux(
    param_set::APS,
    ζ,
    ustar,
    inputs::SurfaceFluxInputs,
)
    κ = SFP.von_karman_const(param_set)
    Δz_eff = effective_height(inputs)
    return -(ustar^3 * ζ) / (κ * Δz_eff)
end

"""
    momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)

Computes the momentum fluxes at the surface.

The momentum fluxes are calculated using the bulk aerodynamic formula:

    ρτxz = -ρ_sfc * Cd * ΔU * Δu_x
    ρτyz = -ρ_sfc * Cd * ΔU * Δu_y

where:
 - `Cd`: Momentum exchange coefficient (drag coefficient)
 - `ΔU`: Magnitude of the wind speed difference (including gustiness)
 - `Δu_x`, `Δu_y`: Components of the wind speed difference
 - `ρ_sfc`: Surface air density

Returns a tuple `(ρτxz, ρτyz)`.

See [`SurfaceFluxInputs`](@ref).

# Arguments
- `Cd`: Drag coefficient.
- `inputs`: [`SurfaceFluxInputs`](@ref) container.
- `ρ_sfc`: Surface air density [kg/m^3].
- `gustiness`: Gustiness velocity scale [m/s].
"""
@inline function momentum_fluxes(
    Cd,
    inputs::SurfaceFluxInputs,
    ρ_sfc,
    gustiness,
)
    Δu = Δu_components(inputs)
    ΔU = windspeed(inputs, gustiness)
    ρτxz = -ρ_sfc * Cd * Δu[1] * ΔU
    ρτyz = -ρ_sfc * Cd * Δu[2] * ΔU
    return (ρτxz, ρτyz)
end

"""
    state_bulk_richardson_number(param_set, inputs, T_sfc, ρ_sfc, ΔU, q_vap_sfc)

Computes the bulk Richardson number from the given state.

# Arguments
- `param_set`: Parameter set.
- `inputs`: [`SurfaceFluxInputs`](@ref) struct.
- `T_sfc`: Surface temperature [K].
- `ρ_sfc`: Surface air density [kg/m³].
- `ΔU`: Wind speed difference [m/s].
- `q_vap_sfc`: Surface vapor specific humidity [kg/kg]. Default: 0.

Returns the bulk Richardson number.
"""
@inline function state_bulk_richardson_number(
    param_set::APS,
    inputs::SurfaceFluxInputs,
    T_sfc,
    ρ_sfc,
    ΔU,
    q_vap_sfc = 0,
)
    FT = eltype(param_set)
    grav = SFP.grav(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)

    q_tot_int = inputs.q_tot_int
    q_liq_int = inputs.q_liq_int
    q_ice_int = inputs.q_ice_int

    # Assume condensate concentration is the same at the surface and in the interior
    q_tot_sfc = q_vap_sfc + q_liq_int + q_ice_int
    theta_v_sfc =
        TD.virtual_pottemp(thermo_params, T_sfc, inputs.ρ_int, q_tot_sfc, q_liq_int, q_ice_int)
    theta_v_int = TD.virtual_pottemp(
        thermo_params,
        inputs.T_int,
        inputs.ρ_int,
        q_tot_int,
        q_liq_int,
        q_ice_int,
    )

    state_sfc = TD.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_tot_sfc)
    state_int = TD.PhaseEquil_ρTq(thermo_params, inputs.ρ_int, inputs.T_int, q_tot_int)
    dse_v_sfc = TD.virtual_dry_static_energy(thermo_params, state_sfc, FT(0))
    dse_v_int = TD.virtual_dry_static_energy(thermo_params, state_int, FT(300))

    Δtheta_v = theta_v_int - theta_v_sfc
    Δdse_v = dse_v_int - dse_v_sfc
    cp_d = SFP.cp_d(param_set)
    dse_v_ref = cp_d * inputs.T_int
    theta_v_ref = theta_v_int
    Δz_eff = effective_height(inputs)

    #Rib_state = (grav * Δz_eff * Δtheta_v) / (theta_v_ref * non_zero(ΔU)^2)
    Rib_state = (grav * Δz_eff * Δdse_v) / (dse_v_ref * non_zero(ΔU)^2)
    return Rib_state
end
