function non_zero(v::FT) where {FT}
    sign_of_v = v == 0 ? 1 : sign(v)
    return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
end


@inline z_in(inputs::SurfaceFluxInputs) = inputs.d + inputs.Δz
@inline z_sfc(inputs::SurfaceFluxInputs) = inputs.d
@inline Δz(inputs::SurfaceFluxInputs) = inputs.Δz
@inline Δz(inputs::Fluxes) = inputs.state_int.z - inputs.state_sfc.z
@inline Δz(inputs::FluxesAndFrictionVelocity) = inputs.state_int.z - inputs.state_sfc.z
@inline Δz(inputs::ValuesOnly) = inputs.state_int.z - inputs.state_sfc.z


@inline function interior_geopotential(param_set::APS, inputs::SurfaceFluxInputs)
    return inputs.Φs + SFP.grav(param_set) * inputs.Δz
end

@inline surface_geopotential(inputs::SurfaceFluxInputs) = inputs.Φs

"""
    surface_density(param_set, T_int, ρ_int, T_sfc, qt_int=0, ql_int=0, qi_int=0)

Estimates the surface air density assuming an adiabatic lapse rate (isentropic process) 
for the temperature ratio between the interior and surface.

# Arguments
- `param_set`: AbstractSurfaceFluxesParameters.
- `T_int`: Interior temperature [K].
- `ρ_int`: Interior density [kg/m^3].
- `T_sfc`: Surface temperature [K].
- `qt_int`: Interior total specific humidity.
- `ql_int`: Interior liquid specific humidity.
- `qi_int`: Interior ice specific humidity.
"""
@inline function surface_density(
    param_set::APS,
    T_int,
    ρ_int,
    T_sfc,
    qt_int = 0,
    ql_int = 0,
    qi_int = 0,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    phase = TD.PhasePartition(qt_int, ql_int, qi_int)
    cv_m = TD.cv_m(thermo_params, phase)
    R_m = TD.gas_constant_air(thermo_params, phase)
    ratio = T_sfc / T_int
    return ρ_int * ratio^(cv_m / R_m)
end

"""
    Δq_vap(q_vap_int, q_vap_sfc)

Computes the difference in vapor specific humidity between the interior and surface.
"""
@inline Δq_vap(q_vap_int, q_vap_sfc) = q_vap_int - q_vap_sfc

"""
    ΔT(T_int, T_sfc)

Computes the difference in temperature between the interior and surface.
"""
@inline ΔT(T_int, T_sfc) = T_int - T_sfc

"""
    θᵥ(param_set, T, ρ, phase)

Computes the virtual potential temperature.

# Arguments
- `param_set`: AbstractSurfaceFluxesParameters.
- `T`: Temperature [K].
- `ρ`: Density [kg/m^3].
- `phase`: Phase partition.
"""
@inline function θᵥ(param_set::APS, T, ρ, phase)
    return TD.virtual_pottemp(SFP.thermodynamics_params(param_set), T, ρ, phase)
end

"""
    virtual_temperature(param_set, T, phase)

Computes the virtual temperature.

# Arguments
- `param_set`: AbstractSurfaceFluxesParameters.
- `T`: Temperature [K].
- `phase`: Phase partition.
"""
@inline function virtual_temperature(param_set::APS, T, phase)
    return TD.virtual_temperature(SFP.thermodynamics_params(param_set), T, phase)
end

"""
    Δθᵥ(param_set, T_int, ρ_int, phase_int, T_sfc, ρ_sfc, phase_sfc)

Computes the difference in virtual potential temperature between the interior and surface.
"""
@inline function Δθᵥ(
    param_set::APS,
    T_int,
    ρ_int,
    phase_int,
    T_sfc,
    ρ_sfc,
    phase_sfc,
)
    return θᵥ(param_set, T_int, ρ_int, phase_int) -
           θᵥ(param_set, T_sfc, ρ_sfc, phase_sfc)
end

