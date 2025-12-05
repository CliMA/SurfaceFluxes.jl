"""
    non_zero(v)

Ensure that `v` is not zero, returning `eps(v)` (preserving sign) if `v` is too small.
"""
function non_zero(v::FT) where {FT}
    sign_of_v = v == 0 ? 1 : sign(v)
    return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
end




@inline function interior_geopotential(param_set::APS, inputs::SurfaceFluxInputs)
    return inputs.Φ_sfc + SFP.grav(param_set) * inputs.Δz
end

@inline surface_geopotential(inputs::SurfaceFluxInputs) = inputs.Φ_sfc

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
    R_m = TD.gas_constant_air(thermo_params, qt_int, ql_int, qi_int)
    cv_m = TD.cv_m(thermo_params, qt_int, ql_int, qi_int)
    ratio = T_sfc / T_int
    return ρ_int * ratio^(cv_m / R_m)
end





