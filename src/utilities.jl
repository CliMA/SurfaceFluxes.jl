
function non_zero(v::FT) where {FT}
    sign_of_v = v == 0 ? 1 : sign(v)
    return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
end

function windspeed(sc::AbstractSurfaceConditions)
    return max(hypot(Δu1(sc), Δu2(sc)), sc.gustiness)
end

### Utilitity functions for calculations of differences between
### atmospheric state properties at the first interior node and

# Thermodynamic States
ts_in(sc::AbstractSurfaceConditions) = sc.state_in.ts
ts_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.ts

# Near-surface layer depth
z_in(sc::AbstractSurfaceConditions) = sc.state_in.z
z_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.z
Δz(sc::AbstractSurfaceConditions) = z_in(sc) - z_sfc(sc)

# Velocity
Δu1(sc::AbstractSurfaceConditions) = sc.state_in.u[1] - sc.state_sfc.u[1]
Δu2(sc::AbstractSurfaceConditions) = sc.state_in.u[2] - sc.state_sfc.u[2]

# Total Specific Humidity
qt_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_in(sc))
qt_sfc(param_set::APS, sc::AbstractSurfaceConditions, args=nothing) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_sfc(sc))
function surface_specific_humidity(param_set::APS, sc::AbstractSurfaceConditions, args=nothing)
    qt_sfc(param_set, sc, args)
end
Δqt(param_set::APS, sc::AbstractSurfaceConditions, args=nothing) =
    qt_in(param_set, sc) - qt_sfc(param_set, sc, args)

# Air temperature
T_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.air_temperature(SFP.thermodynamics_params(param_set), ts_in(sc))
T_sfc(param_set::APS, sc::AbstractSurfaceConditions, args=nothing) =
    TD.air_temperature(SFP.thermodynamics_params(param_set), ts_sfc(sc))
function surface_temperature(param_set::APS, sc::AbstractSurfaceConditions, args=nothing)
    T_sfc(param_set, sc, args)
end
ΔT(param_set::APS, sc::AbstractSurfaceConditions, args=nothing) =
    T_in(param_set, sc) - T_sfc(param_set, sc, args)

# Virtual Potential Temperature
θᵥ_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_pottemp(SFP.thermodynamics_params(param_set), ts_in(sc))
θᵥ_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_pottemp(SFP.thermodynamics_params(param_set), ts_sfc(sc))
Δθᵥ(param_set::APS, sc::AbstractSurfaceConditions) =
    θᵥ_in(param_set, sc) - θᵥ_sfc(param_set, sc)

# Virtual Dry Static Energy
DSEᵥ_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_dry_static_energy(SFP.thermodynamics_params(param_set),
        ts_in(sc),
        SFP.grav(param_set)*z_in(sc))
DSEᵥ_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_dry_static_energy(SFP.thermodynamics_params(param_set),
        ts_sfc(sc),
        SFP.grav(param_set)*z_sfc(sc))
ΔDSEᵥ(param_set::APS, sc::AbstractSurfaceConditions) =
    DSEᵥ_in(param_set, sc) - DSEᵥ_sfc(param_set, sc)
