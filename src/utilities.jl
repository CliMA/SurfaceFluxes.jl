
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
# When args = nothing, qt_sfc and surface_specific_humidity are the same
function surface_specific_humidity(param_set::APS, sc::AbstractSurfaceConditions, args = nothing)
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_sfc(sc))
end
qt_sfc(param_set::APS, sc::AbstractSurfaceConditions, args = nothing) =
    surface_specific_humidity(param_set, sc, args)
# Alternative formulation: when Real arguments are provided, return qt_sfc + sum of all args
function surface_specific_humidity(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    Q1::Real,
    Q2::Real,
)
    q_base = qt_sfc(param_set, sc, nothing)
    return q_base + Q1 + Q2
end
# Method to handle named tuples - unpacks the arguments
function surface_specific_humidity(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    args::NamedTuple,
)
    # Unpack named tuple and pass as varargs
    return surface_specific_humidity(param_set, sc, values(args)...)
end
# Splattable method for alternative formulations - accepts args but ignores them (following obukhov_similarity_solution pattern)
function surface_specific_humidity(param_set::APS, sc::AbstractSurfaceConditions, args...)
    return surface_specific_humidity(param_set, sc, nothing)
end
Δqt(param_set::APS, sc::AbstractSurfaceConditions, args = nothing) =
    qt_in(param_set, sc) - qt_sfc(param_set, sc, args)

# Air temperature
T_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.air_temperature(SFP.thermodynamics_params(param_set), ts_in(sc))
# When args = nothing, T_sfc and surface_temperature are the same
function surface_temperature(param_set::APS, sc::AbstractSurfaceConditions, args = nothing)
    TD.air_temperature(SFP.thermodynamics_params(param_set), ts_sfc(sc))
end
T_sfc(param_set::APS, sc::AbstractSurfaceConditions, args = nothing) =
    surface_temperature(param_set, sc, args)
# Alternative formulation: when Real arguments are provided, return T_sfc + sum of all args
function surface_temperature(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    T1::Real,
    T2::Real,
)
    T_base = T_sfc(param_set, sc, nothing)
    return T_base + T1 + T2
end
# Method to handle named tuples - unpacks the arguments
function surface_temperature(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    args::NamedTuple,
)
    # Unpack named tuple and pass as varargs
    return surface_temperature(param_set, sc, values(args)...)
end
# Splattable method for alternative formulations - accepts args but ignores them (following obukhov_similarity_solution pattern)
function surface_temperature(param_set::APS, sc::AbstractSurfaceConditions, args...)
    return surface_temperature(param_set, sc, nothing)
end
ΔT(param_set::APS, sc::AbstractSurfaceConditions, args...) =
    T_in(param_set, sc) - surface_temperature(param_set, sc, args...)

# Virtual Potential Temperature
θᵥ_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_pottemp(SFP.thermodynamics_params(param_set), ts_in(sc))
# When args = nothing, θᵥ_sfc and surface_virtual_pottemp are the same
function surface_virtual_pottemp(param_set::APS, sc::AbstractSurfaceConditions, args = nothing)
    TD.virtual_pottemp(SFP.thermodynamics_params(param_set), ts_sfc(sc))
end
# Alternative formulation: when Real arguments are provided, return θᵥ_sfc + sum of all args
function surface_virtual_pottemp(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    θ1::Real,
    θ2::Real,
)
    θ_base = surface_virtual_pottemp(param_set, sc, nothing)
    return θ_base + θ1 + θ2
end
# Method to handle named tuples - unpacks the arguments
function surface_virtual_pottemp(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    args::NamedTuple,
)
    # Unpack named tuple and pass as varargs
    return surface_virtual_pottemp(param_set, sc, values(args)...)
end
# Splattable method for alternative formulations - accepts args but ignores them (following obukhov_similarity_solution pattern)
function surface_virtual_pottemp(param_set::APS, sc::AbstractSurfaceConditions, args...)
    return surface_virtual_pottemp(param_set, sc, nothing)
end
θᵥ_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    surface_virtual_pottemp(param_set, sc, nothing)
Δθᵥ(param_set::APS, sc::AbstractSurfaceConditions) =
    θᵥ_in(param_set, sc) - θᵥ_sfc(param_set, sc)

# Virtual Dry Static Energy
DSEᵥ_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_dry_static_energy(SFP.thermodynamics_params(param_set),
        ts_in(sc),
        SFP.grav(param_set) * z_in(sc))
DSEᵥ_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_dry_static_energy(SFP.thermodynamics_params(param_set),
        ts_sfc(sc),
        SFP.grav(param_set) * z_sfc(sc))
ΔDSEᵥ(param_set::APS, sc::AbstractSurfaceConditions) =
    DSEᵥ_in(param_set, sc) - DSEᵥ_sfc(param_set, sc)
