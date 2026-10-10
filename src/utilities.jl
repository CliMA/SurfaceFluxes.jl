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
    float_parameter(FT, x)

Return the model parameter `x` in the floating-point type `FT` of the inputs. Plain
floating-point and integer parameters are converted, so that models constructed with the
default `Float64` parameters (e.g., `ExponentialRSL()`, `ConstantStabilityCap(0.5)`) keep
`Float32` computations in `Float32`. Other numbers (e.g., dual numbers for differentiation
with respect to the parameter) are returned unchanged.
"""
@inline float_parameter(::Type{FT}, x::Union{AbstractFloat, Integer}) where {FT} = FT(x)
@inline float_parameter(::Type{FT}, x) where {FT} = x

"""
    interior_geopotential(param_set, inputs)

Compute the geopotential at the interior (atmospheric) reference level.

# Arguments
- `param_set`: Parameter set containing gravitational constant.
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).

Returns `Φ_sfc + g * Δz` [m²/s²], with `Δz` measured from the surface (see
[`reference_above_surface`](@ref)).
"""
@inline function interior_geopotential(param_set::APS, inputs)
    inputs = reference_above_surface(param_set, inputs)
    return inputs.Φ_sfc + SFP.grav(param_set) * inputs.Δz
end

"""
    surface_geopotential(param_set, inputs)

Compute the geopotential of the surface state, `Φ_sfc + g * d` [m²/s²]. The surface
temperature and humidity apply at the apparent sink of the Monin-Obukhov profiles, which
lies at the displacement height `d` above the surface (the roughness length `z0h` above
it is neglected). Over a canopy, this is the level of the leaves that exchange heat with
the air, so the dry static energy difference to the reference level,
`cp (T_int - T_sfc) + g (Δz - d)`, does not include the air column below the canopy.

# Arguments
- `param_set`: Parameter set containing gravitational constant.
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
"""
@inline function surface_geopotential(param_set::APS, inputs)
    return inputs.Φ_sfc + SFP.grav(param_set) * inputs.d
end

"""
    surface_geopotential(inputs)

Return the geopotential of the ground, `inputs.Φ_sfc` [m²/s²]. Deprecated: the surface
state applies at the displacement height, with the geopotential `Φ_sfc + g d` of
[`surface_geopotential`](@ref)`(param_set, inputs)`; a surface energy balance that uses
this form with a displaced canopy is inconsistent with the sensible heat flux by
`g d / cp`. Kept for callers of the one-argument form; to be removed in the next
breaking release.
"""
@inline surface_geopotential(inputs) = inputs.Φ_sfc

"""
    surface_density(param_set, T_int, ρ_int, T_sfc, Δz, q_tot_int=0, q_liq_int=0, q_ice_int=0, q_vap_sfc=nothing)
    surface_density(param_set, inputs, T_sfc, q_vap_sfc)

Estimates the surface air density assuming hydrostatic balance between the interior and surface.
It effectively extrapolates the interior pressure to the surface using the hydrostatic 
equation with an average virtual temperature, and then computes the surface density using the 
ideal gas law. The form with the inputs container extrapolates over the effective height
`Δz - d`, from the reference level to the displacement height where the surface state
applies (see [`surface_geopotential`](@ref)).

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

    # Using hydrostatic balance: p_sfc = p_int * exp(g * Δz / (R_m_T_avg)) together with
    # ideal gas law ρ = p / (R_m * T), we get:
    ρ_sfc = ρ_int * (R_m_int * T_int) / (R_m_sfc * T_sfc) * exp(grav * Δz / R_m_T_avg)

    # Surface density
    return ρ_sfc
end

@inline function surface_density(param_set::APS, inputs, T_sfc, q_vap_sfc)
    return surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc,
        effective_height(param_set, inputs),
        inputs.q_tot_int,
        inputs.q_liq_int,
        inputs.q_ice_int,
        q_vap_sfc,
    )
end

"""
    effective_height(inputs)
    effective_height(param_set, inputs)

Compute the effective aerodynamic height `z_eff = Δz - d`, the height of the reference
level above the displacement height, which the Monin-Obukhov profiles span and over which
the surface state at `d` (see [`surface_geopotential`](@ref)) is connected to the interior
state. The one-argument form expects inputs under [`ReferenceAboveSurface`](@ref); the
two-argument form converts inputs under [`ReferenceAboveApparentSink`](@ref) first with
[`reference_above_surface`](@ref).

# Arguments
- `param_set`: Parameter set (required when `inputs` may use [`ReferenceAboveApparentSink`](@ref)).
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).

Returns `Δz - d` [m].
"""
@inline function effective_height(inputs)
    FT = typeof(inputs.Δz)
    return max(inputs.Δz - inputs.d, eps(FT))
end

@inline effective_height(param_set::APS, inputs) =
    effective_height(reference_above_surface(param_set, inputs))

"""
    reference_above_surface(param_set, inputs)

Return the inputs with the reference height `Δz` measured from the surface. Under
[`ReferenceAboveSurface`](@ref), the inputs are returned as they are. Under
[`ReferenceAboveApparentSink`](@ref), `Δz` is measured from the apparent sink for momentum
and becomes `Δz + d + z0m`, with the roughness length `z0m` of a roughness model that is
independent of the friction velocity.

# Arguments
- `param_set`: Parameter set.
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
"""
@inline reference_above_surface(param_set::APS, inputs) = reference_above_surface(
    get(inputs, :reference_level, ReferenceAboveSurface()),
    param_set,
    inputs,
)
@inline reference_above_surface(::ReferenceAboveSurface, param_set::APS, inputs) = inputs
@inline function reference_above_surface(
    ::ReferenceAboveApparentSink,
    param_set::APS,
    inputs,
)
    depends_on_ustar(inputs.roughness_model) && throw(
        ArgumentError(
            "ReferenceAboveApparentSink requires a roughness model independent of the friction velocity",
        ),
    )
    # The roughness model is independent of u★, so any value of u★ serves
    z0m = momentum_roughness(
        inputs.roughness_model,
        zero(inputs.Δz),
        param_set,
        inputs.roughness_inputs,
    )
    return merge(
        inputs,
        (; Δz = inputs.Δz + inputs.d + z0m, reference_level = ReferenceAboveSurface()),
    )
end

"""
    interior_vapor_specific_humidity(inputs)

Return the vapor specific humidity of the interior air [kg/kg], the total specific
humidity `q_tot_int` less the condensate `q_liq_int + q_ice_int`.

# Arguments
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
"""
@inline interior_vapor_specific_humidity(inputs) =
    inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int

"""
    reference_height_valid(inputs, z0m, z0h = z0m)

Whether the reference level lies above both roughness lengths, `Δz - d > max(z0m, z0h)`,
so that the Monin-Obukhov profiles of momentum and of scalars between the surface and
the reference level are defined. The scalar roughness length matters when it exceeds the
momentum one, as the COARE 3.0 model gives at low friction velocities.
`Δz` is the height above the surface: inputs under [`ReferenceAboveApparentSink`](@ref)
are converted first with [`reference_above_surface`](@ref), as [`surface_fluxes`](@ref)
does. [`surface_fluxes`](@ref) returns `NaN` fluxes with `converged = false` for inputs
that fail this test, since the solve cannot throw inside a GPU kernel;
[`check_reference_height`](@ref) raises the corresponding error on the host.

# Arguments
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `z0m`: Momentum roughness length [m].
- `z0h`: Scalar roughness length [m]; `z0m` by default.
"""
@inline reference_height_valid(inputs, z0m, z0h = z0m) =
    inputs.Δz - inputs.d > max(z0m, z0h)

"""
    check_reference_height(Δz, d, z0m, z0h = z0m)

Throw an `ArgumentError` unless the reference level lies above both roughness lengths,
`Δz - d > max(z0m, z0h)` (see [`reference_height_valid`](@ref)). For a host-side check
of a model's configuration before fluxes are computed in kernels.

# Arguments
- `Δz`: Height of the reference level above the surface [m].
- `d`: Displacement height [m].
- `z0m`: Momentum roughness length [m].
- `z0h`: Scalar roughness length [m]; `z0m` by default.
"""
function check_reference_height(Δz, d, z0m, z0h = z0m)
    Δz - d > max(z0m, z0h) || throw(
        ArgumentError(
            "The reference height Δz = $Δz m must exceed the displacement height d = $d m plus the larger roughness length max(z0m, z0h) = $(max(z0m, z0h)) m",
        ),
    )
    return nothing
end

# ============================================================================
# Quadrature
# ============================================================================

"""
    gauss_legendre4(f, a, b)

4-point Gauss-Legendre quadrature of callable `f` on `[a, b]`.

Nodes and weights are the exact algebraic 4-point rule on [-1, 1]:
```
t = ±√((3 ∓ 2√(6/5))/7),    w = (18 ± √30)/36
```
mapped affinely via `x = m ± h·t` with `m = (a+b)/2`, `h = (b-a)/2`.
`f` should be a functor (not a capturing closure) for allocation-free /
AD-safe evaluation.
"""
@inline function gauss_legendre4(f::F, a, b) where {F}
    FT = typeof(a)
    s = sqrt(FT(6) / FT(5))
    t_inner = sqrt((FT(3) - FT(2) * s) / FT(7))
    t_outer = sqrt((FT(3) + FT(2) * s) / FT(7))
    w_inner = (FT(18) + sqrt(FT(30))) / FT(36)
    w_outer = (FT(18) - sqrt(FT(30))) / FT(36)

    m = (a + b) / 2
    h = (b - a) / 2
    return h * (
        w_inner * (f(m - h * t_inner) + f(m + h * t_inner)) +
        w_outer * (f(m - h * t_outer) + f(m + h * t_outer))
    )
end
