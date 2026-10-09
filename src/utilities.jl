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

Returns `Φ_sfc + g * Δz` [m²/s²].
"""
@inline function interior_geopotential(param_set::APS, inputs)
    return inputs.Φ_sfc + SFP.grav(param_set) * inputs.Δz
end

"""
    surface_geopotential(inputs)

Return the surface geopotential from the inputs.

# Arguments
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).

Returns `inputs.Φ_sfc` [m²/s²].
"""
@inline surface_geopotential(inputs) = inputs.Φ_sfc

"""
    surface_density(param_set, T_int, ρ_int, T_sfc, Δz, q_tot_int=0, q_liq_int=0, q_ice_int=0, q_vap_sfc=nothing)

Estimates the surface air density assuming hydrostatic balance between the interior and surface.
It effectively extrapolates the interior pressure to the surface using the hydrostatic 
equation with an average virtual temperature, and then computes the surface density using the 
ideal gas law.

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

"""
    effective_height(inputs)

Compute the effective aerodynamic height `z_eff = Δz - d`.

# Arguments
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).

Returns `Δz - d` [m].
"""
@inline function effective_height(inputs)
    FT = typeof(inputs.Δz)
    return max(inputs.Δz - inputs.d, eps(FT))
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
[`surface_fluxes`](@ref) returns `NaN` fluxes with `converged = false` for inputs that
fail this test, since the solve cannot throw inside a GPU kernel;
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
