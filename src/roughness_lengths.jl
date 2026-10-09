# Roughness-length evaluation
#
# These helpers resolve the momentum and scalar roughness length specifications.
# They are intentionally lightweight so they can be inlined inside GPU kernels.

"""
    ConstantRoughnessParams{FT} <: AbstractRoughnessParams

Roughness lengths fixed to constant values.

# Fields
- `z0m`: Momentum roughness length [m].
- `z0s`: Scalar roughness length [m]. Used for both heat (`z0h`) and humidity (`z0q`).

The keyword defaults shown here (`z0m = 2e-4` m, `z0s = 2e-5` m) are also the roughness
lengths used by [`surface_fluxes`](@ref) when `config` is omitted, via
`default_surface_flux_config`. When loading via `ClimaParams`, the values are read from the
TOML file. Most applications pass `z0m`/`z0s` explicitly or load them from `ClimaParams`.
"""
Base.@kwdef struct ConstantRoughnessParams{FT} <: AbstractRoughnessParams
    z0m::FT = 2e-4
    z0s::FT = 2e-5
end

"""
    COARE3RoughnessParams{FT} <: AbstractRoughnessParams

COARE 3.0 roughness parameterization (Fairall et al. 2003).

# References
- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). 
    Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm.
    Journal of Climate, 16, 571–591.
    [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)

The default values specified here are used when constructing the struct manually. When loading
via `ClimaParams`, these values are overwritten by the parameters in the TOML file.

"""
Base.@kwdef struct COARE3RoughnessParams{FT} <: AbstractRoughnessParams
    kinematic_visc::FT = 1.5e-5
    z0m_default::FT = 1e-4
    α_low::FT = 0.011
    α_high::FT = 0.018
    u_low::FT = 10.0
    u_high::FT = 18.0
end

"""
    RaupachRoughnessParams <: AbstractRoughnessParams

Raupach (1994) canopy roughness model: the momentum roughness length `z0m` and the
zero-plane displacement height `d` as functions of the canopy height `h` and the area
index of the canopy (see [`momentum_roughness`](@ref) and [`displacement_height`](@ref)).
The roughness inputs are the canopy height `roughness_inputs.h` and the plant area index
`Λ = roughness_inputs.PAI`, the single-sided area of all canopy elements (leaves, living or
dead, stems, and branches) per unit ground area, which Raupach (1994) calls the canopy
area index: the sum of the leaf and stem area indices, `LAI + SAI`. The field name `LAI`
is accepted as a deprecated alias for `PAI`.

Raupach (1994) writes the drag partition in terms of the frontal area index `λ`, the
frontal area of the canopy elements facing the mean wind per unit ground area, with
`Λ = 2λ` for isotropically oriented elements (see [`frontal_area_index`](@ref) and
[`canopy_area_index`](@ref)).

# Fields
- `C_R`: Drag coefficient of an isolated roughness element (0.3).
- `C_S`: Drag coefficient of the substrate at height `h` (0.003).
- `c_d1`: Constant of the displacement height expression (7.5).
- `stanton_number`: Ratio `z0s / z0m` of the scalar to the momentum roughness length;
  `exp(-kB⁻¹)` for an excess resistance `kB⁻¹ = ln(z0m / z0s)` (0.1, so that
  `kB⁻¹ ≈ 2.3`).
- `frontal_area_ratio`: Frontal area index per unit plant area index,
  `λ = frontal_area_ratio * Λ`; 0.5 for isotropically oriented elements (Raupach 1994).
  It enters the drag partition (Eq. 7) only; the displacement height (Eq. 8) is an
  empirical fit in `Λ`.
- `λ_min`: Floor on the frontal area index (0, so that `PAI` is used as given, following
  Raupach 1994). A positive floor stands in for stems and branches when the input counts
  leaves only, or vanishes for a canopy of nonzero height, so that such a canopy stays
  aerodynamically rough; with a plant area index that includes the stem area, it is
  unnecessary.
- `ustar_Uh_max`: Sheltering limit of `u★ / U(h)` (0.3, Raupach 1994, Eq. 7).
- `c_w`: Ratio `(z_w - d) / (h - d)` of the heights of the roughness-sublayer top `z_w`
  and the canopy top above the displacement height (2, Raupach 1994). It sets the
  roughness-sublayer influence function at the canopy top, `Ψ_h = ln c_w - 1 + 1 / c_w`
  (Eq. 5), which is 0.193 for `c_w = 2`.

# References
- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane displacement 
    as functions of canopy height and area index.
    Boundary-Layer Meteorology, 71, 211–216.
    [DOI: 10.1007/BF00709229](https://doi.org/10.1007/BF00709229)

The default values specified here are used when constructing the struct manually. When loading
via `ClimaParams`, all fields are read from the TOML file: `stanton_number` and the
`raupach_*` parameters.
"""
Base.@kwdef struct RaupachRoughnessParams{FT} <: AbstractRoughnessParams
    C_R::FT = 0.3
    C_S::FT = 0.003
    c_d1::FT = 7.5
    stanton_number::FT = 0.1
    frontal_area_ratio::FT = 0.5
    λ_min::FT = 0
    ustar_Uh_max::FT = 0.3
    c_w::FT = 2
end

"""
    charnock_parameter(mag_u_10, α_low, α_high, u_low, u_high)

Compute the Charnock parameter `α` as a function of the 10-m wind speed `mag_u_10` [m/s] 
estimated from the neutral profile. Piecewise linear interpolation between the lower and 
upper bounds based on COARE 3.0 (Fairall et al. 2003).

# References
- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). 
    Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm.
    Journal of Climate, 16, 571–591.
    [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)
"""
@inline function charnock_parameter(mag_u_10, α_low, α_high, u_low, u_high)
    return ifelse(
        mag_u_10 <= u_low,
        α_low,
        ifelse(
            mag_u_10 >= u_high,
            α_high,
            α_low + (α_high - α_low) * (mag_u_10 - u_low) / (u_high - u_low),
        ),
    )
end

# Accessors, in the floating-point type of the parameter set (see `float_parameter`)
@inline function momentum_roughness(
    spec::ConstantRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    return float_parameter(eltype(sfc_param_set), spec.z0m)
end

@inline function scalar_roughness(
    spec::ConstantRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    return float_parameter(eltype(sfc_param_set), spec.z0s)
end

@inline function momentum_and_scalar_roughness(
    spec::ConstantRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    return (
        momentum_roughness(spec, u★, sfc_param_set, roughness_inputs),
        scalar_roughness(spec, u★, sfc_param_set, roughness_inputs),
    )
end

"""
    momentum_roughness(spec::COARE3RoughnessParams, u★, sfc_param_set, roughness_inputs)

Calculate momentum roughness length using the COARE 3.0 algorithm (Fairall et al. 2003).

# Formulation
The momentum roughness length `z0m` is parameterized as the sum of a smooth flow limit
(Smith 1988) and a rough flow limit (Charnock 1955):
```math
z_{0m} = z_{0m,smooth} + z_{0m,rough}
```
- **Smooth flow**: Dominated by viscous limit, proportional to `ν / u★`.
- **Rough flow**: Dominated by wind stress, proportional to `α * u★^2 / g`.

The Charnock parameter `α` varies with the 10-m wind speed and is interpolated linearly between
lower and upper bounds defined in `spec`.

# Dependencies
- `u★`: Friction velocity [m/s]
- `kinematic_visc`: Kinematic viscosity of air [m^2/s], from `spec`
- `grav`: Gravitational acceleration [m/s^2], from `sfc_param_set`
- `mag_u_10`: 10m wind speed [m/s]. (internally recovered from u★ assuming neutral log profile)

# References
- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). 
    Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm.
    Journal of Climate, 16, 571–591.
    [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)
- Smith, S. D. (1988). Coefficients for sea surface wind stress, heat flux, and wind profiles 
    as a function of wind speed and temperature.
    Journal of Geophysical Research: Oceans, 93, 15467–15472.
    [DOI: 10.1029/JC093iC12p15467](https://doi.org/10.1029/JC093iC12p15467)
- Charnock, H. (1955). Wind stress on a water surface.
    Quarterly Journal of the Royal Meteorological Society, 81, 639–640.
    [DOI: 10.1002/qj.49708135027](https://doi.org/10.1002/qj.49708135027)
"""
@inline function momentum_roughness(
    spec::COARE3RoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    FT = eltype(sfc_param_set)
    grav = SFP.grav(sfc_param_set)
    kinematic_visc = spec.kinematic_visc

    # Recover 10-m wind speed using neutral profile with a proxy roughness length (to avoid 
    # circular dependency)
    z0_proxy = spec.z0m_default
    κ = SFP.von_karman_const(sfc_param_set)
    mag_u_10 = (u★ / κ) * log(FT(10) / z0_proxy)

    α = charnock_parameter(mag_u_10, spec.α_low, spec.α_high, spec.u_low, spec.u_high)

    # Smooth flow limit (Smith 1988)
    u★_safe = max(u★, eps(FT))
    z0_smooth = FT(0.11) * kinematic_visc / u★_safe

    # Rough flow limit (Charnock 1955)
    z0_rough = α * u★_safe^2 / grav

    return z0_smooth + z0_rough
end

"""
    scalar_roughness(spec::COARE3RoughnessParams, u★, sfc_param_set, roughness_inputs)

Calculate scalar roughness length using the COARE 3.0 algorithm (Fairall et al. 2003).

# Formulation
The scalar roughness length `z0s` is parameterized as an empirical fit to COARE and HEXOS data.
It limits `z0s` to a smooth flow limit (`1.1e-4` m) and decreases for rough flow following a 
power law of the roughness Reynolds number (`Re_star`):
```math
z_{0s} = \\min(1.1 \\times 10^{-4}, 5.5 \\times 10^{-5} \\cdot R_{e*}^{-0.6})
```
where ``R_{e*} = z_{0m} u_* / \\nu``.

# Dependencies
- `u★`: Friction velocity [m/s]
- `kinematic_visc`: Kinematic viscosity of air [m^2/s]
- `z0m`: Momentum roughness length (calculated internally)

# References
- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). 
    Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm.
    Journal of Climate, 16, 571–591.
    [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)
"""
@inline function scalar_roughness(
    spec::COARE3RoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    # Forward to combined calculation to avoid code duplication
    _, z0s = momentum_and_scalar_roughness(spec, u★, sfc_param_set, roughness_inputs)
    return z0s
end

@inline function momentum_and_scalar_roughness(
    spec::COARE3RoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    FT = eltype(sfc_param_set)
    z0m = momentum_roughness(spec, u★, sfc_param_set, roughness_inputs)
    kinematic_visc = spec.kinematic_visc
    u★_safe = max(u★, eps(FT))
    Re_star = z0m * u★_safe / kinematic_visc
    z0s = min(FT(1.1e-4), FT(5.5e-5) * Re_star^FT(-0.6))
    return (z0m, z0s)
end


# The plant area index of the roughness inputs, from the field `PAI` or from its
# deprecated alias `LAI`. `hasproperty` on a NamedTuple resolves at compile time, so the
# alias costs nothing in kernels.
@inline _plant_area_index(roughness_inputs) =
    hasproperty(roughness_inputs, :PAI) ? roughness_inputs.PAI : roughness_inputs.LAI

"""
    frontal_area_index(spec::RaupachRoughnessParams, plant_area_index)

Frontal area index `λ = max(frontal_area_ratio * plant_area_index, λ_min)` [m^2/m^2] of
the canopy, which sets the drag partition of the Raupach (1994) roughness model (Eq. 7).
"""
@inline function frontal_area_index(spec::RaupachRoughnessParams, plant_area_index)
    # Integer plant area indices are converted, so that the coefficients take a
    # floating-point type
    PAI = float(plant_area_index)
    FT = typeof(PAI)
    return max(
        float_parameter(FT, spec.frontal_area_ratio) * PAI,
        float_parameter(FT, spec.λ_min),
    )
end

"""
    canopy_area_index(spec::RaupachRoughnessParams, plant_area_index)

Canopy area index `Λ` [m^2/m^2] of the Raupach (1994) displacement height (Eq. 8): the
plant area index, raised to `λ_min / frontal_area_ratio` where the floor on the
[`frontal_area_index`](@ref) is active, so that both indices describe the same canopy.
"""
@inline function canopy_area_index(spec::RaupachRoughnessParams, plant_area_index)
    λ = frontal_area_index(spec, plant_area_index)
    return λ / float_parameter(typeof(λ), spec.frontal_area_ratio)
end

"""
    raupach_displacement_fraction(spec::RaupachRoughnessParams, plant_area_index)

Zero-plane displacement height as a fraction of the canopy height for a plant area index
(Raupach 1994, Eq. 8), with the [`canopy_area_index`](@ref) `Λ`:
```math
d / h = 1 - \\frac{1 - \\exp(-\\sqrt{c_{d1} Λ})}{\\sqrt{c_{d1} Λ}}
```
The fraction tends to zero as `Λ → 0`.
"""
@inline function raupach_displacement_fraction(
    spec::RaupachRoughnessParams,
    plant_area_index,
)
    Λ = canopy_area_index(spec, plant_area_index)
    FT = typeof(Λ)
    c_d1 = float_parameter(FT, spec.c_d1)
    # The floor keeps the fraction finite at Λ = 0, and expm1 keeps it accurate for
    # small Λ
    x = sqrt(max(c_d1 * Λ, eps(FT)))
    return 1 + expm1(-x) / x
end

"""
    raupach_roughness_fraction(spec::RaupachRoughnessParams, κ, plant_area_index)

Momentum roughness length as a fraction of the canopy height for a plant area index and
von Kármán constant `κ` (Raupach 1994, Eqs. 4, 5, 7, and 8):
```math
z_{0m} / h = (1 - d / h) \\exp(-κ U_h / u_* - Ψ_h), \\quad
u_* / U_h = \\min(\\sqrt{C_S + C_R λ}, (u_* / U_h)_{max}),
```
with the [`frontal_area_index`](@ref) `λ`, `d / h` from
[`raupach_displacement_fraction`](@ref), and the roughness-sublayer influence function
`Ψ_h = ln c_w - 1 + 1 / c_w`, the departure of the wind profile immediately above the
canopy from the logarithmic law. The cap on `u_* / U_h` is the sheltering limit, beyond
which `z0m / h` decreases with `λ` (for `λ > 0.29` with the default coefficients).

# Arguments
- `spec`: Coefficients, see [`RaupachRoughnessParams`](@ref).
- `κ`: Von Kármán constant [-].
- `plant_area_index`: Plant area index, the sum of the leaf and stem area indices
  [m^2/m^2].
"""
@inline function raupach_roughness_fraction(
    spec::RaupachRoughnessParams,
    κ,
    plant_area_index,
)
    λ = frontal_area_index(spec, plant_area_index)
    FT = typeof(λ)
    C_S = float_parameter(FT, spec.C_S)
    C_R = float_parameter(FT, spec.C_R)
    ustar_Uh_max = float_parameter(FT, spec.ustar_Uh_max)
    c_w = float_parameter(FT, spec.c_w)
    Ψ_h = log(c_w) - 1 + 1 / c_w
    ustar_over_Uh = min(sqrt(C_S + C_R * λ), ustar_Uh_max)
    return (1 - raupach_displacement_fraction(spec, plant_area_index)) *
           exp(-κ / ustar_over_Uh - Ψ_h)
end

"""
    displacement_height(spec::RaupachRoughnessParams, roughness_inputs)

Zero-plane displacement height [m] of the canopy, `h` times
[`raupach_displacement_fraction`](@ref) at the plant area index `roughness_inputs.PAI`.
Inside [`surface_fluxes`](@ref), the displacement height is the input `d`, which sets the
effective reference height `Δz - d`; callers compute it here from the same canopy inputs
as the roughness length.
"""
@inline function displacement_height(spec::RaupachRoughnessParams, roughness_inputs)
    return roughness_inputs.h *
           raupach_displacement_fraction(spec, _plant_area_index(roughness_inputs))
end

"""
    momentum_roughness(spec::RaupachRoughnessParams, u★, sfc_param_set, roughness_inputs)

Momentum roughness length [m] of the Raupach (1994) canopy roughness model, `h` times
[`raupach_roughness_fraction`](@ref) at the plant area index `roughness_inputs.PAI`, and
at least the fixed roughness length `z0m_fixed` of the parameter set (which also covers a
vanishing canopy height).

# Formulation
The model partitions the surface drag between the substrate (soil) and the roughness
elements (plants), which gives the friction velocity ratio at the canopy top
`u★ / U(h) = min((C_S + C_R λ)^(1/2), (u★ / U(h))_max)` for the frontal area index `λ`
(Eq. 7); the cap is the sheltering limit, beyond which `z0m / h` decreases with `λ` (for
`λ > 0.29`, i.e., `PAI > 0.58` with `frontal_area_ratio = 0.5`). The roughness length
follows from the wind profile at the canopy top with the roughness-sublayer influence
function `Ψ_h` (Eq. 4), set by the roughness-sublayer depth ratio `c_w` (Eq. 5).

# Dependencies
- `roughness_inputs.PAI`: Plant area index `Λ = LAI + SAI` (leaves plus stems) [m^2/m^2],
  from which the frontal area index is `λ = max(frontal_area_ratio * Λ, λ_min)`. `Λ`
  enters the displacement height (Eq. 8) directly and `λ` the drag partition (Eq. 7).
- `roughness_inputs.h`: Canopy height [m].
- `spec`: Coefficients, see [`RaupachRoughnessParams`](@ref).
- `z0m_fixed` and `von_karman_const` from `sfc_param_set`.

# References
- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane 
    displacement as functions of canopy height and area index.
    Boundary-Layer Meteorology, 71, 211–216.
    [DOI: 10.1007/BF00709229](https://doi.org/10.1007/BF00709229)
"""
@inline function momentum_roughness(
    spec::RaupachRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    κ = SFP.von_karman_const(sfc_param_set)
    PAI = _plant_area_index(roughness_inputs)
    z0m = roughness_inputs.h * raupach_roughness_fraction(spec, κ, PAI)
    return max(z0m, SFP.z0m_fixed(sfc_param_set))
end

"""
    scalar_roughness(spec::RaupachRoughnessParams, u★, sfc_param_set, roughness_inputs)

Calculate scalar roughness length from the momentum roughness length using a fixed Stanton number.

# Formulation
Scaled from the momentum roughness length using a fixed Stanton number:
```math
z_{0s} = z_{0m} \\cdot St
```
where ``St`` is `spec.stanton_number`.

# Dependencies
- `spec.stanton_number`: Stanton number specific to the canopy/surface type.
"""
@inline function scalar_roughness(
    spec::RaupachRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    z0m = momentum_roughness(spec, u★, sfc_param_set, roughness_inputs)
    return z0m * spec.stanton_number
end

@inline function momentum_and_scalar_roughness(
    spec::RaupachRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    z0m = momentum_roughness(spec, u★, sfc_param_set, roughness_inputs)
    return (z0m, z0m * spec.stanton_number)
end

"""
    depends_on_ustar(model)

Whether the roughness lengths of a roughness model, or the gustiness of a gustiness
model, depend on the friction velocity `u★`. [`compute_ustar_and_roughness`](@ref) finds
`u★` by root-finding when either model does; otherwise, the roughness lengths, the
gustiness, and `u★` follow directly from the stability parameter. The models independent
of `u★` are [`ConstantRoughnessParams`](@ref), [`RaupachRoughnessParams`](@ref),
[`ConstantGustinessSpec`](@ref), and [`FlooredDeardorffGustinessSpec`](@ref).
"""
@inline depends_on_ustar(::AbstractRoughnessParams) = true
@inline depends_on_ustar(::ConstantRoughnessParams) = false
@inline depends_on_ustar(::RaupachRoughnessParams) = false
@inline depends_on_ustar(::AbstractGustinessSpec) = true
@inline depends_on_ustar(::ConstantGustinessSpec) = false
@inline depends_on_ustar(::FlooredDeardorffGustinessSpec) = false

# =========================================================================================
# Iterative Solver for combined ustar and roughness
# =========================================================================================

struct UstarResidual{FT, PS, I, S} <: Function
    param_set::PS
    inputs::I
    scheme::S
    ζ::FT
end

function (ur::UstarResidual)(ustar)
    param_set = ur.param_set
    FT = eltype(param_set)
    inputs = ur.inputs
    scheme = ur.scheme
    ζ = ur.ζ

    # Ensure ustar is positive for physical consistency
    ustar_safe = max(ustar, FT(0))
    z0m, z0s = momentum_and_scalar_roughness(
        inputs.roughness_model,
        ustar_safe,
        param_set,
        inputs.roughness_inputs,
    )

    gustiness_val = gustiness_value(
        inputs.gustiness_model,
        param_set,
        ζ,
        ustar_safe,
        inputs,
        scheme,
    )
    ustar_calc = compute_ustar(param_set, ζ, z0m, inputs, scheme, gustiness_val)

    return ustar - ustar_calc
end

"""
    compute_ustar_and_roughness(param_set, ζ, inputs, scheme)

Computes friction velocity `ustar` and roughness lengths `z0m`, `z0h` for a given stability `ζ`.

- If `inputs.ustar` is prescribed, it is returned directly.
- If the roughness and gustiness models are independent of `ustar` (see
  [`depends_on_ustar`](@ref)), `z0m`, `z0h`, and `ustar` follow directly from `ζ`.
- Otherwise, three iterations of Brent's method on `ustar ∈ [1e-4, 4]` m/s find the
  `ustar` consistent with the roughness and gustiness models; the result lies in this
  bracket. If no consistent `ustar` lies in the bracket, the endpoint on the side of the
  root is returned: `4` m/s when the friction velocity implied by `ζ` and the gustiness it
  generates exceeds the bracket for every `ustar`, and `1e-4` m/s in calm conditions.
  With [`DeardorffGustinessSpec`](@ref), the gustiness at fixed `ζ` is proportional to
  `ustar`, and no consistent `ustar` exists for `ζ` more unstable than the free-convection
  limit; the solve for `ζ` then settles where a consistent `ustar` exists.
"""
function compute_ustar_and_roughness(
    param_set::APS,
    ζ,
    inputs,
    scheme,
)
    FT = eltype(param_set)

    if inputs.ustar !== nothing
        ustar = inputs.ustar
        z0m, z0s = momentum_and_scalar_roughness(
            inputs.roughness_model,
            ustar,
            param_set,
            inputs.roughness_inputs,
        )
        return ustar, z0m, z0s
    end

    # For roughness and gustiness independent of ustar, the residual below is linear in
    # ustar, and its root is the friction velocity implied by ζ. The branch condition is
    # known at compile time from the model types.
    if !(
        depends_on_ustar(inputs.roughness_model) ||
        depends_on_ustar(inputs.gustiness_model)
    )
        z0m, z0s = momentum_and_scalar_roughness(
            inputs.roughness_model,
            zero(FT),
            param_set,
            inputs.roughness_inputs,
        )
        gustiness_val = gustiness_value(
            inputs.gustiness_model,
            param_set,
            ζ,
            zero(FT),
            inputs,
            scheme,
        )
        ustar = compute_ustar(param_set, ζ, z0m, inputs, scheme, gustiness_val)
        return ustar, z0m, z0s
    end

    rf = UstarResidual(param_set, inputs, scheme, ζ)

    ustar_min = FT(1e-4) # Sufficient for very calm conditions
    ustar_max = FT(4.0)  # Sufficient for hurricane-force winds
    maxiter = 3          # Relatively small maxiter for performance

    # Set tolerance to 0 to force the solver to use exactly `maxiter` iterations.
    # This avoids branch divergence on GPUs, improving performance.
    rtol = FT(0)

    sol = RS.find_zero(
        rf,
        RS.BrentsMethod(ustar_min, ustar_max),
        RS.TwoPointSolution(),
        RS.RelativeSolutionTolerance(rtol),
        maxiter,
    )
    # Without a sign change in the bracket, the residual `ustar - ustar_calc` has one sign
    # throughout: negative when the friction velocity implied by ζ (with the gustiness it
    # generates) exceeds `ustar_max` for every `ustar`, positive when it stays below
    # `ustar_min` (calm). The endpoint on the side of the root is returned, so that the
    # effective wind seen by the ζ solve varies continuously. With Deardorff gustiness at
    # fixed ζ, the gustiness is proportional to ustar, and the first case occurs for ζ
    # more unstable than the free-convection limit, where the consistent ustar diverges;
    # the large ustar then gives a small state Richardson number and a negative ζ
    # residual, which steers the ζ solve back to the region with a consistent ustar.
    # `sol.root` alone is the endpoint of smaller residual, `ustar_min`, which lets the ζ
    # solve converge to roots with a vanishing friction velocity over rough surfaces.
    # The candidates are promoted to the type of the residual, which carries the
    # derivative information under AD; the endpoints, and `sol.root` without a sign
    # change, are plain floats. The type comes from inference of the residual (the
    # solution's field types are a union over the two cases).
    RT = Base.promote_op(rf, FT)
    no_sign_change = sol.y0 * sol.y1 > 0
    ustar = ifelse(
        no_sign_change,
        ifelse(sol.y1 < 0, convert(RT, ustar_max), convert(RT, ustar_min)),
        convert(RT, sol.root),
    )

    z0m, z0s = momentum_and_scalar_roughness(
        inputs.roughness_model,
        ustar,
        param_set,
        inputs.roughness_inputs,
    )

    return ustar, z0m, z0s
end
