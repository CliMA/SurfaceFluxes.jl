# Roughness sublayer (RSL) corrections for surface flux calculations.
#
# The roughness sublayer is the layer immediately above tall roughness elements
# (plant canopies, urban canopies) in which turbulent mixing is enhanced by
# organized eddies shed from the roughness elements, so that the dimensionless
# gradients are smaller than the Monin-Obukhov similarity theory (MOST) predicts:
#
#     φ̂(z) = φ(z / L) μ(z),     μ_min ≤ μ(z) ≤ 1,   μ(z) = 1 for z ≥ z_RSL,
#
# where z is the height above the displacement height d and z_RSL is the depth of
# the RSL above d. The RSL-corrected dimensionless profile is
#
#     F̂ = F + P,
#
# where F is the MOST profile and P ≥ 0 the RSL correction, the integral of
# φ (1 - μ) / z from z to the top of the RSL, evaluated here by Gauss-Legendre
# quadrature in ln z. Because it weights the RSL factor with φ(z/L), the correction
# is consistent with the stability dependence of the MOST profiles and with the
# neutral Prandtl number (φ_h(0) = Pr_0) for scalars.
#
# The corrected profile is anchored at the top of the RSL: it coincides with the
# MOST profile above the RSL, so z0 and d are the *apparent* roughness length and
# displacement height, as obtained from standard canopy relations or by fitting
# MOST profiles above the RSL (Physick & Garratt 1995, Eqs. 7 and 9; Harman &
# Finnigan 2007, 2008; Bonan 2019, Fig. 6.8).

Base.broadcastable(m::AbstractRoughnessSubLayerModel) = tuple(m)

"""
    NoRoughnessSubLayer <: AbstractRoughnessSubLayerModel

No roughness sublayer correction. Standard Monin-Obukhov similarity theory (MOST)
is applied without modification. This is the default when no RSL model is specified.
"""
struct NoRoughnessSubLayer <: AbstractRoughnessSubLayerModel end

@inline function _check_rsl_parameters(name, c_m, c_h, z_RSL, c_max)
    (0 <= c_m < c_max && 0 <= c_h < c_max) || throw(
        ArgumentError(
            "$name: RSL coefficients must satisfy 0 ≤ c < $c_max (got c_m = $c_m, c_h = $c_h)",
        ),
    )
    z_RSL >= 0 || throw(ArgumentError("$name: z_RSL must be ≥ 0 (got $z_RSL)"))
    return nothing
end

"""
    LinearRSL{FT} <: AbstractRoughnessSubLayerModel
    LinearRSL(; c_m = 0.4, c_h = 0.4, z_RSL = 10.0)
    LinearRSL(FT; kwargs...)

Roughness sublayer model with a linear RSL factor,

```
μ(z) = 1 - c (1 - z / z_RSL)   for z < z_RSL,     μ(z) = 1   for z ≥ z_RSL,
```

where `z` is the height above the displacement height `d`, so that the dimensionless
gradients are `φ̂(z) = φ(z / L) μ(z)`. The factor increases linearly from `1 - c` at
the displacement height to 1 at the top of the RSL. It is the first-order (small-`c`)
approximation of [`ExponentialRSL`](@ref).

The corrected profile coincides with the MOST profile above the RSL, so the roughness
length `z0` and displacement height `d` are the *apparent* values, as obtained from
standard canopy relations (e.g., `z0 ≈ 0.1 h`, `d ≈ 0.67 h`, or
[`RaupachRoughnessParams`](@ref)) or by fitting MOST profiles above the RSL. Within the
RSL, the wind speed and scalar differences are larger, and the exchange coefficients
smaller, than those of the MOST profile extrapolated downward with the apparent `z0`
and `d`.

# Fields
- `c_m`: RSL strength for momentum, `μ(0) = 1 - c_m`, with `0 ≤ c_m < 1` [-].
- `c_h`: RSL strength for scalars (heat, moisture), with `0 ≤ c_h < 1` [-].
- `z_RSL`: RSL depth above the displacement height, `≥ 0` [m].

The RSL top is typically at 2–3 canopy heights `h` above the ground (Garratt 1980;
Raupach et al. 1991), i.e., `z_RSL ≈ 1.3–2.3 h` for `d ≈ 0.67 h`; Physick & Garratt
(1995) use a depth of `50 z0`. `z_RSL` is a fixed depth, to be scaled with the canopy
height when the model is constructed.

`FT` is the floating-point type of the parameters (default `Float64`). In the flux
computation, the parameters are converted to the floating-point type of the inputs, so
`Float32` computations remain in `Float32` with the default constructors;
`LinearRSL(Float32; ...)` constructs a model with `Float32` parameters.

# Examples
```julia
h = 30.0  # canopy height [m]
rsl = LinearRSL(c_m = 0.4, c_h = 0.4, z_RSL = 2h - 0.67h)
rsl32 = LinearRSL(Float32; z_RSL = 20.0)
```

# References
- Garratt, J. R. (1980). Surface influence upon vertical profiles in the atmospheric
    near-surface layer. Quarterly Journal of the Royal Meteorological Society, 106, 803–819.
- Physick, W. L., & Garratt, J. R. (1995). Incorporation of a high-roughness lower boundary
    into a mesoscale model for studies of dry deposition over complex terrain.
    Boundary-Layer Meteorology, 74, 55–71.
- Raupach, M. R., Antonia, R. A., & Rajagopalan, S. (1991). Rough-wall turbulent boundary
    layers. Applied Mechanics Reviews, 44, 1–25.
"""
struct LinearRSL{FT} <: AbstractRoughnessSubLayerModel
    c_m::FT
    c_h::FT
    z_RSL::FT
    function LinearRSL{FT}(c_m, c_h, z_RSL) where {FT}
        _check_rsl_parameters("LinearRSL", c_m, c_h, z_RSL, 1)
        return new{FT}(c_m, c_h, z_RSL)
    end
end

LinearRSL(c_m::FT, c_h::FT, z_RSL::FT) where {FT} = LinearRSL{FT}(c_m, c_h, z_RSL)
LinearRSL(; c_m = 0.4, c_h = 0.4, z_RSL = 10.0) = LinearRSL(c_m, c_h, z_RSL)
LinearRSL(::Type{FT}; c_m = 0.4, c_h = 0.4, z_RSL = 10.0) where {FT} =
    LinearRSL(FT(c_m), FT(c_h), FT(z_RSL))

"""
    ExponentialRSL{FT} <: AbstractRoughnessSubLayerModel
    ExponentialRSL(; c_m = 0.7, c_h = 0.7, z_RSL = 10.0)
    ExponentialRSL(FT; kwargs...)

Roughness sublayer model with an exponential RSL factor,

```
μ(z) = exp(-c (1 - z / z_RSL))   for z < z_RSL,     μ(z) = 1   for z ≥ z_RSL,
```

where `z` is the height above the displacement height `d`, so that the dimensionless
gradients are `φ̂(z) = φ(z / L) μ(z)`. The factor increases from `exp(-c)` at the
displacement height to 1 at the top of the RSL. This is the form of Garratt (1980, 1983),
used by Physick & Garratt (1995) with `c = 0.7` (their `0.5 exp(0.7 z / z_RSL)`, as
`ln 2 ≈ 0.7`) and shown in Bonan (2019, Fig. 6.8). Physick & Garratt (1995) also anchor
the corrected profiles at the RSL top (their Eqs. 7 and 9), as done here.

The RSL factor `μ` is a prescribed function of height; the stability dependence of the
corrected profiles enters through `φ(z / L)`. As for [`LinearRSL`](@ref), the corrected
profile coincides with MOST above the RSL, so `z0` and `d` are the apparent roughness
length and displacement height.

# Fields
- `c_m`: RSL exponent for momentum, `μ(0) = exp(-c_m)`, with `c_m ≥ 0` [-].
- `c_h`: RSL exponent for scalars (heat, moisture), with `c_h ≥ 0` [-].
- `z_RSL`: RSL depth above the displacement height, `≥ 0` [m].

See [`LinearRSL`](@ref) for typical RSL depths and the floating-point type.

# Examples
```julia
h = 30.0  # canopy height [m]
rsl = ExponentialRSL(c_m = 0.7, c_h = 0.7, z_RSL = 2h - 0.67h)
```

# References
- Garratt, J. R. (1980). Surface influence upon vertical profiles in the atmospheric
    near-surface layer. Quarterly Journal of the Royal Meteorological Society, 106, 803–819.
- Physick, W. L., & Garratt, J. R. (1995). Incorporation of a high-roughness lower boundary
    into a mesoscale model for studies of dry deposition over complex terrain.
    Boundary-Layer Meteorology, 74, 55–71.
- Harman, I. N., & Finnigan, J. J. (2007). A simple unified theory for flow in the canopy
    and roughness sublayer. Boundary-Layer Meteorology, 123, 339–363.
- Bonan, G. (2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge University Press.
"""
struct ExponentialRSL{FT} <: AbstractRoughnessSubLayerModel
    c_m::FT
    c_h::FT
    z_RSL::FT
    function ExponentialRSL{FT}(c_m, c_h, z_RSL) where {FT}
        _check_rsl_parameters("ExponentialRSL", c_m, c_h, z_RSL, Inf)
        return new{FT}(c_m, c_h, z_RSL)
    end
end

ExponentialRSL(c_m::FT, c_h::FT, z_RSL::FT) where {FT} =
    ExponentialRSL{FT}(c_m, c_h, z_RSL)
ExponentialRSL(; c_m = 0.7, c_h = 0.7, z_RSL = 10.0) = ExponentialRSL(c_m, c_h, z_RSL)
ExponentialRSL(::Type{FT}; c_m = 0.7, c_h = 0.7, z_RSL = 10.0) where {FT} =
    ExponentialRSL(FT(c_m), FT(c_h), FT(z_RSL))

const ParametricRSL = Union{LinearRSL, ExponentialRSL}

# RSL coefficient for the transport type
@inline rsl_coefficient(m::ParametricRSL, ::UF.MomentumTransport) = m.c_m
@inline rsl_coefficient(m::ParametricRSL, ::UF.HeatTransport) = m.c_h

# 1 - μ as a function of x = max(1 - z / z_RSL, 0) ∈ [0, 1]
@inline rsl_one_minus_mu(::LinearRSL, c, x) = c * x
@inline rsl_one_minus_mu(::ExponentialRSL, c, x) = -expm1(-c * x)

"""
    _RSLIntegrand{LOG}(model, uf_params, transport, c, inv_z_RSL, inv_L)

Integrand `φ(z / L) (1 - μ(z))` of the RSL correction, as a function of `u = ln z`
(`LOG = true`, for integrals with respect to `dz / z = du`) or of `z` (`LOG = false`,
for integrals with respect to `dz`). Called from [`rsl_profile_correction`](@ref).
"""
struct _RSLIntegrand{LOG, M, UFP, TR, C, IZ, IL}
    model::M
    uf_params::UFP
    transport::TR
    c::C
    inv_z_RSL::IZ
    inv_L::IL
end
# The numeric fields keep their own types: the RSL coefficient and depth are in the type
# of the inputs, or dual numbers when differentiating with respect to the model
# parameters (see `float_parameter`), and `inv_L` is in the type of the inputs (a dual
# number under AD).
_RSLIntegrand{LOG}(model::M, uf_params::UFP, transport::TR, c::C, inv_z_RSL::IZ,
    inv_L::IL) where {LOG, M, UFP, TR, C, IZ, IL} =
    _RSLIntegrand{LOG, M, UFP, TR, C, IZ, IL}(model, uf_params, transport, c,
        inv_z_RSL, inv_L)
@inline function (f::_RSLIntegrand{LOG})(v) where {LOG}
    z = LOG ? exp(v) : v
    x = max(1 - z * f.inv_z_RSL, zero(z))
    return UF.phi(f.uf_params, z * f.inv_L, f.transport) * rsl_one_minus_mu(f.model, f.c, x)
end

# Two-panel 4-point Gauss-Legendre quadrature on [a, b]
@inline function _gl4_2panel(f::F, a, b) where {F}
    m = (a + b) / 2
    return gauss_legendre4(f, a, m) + gauss_legendre4(f, m, b)
end

# ∫_{z_lo}^{z_hi} φ (1 - μ) dz / z  (≥ 0), by quadrature in ln z
@inline function _rsl_integral(m, uf_params, transport, c, inv_z_RSL, inv_L, z_lo, z_hi)
    f = _RSLIntegrand{true}(m, uf_params, transport, c, inv_z_RSL, inv_L)
    return _gl4_2panel(f, log(z_lo), log(z_hi))
end

# ∫_{z_lo}^{z_hi} φ (1 - μ) dz  (≥ 0), by quadrature in z
@inline function _rsl_integral_linear(
    m,
    uf_params,
    transport,
    c,
    inv_z_RSL,
    inv_L,
    z_lo,
    z_hi,
)
    f = _RSLIntegrand{false}(m, uf_params, transport, c, inv_z_RSL, inv_L)
    return _gl4_2panel(f, z_lo, z_hi)
end

"""
    rsl_profile_correction(uf_params, rsl_model, Δz_eff, ζ, z0, transport, scheme = PointValueScheme())

Return the roughness sublayer correction `P = F̂ - F ≥ 0` to the MOST dimensionless
profile `F` for the given transport type, stability parameter `ζ = Δz_eff / L`, and
discretization scheme, where `F̂` is the RSL-corrected profile
(see [`rsl_corrected_profile`](@ref)).

For point values, with `z_c = min(Δz_eff, z_RSL)` (limited to `[z0, max(z_RSL, z0)]`),
```
P = ∫_{z_c}^{z_RSL} φ(z/L) (1 - μ(z)) dz/z,
```
which vanishes above the RSL. For layer averages ([`LayerAverageScheme`](@ref)), `P` is
the layer average `(1/Δz_eff) ∫_{z0}^{Δz_eff} P_point(z) dz` of the point-value
correction, consistent with the layer-averaged MOST profile. The integrals are evaluated
by Gauss-Legendre quadrature in `ln z` (in `z` for the part of the layer average that is
an integral with respect to `z`).

# Arguments
- `uf_params`: Universal function parameters.
- `rsl_model`: RSL model ([`LinearRSL`](@ref), [`ExponentialRSL`](@ref), or
  [`NoRoughnessSubLayer`](@ref)).
- `Δz_eff`: Effective height `Δz - d` above the displacement height [m].
- `ζ`: Stability parameter `Δz_eff / L` [-].
- `z0`: Roughness length for the transported quantity [m].
- `transport`: [`UF.MomentumTransport`](@ref) or [`UF.HeatTransport`](@ref).
- `scheme`: [`PointValueScheme`](@ref) (default) or [`LayerAverageScheme`](@ref).

# Returns
The correction `P` [-], zero for [`NoRoughnessSubLayer`](@ref) and above the RSL.
"""
@inline rsl_profile_correction(
    uf_params,
    ::NoRoughnessSubLayer,
    Δz_eff,
    ζ,
    z0,
    transport,
    scheme = UF.PointValueScheme(),
) = zero(ζ)

@inline function rsl_profile_correction(
    uf_params,
    m::ParametricRSL,
    Δz_eff,
    ζ,
    z0,
    transport,
    scheme = UF.PointValueScheme(),
)
    FT = eltype(ζ)
    # Model parameters in the type of the inputs; dual parameters are kept
    c = float_parameter(FT, rsl_coefficient(m, transport))
    z_RSL = float_parameter(FT, m.z_RSL)
    z0_safe = max(z0, eps(FT))
    z_top = max(z_RSL, z0_safe)                      # RSL top (≥ z0)
    z_c = max(min(Δz_eff, z_top), z0_safe)           # z_c ∈ [z0, z_top]
    inv_z_RSL = 1 / z_RSL                            # Inf for z_RSL = 0 (then μ ≡ 1)
    inv_L = ζ / Δz_eff
    P = _rsl_correction(scheme, m, uf_params, transport, c, inv_z_RSL, inv_L,
        Δz_eff, z0_safe, z_c, z_top)
    # The exact correction is ≥ 0 (μ ≤ 1); the bound guards against quadrature error
    return max(P, zero(P))
end

"""
    rsl_corrected_profile(uf_params, rsl_model, Δz_eff, ζ, z0, transport, scheme = PointValueScheme())

Return the RSL-corrected dimensionless profile `F̂ = F + P`, where `F` is the MOST
profile `UF.dimensionless_profile` and `P` the roughness sublayer correction
(see [`rsl_profile_correction`](@ref), also for the arguments). All exchange
coefficients, similarity scales, the bulk Richardson number, and profile recovery use
this function, so that the RSL correction is applied consistently.

Since `φ̂ = φ μ` with `μ ≤ 1`, the exact corrected profile satisfies `F̂ ≥ F`; this bound
is enforced (`P ≥ 0`) to guard against quadrature error, so `F̂ > 0` whenever `F > 0`.
"""
@inline rsl_corrected_profile(
    uf_params,
    ::NoRoughnessSubLayer,
    Δz_eff,
    ζ,
    z0,
    transport,
    scheme = UF.PointValueScheme(),
) = UF.dimensionless_profile(uf_params, Δz_eff, ζ, z0, transport, scheme)

@inline rsl_corrected_profile(
    uf_params,
    m::ParametricRSL,
    Δz_eff,
    ζ,
    z0,
    transport,
    scheme = UF.PointValueScheme(),
) =
    UF.dimensionless_profile(uf_params, Δz_eff, ζ, z0, transport, scheme) +
    rsl_profile_correction(uf_params, m, Δz_eff, ζ, z0, transport, scheme)

# Point values: P = ∫_{z_c}^{z_top} φ (1 - μ) dz/z
@inline function _rsl_correction(::UF.PointValueScheme, m,
    uf_params, transport, c, inv_z_RSL, inv_L, Δz_eff, z0, z_c, z_top)
    return _rsl_integral(m, uf_params, transport, c, inv_z_RSL, inv_L, z_c, z_top)
end

# Layer averages:
# P = (1/Δz) ∫_{z0}^{Δz} ∫_{min(z, z_top)}^{z_top} φ (1 - μ) dz'/z' dz
#   = ∫_{z0}^{z_c} φ (1 - μ) (z' - z0)/Δz dz'/z'
#     + ∫_{z_c}^{z_top} φ (1 - μ) (Δz - z0)/Δz dz'/z'
# The weight (z' - z0)/z' of the first integral varies like exp(u) in u = ln z', which
# the quadrature in ln z' resolves poorly when z_c/z0 is large. It is therefore split
# into (1/Δz) ∫ φ (1 - μ) dz' (quadrature in z') and -(z0/Δz) ∫ φ (1 - μ) dz'/z'
# (quadrature in ln z'), whose integrands are smooth in their quadrature variables.
@inline function _rsl_correction(::UF.LayerAverageScheme, m,
    uf_params, transport, c, inv_z_RSL, inv_L, Δz_eff, z0, z_c, z_top)
    inv_Δz = 1 / Δz_eff
    lower_lin = _rsl_integral_linear(m, uf_params, transport, c, inv_z_RSL, inv_L, z0, z_c)
    lower_log = _rsl_integral(m, uf_params, transport, c, inv_z_RSL, inv_L, z0, z_c)
    upper = _rsl_integral(m, uf_params, transport, c, inv_z_RSL, inv_L, z_c, z_top)
    return (lower_lin - z0 * lower_log) * inv_Δz +
           max(Δz_eff - z0, zero(z0)) * inv_Δz * upper
end
