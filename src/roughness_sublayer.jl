# Roughness sublayer (RSL) corrections for surface flux calculations.
#
# The roughness sublayer (RSL) is the region immediately above rough surfaces
# (canopies, urban elements) where standard Monin-Obukhov similarity theory (MOST)
# overestimates the wind and scalar gradients because turbulent mixing is enhanced
# by organized eddies shed from the roughness elements.
#
# These models correct the dimensionless profile F (as used in MOST) by adding a
# negative RSL correction P, so that F̂ = F + P < F. The reduced profile gives
# larger friction velocity u* = κ ΔU / F̂ and hence larger drag and scalar exchange.

abstract type AbstractRoughnessSubLayerModel end

Base.broadcastable(m::AbstractRoughnessSubLayerModel) = tuple(m)

"""
    NoRoughnessSubLayer <: AbstractRoughnessSubLayerModel

No roughness sublayer correction. Standard Monin-Obukhov similarity theory (MOST)
is applied without modification. This is the default when no RSL model is specified.
"""
struct NoRoughnessSubLayer <: AbstractRoughnessSubLayerModel end

"""
    PhysickGarrattRSL{FT} <: AbstractRoughnessSubLayerModel

Roughness sublayer correction following Physick & Garratt (1995), based on the
Raupach et al. (1991) mixing-layer analogy.

## Physics

Standard MOST assumes homogeneous turbulence throughout the surface layer, but above
rough surfaces (e.g., forests, urban canopies) a roughness sublayer (RSL) of depth
`z_RSL` develops where the turbulent diffusivity is enhanced relative to MOST predictions.
Within the RSL, the non-dimensional gradient function is reduced:

```
ϕ̂(z) = ϕ_MOST(ζ) · μ(z)   with μ(z) < 1 for z < d + z_RSL
```

This model uses a **linear-ramp enhancement function**:
```
μ(z) = 1 - c · (1 - (z - d) / z_RSL)   for z - d < z_RSL
μ(z) = 1                                  for z - d ≥ z_RSL
```

Integrating from `z₀` to `min(Δz_eff, z_RSL)` yields the closed-form RSL correction
(negative, reducing the profile):
```
P = -c · [ln(z_clip / z₀) - (z_clip - z₀) / z_RSL]
```
where `z_clip = min(Δz_eff, z_RSL)`. Separate coefficients `c_m` and `c_h` are used
for momentum and scalar (heat/humidity) transport respectively.

## Parameters

| Field    | Typical range | Meaning |
|----------|--------------|---------|
| `c_m`    | 0.1–0.6      | RSL correction strength for momentum |
| `c_h`    | 0.1–0.6      | RSL correction strength for scalars  |
| `z_RSL`  | 2–5 × h_c   | RSL depth above displacement height [m] |

## Effect on fluxes

The RSL correction increases the effective drag coefficient `Cd` and heat exchange
coefficient `Ch` relative to standard MOST, which is consistent with observations over
forests and urban surfaces where turbulent exchange is enhanced beyond the MOST prediction.

## References

- Physick, W. L., & Garratt, J. R. (1995). Incorporation of a high-roughness lower boundary
    into a mesoscale model for studies of dry deposition over complex terrain.
    Boundary-Layer Meteorology, 74, 55–71.
    [DOI: 10.1007/BF00715710](https://doi.org/10.1007/BF00715710)
- Raupach, M. R., Finnigan, J. J., & Brunet, Y. (1991). Coherent eddies and turbulence
    in vegetation canopies: the mixing-layer analogy. Boundary-Layer Meteorology, 60, 375–395.
    [DOI: 10.1007/BF00155877](https://doi.org/10.1007/BF00155877)
- Garratt, J. R. (1992). The Atmospheric Boundary Layer. Cambridge University Press.
- Harman, I. N., & Finnigan, J. J. (2007). A simple unified theory for flow in the canopy
    and roughness sublayer. Boundary-Layer Meteorology, 123, 339–363.
    [DOI: 10.1007/s10546-006-9145-6](https://doi.org/10.1007/s10546-006-9145-6)
"""
Base.@kwdef struct PhysickGarrattRSL{FT} <: AbstractRoughnessSubLayerModel
    c_m::FT = 0.4
    c_h::FT = 0.4
    z_RSL::FT = 10.0
end

"""
    HarmanFinniganRSL{FT} <: AbstractRoughnessSubLayerModel

Roughness sublayer correction following Harman & Finnigan (2007), based on the
mixing-layer analogy with an exponential modification to the dimensionless gradient.

## Physics

Within the roughness sublayer (depth `z_RSL` above the displacement height `d`), the
local non-dimensional gradient function is multiplied by an exponential correction:

```
ϕ̂(z) = ϕ_MOST(ζ) · exp(-c₁ · (1 - (z - d) / z_RSL))   for z - d ≤ z_RSL
ϕ̂(z) = ϕ_MOST(ζ)                                          for z - d > z_RSL
```

This gives an RSL correction factor that is exp(−c₁) at the surface and increases to 1
at the top of the RSL (z = d + z_RSL), consistent with enhanced turbulent diffusion near
the canopy and a smooth transition to standard MOST above.

Integrating from `z₀` to `z_clip = min(Δz_eff, z_RSL)` gives the RSL profile correction:

```
P = ∫_{z₀}^{z_clip} expm1(-c₁ + c₁·z/z_RSL) / z  dz   (≤ 0)
```

This integral has no elementary closed form (it involves the exponential integral Ei).
A change of variables `u = ln z` removes the `1/z` factor,
```
P = ∫_{ln z₀}^{ln z_clip} expm1(-c₁ + c₁·eᵘ/z_RSL) du,
```
and the result is evaluated by **4-point Gauss-Legendre quadrature** with the exact
algebraic nodes and weights (nested radicals). No tabulated constants or extra
dependencies are required; the evaluation is GPU-compatible and branch-free.

For small `c₁` the exponential reduces to a linear ramp (`exp(-x) ≈ 1 - x`), and the HF
correction converges to the Physick-Garratt (1995) correction with the same coefficient.
For any `c₁ > 0`, the HF correction is weaker than PG (less negative P) for the same
parameter value, because `exp(-x) > 1 - x` for `x > 0`.

Separate parameters `c1_m` and `c1_h` allow independent tuning for momentum and scalar
(heat/humidity) transport.

## Parameters

| Field    | Typical range  | Meaning |
|----------|----------------|---------|
| `c1_m`   | 0.3 – 0.7      | Exponential correction exponent for momentum |
| `c1_h`   | 0.3 – 0.7      | Exponential correction exponent for scalars  |
| `z_RSL`  | 0.5 – 1.5 h_c  | RSL depth above displacement height [m]      |

## References

- Harman, I. N., & Finnigan, J. J. (2007). A simple unified theory for flow in the canopy
    and roughness sublayer. Boundary-Layer Meteorology, 123, 339–363.
    [DOI: 10.1007/s10546-006-9145-6](https://doi.org/10.1007/s10546-006-9145-6)
- Raupach, M. R., Finnigan, J. J., & Brunet, Y. (1996). Coherent eddies and turbulence in
    vegetation canopies: the mixing-layer analogy. Boundary-Layer Meteorology, 78, 351–382.
    [DOI: 10.1007/BF00120941](https://doi.org/10.1007/BF00120941)
"""
Base.@kwdef struct HarmanFinniganRSL{FT} <: AbstractRoughnessSubLayerModel
    c1_m::FT = 0.5
    c1_h::FT = 0.5
    z_RSL::FT = 10.0
end

"""
    rsl_profile_correction(rsl_model, Δz_eff, z0, transport) -> P

Compute the roughness sublayer correction `P` (≤ 0) to the MOST dimensionless profile
`F` for the given transport type. The corrected profile is `F̂ = F + P`.

The correction is negative (enhances exchange) and saturates above the RSL height.

# Arguments
- `rsl_model`: RSL model (e.g. [`PhysickGarrattRSL`](@ref) or [`NoRoughnessSubLayer`](@ref)).
- `Δz_eff`: Effective height `z - d` above displacement height [m].
- `z0`: Roughness length for the transported quantity [m].
- `transport`: [`UF.MomentumTransport`](@ref) or [`UF.HeatTransport`](@ref).
"""
@inline function rsl_profile_correction(
    ::NoRoughnessSubLayer,
    Δz_eff,
    z0,
    transport,
)
    return zero(Δz_eff)
end

@inline function rsl_profile_correction(
    rsl_model::PhysickGarrattRSL,
    Δz_eff,
    z0,
    ::UF.MomentumTransport,
)
    return pg_rsl_correction(rsl_model.c_m, Δz_eff, z0, rsl_model.z_RSL)
end

@inline function rsl_profile_correction(
    rsl_model::PhysickGarrattRSL,
    Δz_eff,
    z0,
    ::UF.HeatTransport,
)
    return pg_rsl_correction(rsl_model.c_h, Δz_eff, z0, rsl_model.z_RSL)
end

@inline function rsl_profile_correction(
    rsl_model::HarmanFinniganRSL,
    Δz_eff,
    z0,
    ::UF.MomentumTransport,
)
    return hf_rsl_correction(rsl_model.c1_m, Δz_eff, z0, rsl_model.z_RSL)
end

@inline function rsl_profile_correction(
    rsl_model::HarmanFinniganRSL,
    Δz_eff,
    z0,
    ::UF.HeatTransport,
)
    return hf_rsl_correction(rsl_model.c1_h, Δz_eff, z0, rsl_model.z_RSL)
end

"""
    pg_rsl_correction(c, Δz_eff, z0, z_RSL) -> P

Physick-Garratt (1995) RSL correction (shared kernel for momentum and scalar).

For the linear-ramp RSL enhancement μ(z) = 1 - c(1 - z/z_RSL), integrating from
z₀ to z_clip = min(Δz_eff, z_RSL) gives:
```
P = -c · [ln(z_clip / z₀) - (z_clip - z₀) / z_RSL]
```

Properties:
- P = 0 at Δz_eff = z₀ (no correction at the roughness length)
- P is negative (reduces F̂ below F, enhancing exchange)
- P saturates at z_clip = z_RSL for Δz_eff > z_RSL
"""
@inline function pg_rsl_correction(c, Δz_eff, z0, z_RSL)
    FT = typeof(c)
    z0_safe = max(z0, eps(FT))
    z_clip = min(Δz_eff, z_RSL)
    # Clamp ratio to ≥ 1 so the logarithm is non-negative (physical bound)
    z_ratio = max(z_clip / z0_safe, FT(1))
    lin_term = max(z_clip - z0_safe, FT(0)) / z_RSL
    return -c * (log(z_ratio) - lin_term)
end

"""
Functor for the log-mapped Harman-Finnigan RSL integrand.

With `u = ln z`, the profile correction becomes
`∫ expm1(-c₁ + c₁·eᵘ/z_RSL) du`, so the `1/z` factor is absorbed into `du`.
"""
struct _HFLogIntegrand{FT}
    c1::FT
    c1_over_z_RSL::FT
end
@inline (∫::_HFLogIntegrand)(u) = expm1(-∫.c1 + ∫.c1_over_z_RSL * exp(u))

"""
    hf_rsl_correction(c1, Δz_eff, z0, z_RSL) -> P

Harman-Finnigan (2007) RSL correction via log-mapped 4-point Gauss-Legendre quadrature.

Computes
```
P = ∫_{z₀}^{z_clip} expm1(-c₁ + c₁·z/z_RSL) / z  dz
  = ∫_{ln z₀}^{ln z_clip} expm1(-c₁ + c₁·eᵘ/z_RSL) du   (≤ 0)
```
where `z_clip = min(Δz_eff, z_RSL)`. The log map removes the `1/z` singularity
structure; `expm1` avoids cancellation near the RSL top. Uses
[`gauss_legendre4`](@ref).

Properties:
- P = 0 when c₁ = 0 (no correction)
- P ≤ 0 for c₁ > 0 (enhances exchange, reduces F̂)
- P saturates when Δz_eff ≥ z_RSL (integrand is zero above RSL top)
- |P| < |P_PG(c₁)| for the same parameter value (since exp(-x) > 1-x for x > 0)
"""
@inline function hf_rsl_correction(c1, Δz_eff, z0, z_RSL)
    FT = typeof(c1)
    z0_safe = max(z0, eps(FT))
    z_clip = max(min(Δz_eff, z_RSL), z0_safe)
    return gauss_legendre4(
        _HFLogIntegrand(c1, c1 / z_RSL),
        log(z0_safe),
        log(z_clip),
    )
end
