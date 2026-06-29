"""
    UniversalFunctions

Universal stability and stability correction functions for the `SurfaceFluxes` module.

Supports the following flux-profile (`ϕ`, `ψ`, `Ψ`) parameterizations:
 - `Businger`: Businger et al. (1971), Dyer (1974)
 - `Gryanik`: Gryanik et al. (2020)
 - `Grachev`: Grachev et al. (2007)

Both standard finite-difference (point-value) and finite-volume (layer-averaged) schemes are
supported; the finite-volume scheme follows Nishizawa & Kitamura (2018).

The module also provides variance/TKE similarity functions (`MomentumVariance`, `HeatVariance`)
from Panofsky et al. (1977), Wyngaard et al. (1971), and Tan et al. (2018). These are empirical
surface-layer closures that are **independent of the flux-profile parameterization** above (see
the `phi(..., MomentumVariance())` / `phi(..., HeatVariance())` docstrings for their range of
validity).

# References
 - Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971). Flux-profile
   relationships in the atmospheric surface layer. Journal of the Atmospheric Sciences, 28,
   181–189. [DOI: 10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2](https://doi.org/10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2)
 - Dyer, A. J. (1974). A review of flux-profile relationships. Boundary-Layer Meteorology, 7,
   363–372. [DOI: 10.1007/BF00240838](https://doi.org/10.1007/BF00240838)
 - Gryanik, V. M., Lüpkes, C., Grachev, A., & Sidorenko, D. (2020). New modified and extended
   stability functions for the stable boundary layer based on SHEBA and parametrizations of
   bulk transfer coefficients for climate models. Journal of the Atmospheric Sciences, 77,
   2687–2716. [DOI: 10.1175/JAS-D-19-0255.1](https://doi.org/10.1175/JAS-D-19-0255.1)
 - Grachev, A. A., Andreas, E. L., Fairall, C. W., Guest, P. S., & Persson, P. O. G. (2007).
   SHEBA flux–profile relationships in the stable atmospheric boundary layer. Boundary-Layer
   Meteorology, 124, 315–333. [DOI: 10.1007/s10546-007-9177-6](https://doi.org/10.1007/s10546-007-9177-6)
 - Nishizawa, S., & Kitamura, Y. (2018). A surface flux scheme based on the Monin-Obukhov
   similarity for finite volume models. Journal of Advances in Modeling Earth Systems, 10,
   3159–3175. [DOI: 10.1029/2018MS001534](https://doi.org/10.1029/2018MS001534)
 - Panofsky, H. A., Tennekes, H., Lenschow, D. H., & Wyngaard, J. C. (1977). The characteristics
   of turbulent velocity components in the surface layer under convective conditions.
   Boundary-Layer Meteorology, 11, 355–361. [DOI: 10.1007/BF02186086](https://doi.org/10.1007/BF02186086)
 - Wyngaard, J. C., Coté, O. R., & Izumi, Y. (1971). Local free convection, similarity, and the
   budgets of shear stress and heat flux. Journal of the Atmospheric Sciences, 28, 1171–1182.
   [DOI: 10.1175/1520-0469(1971)028<1171:LFCSAT>2.0.CO;2](https://doi.org/10.1175/1520-0469(1971)028<1171:LFCSAT>2.0.CO;2)
 - Panofsky, H. A., & Dutton, J. A. (1984). Atmospheric Turbulence: Models and Methods for
   Engineering Applications. Wiley, New York, 397 pp.
 - Tan, Z., Kaul, C. M., Pressel, K. G., Cohen, Y., Schneider, T., & Teixeira, J. (2018). An
   extended eddy-diffusivity mass-flux scheme for unified representation of subgrid-scale
   turbulence and convection. Journal of Advances in Modeling Earth Systems, 10, 770–800.
   [DOI: 10.1002/2017MS001162](https://doi.org/10.1002/2017MS001162)
"""
module UniversalFunctions

abstract type AbstractUniversalFunctionParameters{FT <: Real} end
const AUFP = AbstractUniversalFunctionParameters

#####
##### Solver Schemes
#####

"""
    SolverScheme

Abstract type for surface flux solver schemes.
"""
abstract type SolverScheme end

"""
    LayerAverageScheme <: SolverScheme

Finite volume approximation scheme following Nishizawa & Kitamura (2018).

Supported with `BusingerParams` and `GryanikParams`. Not supported with
`GrachevParams` (use `PointValueScheme` instead).
"""
struct LayerAverageScheme <: SolverScheme end

"""
    PointValueScheme <: SolverScheme

Standard finite difference scheme using point values.
"""
struct PointValueScheme <: SolverScheme end


#####
##### Interface
#####

abstract type AbstractTransportType end

"""
    MomentumTransport

Type selecting momentum-transfer stability functions (ϕₘ, ψₘ, Ψₘ).
"""
struct MomentumTransport <: AbstractTransportType end

"""
    HeatTransport

Type selecting heat/scalar-transfer stability functions (ϕₕ, ψₕ, Ψₕ).
"""
struct HeatTransport <: AbstractTransportType end

"""
    MomentumVariance

Type selecting momentum variance stability functions (ϕ_σu).
"""
struct MomentumVariance <: AbstractTransportType end

"""
    HeatVariance

Type selecting heat/scalar variance stability functions (ϕ_σθ).
"""
struct HeatVariance <: AbstractTransportType end

Base.broadcastable(tt::AbstractTransportType) = tuple(tt)
Base.broadcastable(p::AbstractUniversalFunctionParameters) = tuple(p)

"""
    phi(p, ζ, transport)

Universal (similarity) function: the non-dimensional vertical gradient of wind
shear (`ϕ_m`, [`MomentumTransport`](@ref)) or of temperature/scalars (`ϕ_h`,
[`HeatTransport`](@ref)) at stability parameter `ζ`.

Dispatches on the parameterization type of `p` ([`BusingerParams`](@ref),
[`GryanikParams`](@ref), [`GrachevParams`](@ref)) and on `transport`.
"""
function phi end

"""
    phi(p, ζ, u_star, w_star, transport)

Extended similarity function allowing dependence on velocity scales (u_*, w_*).
Default implementation falls back to `phi(p, ζ, transport)`.
"""
@inline phi(p, ζ, u_star, w_star, transport) = phi(p, ζ, transport)

"""
    psi(p, ζ, transport_type)

The standard integrated stability correction function `ψ(ζ)`.
Defined as:
    ψ(ζ) = ∫`[0 to ζ]` (ϕ(0) - ϕ(x)) / x dx

This is the standard correction used in point-based Monin-Obukhov
similarity theory. 
"""
function psi end

"""
    Psi(p, ζ, transport_type)

The volume-averaged stability correction function `Ψ(ζ)`.
Mathematically, this is defined as:
    Ψ(ζ) = (1/ζ) ∫`[0 to ζ]` ψ(x) dx

This function is required for finite-volume models where fluxes are
calculated using cell-averaged values rather than point values at the
cell center.

See Nishizawa & Kitamura (2018), Eqs. 14 & 15.
"""
function Psi end

#####
##### Forwarding methods for free parameters
#####

# Parameter-struct accessors (allow calling phi/psi/Psi with params directly)
Pr_0(p::AUFP) = p.Pr_0
a_m(p::AUFP) = p.a_m
a_h(p::AUFP) = p.a_h
b_m(p::AUFP) = p.b_m
b_h(p::AUFP) = p.b_h
c_h(p::AUFP) = p.c_h
ζ_a(p::AUFP) = p.ζ_a
γ(p::AUFP) = p.γ


#####
##### Private Helpers (Unstable Businger Logic)
#####

# Recycled logic for Businger unstable regimes to avoid code duplication
# and ensure consistency across Businger, Gryanik, and Grachev formulations.

@inline function _phi_m_unstable(ζ, γ)
    FT = eltype(ζ)
    return FT(1) / sqrt(sqrt(FT(1) - γ * ζ))
end

@inline function _phi_h_unstable(ζ, γ)
    FT = eltype(ζ)
    return FT(1) / sqrt(FT(1) - γ * ζ)
end

@inline function _psi_m_unstable(ζ, γ)
    FT = eltype(ζ)
    x = sqrt(sqrt(FT(1) - γ * ζ))
    log_term = log((FT(1) + x)^2 * (FT(1) + x^2) / FT(8))
    return log_term - FT(2) * atan(x) + FT(π) / FT(2)
end

@inline function _psi_h_unstable(ζ, γ)
    FT = eltype(ζ)
    y = sqrt(FT(1) - γ * ζ)
    return FT(2) * log((FT(1) + y) / FT(2))
end

"""
    _Psi_m_unstable(ζ, γ)

Finite-volume integral for momentum transfer in unstable conditions.
Derivation follows Nishizawa & Kitamura (2018, Eq. A5), but corrects the 
coefficient in the cubic term. 

Eq. A5 in the paper uses a hardcoded denominator of `12ζ`, which is only 
valid for γ = 16. For the general case, the denominator is `3γζ/4`.

Note: `γ` here refers to the **coefficient inside the sqrt/cbrt** 
(e.g., 15 in `(1 - 15ζ)^(-1/4)`), not the linear coefficient `a_m`.
"""
@inline function _Psi_m_unstable(ζ, γ)
    FT = eltype(ζ)
    # Small-ζ limit (Nishizawa & Kitamura 2018, Eq. A13)
    val_small = -γ * ζ / FT(8)

    # Safe ζ for large-ζ limit to avoid division by zero in cubic_term
    ζ_safe = min(ζ, -eps(FT))

    # Full computation (Nishizawa & Kitamura 2018, Eq. A5)
    # Ψ_m = log[(1+x)²(1+x²)/8] - 2tan⁻¹(x) + π/2 - 1 + (1-x³)/(3γζ/4)

    # 1. Compute x = (1 - γζ)^(1/4)
    x = sqrt(sqrt(FT(1) - γ * ζ_safe))

    # 2. Logarithmic term: log[(1+x)²(1+x²)/8]
    # (Combination of 2*log((1+x)/2) + log((1+x^2)/2))
    log_term = log((FT(1) + x)^2 * (FT(1) + x^2) / FT(8))

    # 3. Angular terms: -2tan⁻¹(x) + π/2
    π_term = FT(π) / FT(2)
    tan_term = FT(2) * atan(x)

    # 4. Cubic term: (1 - x³) / (3γζ/4)
    # 
    # Optimization: We use expm1/log1p for precision near ζ=0.
    # Naive calculation (1 - x^3) suffers from catastrophic cancellation when x ≈ 1.
    # Identity: 1 - x³ = 1 - (1 - γζ)^(3/4) 
    #                  = 1 - exp(0.75 * log(1 - γζ))
    #                  = -expm1(0.75 * log1p(-γ * ζ_safe))
    cubic_num = -expm1(FT(0.75) * log1p(-γ * ζ_safe))

    # Denominator Correction: 
    # The paper (Eq. A5) hardcodes the denominator as 12ζ (assuming γ=16).
    # The general form for arbitrary γ is 3γζ/4.
    cubic_denom = (FT(3) * γ / FT(4)) * ζ_safe
    cubic_term = cubic_num / cubic_denom

    # Final Sum: log_term - tan_term + π_term - 1 + cubic_term
    val_large = log_term - tan_term + π_term - FT(1) + cubic_term

    return ifelse(abs(ζ) < eps(FT), val_small, val_large)
end

"""
    _Psi_h_unstable(ζ, γ)

Finite-volume integral for unstable heat/scalar.
Matches Nishizawa & Kitamura (2018, Eq. A6).

Note: `γ` here refers to the **coefficient inside the sqrt**
(e.g., 9 in `(1 - 9ζ)^(-1/2)`), not the linear coefficient `a_h`.
"""
@inline function _Psi_h_unstable(ζ, γ)
    FT = eltype(ζ)
    # Small-ζ limit (Nishizawa & Kitamura 2018, Eq. A14)
    val_small = -γ * ζ / FT(4)

    # Safe ζ for large-ζ limit to avoid division by zero
    ζ_safe = min(ζ, -eps(FT))

    # Full computation (Nishizawa & Kitamura 2018, Eq. A6)
    # Using expm1 for precision near ζ=0
    # y = (1 - γζ)^(1/2)
    # 1 - y = 1 - (1 - γζ)^(1/2) = -expm1(0.5 * log1p(-γζ))

    y = sqrt(FT(1) - γ * ζ_safe)
    log_term = FT(2) * log((FT(1) + y) / FT(2))

    # Optimized linear term
    lin_num = -expm1(FT(0.5) * log1p(-γ * ζ_safe))
    lin_term = FT(2) * lin_num / (γ * ζ_safe)

    val_large = log_term + lin_term - FT(1)

    return ifelse(abs(ζ) < eps(FT), val_small, val_large)
end

#####
##### Businger
#####
"""
    BusingerParams{FT}

Parameter bundle for the Businger (1971) similarity relations.
Mappings to Nishizawa & Kitamura (2018) coefficients:
 - `a_m`, `a_h`: The linear coefficients for stable conditions (β in some texts).
 - `b_m`, `b_h`: The coefficients γ inside the unstable sqrt/cbrt terms (e.g., (1 - γζ)).
 - `Pr_0`: The neutral Prandtl number.

See Businger et al. (1971) and Nishizawa & Kitamura (2018).
"""
Base.@kwdef struct BusingerParams{FT} <: AbstractUniversalFunctionParameters{FT}
    Pr_0::FT
    a_m::FT
    a_h::FT
    b_m::FT
    b_h::FT
    ζ_a::FT
    γ::FT
end

"""
    phi(p::BusingerParams, ζ, ::MomentumTransport)

Businger momentum similarity `ϕ_m`.

# References
 - Stable (ζ >= 0): Eq. A1 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A1 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function phi(p::BusingerParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    # Ensure invalid branch inputs (ζ >= 0 for unstable, ζ < 0 for stable) don't error
    # even though their result is discarded by ifelse.
    return ifelse(
        ζ < 0,
        _phi_m_unstable(min(ζ, FT(0)), b_m(p)),
        FT(a_m(p)) * max(ζ, FT(0)) + FT(1),
    )
end

"""
    phi(p::BusingerParams, ζ, ::HeatTransport)

Businger heat/scalar-gradient similarity `ϕ_h`.

Returns the non-dimensional gradient such that `ϕ_h(0) = Pr_0` in the neutral limit.

# References
 - Stable (ζ >= 0): Eq. A2 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A2 (L < 0) in Nishizawa & Kitamura (2018), scaled by `Pr_0`.
"""
@inline function phi(p::BusingerParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    _Pr_0 = FT(Pr_0(p))
    _a_h = FT(a_h(p))
    return ifelse(
        ζ < 0,
        _Pr_0 * _phi_h_unstable(min(ζ, FT(0)), b_h(p)),
        _Pr_0 + _a_h * max(ζ, FT(0)),
    )
end

"""
    psi(p::BusingerParams, ζ, ::MomentumTransport)

Businger momentum stability correction `ψ_m`.

# References
 - Stable (ζ >= 0): Eq. A3 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A3 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::BusingerParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    return ifelse(
        ζ < 0,
        _psi_m_unstable(min(ζ, FT(0)), b_m(p)),
        -FT(a_m(p)) * max(ζ, FT(0)),
    )
end

"""
    psi(p::BusingerParams, ζ, ::HeatTransport)

Businger heat/scalar stability correction `ψ_h`.

# References
 - Stable (ζ >= 0): Eq. A4 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A4 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::BusingerParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    _Pr_0 = FT(Pr_0(p))
    _a_h = FT(a_h(p))
    return ifelse(
        ζ < 0,
        _Pr_0 * _psi_h_unstable(min(ζ, FT(0)), b_h(p)),
        -_a_h * max(ζ, FT(0)),
    )
end

"""
    Psi(p::BusingerParams, ζ, ::MomentumTransport)

Volume-averaged Businger momentum stability correction `Ψ_m`.

# References
 - Stable (ζ >= 0): Eqs. A5 and A13 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A5 (L < 0) in Nishizawa & Kitamura (2018).
 - Small ζ limit: Eq. A13 (ζ < 0) in Nishizawa & Kitamura (2018).
"""
@inline function Psi(p::BusingerParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    return ifelse(
        ζ >= 0,
        -FT(a_m(p)) * max(ζ, FT(0)) / FT(2),
        _Psi_m_unstable(min(ζ, FT(0)), b_m(p)),
    )
end

"""
    Psi(p::BusingerParams, ζ, ::HeatTransport)

Volume-averaged Businger heat/scalar stability correction `Ψ_h`.

# References
 - Stable (ζ >= 0): Eqs. A6 and A14 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A6 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function Psi(p::BusingerParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    _a_h = FT(a_h(p))
    _Pr_0 = FT(Pr_0(p))
    return ifelse(
        ζ >= 0,
        -_a_h * max(ζ, FT(0)) / FT(2),
        _Pr_0 * _Psi_h_unstable(min(ζ, FT(0)), b_h(p)),
    )
end

# --- Variance Functions (parameterization-independent) ---
#
# The variance/TKE similarity functions use empirical constants from the literature
# (Panofsky et al. 1977, Wyngaard et al. 1971, Tan et al. 2018) that do not depend on
# the stability-function parameterization. They are therefore defined once on the abstract
# parameter type `AUFP` and shared by Businger, Gryanik, and Grachev.
#
# IMPORTANT (range of validity): these are *convective* surface-layer closures with constant
# neutral values on the stable side. Neither Grachev et al. (2007) nor Gryanik et al. (2020)
# defines variance functions, so selecting those parameterizations does NOT change the
# variances returned here. In the stable boundary layer the constant values are a crude
# approximation: observed σ_u/u_* increases with stability rather than staying constant, and
# Monin-Obukhov scaling of the horizontal-velocity variances breaks down. Treat the stable-side
# variances as rough estimates, not validated stable-boundary-layer similarity.

"""
    phi(p::AUFP, ζ, ::MomentumVariance)

Streamwise velocity variance similarity `ϕ_σu = σ_u / u_*`.

This is a parameterization-independent convective closure; see the note below. The exported
[`surface_tke`](@ref SurfaceFluxes.surface_tke) uses the TKE form (the 5-argument method) instead.

# References
 - Unstable (ζ < 0): Panofsky et al. (1977). In the original, the argument is the mixed-layer
   stability `z_i / L` (the horizontal variances scale with boundary-layer depth, not local
   height); this method is evaluated at the local `ζ` passed by the caller.
 - Stable (ζ >= 0): neutral-limit constant (2.3; Panofsky & Dutton 1984). Observations show
   `σ_u/u_*` actually *increases* with stability, so this constant is only a rough estimate.

See the module docstring for full references.
"""
@inline function phi(p::AUFP, ζ, ::MomentumVariance)
    FT = eltype(ζ)
    # Panofsky et al. (1977) for unstable: (12 - 0.5 * ζ)^(1/3)
    # Neutral limit for stable: 2.3
    # Safe input for unstable: min(ζ, 0)
    return ifelse(
        ζ < 0,
        cbrt(FT(12) - FT(0.5) * min(ζ, FT(0))),
        FT(2.3),
    )
end

"""
    phi(p::AUFP, ζ, u_star, w_star, ::MomentumVariance)

Turbulent-kinetic-energy similarity `sqrt(TKE) / u_*`, the EDMF surface-layer TKE closure of
Tan et al. (2018), Eq. 22: `TKE = 3.75 u_*^2 + 0.2 w_*^2 + u_*^2 (-ζ)^{2/3}` for `ζ < 0`,
reducing to `3.75 u_*^2` for `ζ >= 0`. This is the form used by [`surface_tke`](@ref SurfaceFluxes.surface_tke).

Parameterization-independent (see the variance-section note above): it is a surface boundary
condition, not a stable-boundary-layer variance similarity. The stable-side value is a constant.

# References
 - Tan et al. (2018). See the module docstring for the full reference.
"""
@inline function phi(p::AUFP, ζ, u_star, w_star, ::MomentumVariance)
    FT = eltype(ζ)
    # Tan et al. (2018) Eq. 22 for unstable
    # TKE = 3.75 u_*^2 + 0.2 w_*^2 + u_*^2 * (-ζ)^(2/3)
    w_ratio = w_star / u_star
    tke_norm = FT(3.75) + FT(0.2) * w_ratio^2 + cbrt(-min(ζ, FT(0)))^2

    return ifelse(
        ζ < 0,
        sqrt(tke_norm),
        sqrt(FT(3.75)),
    )
end

"""
    phi(p::AUFP, ζ, ::HeatVariance)

Temperature variance similarity `ϕ_σθ = σ_θ / |θ_*|`.

Parameterization-independent (see the variance-section note above).

# References
 - Unstable (ζ < 0): free-convection form of Wyngaard et al. (1971), with constants from
   Tan et al. (2018): `2 (1 - 8.3 ζ)^{-1/3}`.
 - Stable (ζ >= 0): constant (2.0). Of the variance functions this is the best supported on the
   stable side: the non-dimensional temperature standard deviation approaches a constant in the
   very stable (z-less) limit, though the value here is not calibrated to stable-boundary-layer
   data. See the module docstring for full references.
"""
@inline function phi(p::AUFP, ζ, ::HeatVariance)
    FT = eltype(ζ)
    # Tan et al. (2018) for unstable: 2 * (1 - 8.3ζ)^(-1/3)
    # Stable: 2.0
    return ifelse(
        ζ < 0,
        FT(2) * cbrt(FT(1) / (FT(1) - FT(8.3) * min(ζ, FT(0)))),
        FT(2.0),
    )
end

#####
##### Gryanik
#####

"""
    GryanikParams{FT}

Parameter bundle for the Gryanik et al. (2020) similarity relations.
These functions are designed to be valid across the entire stability range,
including very stable conditions.

 - `a_m`, `b_m`: Coefficients for momentum stability function (Eq. 32).
 - `a_h`, `b_h`: Coefficients for heat stability function (Eq. 33).
 - `Pr_0`: Neutral Prandtl number. (The paper recommends Pr_0 ≈ 0.98.)

Reference: Gryanik et al. (2020).
"""
Base.@kwdef struct GryanikParams{FT} <: AbstractUniversalFunctionParameters{FT}
    Pr_0::FT
    a_m::FT
    a_h::FT
    b_m::FT
    b_h::FT
    b_m_unstable::FT
    b_h_unstable::FT
    ζ_a::FT
    γ::FT
end

"""
    phi(p::GryanikParams, ζ, ::MomentumTransport)

Gryanik momentum similarity `ϕ_m`.

# References
 - Stable (ζ > 0): Eq. 32 in Gryanik et al. (2020).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A1 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function phi(p::GryanikParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    # Gryanik stable: (1 + b_m*ζ)^(2/3). Safe for ζ > 0.
    _a_m = FT(a_m(p))
    _b_m = FT(b_m(p))
    denom = cbrt((FT(1) + _b_m * max(ζ, FT(0)))^2)
    phi_stable = FT(1) + (_a_m * max(ζ, FT(0))) / denom

    return ifelse(
        ζ > 0,
        phi_stable,
        _phi_m_unstable(min(ζ, FT(0)), FT(b_m_unstable(p))), # Businger fallback
    )
end

"""
    phi(p::GryanikParams, ζ, ::HeatTransport)

Gryanik heat/scalar-gradient similarity `ϕ_h`.

# References
 - Stable (ζ >= 0): Eq. 33 in Gryanik et al. (2020). Matches the neutral limit `ϕ_h(0) = Pr_0`.
 - Unstable (ζ < 0): Scaled Businger form to ensure continuity at ζ=0.
"""
@inline function phi(p::GryanikParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    _Pr_0 = FT(Pr_0(p))
    _a_h = FT(a_h(p))
    _b_h = FT(b_h(p))

    # Stable: Gryanik et al. (2020), Eq. 33
    phi_stable = _Pr_0 * (FT(1) + (_a_h * max(ζ, FT(0))) / (FT(1) + _b_h * max(ζ, FT(0))))
    return ifelse(
        ζ >= 0,
        phi_stable,
        _Pr_0 * _phi_h_unstable(min(ζ, FT(0)), FT(b_h_unstable(p))), # Businger fallback (scaled by Pr_0)
    )
end

"""
    psi(p::GryanikParams, ζ, ::MomentumTransport)

Gryanik momentum stability correction `ψ_m`.

# References
 - Stable (ζ > 0): Eq. 34 in Gryanik et al. (2020).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A3 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::GryanikParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    _a_m = FT(a_m(p))
    _b_m = FT(b_m(p))
    # Stable: Gryanik et al. (2020), Eq. 34
    # (1 + b_m*ζ)^(1/3) is safe for ζ > 0.
    psi_stable = -FT(3) * (_a_m / _b_m) * (cbrt(FT(1) + _b_m * max(ζ, FT(0))) - FT(1))

    return ifelse(
        ζ > 0,
        psi_stable,
        _psi_m_unstable(min(ζ, FT(0)), FT(b_m_unstable(p))), # Businger fallback
    )
end

"""
    psi(p::GryanikParams, ζ, ::HeatTransport)

Gryanik heat/scalar stability correction `ψ_h`.

# References
 - Stable (ζ > 0): Eq. 35 in Gryanik et al. (2020).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A4 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::GryanikParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    _Pr_0 = FT(Pr_0(p))
    _a_h = FT(a_h(p))
    _b_h = FT(b_h(p))

    # Stable: Gryanik et al. (2020), Eq. 35
    # log1p(b_h * ζ) (safe for ζ > 0)
    psi_stable = -_Pr_0 * (_a_h / _b_h) * log1p(_b_h * max(ζ, FT(0)))

    return ifelse(
        ζ > 0,
        psi_stable,
        _Pr_0 * _psi_h_unstable(min(ζ, FT(0)), FT(b_h_unstable(p))), # Businger fallback (scaled by Pr_0)
    )
end

"""
    Psi(p::GryanikParams, ζ, ::MomentumTransport)

Volume-averaged Gryanik momentum stability correction `Ψ_m`.

# References
 - Stable (ζ >= 0): Analytically derived from Eq. 34 in Gryanik et al. (2020).
 - Unstable (ζ < 0): Falls back to Businger form, Eq. A5 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function Psi(p::GryanikParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    _a_m = FT(a_m(p))
    _b_m = FT(b_m(p))

    # Analytical Integral of Gryanik et al. (2020), Eq. 34:
    # The volume-averaged stability function Ψ_m is defined as (1/ζ) ∫ ψ_m(x) dx.
    # Integrating Eq. 34 yields:
    # Ψ_m(ζ) = 3(a_m/b_m) - [9 a_m / (4 b_m² ζ)] * ((1 + b_m ζ)^(4/3) - 1)

    # Numerical Optimization:
    # The term ((1 + b_m ζ)^(4/3) - 1) suffers from catastrophic cancellation 
    # when ζ is small (the result approaches 0).
    # We use the identity: u^p - 1 = exp(p * ln(u)) - 1 = expm1(p * log(u))
    # substituting u = (1 + b_m ζ) and p = 4/3.
    # We further use log1p(x) for log(1+x) to maintain precision.

    ζ_safe = max(ζ, FT(0))

    # Stable large:
    base_val = _b_m * ζ_safe
    term_diff = expm1(FT(4) / FT(3) * log1p(base_val))
    numerator = FT(9) * _a_m * term_diff

    # Avoid div by zero if ζ_safe is 0 in inactive branch
    denominator = FT(4) * max(ζ_safe, eps(FT)) * _b_m^2

    val_stable_large = FT(3) * (_a_m / _b_m) - numerator / denominator

    # Stable small limit (ζ -> 0): -_a_m * ζ / 2
    val_stable_small = -_a_m * ζ_safe / FT(2)

    is_small = abs(ζ) < eps(FT)
    val_stable = ifelse(is_small, val_stable_small, val_stable_large)

    return ifelse(
        ζ >= 0,
        val_stable,
        _Psi_m_unstable(min(ζ, FT(0)), FT(b_m_unstable(p))),
    )
end

"""
    Psi(p::GryanikParams, ζ, ::HeatTransport)

Volume-averaged Gryanik heat/scalar stability correction `Ψ_h`.

# References
 - Stable (ζ >= 0): Analytically derived from Eq. 35 in Gryanik et al. (2020).
 - Unstable (ζ < 0): Falls back to Businger form, Eq. A6 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function Psi(p::GryanikParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    _a_h = FT(a_h(p))
    _b_h = FT(b_h(p))
    _Pr_0 = FT(Pr_0(p))

    ζ_safe = max(ζ, FT(0))

    # Analytical Integral of Gryanik et al. (2020), Eq. 35:
    # The volume-averaged stability function Ψ_h is defined as (1/ζ) ∫ ψ_h(x) dx.
    # Eq. 35 gives: ψ_h(ζ) = -Pr_0 * (a_h/b_h) * ln(1 + b_h * ζ).
    #
    # Integrating yields:
    # Ψ_h(ζ) = -[Pr_0 * a_h / (b_h^2 * ζ)] * ((1 + b_h * ζ) * ln(1 + b_h * ζ) - b_h * ζ)
    #
    # The implementation factors out (1/b_h) inside the bracket for structure:
    # Ψ_h(ζ) = -[Pr_0 * a_h / (b_h * ζ)] * ((1/b_h + ζ) * ln(1 + b_h * ζ) - ζ)

    # Factor 1/ζ creates singularity at 0.
    term_paren = (FT(1) / _b_h + ζ_safe) * log1p(_b_h * ζ_safe) - ζ_safe
    val_stable_large = -_a_h / _b_h / max(ζ_safe, eps(FT)) * _Pr_0 * term_paren

    # Stable small limit (linear approximation for very small ζ)
    val_stable_small = -_a_h * ζ_safe * _Pr_0 / FT(2)

    is_small = abs(ζ) < eps(FT)
    val_stable = ifelse(is_small, val_stable_small, val_stable_large)

    return ifelse(
        ζ >= 0,
        val_stable,
        _Pr_0 * _Psi_h_unstable(min(ζ, FT(0)), FT(b_h_unstable(p))),
    )
end

#####
##### Grachev
#####

"""
    GrachevParams{FT}

Parameter bundle for the Grachev et al. (2007) similarity relations,
based on SHEBA data.

 - `a_h`, `b_h`, `c_h`: Coefficients for heat/scalar stability function (Eq. 9b).
   Note: `c_h` is the coefficient for the linear ζ term in the denominator.
 - `Pr_0`: Neutral Prandtl number. In the original Grachev et al. (2007) derivation,
    this is 1.0. It is included here for structural consistency with other
    parameterizations.

!!! note "LayerAverageScheme not supported"
    The volume-averaged stability corrections `Ψ` are not implemented for `GrachevParams`
    because the analytical integrals of Eqs. 12–13 in Grachev et al. (2007) do not admit
    closed-form expressions suitable for efficient evaluation. Use `PointValueScheme`, or
    use `BusingerParams`/`GryanikParams` if layer averaging is required.

Reference: Grachev et al. (2007).
"""
Base.@kwdef struct GrachevParams{FT} <: AbstractUniversalFunctionParameters{FT}
    Pr_0::FT
    a_m::FT
    a_h::FT
    b_m::FT
    b_h::FT
    c_h::FT
    b_m_unstable::FT
    b_h_unstable::FT
    ζ_a::FT
    γ::FT
end

# Accessor methods for unstable coefficients (defined after structs)
b_m_unstable(p::GryanikParams) = p.b_m_unstable
b_h_unstable(p::GryanikParams) = p.b_h_unstable
b_m_unstable(p::GrachevParams) = p.b_m_unstable
b_h_unstable(p::GrachevParams) = p.b_h_unstable

"""
    phi(p::GrachevParams, ζ, ::MomentumTransport)

Grachev momentum similarity `ϕ_m`.

# References
 - Stable (ζ > 0): Eq. 9a in Grachev et al. (2007).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A1 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function phi(p::GrachevParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    _a_m = FT(a_m(p))
    _b_m = FT(b_m(p))

    # Stable: Grachev et al. (2007) Eq. 9a:
    # ϕ_m(ζ) = 1 + (a_m * ζ * (1 + ζ)^(1/3)) / (1 + b_m * ζ)
    #
    # We use cbrt(1 + ζ) instead of ^(1/3) for better performance on GPUs.
    ζ_safe = max(ζ, FT(0))
    phi_stable = FT(1) + _a_m * ζ_safe * cbrt(FT(1) + ζ_safe) / (FT(1) + _b_m * ζ_safe)

    return ifelse(
        ζ > 0,
        phi_stable,
        _phi_m_unstable(min(ζ, FT(0)), FT(b_m_unstable(p))),  # Businger fallback
    )
end

"""
    phi(p::GrachevParams, ζ, ::HeatTransport)

Grachev heat/scalar-gradient similarity `ϕ_h`.

Calculates the gradient function scaled by `Pr_0`. While Grachev et al. (2007)
assume `ϕ_h(0) = 1`, this implementation returns `ϕ_h(ζ) = Pr_0 * (...)` 
to maintain the invariant `ϕ_h(0) = Pr_0` used by the solver.

# References
 - Stable (ζ > 0): Eq. 9b in Grachev et al. (2007).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A2 (L < 0) in Nishizawa & Kitamura (2018).
"""

@inline function phi(p::GrachevParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    # Stable: Grachev et al. (2007) Eq. 9b:
    # ϕ_h(ζ) = 1 + (a_h * ζ + b_h * ζ^2) / (1 + c_h * ζ + ζ^2)
    # Scaled by Pr_0 (which is 1.0 for Grachev) for consistency.
    _a_h = FT(a_h(p))
    _b_h = FT(b_h(p))
    _c_h = FT(c_h(p))
    _Pr_0 = FT(Pr_0(p))
    ζ_safe = max(ζ, FT(0))
    phi_stable =
        _Pr_0 *
        (FT(1) + (_a_h * ζ_safe + _b_h * ζ_safe^2) / (FT(1) + _c_h * ζ_safe + ζ_safe^2))

    return ifelse(
        ζ > 0,
        phi_stable,
        _Pr_0 * _phi_h_unstable(min(ζ, FT(0)), FT(b_h_unstable(p))),  # Businger fallback (scaled by Pr_0)
    )
end

"""
    psi(p::GrachevParams, ζ, ::MomentumTransport)

Grachev momentum stability correction `ψ_m`.

# References
 - Stable (ζ > 0): Eq. 12 in Grachev et al. (2007).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A3 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::GrachevParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    _a_m = FT(a_m(p))
    _b_m = FT(b_m(p))

    # Grachev et al. (2007) Eq. 12:
    # ψ_m(ζ) = -3(a_m/b_m)(x - 1) + (a_m B_m / (2 b_m)) * [ ln_terms + atan_terms ]
    #
    # Definitions: x = (1 + ζ)^(1/3), B_m = ((1 - b_m)/b_m)^(1/3)

    ζ_safe = max(ζ, FT(0))

    # Auxiliary variables
    x = cbrt(FT(1) + ζ_safe)
    B_m = cbrt(FT(1) / _b_m - FT(1))
    sqrt3 = sqrt(FT(3))

    # 1. Linear term: -3(a_m/b_m)(x - 1)
    linear_term = -FT(3) * (_a_m / _b_m) * (x - FT(1))

    # 2. Logarithmic terms inside the bracket:
    # 2 * ln((x + B_m) / (1 + B_m)) - ln((x^2 - x*B_m + B_m^2) / (1 - B_m + B_m^2))
    log_1 = FT(2) * log((x + B_m) / (FT(1) + B_m))
    log_2 = -log((x^2 - x * B_m + B_m^2) / (FT(1) - B_m + B_m^2))

    # 3. Arctangent terms inside the bracket:
    # 2√3 * (atan((2x - B_m)/(√3 B_m)) - atan((2 - B_m)/(√3 B_m)))
    atan_1 = atan((FT(2) * x - B_m) / (sqrt3 * B_m))
    atan_2 = atan((FT(2) - B_m) / (sqrt3 * B_m))
    atan_part = FT(2) * sqrt3 * (atan_1 - atan_2)

    # Combine all parts
    bracket_sum = log_1 + log_2 + atan_part
    psi_stable = linear_term + (_a_m * B_m) / (FT(2) * _b_m) * bracket_sum

    return ifelse(
        ζ > 0,
        psi_stable,
        _psi_m_unstable(min(ζ, FT(0)), FT(b_m_unstable(p))),  # Businger fallback
    )
end

"""
    psi(p::GrachevParams, ζ, ::HeatTransport)

Grachev heat/scalar stability correction `ψ_h`.

# References
 - Stable (ζ > 0): Eq. 13 in Grachev et al. (2007).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A4 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::GrachevParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    # Stable: Grachev et al. (2007) Eq. 13:
    # ψ_h(ζ) = - (b_h/2) * ln(1 + c_h*ζ + ζ^2) 
    #          + [ -a_h/B_h + (b_h*c_h)/(2*B_h) ] * [ ln(...) - ln(...) ]
    #
    # where B_h = sqrt(c_h^2 - 4)

    _a_h = FT(a_h(p))
    _b_h = FT(b_h(p))
    _c_h = FT(c_h(p))
    ζ_safe = max(ζ, FT(0))

    # 1. Auxiliary constant B_h
    B_h = sqrt(_c_h^2 - FT(4))

    # 2. Coefficient for the fractional log terms
    # Corresponds to: [ -a_h/B_h + (b_h*c_h)/(2*B_h) ]
    # We compute the positive magnitude and negate it in the final return.
    # coeff = a_h/B_h - (b_h*c_h)/(2*B_h)
    coeff = _a_h / B_h - (_b_h * _c_h) / (FT(2) * B_h)

    # 3. Fractional logarithmic terms
    # ln((2ζ + c_h - B_h) / (2ζ + c_h + B_h)) - ln((c_h - B_h) / (c_h + B_h))
    # The second log term ensures ψ_h(0) = 0.
    log_num = FT(2) * ζ_safe + _c_h - B_h
    log_den = FT(2) * ζ_safe + _c_h + B_h
    log_0_num = _c_h - B_h
    log_0_den = _c_h + B_h

    fractional_logs = log(log_num / log_den) - log(log_0_num / log_0_den)

    # 4. Quadratic Logarithmic Term: (b_h/2) * ln(1 + c_h*ζ + ζ^2)
    quadratic_log = (_b_h / FT(2)) * log1p(_c_h * ζ_safe + ζ_safe^2)

    # Final Sum: -coeff * fractional_logs - quadratic_log
    psi_stable = -coeff * fractional_logs - quadratic_log
    _Pr_0 = FT(Pr_0(p))

    return ifelse(
        ζ > 0,
        _Pr_0 * psi_stable,
        _Pr_0 * _psi_h_unstable(min(ζ, FT(0)), FT(b_h_unstable(p))),
    )
end

"""
    bulk_richardson_number(uf_params, Δz_eff, ζ, z0m, z0h)

Compute the bulk Richardson number at a given stability parameter ζ,
defined as:

    Ri_b(ζ) = ζ * F_h(ζ) / F_m(ζ)^2

where F_m and F_h are the dimensionless profiles for momentum and heat/scalars.
"""
function bulk_richardson_number(uf_params, Δz_eff, ζ, z0m, z0h, scheme)
    F_m = dimensionless_profile(uf_params, Δz_eff, ζ, z0m, MomentumTransport(), scheme)
    F_h = dimensionless_profile(uf_params, Δz_eff, ζ, z0h, HeatTransport(), scheme)
    return ζ * F_h / F_m^2
end

# Default to PointValueScheme
function bulk_richardson_number(uf_params, Δz_eff, ζ, z0m, z0h)
    return bulk_richardson_number(
        uf_params,
        Δz_eff,
        ζ,
        z0m,
        z0h,
        PointValueScheme(),
    )
end


"""
    dimensionless_profile(uf_params, Δz_eff, ζ, z0, transport, scheme)

The dimensionless vertical profile of the variable (momentum or scalar).
"""
function dimensionless_profile end

"""
    dimensionless_profile(uf_params, Δz_eff, ζ, z0, transport, ::PointValueScheme)

The dimensionless vertical profile of the variable (momentum or scalar) using
point values (standard Monin-Obukhov Similarity Theory).

Defined as

    F(Δz_eff) = ϕ(0) * ln(Δz_eff/z0) - ψ(ζ) + ψ(ζ * z0/Δz_eff),

This represents the integral of the dimensionless gradient function ϕ(ζ)/z
from the roughness length z0 to the given height `Δz_eff`. 

The neutral limit `ϕ(0)` determines the slope of the logarithmic profile:
- For `MomentumTransport`: `ϕ_m(0) = 1`.
- For `HeatTransport`: `ϕ_h(0) = Pr_0` (Neutral Prandtl number).

The resulting profile is:
    F(z) = ϕ(0) * ln(z/z0) - ψ(ζ) + ψ(ζ_0)

Note that ϕ(0) corresponds
to the neutral dimensionless gradient (slope), which depends on the parameterization 
(e.g., `Pr_0` for Businger/Gryanik, 1 for Grachev) and transport type.
"""
@inline function dimensionless_profile(
    uf_params,
    Δz_eff,
    ζ,
    z0,
    transport,
    ::PointValueScheme,
)
    FT = eltype(ζ)
    slope = phi(uf_params, FT(0), transport)
    return slope * log(Δz_eff / z0) - psi(uf_params, ζ, transport) +
           psi(uf_params, z0 * ζ / Δz_eff, transport)
end

# Default to PointValueScheme
@inline function dimensionless_profile(uf_params, Δz_eff, ζ, z0, transport)
    return dimensionless_profile(
        uf_params,
        Δz_eff,
        ζ,
        z0,
        transport,
        PointValueScheme(),
    )
end

"""
    dimensionless_profile(uf_params, Δz_eff, ζ, z0, transport, ::LayerAverageScheme)

The dimensionless vertical profile of the variable (momentum or scalar) using
layer-averaged values (finite volume formulation).

Derivation follows Nishizawa & Kitamura (2018), adapted for generalized neutral limits.

    F_ave(Δz_eff) = Slope * (ln(Δz_eff/z0) - R_z0) - Ψ(ζ) + (z0/Δz_eff) * Ψ(ζ * z0/Δz_eff) + R_z0 * ψ(ζ * z0/Δz_eff)

where:
 - Slope = ϕ(0). This is `1` for momentum and `Pr_0` for heat/scalars.
 - R_z0 = 1 - z0/Δz_eff (Approximation of geometric factor)
"""
@inline function dimensionless_profile(
    uf_params,
    Δz_eff,
    ζ,
    z0,
    transport,
    ::LayerAverageScheme,
)
    FT = eltype(ζ)
    slope = phi(uf_params, FT(0), transport)
    ζ_z0 = ζ * z0 / Δz_eff
    R_z0 = FT(1) - z0 / Δz_eff

    # Note: Nishizawa & Kitamura (2018) Eq 24:
    # F = ln(Δz/z0) - Ψ_M(ζ) + (z0/Δz)Ψ_M(ζ_z0) + ε_M(ζ_z0)
    # where ε_M(ζ) = (1 - z0/Δz)(ψ_M(ζ) - 1)
    #
    # Generalized for slope != 1:
    # The logarithmic term ln(z/z0) becomes Slope * ln(z/z0).
    # The linear term -1 in (ψ - 1) comes from the integral of 1/z -> ln(z).
    # If slope is phi(0), then the integral is slope * ln(z).
    # So the term (ψ - 1) becomes (ψ - slope).

    return slope * log(Δz_eff / z0) - Psi(uf_params, ζ, transport) +
           (z0 / Δz_eff) * Psi(uf_params, ζ_z0, transport) +
           R_z0 * (psi(uf_params, ζ_z0, transport) - slope)
end

end # module
