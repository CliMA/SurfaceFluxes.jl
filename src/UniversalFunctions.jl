"""
    UniversalFunctions

Universal stability and stability correction functions for `SurfaceFluxes` module. 
Supports the following universal functions:
 - `Businger`: Businger et al. (1971), Dyer (1974)
 - `Gryanik`: Gryanik et al. (2020)
 - `Grachev`: Grachev et al. (2007)
"""
module UniversalFunctions

abstract type AbstractUniversalFunctionParameters{FT <: Real} end
const AUFP = AbstractUniversalFunctionParameters

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

Type selecting heat-transfer stability functions (ϕₕ, ψₕ, Ψₕ).
"""
struct HeatTransport <: AbstractTransportType end

Base.broadcastable(tt::AbstractTransportType) = tuple(tt)
Base.broadcastable(p::AbstractUniversalFunctionParameters) = tuple(p)

"""
    phi

Universal stability function for wind shear (`ϕ_m`) and 
temperature gradient (`ϕ_h`)
"""
function phi end

"""
    psi(p, ζ, transport_type)

The standard integrated stability correction function `ψ(ζ)`.
Defined as:
    ψ(ζ) = ∫[0 to ζ] (ϕ(0) - ϕ(x)) / x dx

This is the standard correction used in point-based Monin-Obukhov
similarity theory. Note that while ϕ(0) is typically 1, it may differ 
(e.g. ϕ_h(0) = Pr_0 in Gryanik et al., 2020).
"""
function psi end

"""
    Psi(p, ζ, transport_type)

The volume-averaged stability correction function `Ψ(ζ)`.
Mathematically, this is defined as:
    Ψ(ζ) = (1/ζ) ∫[0 to ζ] ψ(x) dx

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
c_m(p::AUFP) = p.c_m
d_h(p::AUFP) = p.d_h
d_m(p::AUFP) = p.d_m
ζ_a(p::AUFP) = p.ζ_a
γ(p::AUFP) = p.γ

# Parameter-struct π-group (avoid constructing UF just to get scalar π)
π_group(p::AUFP, ::HeatTransport) = Pr_0(p)
π_group(::AUFP, ::MomentumTransport) = 1

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
    if abs(ζ) < eps(FT)
        return -γ * ζ / FT(8)
    end

    # Full computation (Nishizawa & Kitamura 2018, Eq. A5)
    # Ψ_m = log[(1+x)²(1+x²)/8] - 2tan⁻¹(x) + π/2 - 1 + (1-x³)/(3γζ/4)

    # 1. Compute x = (1 - γζ)^(1/4)
    x = sqrt(sqrt(FT(1) - γ * ζ))

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
    #                  = -expm1(0.75 * log1p(-γζ))
    cubic_num = -expm1(FT(0.75) * log1p(-γ * ζ))

    # Denominator Correction: 
    # The paper (Eq. A5) hardcodes the denominator as 12ζ (assuming γ=16).
    # The general form for arbitrary γ is 3γζ/4.
    cubic_denom = (FT(3) * γ / FT(4)) * ζ
    cubic_term = cubic_num / cubic_denom

    # Final Sum: log_term - tan_term + π_term - 1 + cubic_term
    return log_term - tan_term + π_term - FT(1) + cubic_term
end

"""
    _Psi_h_unstable(ζ, γ)

Finite-volume integral for unstable heat.
Matches Nishizawa & Kitamura (2018, Eq. A6).

Note: `γ` here refers to the **coefficient inside the sqrt**
(e.g., 9 in `(1 - 9ζ)^(-1/2)`), not the linear coefficient `a_h`.
"""
@inline function _Psi_h_unstable(ζ, γ)
    FT = eltype(ζ)
    # Small-ζ limit (Nishizawa & Kitamura 2018, Eq. A14)
    if abs(ζ) < eps(FT)
        return -γ * ζ / FT(4)
    end

    # Full computation (Nishizawa & Kitamura 2018, Eq. A6)
    # Using expm1 for precision near ζ=0
    # y = (1 - γζ)^(1/2)
    # 1 - y = 1 - (1 - γζ)^(1/2) = -expm1(0.5 * log1p(-γζ))

    y = sqrt(FT(1) - γ * ζ)
    log_term = FT(2) * log((FT(1) + y) / FT(2))

    # Optimized linear term
    lin_num = -expm1(FT(0.5) * log1p(-γ * ζ))
    lin_term = FT(2) * lin_num / (γ * ζ)

    return log_term + lin_term - FT(1)
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
    if ζ < 0
        return _phi_m_unstable(ζ, b_m(p))
    else
        return FT(a_m(p)) * ζ + FT(1)
    end
end

"""
    phi(p::BusingerParams, ζ, ::HeatTransport)

Businger heat-gradient similarity `ϕ_h`.

# References
 - Stable (ζ >= 0): Eq. A2 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A2 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function phi(p::BusingerParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if ζ < 0
        return _phi_h_unstable(ζ, b_h(p))
    else
        _a_h = FT(a_h(p))
        _π_group = FT(π_group(p, tt))
        return _a_h * ζ / _π_group + FT(1)
    end
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
    if ζ < 0
        return _psi_m_unstable(ζ, b_m(p))
    else
        return -FT(a_m(p)) * ζ
    end
end

"""
    psi(p::BusingerParams, ζ, ::HeatTransport)

Businger heat stability correction `ψ_h`.

# References
 - Stable (ζ >= 0): Eq. A4 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A4 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::BusingerParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if ζ < 0
        return _psi_h_unstable(ζ, b_h(p))
    else
        _a_h = FT(a_h(p))
        _π_group = FT(π_group(p, tt))
        return -_a_h * ζ / _π_group
    end
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
    if ζ >= 0
        return -FT(a_m(p)) * ζ / FT(2)
    else
        return _Psi_m_unstable(ζ, b_m(p))
    end
end

"""
    Psi(p::BusingerParams, ζ, ::HeatTransport)

Volume-averaged Businger heat stability correction `Ψ_h`.

# References
 - Stable (ζ >= 0): Eqs. A6 and A14 (L >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A6 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function Psi(p::BusingerParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if ζ >= 0
        _a_h = FT(a_h(p))
        _π_group = FT(π_group(p, tt))
        return -_a_h * ζ / (FT(2) * _π_group)
    else
        return _Psi_h_unstable(ζ, b_h(p))
    end
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
    if ζ > 0
        # Gryanik et al. (2020), Eq. 32
        _a_m = FT(a_m(p))
        _b_m = FT(b_m(p))
        # Optimization: (1+b_m*ζ)^(2/3) -> cbrt((1+b_m*ζ)^2)
        denom = cbrt((FT(1) + _b_m * ζ)^2)
        return FT(1) + (_a_m * ζ) / denom
    else
        # Fallback to Businger form
        return _phi_m_unstable(ζ, FT(b_m_unstable(p)))
    end
end

"""
    phi(p::GryanikParams, ζ, ::HeatTransport)

Gryanik heat-gradient similarity `ϕ_h`.

# References
 - Stable (ζ >= 0): Eq. 33 in Gryanik et al. (2020). Matches the neutral limit `ϕ_h(0) = Pr_0`.
 - Unstable (ζ < 0): Scaled Businger form to ensure continuity at ζ=0.
"""
@inline function phi(p::GryanikParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    _Pr_0 = FT(Pr_0(p))
    if ζ >= 0
        # Gryanik et al. (2020), Eq. 33
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        return _Pr_0 * (FT(1) + (_a_h * ζ) / (FT(1) + _b_h * ζ))
    else
        # Fallback to Businger form but scale unstable branch by Pr_0 to ensure continuity at 0
        return _Pr_0 * _phi_h_unstable(ζ, FT(b_h_unstable(p)))
    end
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
    if ζ > 0
        # Gryanik et al. (2020), Eq. 34
        _a_m = FT(a_m(p))
        _b_m = FT(b_m(p))
        # Optimization: (1 + b_m*ζ)^(1/3) -> cbrt(...)
        return -FT(3) * (_a_m / _b_m) * (cbrt(FT(1) + _b_m * ζ) - FT(1))
    else
        # Fallback to Businger form
        return _psi_m_unstable(ζ, FT(b_m_unstable(p)))
    end
end

"""
    psi(p::GryanikParams, ζ, ::HeatTransport)

Gryanik heat stability correction `ψ_h`.

# References
 - Stable (ζ > 0): Eq. 35 in Gryanik et al. (2020).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A4 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::GryanikParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    _Pr_0 = FT(Pr_0(p))
    if ζ > 0
        # Gryanik et al. (2020), Eq. 35
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        return -_Pr_0 * (_a_h / _b_h) * log1p(_b_h * ζ)
    else
        # Fallback to Businger form but scale unstable branch by Pr_0 to ensure continuity at 0
        return _Pr_0 * _psi_h_unstable(ζ, FT(b_h_unstable(p)))
    end
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
    if ζ >= 0
        # Limit at ζ -> 0 is 0. Required to avoid NaN (0/0) at exactly ζ = 0.
        if abs(ζ) < eps(FT)
            return -_a_m * ζ / FT(2)
        end

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

        base_val = _b_m * ζ
        term_diff = expm1(FT(4) / FT(3) * log1p(base_val)) # Equivalent to ((1 + bζ)^(4/3) - 1)

        numerator = FT(9) * _a_m * term_diff
        denominator = FT(4) * ζ * _b_m^2

        return FT(3) * (_a_m / _b_m) - numerator / denominator
    else
        # Fallback to Businger form
        return _Psi_m_unstable(ζ, FT(b_m_unstable(p)))
    end
end

"""
    Psi(p::GryanikParams, ζ, ::HeatTransport)

Volume-averaged Gryanik heat stability correction `Ψ_h`.

# References
 - Stable (ζ >= 0): Analytically derived from Eq. 35 in Gryanik et al. (2020).
 - Unstable (ζ < 0): Falls back to Businger form, Eq. A6 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function Psi(p::GryanikParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    _a_h = FT(a_h(p))
    _b_h = FT(b_h(p))
    _Pr_0 = FT(Pr_0(p))
    if ζ >= 0
        # Limit at ζ -> 0 is 0
        if abs(ζ) < eps(FT)
            return -_a_h * ζ * _Pr_0 / FT(2) # Linear approximation for very small ζ
        end

        # Analytical Integral of Gryanik et al. (2020), Eq. 35:
        # The volume-averaged stability function Ψ_h is defined as (1/ζ) ∫ ψ_h(x) dx.
        # Eq. 35 gives: ψ_h(ζ) = -Pr_0 * (a_h/b_h) * ln(1 + b_h * ζ).
        #
        # Integrating yields:
        # Ψ_h(ζ) = -[Pr_0 * a_h / (b_h^2 * ζ)] * ((1 + b_h * ζ) * ln(1 + b_h * ζ) - b_h * ζ)
        #
        # The implementation factors out (1/b_h) inside the bracket for structure:
        # Ψ_h(ζ) = -[Pr_0 * a_h / (b_h * ζ)] * ((1/b_h + ζ) * ln(1 + b_h * ζ) - ζ)
        return -_a_h / _b_h / ζ * _Pr_0 * ((FT(1) / _b_h + ζ) * log1p(_b_h * ζ) - ζ)
    else
        # Fallback to Businger form but scale unstable branch by Pr_0 to ensure continuity at 0
        return _Pr_0 * _Psi_h_unstable(ζ, FT(b_h_unstable(p)))
    end
end

#####
##### Grachev
#####

"""
    GrachevParams{FT}

Parameter bundle for the Grachev et al. (2007) similarity relations,
based on SHEBA data.

 - `a_m`, `b_m`: Coefficients for momentum stability function (Eq. 9a).
 - `a_h`, `b_h`, `c_h`: Coefficients for heat stability function (Eq. 9b).
   Note: `c_h` is the coefficient for the linear ζ term in the denominator.

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
    if ζ > 0
        # Grachev et al. (2007) Eq. 9a:
        # ϕ_m(ζ) = 1 + (a_m * ζ * (1 + ζ)^(1/3)) / (1 + b_m * ζ)
        #
        # We use cbrt(1 + ζ) instead of ^(1/3) for better performance on GPUs.
        _a_m = FT(a_m(p))
        _b_m = FT(b_m(p))
        return FT(1) + _a_m * ζ * cbrt(FT(1) + ζ) / (FT(1) + _b_m * ζ)
    else
        # Fallback to Businger form
        return _phi_m_unstable(ζ, FT(b_m_unstable(p)))
    end
end

"""
    phi(p::GrachevParams, ζ, ::HeatTransport)

Grachev heat-gradient similarity `ϕ_h`.

# References
 - Stable (ζ > 0): Eq. 9b in Grachev et al. (2007).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A2 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function phi(p::GrachevParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    if ζ > 0
        # Grachev et al. (2007) Eq. 9b:
        # ϕ_h(ζ) = 1 + (a_h * ζ + b_h * ζ^2) / (1 + c_h * ζ + ζ^2)
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        _c_h = FT(c_h(p))
        return FT(1) + (_a_h * ζ + _b_h * ζ^2) / (FT(1) + _c_h * ζ + ζ^2)
    else
        # Fallback to Businger form
        return _phi_h_unstable(ζ, FT(b_h_unstable(p)))
    end
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
    if ζ > 0
        # Grachev et al. (2007) Eq. 12:
        # ψ_m(ζ) = -3(a_m/b_m)(x - 1) + (a_m B_m / (2 b_m)) * [ ln_terms + atan_terms ]
        #
        # Definitions: x = (1 + ζ)^(1/3), B_m = ((1 - b_m)/b_m)^(1/3)

        _a_m = FT(a_m(p))
        _b_m = FT(b_m(p))

        # Auxiliary variables
        x = cbrt(FT(1) + ζ)
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
        return linear_term + (_a_m * B_m) / (FT(2) * _b_m) * bracket_sum
    else
        # Fallback to Businger form
        return _psi_m_unstable(ζ, FT(b_m_unstable(p)))
    end
end

"""
    psi(p::GrachevParams, ζ, ::HeatTransport)

Grachev heat stability correction `ψ_h`.

# References
 - Stable (ζ > 0): Eq. 13 in Grachev et al. (2007).
 - Unstable (ζ <= 0): Falls back to Businger form, Eq. A4 (L < 0) in Nishizawa & Kitamura (2018).
"""
@inline function psi(p::GrachevParams, ζ, ::HeatTransport)
    FT = eltype(ζ)
    if ζ > 0
        # Grachev et al. (2007) Eq. 13:
        # ψ_h(ζ) = - (b_h/2) * ln(1 + c_h*ζ + ζ^2) 
        #          + [ -a_h/B_h + (b_h*c_h)/(2*B_h) ] * [ ln(...) - ln(...) ]
        #
        # where B_h = sqrt(c_h^2 - 4)

        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        _c_h = FT(c_h(p))

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
        log_num = FT(2) * ζ + _c_h - B_h
        log_den = FT(2) * ζ + _c_h + B_h
        log_0_num = _c_h - B_h
        log_0_den = _c_h + B_h

        fractional_logs = log(log_num / log_den) - log(log_0_num / log_0_den)

        # 4. Quadratic Logarithmic Term: (b_h/2) * ln(1 + c_h*ζ + ζ^2)
        quadratic_log = (_b_h / FT(2)) * log1p(_c_h * ζ + ζ^2)

        # Final Sum: -coeff * fractional_logs - quadratic_log
        return -coeff * fractional_logs - quadratic_log
    else
        # Fallback to Businger form
        return _psi_h_unstable(ζ, FT(b_h_unstable(p)))
    end
end

"""
    richardson_number(uf_params, ζ)

Compute the Gradient Richardson number at a given stability parameter ζ.
Defined as:
    Ri(ζ) = ζ * ϕ_h(ζ) / ϕ_m(ζ)^2
"""
function richardson_number(uf_params, ζ)
    ϕ_m = phi(uf_params, ζ, MomentumTransport())
    ϕ_h = phi(uf_params, ζ, HeatTransport())
    return ζ * ϕ_h / ϕ_m^2
end

"""
    dimensionless_profile(uf_params, Δz, ζ, z0, transport)

The dimensionless vertical profile of the variable (momentum or scalar).
Defined as:

    F(z) = ln(z/z0) - ψ(ζ) + ψ(ζ * z0/z)

This represents the integral of the dimensionless gradient function ϕ(ζ)/z
from roughness length z0 to the given height z.
"""
@inline function dimensionless_profile(uf_params, Δz, ζ, z0, transport)
    return log(Δz / z0) -
           psi(uf_params, ζ, transport) +
           psi(uf_params, z0 * ζ / Δz, transport)
end

end # module
