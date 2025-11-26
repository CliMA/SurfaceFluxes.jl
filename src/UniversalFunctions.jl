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
(e.g. ϕ_h(0) [cite_start]= Pr_0 in Gryanik et al., 2020)[cite: 1418].
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
# and ensure consistency across Businger, Gryanik, and Grachev.

@inline _phi_m_unstable(ζ, γ) = 1 / sqrt(sqrt(1 - γ * ζ))

@inline _phi_h_unstable(ζ, γ) = 1 / sqrt(1 - γ * ζ)

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

Finite-volume integral for unstable momentum.
Matches Nishizawa & Kitamura (2018) Eq. A5.

Note: `γ` here refers to the **coefficient inside the sqrt/cbrt** (e.g., 15 in `(1 - 15ζ)^(-1/4)`), not the linear coefficient `a_m`.
"""
function _Psi_m_unstable(ζ, γ)
    FT = eltype(ζ)
    # Small-ζ limit (Nishizawa2018 Eq. A13)
    limit_val = -γ * ζ / FT(8)
    # Full computation (Nishizawa2018 Eq. A5) with safe division to avoid zero
    ζ_safe = ifelse(abs(ζ) < eps(FT), copysign(eps(FT), ζ), ζ)
    x = sqrt(sqrt(FT(1) - γ * ζ))
    log_term = log((FT(1) + x)^2 * (FT(1) + x^2) / FT(8))
    π_term = FT(π) / FT(2)
    tan_term = FT(2) * atan(x)
    cubic_term = (FT(1) - x^3) / (FT(12) * ζ_safe)
    full_val = log_term - tan_term + π_term - FT(1) + cubic_term
    return ifelse(abs(ζ) < eps(FT), limit_val, full_val)
end

"""
    _Psi_h_unstable(ζ, γ)

Finite-volume integral for unstable heat.
Matches Nishizawa & Kitamura (2018) Eq. A6.

Note: `γ` here refers to the **coefficient inside the sqrt**
(e.g., 9 in `(1 - 9ζ)^(-1/2)`), not the linear coefficient `a_h`.
"""
function _Psi_h_unstable(ζ, γ)
    FT = eltype(ζ)
    # Small-ζ limit (Nishizawa2018 Eq. A14)
    limit_val = -γ * ζ / FT(4)
    # Full computation (Nishizawa2018 Eq. A6) with safe division to avoid zero
    ζ_safe = ifelse(abs(ζ) < eps(FT), copysign(eps(FT), ζ), ζ)
    y = sqrt(FT(1) - γ * ζ)
    log_term = FT(2) * log((FT(1) + y) / FT(2))
    full_val = log_term + FT(2) * (FT(1) - y) / (γ * ζ_safe) - FT(1)
    return ifelse(abs(ζ) < eps(FT), limit_val, full_val)
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
 - Stable (ζ >= 0): Derived from linear form in Eq. A13 (ζ >= 0) in Nishizawa & Kitamura (2018).
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
 - Stable (ζ >= 0): Derived from linear form in Eq. A14 (ζ >= 0) in Nishizawa & Kitamura (2018).
 - Unstable (ζ < 0): Eq. A6 (L < 0) in Nishizawa & Kitamura (2018).
 - Small ζ limit: Eq. A14 (ζ < 0) in Nishizawa & Kitamura (2018).
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
 - `Pr_0`: Neutral Prandtl number. The paper recommends Pr_0 ≈ 0.98.

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
        # Fallback to Businger form with γ = 15
        return _phi_m_unstable(ζ, FT(15))
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
        # Gryanik2020 Eq. 33
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        return _Pr_0 * (FT(1) + (_a_h * ζ) / (FT(1) + _b_h * ζ))
    else
        # Fallback to Businger form but scale unstable branch by Pr_0 to ensure continuity at 0
        return _Pr_0 * _phi_h_unstable(ζ, FT(9))
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
        # Fallback to Businger form with γ = 15
        return _psi_m_unstable(ζ, FT(15))
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
        # Gryanik2020 Eq. 35
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        return -_Pr_0 * (_a_h / _b_h) * log1p(_b_h * ζ)
    else
        # Fallback to Businger form but scale unstable branch by Pr_0 to ensure continuity at 0
        return _Pr_0 * _psi_h_unstable(ζ, FT(9))
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

        # Optimization: use expm1/log1p to preserve precision for small ζ
        # (1 + b_m * ζ)^(4/3) - 1  == expm1(4/3 * log1p(b_m * ζ))
        term_diff = expm1(FT(4) / FT(3) * log1p(_b_m * ζ))
        numerator = FT(9) * _a_m * term_diff
        denominator = FT(4) * ζ * _b_m^2

        return FT(3) * (_a_m / _b_m) - numerator / denominator
    else
        # Fallback to Businger form with γ = 15
        return _Psi_m_unstable(ζ, FT(15))
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

        # Gryanik2020 Eq. 35 integral
        return -_a_h / _b_h / ζ * _Pr_0 * ((FT(1) / _b_h + ζ) * log1p(_b_h * ζ) - ζ)
    else
        # Fallback to Businger form but scale unstable branch by Pr_0 to ensure continuity at 0
        return _Pr_0 * _Psi_h_unstable(ζ, FT(9))
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
        # Optimization: (1+ζ)^(1/3) -> cbrt(1+ζ)
        return FT(1) + _a_m * ζ * cbrt(FT(1) + ζ) / (FT(1) + _b_m * ζ)
    else
        # Fallback to Businger form with γ = 15
        return _phi_m_unstable(ζ, FT(15))
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
        # Fallback to Businger form with γ = 9
        return _phi_h_unstable(ζ, FT(9))
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
        # Optimization: cbrt
        B_m = cbrt(FT(1) / _b_m - FT(1))
        x = cbrt(FT(1) + ζ)
        sqrt3 = sqrt(FT(3))
        linear_term = -FT(3) * (_a_m / _b_m) * (x - FT(1))
        log_term_1 = FT(2) * log((x + B_m) / (FT(1) + B_m))
        log_term_2 = log((x^2 - x * B_m + B_m^2) / (FT(1) - B_m + B_m^2))
        atan_term_1 = atan((FT(2) * x - B_m) / (sqrt3 * B_m))
        atan_term_2 = atan((FT(2) - B_m) / (sqrt3 * B_m))
        atan_terms = atan_term_1 - atan_term_2
        bracket_term = log_term_1 - log_term_2 + FT(2) * sqrt3 * atan_terms
        return linear_term + _a_m * B_m / (FT(2) * _b_m) * bracket_term
    else
        # Fallback to Businger form with γ = 15
        return _psi_m_unstable(ζ, FT(15))
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
        # Grachev2007 Eq. 13
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        _c_h = FT(c_h(p))
        # Constants computed on the fly 
        B_h = sqrt(_c_h^2 - FT(4))
        coeff = _a_h / B_h - _b_h * _c_h / (FT(2) * B_h)
        log_term_1 = log((FT(2) * ζ + _c_h - B_h) / (FT(2) * ζ + _c_h + B_h))
        log_term_2 = log((_c_h - B_h) / (_c_h + B_h))
        log_terms = log_term_1 - log_term_2
        term_2 = _b_h / FT(2) * log1p(_c_h * ζ + ζ^2)
        return -coeff * log_terms - term_2
    else
        # Fallback to Businger form with γ = 9
        return _psi_h_unstable(ζ, FT(9))
    end
end

end # module
