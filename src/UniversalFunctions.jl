"""
    UniversalFunctions

Universal stability and stability correction
functions for `SurfaceFluxes` module. Supports
universal functions:
 - `Businger`
 - `Gryanik`
 - `Grachev`
"""
module UniversalFunctions

import DocStringExtensions
const DSE = DocStringExtensions

abstract type AbstractUniversalFunctionParameters{FT <: Real} end
const AUFP = AbstractUniversalFunctionParameters

#####
##### Interface
#####

abstract type AbstractTransportType end
struct MomentumTransport <: AbstractTransportType end
struct HeatTransport <: AbstractTransportType end

Base.broadcastable(tt::AbstractTransportType) = tuple(tt)
Base.broadcastable(p::AbstractUniversalFunctionParameters) = tuple(p)

"""
    phi

Universal stability function for wind shear
(`ϕ_m`) and temperature gradient (`ϕ_h`)
"""
function phi end

"""
    psi

Universal stability correction function for
momentum (`ψ_m`) and heat (`ψ_h`)
"""
function psi end

"""
    Psi

Integral of universal stability correction
function for momentum (`ψ_m`) and heat (`ψ_h`)
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
##### Businger
#####
"""
    BusingerParams{FT} <: AbstractUniversalFunctionParameters{FT}

Free parameters for the Businger universal stability and stability correction
functions.

# Fields

$(DSE.FIELDS)
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

# Nishizawa2018 Eq. A7
f_momentum(p::BusingerParams, ζ) = sqrt(sqrt(1 - typeof(ζ)(b_m(p)) * ζ))

# Nishizawa2018 Eq. A8
f_heat(p::BusingerParams, ζ) = sqrt(1 - typeof(ζ)(b_h(p)) * ζ)

function phi(p::BusingerParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    if ζ < 0
        # Businger1971 Eq. A1 (ζ < 0)
        f_m = f_momentum(p, ζ)
        return FT(1) / FT(f_m)
    else
        # Businger1971 Eq. A1 (ζ >= 0)
        _a_m = FT(a_m(p))
        return _a_m * ζ + FT(1)
    end
end

function phi(p::BusingerParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if ζ < 0
        # Businger1971 Eq. A2 (ζ < 0)
        f_h = f_heat(p, ζ)
        return FT(1) / FT(f_h)
    else
        # Businger1971 Eq. A2 (ζ >= 0)
        _a_h = FT(a_h(p))
        _π_group = FT(π_group(p, tt))
        return _a_h * ζ / _π_group + FT(1)
    end
end

function psi(p::BusingerParams, ζ, ::MomentumTransport)
    FT = eltype(ζ)
    if abs(ζ) < eps(FT)
        return FT(0)
    end
    if ζ < 0
        # Businger1971 Eq. A3 (ζ < 0)
        f_m = f_momentum(p, ζ)
        log_term = log((1 + f_m)^2 * (1 + f_m^2) / 8)
        return FT(log_term - 2 * atan(f_m) + π / 2)
    else
        # Businger1971 Eq. A3 (ζ >= 0)
        _a_m = FT(a_m(p))
        return -_a_m * ζ
    end
end

function psi(p::BusingerParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if abs(ζ) < eps(FT)
        return FT(0)
    end
    if ζ < 0
        # Businger1971 Eq. A4 (ζ < 0)
        f_h = f_heat(p, ζ)
        return FT(2 * log((1 + f_h) / 2))
    else
        # Businger1971 Eq. A4 (ζ >= 0)
        _a_h = FT(a_h(p))
        _π_group = FT(π_group(p, tt))
        return -_a_h * ζ / _π_group
    end
end

function Psi(p::BusingerParams, ζ, tt::MomentumTransport)
    FT = eltype(ζ)
    if ζ >= 0
        # Nishizawa2018 Eq. A5 and A13 (ζ >= 0)
        _a_m = FT(a_m(p))
        return -_a_m * ζ / 2
    else
        if abs(ζ) < eps(FT)
            # Nishizawa2018 Eq. A13 (ζ < 0)
            return -FT(b_m(p)) * ζ / FT(8)
        else
            # Nishizawa2018 Eq. A5 (ζ < 0)
            f_m = f_momentum(p, ζ)
            log_term = log((1 + f_m)^2 * (1 + f_m^2) / 8)
            π_term = FT(π) / 2
            tan_term = 2 * atan(f_m)
            cubic_term = (1 - f_m^3) / (12 * ζ)
            return FT(log_term - tan_term + π_term - 1 + cubic_term)
        end
    end
end

function Psi(p::BusingerParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    _a_h = FT(a_h(p))
    if ζ >= 0
        # Nishizawa2018 Eq. A6 and A14 (ζ >= 0)
        _π_group = FT(π_group(p, tt))
        return -_a_h * ζ / (2 * _π_group)
    else
        if abs(ζ) < eps(FT)
            # Nishizawa2018 Eq. A14 (ζ < 0)
            return -FT(b_h(p)) * ζ / 4
        else
            # Nishizawa2018 Eq. A6 (ζ < 0)
            f_h = f_heat(p, ζ)
            log_term = 2 * log((1 + f_h) / 2)
            return FT(log_term + 2 * (1 - f_h) / (FT(b_h(p)) * ζ) - 1)
        end
    end
end

#####
##### Gryanik
#####

"""
    GryanikParams{FT} <: AbstractUniversalFunctionParameters{FT}

Free parameters for the Gryanik universal stability and stability correction
functions.

# Fields

$(DSE.FIELDS)
"""
Base.@kwdef struct GryanikParams{FT} <: AbstractUniversalFunctionParameters{FT}
    Pr_0::FT
    a_m::FT
    a_h::FT
    b_m::FT
    b_h::FT
    ζ_a::FT
    γ::FT
end

# Nishizawa2018 Eq. A7
f_momentum(::GryanikParams, ζ) = sqrt(sqrt(1 - 15 * ζ))

# Nishizawa2018 Eq. A8
f_heat(::GryanikParams, ζ) = sqrt(1 - 9 * ζ)

function phi(p::GryanikParams, ζ, tt::MomentumTransport)
    FT = eltype(ζ)
    if ζ > 0
        # Gryanik2020 Eq. 32
        _a_m = FT(a_m(p))
        _b_m = FT(b_m(p))
        return FT(1) + (_a_m * ζ) / (1 + _b_m * ζ)^(FT(2 / 3))
    else
        # Nishizawa2018 Eq. A1 (ζ <= 0)
        f_m = f_momentum(p, ζ)
        return FT(1) / FT(f_m)
    end
end

function phi(p::GryanikParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if ζ > 0
        # Gryanik2020 Eq. 33
        _Pr_0 = FT(Pr_0(p))
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        return FT(1) + (ζ * _Pr_0 * _a_h) / (1 + _b_h * ζ)
    else
        # Nishizawa2018 Eq. A2 (ζ <= 0)
        f_h = f_heat(p, ζ)
        return FT(1) / FT(f_h)
    end
end

function psi(p::GryanikParams, ζ, tt::MomentumTransport)
    FT = eltype(ζ)
    if abs(ζ) < eps(FT)
        return FT(0)
    end
    if ζ > 0
        # Gryanik2020 Eq. 34
        _a_m = FT(a_m(p))
        _b_m = FT(b_m(p))
        return -3 * (_a_m / _b_m) * ((1 + _b_m * ζ)^(FT(1 / 3)) - 1)
    else
        # Nishizawa2018 Eq. A3 (ζ <= 0)
        f_m = f_momentum(p, ζ)
        log_term = log((1 + f_m)^2 * (1 + f_m^2) / 8)
        return FT(log_term - 2 * atan(f_m) + π / 2)
    end
end

function psi(p::GryanikParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if abs(ζ) < eps(FT)
        return FT(0)
    end
    if ζ > 0
        # Gryanik2020 Eq. 35
        _Pr_0 = FT(Pr_0(p))
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        return -_Pr_0 * (_a_h / _b_h) * log1p(_b_h * ζ)
    else
        # Nishizawa2018 Eq. A4 (ζ <= 0)
        f_h = f_heat(p, ζ)
        return FT(2 * log((1 + f_h) / 2))
    end
end

function Psi(p::GryanikParams, ζ, tt::MomentumTransport)
    FT = eltype(ζ)
    _a_m = FT(a_m(p))
    _b_m = FT(b_m(p))
    if ζ >= 0
        return 3 * (_a_m / _b_m) - FT(9) * _a_m * ((_b_m * ζ + FT(1))^(FT(4 / 3)) - 1) / ζ / FT(4) / _b_m^FT(2)
    else
        if abs(ζ) < eps(FT)
            return -FT(15) * ζ / FT(8)
        else
            f_m = f_momentum(p, ζ)
            log_term = log((1 + f_m)^2 * (1 + f_m^2) / 8)
            π_term = FT(π) / 2
            tan_term = 2 * atan(f_m)
            cubic_term = (1 - f_m^3) / (12 * ζ)
            return FT(log_term - tan_term + π_term - 1 + cubic_term)
        end
    end
end

function Psi(p::GryanikParams, ζ, tt::HeatTransport)
    FT = typeof(ζ)
    _a_h = FT(a_h(p))
    _b_h = FT(b_h(p))
    Pr0 = FT(Pr_0(p))
    if ζ >= 0
        return -_a_h / _b_h / ζ * Pr0 * ((1 / _b_h + ζ) * log1p(_b_h * ζ) - ζ)
    else
        if abs(ζ) < eps(FT)
            return -FT(9) * ζ / 4
        else
            f_h = f_heat(p, ζ)
            log_term = 2 * log((1 + f_h) / 2)
            return FT(log_term + 2 * (1 - f_h) / (9 * ζ) - FT(1))
        end
    end
end

#####
##### Grachev
#####

"""
    GrachevParams{FT} <: AbstractUniversalFunctionParameters{FT}

Free parameters for the Grachev universal stability and stability correction
functions.

# Fields

$(DSE.FIELDS)
"""
Base.@kwdef struct GrachevParams{FT} <: AbstractUniversalFunctionParameters{FT}
    Pr_0::FT
    a_m::FT
    a_h::FT
    b_m::FT
    b_h::FT
    c_h::FT
    ζ_a::FT
    γ::FT
end

# Nishizawa2018 Eq. A7
f_momentum(::GrachevParams, ζ) = sqrt(sqrt(1 - 15 * ζ))

# Nishizawa2018 Eq. A8
f_heat(::GrachevParams, ζ) = sqrt(1 - 9 * ζ)

function phi(p::GrachevParams, ζ, tt::MomentumTransport)
    FT = eltype(ζ)
    if ζ > 0
        # Grachev2007 Eq. 9a
        _a_m = FT(a_m(p))
        _b_m = FT(b_m(p))
        return FT(1) + _a_m * ζ * (1 + ζ)^FT(1 / 3) / (1 + _b_m * ζ)
    else
        # Nishizawa2018 Eq. A1 (ζ < 0)
        f_m = f_momentum(p, ζ)
        return FT(1) / FT(f_m)
    end
end

function phi(p::GrachevParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if ζ > 0
        # Grachev2007 Eq. 9b
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        _c_h = FT(c_h(p))
        return FT(1) + (_a_h * ζ + _b_h * ζ^2) / (1 + _c_h * ζ + ζ^2)
    else
        # Nishizawa2018 Eq. A2 (ζ < 0)
        f_h = f_heat(p, ζ)
        return FT(1) / FT(f_h)
    end
end

function psi(p::GrachevParams, ζ, tt::MomentumTransport)
    FT = eltype(ζ)
    if abs(ζ) < eps(FT)
        return FT(0)
    end
    if ζ > 0
        # Grachev2007 Eq. 12
        _a_m = FT(a_m(p))
        _b_m = FT(b_m(p))
        B_m = cbrt(1 / _b_m - 1)
        x = cbrt(1 + ζ)
        sqrt3 = FT(sqrt(3))
        linear_term = -3 * (_a_m / _b_m) * (x - 1)
        log_term_1 = 2 * log((x + B_m) / (1 + B_m))
        log_term_2 = log((x^2 - x * B_m + B_m^2) / (1 - B_m + B_m^2))
        atan_term_1 = atan((2 * x - B_m) / (sqrt3 * B_m))
        atan_term_2 = atan((2 - B_m) / (sqrt3 * B_m))
        atan_terms = atan_term_1 - atan_term_2
        bracket_term = log_term_1 - log_term_2 + 2 * sqrt3 * atan_terms
        return linear_term + _a_m * B_m / (2 * _b_m) * bracket_term
    else
        # Nishizawa2018 Eq. A3 (ζ < 0)
        f_m = f_momentum(p, ζ)
        log_term = log((1 + f_m)^2 * (1 + f_m^2) / 8)
        return FT(log_term - 2 * atan(f_m) + π / 2)
    end
end

function psi(p::GrachevParams, ζ, tt::HeatTransport)
    FT = eltype(ζ)
    if abs(ζ) < eps(FT)
        return FT(0)
    end
    if ζ > 0
        # Grachev2007 Eq. 13
        _Pr_0 = FT(Pr_0(p))
        _a_h = FT(a_h(p))
        _b_h = FT(b_h(p))
        _c_h = FT(c_h(p))
        B_h = sqrt(_c_h^2 - 4)
        coeff = _a_h / B_h - _b_h * _c_h / (2 * B_h)
        log_term_1 = log((2 * ζ + _c_h - B_h) / (2 * ζ + _c_h + B_h))
        log_term_2 = log((_c_h - B_h) / (_c_h + B_h))
        log_terms = log_term_1 - log_term_2
        term_2 = _b_h / 2 * log1p(_c_h * ζ + ζ^2)
        return -coeff * log_terms - term_2
    else
        # Nishizawa2018 Eq. A4 (ζ < 0)
        f_h = f_heat(p, ζ)
        return FT(2 * log((1 + f_h) / 2))
    end
end

end # module
