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

const FTypes = Union{Real, AbstractArray}

abstract type AbstractUniversalFunction{FT <: FTypes} end
const AUF = AbstractUniversalFunction

#=
    AbstractUniversalFunctionType

Internal abstract type. Subtypes mirror `AbstractUniversalFunction`s.
These mirrored subtypes are needed due to several constraints:
 - Types we pass in concrete types to avoid UnionAll types (which incur allocations)
 - We cannot pass in concrete types, e.g. `Businger{FT, typeof(param_set)}`
   because we must be able to use ForwardDiff, which requires changing `FT`.
=#
abstract type AbstractUniversalFunctionType end
const AUFT = AbstractUniversalFunctionType

abstract type AbstractUniversalFunctionParameters{FT <: Real} end
const AUFP = AbstractUniversalFunctionParameters

Base.eltype(uf::AbstractUniversalFunction{FT}) where {FT} = FT

#####
##### Interface
#####

abstract type AbstractTransportType end
struct MomentumTransport <: AbstractTransportType end
struct HeatTransport <: AbstractTransportType end

Base.broadcastable(tt::AbstractUniversalFunction) = Ref(tt)
Base.broadcastable(tt::AbstractTransportType) = Ref(tt)

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

Pr_0(uf::AUF) = uf.params.Pr_0
a_m(uf::AUF) = uf.params.a_m
a_h(uf::AUF) = uf.params.a_h
b_m(uf::AUF) = uf.params.b_m
b_h(uf::AUF) = uf.params.b_h
c_h(uf::AUF) = uf.params.c_h

π_group(uf::AUF, ::HeatTransport) = Pr_0(uf)
π_group(::AUF, ::MomentumTransport) = 1

#####
##### Businger
#####

Base.@kwdef struct BusingerParams{FT} <: AbstractUniversalFunctionParameters{FT}
    Pr_0::FT
    a_m::FT
    a_h::FT
end

"""
    Businger

# Reference

 - [Nishizawa2018](@cite)

# Original research

 - [Businger1971](@cite)

# Equations in reference:

    `ϕ_m`: Eq. A1
    `ϕ_h`: Eq. A2
    `ψ_m`: Eq. A3
    `ψ_h`: Eq. A4

# Fields

$(DSE.FIELDS)
"""
struct Businger{FT, PS <: BusingerParams} <: AbstractUniversalFunction{FT}
    "Monin-Obhukov Length"
    L::FT
    params::PS
end

struct BusingerType <: AbstractUniversalFunctionType end
Businger() = BusingerType()

f_momentum(uf::Businger, ζ) = sqrt(sqrt(1 - 15 * ζ))

f_heat(uf::Businger, ζ) = sqrt(1 - 9 * ζ)

function phi(uf::Businger, ζ, ::MomentumTransport)
    if ζ < 0
        f_m = f_momentum(uf, ζ)
        return 1 / f_m
    else
        FT = eltype(uf)
        _a_m = FT(a_m(uf))
        return _a_m * ζ + 1
    end
end

function phi(uf::Businger, ζ, tt::HeatTransport)
    if ζ < 0
        f_h = f_heat(uf, ζ)
        return 1 / f_h
    else
        FT = eltype(uf)
        _a_h = FT(a_h(uf))
        _π_group = FT(π_group(uf, tt))
        return _a_h * ζ / _π_group + 1
    end
end

function psi(uf::Businger, ζ, ::MomentumTransport)
    FT = eltype(uf.L)
    if ζ < 0
        f_m = f_momentum(uf, ζ)
        log_term = log((1 + f_m)^2 * (1 + f_m^2) / 8)
        return log_term - 2 * atan(f_m) + FT(π) / 2
    else
        _a_m = FT(a_m(uf))
        return -_a_m * ζ
    end
end

function psi(uf::Businger, ζ, tt::HeatTransport)
    if ζ < 0
        f_h = f_heat(uf, ζ)
        return 2 * log((1 + f_h) / 2)
    else
        FT = eltype(uf)
        _a_h = FT(a_h(uf))
        _π_group = FT(π_group(uf, tt))
        return -_a_h * ζ / _π_group
    end
end

function Psi(uf::Businger, ζ, tt::MomentumTransport)
    FT = eltype(uf)
    if abs(ζ) < eps(FT)
        # Psi_m in Eq. A13
        if ζ >= 0
            _a_m = FT(a_m(uf))
            return -_a_m * ζ / 2
        else
            return -FT(15) * ζ / FT(8)
        end
    else
        # Note that "1-f^3" in is a typo, it is
        # supposed to be "1-f_m^3". This was
        # confirmed by communication with the
        # author.
        if ζ < 0
            f_m = f_momentum(uf, ζ)
            log_term = log((1 + f_m)^2 * (1 + f_m^2) / 8)
            π_term = FT(π) / 2
            tan_term = 2 * atan(f_m)
            cubic_term = (1 - f_m^3) / (12 * ζ)
            return log_term - tan_term + π_term - 1 + cubic_term
        else
            _a_m = FT(a_m(uf))
            return -_a_m * ζ / 2
        end
    end
end

function Psi(uf::Businger, ζ, tt::HeatTransport)
    FT = eltype(uf)
    _a_h = FT(a_h(uf))
    if abs(ζ) < eps(typeof(uf.L))
        # Psi_h in Eq. A14
        if ζ >= 0
            _π_group = FT(π_group(uf, tt))
            return -_a_h * ζ / (2 * _π_group)
        else
            return -9 * ζ / 4
        end
    else
        if ζ >= 0
            _π_group = FT(π_group(uf, tt))
            return -_a_h * ζ / (2 * _π_group)
        else
            f_h = f_heat(uf, ζ)
            log_term = 2 * log((1 + f_h) / 2)
            return log_term + 2 * (1 - f_h) / (9 * ζ) - 1
        end
    end
end

#####
##### Gryanik
#####

Base.@kwdef struct GryanikParams{FT} <: AbstractUniversalFunctionParameters{FT}
    Pr_0::FT
    a_m::FT
    a_h::FT
    b_m::FT
    b_h::FT
end

"""
    Gryanik <: AbstractUniversalFunction{FT}

# References
 - [Gryanik2020](@cite)

# Equations in reference:

    `ϕ_m`: Eq. 13
    `ϕ_h`: Eq. 13
    `ψ_m`: Eq. 14
    `ψ_h`: Eq. 14

# Fields

$(DSE.FIELDS)
"""
struct Gryanik{FT, PS <: GryanikParams} <: AbstractUniversalFunction{FT}
    "Monin-Obhukov Length"
    L::FT
    params::PS
end

struct GryanikType <: AbstractUniversalFunctionType end
Gryanik() = GryanikType()

function phi(uf::Gryanik, ζ, tt::MomentumTransport)
    FT = eltype(uf)
    if 0 < ζ
        _a_m = FT(a_m(uf))
        _b_m = FT(b_m(uf))
        return 1 + (_a_m * ζ) / (1 + _b_m * ζ)^(FT(2 / 3))
    else
        _a_m = FT(a_m(uf))
        return _a_m * ζ + 1
    end
end

function phi(uf::Gryanik, ζ, tt::HeatTransport)
    FT = eltype(uf)
    if 0 < ζ
        _Pr_0 = FT(Pr_0(uf))
        _a_h = FT(a_h(uf))
        _b_h = FT(b_h(uf))
        return _Pr_0 * (1 + (_a_h * ζ) / (1 + _b_h * ζ))
    else
        _a_h = FT(a_h(uf))
        _π_group = FT(π_group(uf, tt))
        return _a_h * ζ / _π_group + 1
    end
end

function psi(uf::Gryanik, ζ, tt::MomentumTransport)
    FT = eltype(uf)
    if 0 < ζ
        _a_m = FT(a_m(uf))
        _b_m = FT(b_m(uf))
        return -3 * (_a_m / _b_m) * ((1 + _b_m * ζ)^(FT(1 / 3)) - 1)
    else
        _a_m = FT(a_m(uf))
        return -_a_m * ζ
    end
end

function psi(uf::Gryanik, ζ, tt::HeatTransport)
    FT = eltype(uf)
    if 0 < ζ
        _Pr_0 = FT(Pr_0(uf))
        _a_h = FT(a_h(uf))
        _b_h = FT(b_h(uf))
        return -_Pr_0 * (_a_h / _b_h) * log1p(_b_h * ζ)
    else
        _a_h = FT(a_h(uf))
        _π_group = FT(π_group(uf, tt))
        return -_a_h * ζ / _π_group
    end
end


#####
##### Grachev
#####

Base.@kwdef struct GrachevParams{FT} <: AbstractUniversalFunctionParameters{FT}
    Pr_0::FT
    a_m::FT
    a_h::FT
    b_m::FT
    b_h::FT
    c_h::FT
end

"""
    Grachev <: AbstractUniversalFunction{FT}

# References
 - [Grachev2007](@cite)

Equations in reference:

    `ϕ_m`: Eq. 13
    `ϕ_h`: Eq. 13
    `ψ_m`: Eq. 14
    `ψ_h`: Eq. 14

# Fields

$(DSE.FIELDS)
"""
struct Grachev{FT, PS <: GrachevParams} <: AbstractUniversalFunction{FT}
    "Monin-Obhukov Length"
    L::FT
    params::PS
end

struct GrachevType <: AbstractUniversalFunctionType end
Grachev() = GrachevType()

function phi(uf::Grachev, ζ, tt::MomentumTransport)
    if 0 < ζ
        FT = eltype(uf)
        _a_m = FT(a_m(uf))
        _b_m = FT(b_m(uf))
        return 1 + _a_m * ζ * (1 + ζ)^FT(1 / 3) / (1 + _b_m * ζ)
    else
        FT = eltype(uf)
        _a_m = FT(a_m(uf))
        return _a_m * ζ + 1
    end
end

function phi(uf::Grachev, ζ, tt::HeatTransport)
    if 0 < ζ
        FT = eltype(uf)
        _a_h = FT(a_h(uf))
        _b_h = FT(b_h(uf))
        _c_h = FT(c_h(uf))
        return 1 + (_a_h * ζ + _b_h * ζ^2) / (1 + _c_h * ζ + ζ^2)
    else
        FT = eltype(uf)
        _a_h = FT(a_h(uf))
        _π_group = FT(π_group(uf, tt))
        return _a_h * ζ / _π_group + 1
    end
end

function psi(uf::Grachev, ζ, tt::MomentumTransport)
    FT = eltype(uf)
    if 0 < ζ
        _a_m = FT(a_m(uf))
        _b_m = FT(b_m(uf))
        B_m = cbrt(1 / _b_m - 1)
        x = cbrt(1 + ζ)
        sqrt3 = FT(sqrt(3))
        # Note: there is a mismatch between
        # Gryanik, Eq. 26, and Grachev Eq. 12.
        # We use the Grachev Eq. 12.
        linear_term = -3 * (_a_m / _b_m) * (x - 1)
        log_term_1 = 2 * log((x + B_m) / (1 + B_m))
        log_term_2 = log((x^2 - x * B_m + B_m^2) / (1 - B_m + B_m^2))
        atan_term_1 = atan((2 * x - B_m) / (sqrt3 * B_m))
        atan_term_2 = atan((2 - B_m) / (sqrt3 * B_m))
        atan_terms = atan_term_1 - atan_term_2
        bracket_term = log_term_1 - log_term_2 + 2 * sqrt3 * atan_terms
        return linear_term + _a_m * B_m / (2 * _b_m) * bracket_term
    else
        _a_m = FT(a_m(uf))
        return -_a_m * ζ
    end
end

function psi(uf::Grachev, ζ, tt::HeatTransport)
    if 0 < ζ
        FT = eltype(uf)
        _Pr_0 = FT(Pr_0(uf))
        _a_h = FT(a_h(uf))
        _b_h = FT(b_h(uf))
        _c_h = FT(c_h(uf))
        B_h = sqrt(_c_h^2 - 4)
        coeff = _a_h / B_h - _b_h * _c_h / (2 * B_h)
        log_term_1 = log((2 * ζ + _c_h - B_h) / (2 * ζ + _c_h + B_h))
        log_term_2 = log((_c_h - B_h) / (_c_h + B_h))
        log_terms = log_term_1 - log_term_2
        term_2 = _b_h / 2 * log1p(_c_h * ζ + ζ^2)
        return -coeff * log_terms - term_2
    else
        FT = eltype(uf)
        _a_h = FT(a_h(uf))
        _π_group = FT(π_group(uf, tt))
        return -_a_h * ζ / _π_group
    end
end

universal_func(::BusingerType, L_MO::Real, params::BusingerParams) = Businger(L_MO, params)
universal_func(::GryanikType, L_MO::Real, params::GryanikParams) = Gryanik(L_MO, params)
universal_func(::GrachevType, L_MO::Real, params::GrachevParams) = Grachev(L_MO, params)

universal_func_type(::Type{T}) where {T <: BusingerParams} = BusingerType()
universal_func_type(::Type{T}) where {T <: GryanikParams} = GryanikType()
universal_func_type(::Type{T}) where {T <: GrachevParams} = GrachevType()

end # module
