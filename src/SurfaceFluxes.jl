"""
    SurfaceFluxes

## Interface
  - [`surface_conditions`](@ref) computes
    - Monin-Obukhov length
    - Potential temperature flux (if not given) using Monin-Obukhov theory
    - transport fluxes using Monin-Obukhov theory
    - friction velocity/temperature scale/tracer scales
    - exchange coefficients

## References
 - [Nishizawa2018](@cite)

"""
module SurfaceFluxes

import NonlinearSolvers
const NS = NonlinearSolvers

import KernelAbstractions
const KA = KernelAbstractions

using DocStringExtensions
const DSE = DocStringExtensions

using Thermodynamics
const TD = Thermodynamics

import CLIMAParameters
const CP = CLIMAParameters
const CPP = CP.Planet
const APS = CP.AbstractEarthParameterSet
const CPSGS = CP.SubgridScale

import StaticArrays
const SA = StaticArrays

include("UniversalFunctions.jl")
import .UniversalFunctions
const UF = UniversalFunctions

abstract type SurfaceFluxesModel end

struct FVScheme end
struct DGScheme end

export b_star_from_ts, surface_conditions,
    exchange_coefficients, recover_profile, monin_obukhov_length, 
    monin_obukhov_length_from_b

"""
    SurfaceFluxConditions{FT}

Surface flux conditions, returned from `surface_conditions`.

# Fields

$(DSE.FIELDS)
"""
struct SurfaceFluxConditions{FT, VFT}
    L_MO::FT
    wb_flux_star::FT
    flux::VFT
    x_star::VFT
    C_exchange::VFT
end

function Base.show(io::IO, sfc::SurfaceFluxConditions)
    println(io, "----------------------- SurfaceFluxConditions")
    println(io, "L_MO           = ", sfc.L_MO)
    println(io, "wb_flux_star = ", sfc.wb_flux_star)
    println(io, "flux           = ", sfc.flux)
    println(io, "x_star         = ", sfc.x_star)
    println(io, "C_exchange     = ", sfc.C_exchange)
    println(io, "-----------------------")
end

function surface_fluxes_f!(F, x, nt)
    param_set = nt.param_set
    universal_func = nt.universal_func
    z_in = nt.z_in
    z_0 = length(nt.z_0) > 1 ? Tuple(nt.z_0) : nt.z_0
    n_vars = nt.n_vars
    scheme = nt.scheme
    u_in = nt.u_in
    u_sfc = nt.u_sfc
    ts_sfc = nt.ts_sfc
    ts_in = nt.ts_in
    LMO_init = nt.L_MO_init

    x_tup = Tuple(x)

    L_MO = x_tup[1]
    uf = universal_func(param_set, L_MO)
    Δz = z_in
    u_star = u_in * compute_physical_scale_coeff(
                                           uf,
                                           z_in, 
                                           z_0[1],
                                           L_MO,
                                           UF.MomentumTransport(),
                                           scheme)
    b_star = b_star_from_ts(param_set, 
                                    z_in, z_0[3], 
                                    ts_sfc, ts_in, 
                                    L_MO)
    L_MO = monin_obukhov_length_from_b(param_set, u_star, b_star)
    u_star = u_in * compute_physical_scale_coeff(
                                           uf,
                                           z_in, 
                                           z_0[1],
                                           L_MO,
                                           UF.MomentumTransport(),
                                           scheme
                                          )
    F_nt =
        x_tup[1] - monin_obukhov_length_from_b(
            param_set,
            u_star,
            b_star
        )
    F .= F_nt
end

function surface_conditions(
    param_set::APS,
    L_MO_init::FT,
    ts_in::TD.ThermodynamicState,
    ts_sfc::TD.ThermodynamicState,
    u_in, u_sfc,
    z_0::Union{AbstractVector, FT},
    z_in::FT,
    scheme,
    universal_func::Union{Nothing, F} = UF.Businger,
    sol_type::NS.SolutionType = NS.CompactSolution(),
    tol::NS.AbstractTolerance = NS.ResidualTolerance{FT}(sqrt(eps(FT))),
    maxiter::Int = 10,
) where {FT <: AbstractFloat, APS, F}

    # Nonlinear Solver Solution in Local Scope
    local sol

    uf = universal_func(param_set, L_MO_init)
    # Unpack thermodynamic parameters and constants
    von_karman_const::FT = CPSGS.von_karman_const(uf.param_set)

    # In this implementation we expect a user to provide only the 
    # initial guess for Monin Obukhov length (initial guesses for 
    # tracer scale quantities are then constrained by this L_MO_init)
    n_vars = 2
    z_0b = z_0[3]

    local sol

    args = (;
        param_set,
        L_MO_init,
        u_in,  u_sfc,
        ts_in, ts_sfc, 
        z_in,
        z_0,
        n_vars,
        scheme,
        universal_func,
    )
    
    MO_param_guess =
        Array(FT[L_MO_init])

    # Define closure over args
    f!(F, x_all) = surface_fluxes_f!(F, x_all, args)

    nls = NS.NewtonsMethodAD(f!, MO_param_guess)
    sol = NS.solve!(nls, sol_type, tol, maxiter)

    root_tup = Tuple(sol.root)
    if sol.converged
        L_MO = root_tup[1] 
    else
        KA.@print("Warning: Unconverged Surface Fluxes\n")
        L_MO = root_tup[1]
    end

    @show L_MO

    # TODO REVAMP get_flux_coefficients()
    #C_exchange = get_flux_coefficients(
    #    param_set,
    #    z_in,
    #    x_star,
    #    x_s,
    #    L_MO,
    #    z_0,
    #    scheme,
    #    universal_func,
    #)
    #VFT = typeof(flux)
    #return SurfaceFluxConditions{FT, VFT}(
    #    L_MO,
    #    wb_flux_star,
    #    flux,
    #    x_star,
    #    C_exchange,
    #)
end

"""
    get_energy_flux(surf_conds::SurfaceFluxConditions{FT, VFT})

Returns the potential temperature and q_tot fluxes needed to
compute the energy flux.
"""
function get_energy_flux(surf_conds::SurfaceFluxConditions)
    b_flux = surf_conds.flux[2]
    q_tot_flux = surf_conds.flux[3]
    return b_flux, q_tot_flux
end

"""
    compute_physical_scale(uf, z_in, z_0, x_in, x_s, transport, ::FVScheme)

Returns u_star, b_star, ... given equations (17), (18) for FV and (8), (9) for DG
from Nishizawa & Kitamura (2018). Consistent with computation of fluxes using 
buoyancy scale in FMS GCM code. Allows extension to multiple tracers by including
their influence on the buoyancy scale parameter.

## Arguments
 - `uf` universal function family
 - `z_in` Input height for the similarity functions. It is Δz for FV and the height
    of the first nodal point for DG
 - `z_0` roughness lengths for state variable array `x`
 - `x_in` Inner values for state variable array `x`. They correspond to volume averages in FV and
    nodal point values in DG
 - `x_s` surface values for state variable array `x`
 - `transport` Type of transport (momentum or heat)

"""
function compute_physical_scale(
    uf::UF.AbstractUniversalFunction{FT},
    z_in,
    z_0,
    x_in,
    x_s,
    transport,
    ::FVScheme,
) where {FT}
    von_karman_const::FT = CPSGS.von_karman_const(uf.param_set)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    R_z0 = 1 - z_0 / z_in
    temp1 = log(z_in / z_0)
    temp2 = -UF.Psi(uf, z_in / uf.L, transport)
    temp3 = z_0 / z_in * UF.Psi(uf, z_0 / uf.L, transport)
    temp4 = R_z0 * (UF.psi(uf, z_0 / uf.L, transport) - 1)
    Σterms = temp1 + temp2 + temp3 + temp4
    return _π_group⁻¹ * von_karman_const / Σterms * (x_in-x_s)
end

function compute_physical_scale(
    uf::UF.AbstractUniversalFunction{FT},
    z_in,
    z_0,
    x_in,
    x_s,
    transport,
    ::DGScheme,
) where {FT}
    von_karman_const::FT = CPSGS.von_karman_const(uf.param_set)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    temp1 = log(z_in / z_0)
    temp2 = -UF.psi(uf, z_in / uf.L, transport)
    temp3 = UF.psi(uf, z_0 / uf.L, transport)
    Σterms = temp1 + temp2 + temp3
    return _π_group⁻¹ * von_karman_const / Σterms * (x_in - x_s)
end

function compute_physical_scale_coeff(
    uf::UF.AbstractUniversalFunction{FT},
    z_in,
    z_0,
    L_MO,
    transport, 
    ::FVScheme,
) where {FT}
    von_karman_const::FT = CPSGS.von_karman_const(uf.param_set)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    R_z0 = 1 - z_0 / z_in
    temp1 = log(z_in / z_0)
    temp2 = -UF.Psi(uf, z_in / uf.L, transport)
    temp3 = z_0 / z_in * UF.Psi(uf, z_0 / uf.L, transport)
    temp4 = R_z0 * (UF.psi(uf, z_0 / uf.L, transport) - 1)
    Σterms = temp1 + temp2 + temp3 + temp4
    return _π_group⁻¹ * von_karman_const / Σterms
end

function compute_physical_scale_coeff(
    uf::UF.AbstractUniversalFunction{FT},
    z_in,
    z_0,
    L_MO, 
    transport,
    ::DGScheme,
) where {FT}
    von_karman_const::FT = CPSGS.von_karman_const(uf.param_set)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    temp1 = log(z_in / z_0)
    temp2 = -UF.psi(uf, z_in / uf.L, transport)
    temp3 = UF.psi(uf, z_0 / uf.L, transport)
    Σterms = temp1 + temp2 + temp3
    return _π_group⁻¹ * von_karman_const / Σterms
end

"""
    recover_profile(z, x_star, x_s, z_0, L_MO, transport, ::DGScheme, uf)

Recover vertical profiles u(z), θ(z), ... using equations (4) and (5)
for DG and (12), (13) for FV from Nishizawa & Kitamura (2018).

## Arguments
 - `z` Input height for evaluation
 - `x_star` Surface layer scales for state variable array `x`.
 - `x_s` surface values for state variable array `x`
 - `z_0` roughness lengths for state variable array `x`
 - `L_MO` Obukhov length
 - `transport` Type of transport (momentum or heat)
 - `uf` universal function family

"""
function recover_profile(
    param_set::APS,
    z,
    x_star,
    x_s,
    z_0::Union{AbstractVector, FT},
    L_MO,
    transport,
    ::DGScheme,
    universal_func = UF.Businger,
) where {FT}
    uf = universal_func(param_set, L_MO)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    _π_group = FT(UF.π_group(uf, transport))
    temp1 = log(z / z_0)
    temp2 = -UF.psi(uf, z / uf.L, transport)
    temp3 = UF.psi(uf, z_0 / uf.L, transport)
    Σterms = temp1 + temp2 + temp3
    return _π_group * x_star * Σterms / von_karman_const + x_s
end

function recover_profile(
    param_set::APS,
    z,
    x_star,
    x_s,
    z_0::Union{AbstractVector, FT},
    L_MO,
    transport,
    ::FVScheme,
    universal_func = UF.Businger,
) where {FT}
    uf = universal_func(param_set, L_MO)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    _π_group = FT(UF.π_group(uf, transport))
    R_z0 = 1 - z_0 / z
    temp1 = log(z / z_0)
    temp2 = -UF.Psi(uf, z / uf.L, transport)
    temp3 = z_0 / z * UF.Psi(uf, z_0 / uf.L, transport)
    temp4 = R_z0 * (UF.psi(uf, z_0 / uf.L, transport) - 1)
    Σterms = temp1 + temp2 + temp3 + temp4
    return _π_group * x_star * Σterms / von_karman_const + x_s
end

### Thermodynamic variables and scale parameters

function b_star_from_ρ(param_set, 
                       ts_sfc, ts_in, 
                       ::FVScheme
                       )
    # Unpack air densities from thermodynamic states
    ρ_in = TD.air_density(ts_in)
    ρ_sfc = TD.air_density(ts_sfc)
    # Get buoyancy scale parameter from surface and interior densities
    ρ_star = compute_physical_scale_coeff(uf, Δz, z_0, L_MO)*(ρ_in - ρ_sfc)
    b_star = -grav * ρ_star/ ρ_in
    return b_star
end

function b_star_from_T(param_set, 
                       ts_sfc, ts_in,
                       ::DGScheme
                       )
    # Unpack air densities from thermodynamic states
    ρ_in = TD.air_density(ts_in)
    ρ_sfc = TD.air_density(ts_sfc)
    # Get buoyancy scale parameter from surface and interior densities
    ρ_star = compute_physical_scale_coeff(uf, Δz, z_0, L_MO)*(ρ_in - ρ_sfc)
    b_star = -grav * ρ_star/ ρ_in
    return b_star
end

"""
    b_star_from_θ(param_set, θ_scale, qt_scale, qt_star)
Returns b⋆ (calculated from potential temperature scale)
If θ is θᵥ then moisture is accounted for. 
"""

function b_star_from_θ(param_set, 
                               u_star, 
                               θ_star, θ_scale, 
                               qt_scale, qt_star)
    FT = typeof(θ_scale)
    grav::FT = CPP.grav(param_set)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    ϵ::FT = CPP.molmass_ratio(param_set)
    return (1 + (ϵ-1)*qt_scale)*θ_star + (ϵ-1)*θ_scale*qt_star
end

"""
    b_star_from_T(param_set, 
                               p_sfc, ρ_sfc, 
                               T_sfc, T_star, 
                               R_d, R_m)
Returns b⋆ (calculated from potential temperature scale) 
"""
function b_star_from_T(param_set, 
                               p_sfc, ρ_sfc, 
                               T_sfc, T_in,
                               R_d, R_m, 
                               )
    FT = typeof(T_star)
    grav::FT = CPP.grav(param_set)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    ϵ::FT = CPP.molmass_ratio(param_set)
    T_star = compute_physical_scale_coeff(uf, Δz, z_0, L_MO)*(T_in - T_sfc)
    qt_star = compute_physical_scale_coeff(uf, Δz, z_0, L_MO)*(qt_in - qt_sfc)
    qc_star = compute_physical_scale_coeff(uf, Δz, z_0, L_MO)*(qc_in - qc_sfc)
    b_star = -grav * p_in/ρ_in/R_m/T_in * (T_star/T_in - R_d/R_m*(ϵ-1)*qt_star + R_d/R_m *ϵ*qc_star) 
    return b_star
end

"""
    b_star_from_ts(param_set, Δz, z0,
                               ts_sfc, 
                               ts_in, 
                               L_MO)
Returns b⋆ (calculated from thermodynamic states and scales for moisture variables)
ts_sfc is the thermodynamic state at the surface
ts_in is the thermodynamic state at the first input state.
"""

function b_star_from_ts(param_set, 
                               Δz, z_0,
                               ts_sfc, 
                               ts_in,
                               L_MO)
    FT = typeof(L_MO)
    universal_func = UF.Businger
    uf = universal_func(param_set, L_MO)
    grav::FT = CPP.grav(param_set)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    ϵ::FT = CPP.molmass_ratio(param_set)
    
    # Unpack Thermodynamic States
    p_in = TD.air_pressure(ts_in)
    qt_in = TD.total_specific_humidity(ts_in)
    T_in = TD.air_temperature(ts_in)
    T_sfc = TD.air_temperature(ts_sfc)
    qt_sfc = TD.total_specific_humidity(ts_sfc)
    ρ_in = TD.air_density(ts_in)
    ρ_sfc = TD.air_density(ts_sfc)
    R_d = CPP.R_d(param_set)
    R_m = TD.gas_constant_air(ts_in)
    
    T_star = compute_physical_scale_coeff(uf, Δz, z_0, L_MO, UF.HeatTransport(), DGScheme())*(T_in - T_sfc)
    qt_star = compute_physical_scale_coeff(uf, Δz, z_0, L_MO, UF.HeatTransport(), DGScheme())*(qt_in - qt_sfc)
   #qc_star = compute_physical_scale_coeff(uf, Δz, z_0, L_MO)*(qc_in - qc_sfc)
    
    return grav * p_in/ρ_in/R_m/T_in * (T_star/T_in - R_d/R_m*(ϵ-1)*qt_star)# + R_d/R_m *ϵ*qc_star) 
end


###
### Monin Obukhov Length Calculation
###
"""
    monin_obukhov_length(param_set, u_star, θ_star, θ_scale, qt_scale, qt_star)

Returns the Monin-Obukhov length.
"""
function monin_obukhov_length(param_set::APS, 
                             u_star, 
                             θ_star, θ_scale, 
                             qt_scale, qt_star)

    FT = typeof(u_star)
    grav::FT = CPP.grav(param_set)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    b_star = b_star_from_θ(param_set::APS,
                                           u_star, 
                                           θ_star, θ_scale, 
                                           qt_scale, qt_star)
    return u_star^2 / (von_karman_const * b_star + eps(FT))
end

"""
    monin_obukhov_length_from_b(param_set, u_star, b_star)

Returns the Monin-Obukhov length as a function of b_star
"""
function monin_obukhov_length_from_b(param_set::APS, 
                                     u_star, 
                                     b_star)
    FT = typeof(u_star)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    return u_star^2 / (von_karman_const * b_star + eps(FT))
end

monin_obukhov_length(sfc::SurfaceFluxConditions) = sfc.L_MO

"""
    compute_u_star(Δu, Δz, z_0, L_MO)

Δu = ⟨u⟩ - u_sfc
Δz = z - z_sfc
z_0 is the aerodynamic roughness length corresponding to momentum transport
Given the Monin-Obukhov lengthscale, compute u_star, 
the velocity scale. 
"""
function compute_u_star(param_set, 
                        Δu, Δz, z_0, 
                        L_MO) 
    uf = universal_func(param_set, L_MO)
    return Δu * compute_physical_scale_coeff(uf, Δz, z_0, L_MO)
end

"""
    compute_X_star(ΔX, Δz, z_0, L_MO)

ΔX = ⟨X⟩ - X_sfc
X is an arbitrary tracer variable, assumes Monin-Obukhov length 
is known. 
z_0 is the aerodynamic roughness length corresponding to variable X.
"""
function compute_X_star(param_set, 
                        ΔX, Δz, z_0, 
                        L_MO) 
    uf = universal_func(param_set, L_MO)
    return ΔX / Pr_t * compute_physical_scale_coeff(uf, Δz, z_0, L_MO)
end

"""
    exchange_coefficients(
        param_set,
        z,
        F_exchange,
        x_star::VFT,
        L_MO,
        universal_func,
    )

Computes exchange transfer coefficients
  - `K_exchange` exchange coefficients
"""
function exchange_coefficients(
    param_set,
    z,
    F_exchange,
    x_star::VFT,
    L_MO,
    universal_func,
) where {VFT}
    N = length(F_exchange)
    FT = typeof(z)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    uf = universal_func(param_set, L_MO)
    x_star_tup = Tuple(x_star)
    K_exchange = similar(x_star)
    F_exchange_tup = Tuple(F_exchange)
    K_exchange .= ntuple(Val(length(x_star))) do i
        transport = i == 1 ? UF.MomentumTransport() : UF.HeatTransport()
        phi_t = phi(uf, z / L_MO, transport)
        _π_group = FT(UF.π_group(uf, transport))
        num = -F_exchange_tup[i] * von_karman_const * z
        den = _π_group * (x_star_tup[i] * phi_t)
        K_exch = num / den # Eq. 19 in
    end
    return K_exchange
end

"""
    get_flux_coefficients(
        param_set,
        z,
        x_star::VFT,
        L_MO,
        z0,
        universal_func,
    )

Returns the exchange coefficients for bulk transfer formulas associated
with the frictional parameters `x_star` and the surface layer similarity
profiles. Taken from equations (36) and (37) of Byun (1990).
"""
function get_flux_coefficients(
    param_set,
    z_in,
    x_star::VFT,
    x_s,
    L_MO,
    z0::Union{AbstractVector, FT},
    scheme,
    universal_func,
) where {VFT, FT}
    N = length(x_star)
    z0_tup = Tuple(z0)
    x_star_tup = Tuple(x_star)
    x_s_tup = Tuple(x_s)
    u_in = recover_profile(
        param_set,
        z_in,
        x_star_tup[1],
        x_s_tup[1],
        (length(z0) > 1 ? z0_tup[1] : z0),
        L_MO,
        UF.MomentumTransport(),
        scheme,
        universal_func,
    )
    C = similar(x_star)
    C .= ntuple(Val(length(x_star))) do i
        if i == 1
            C_i = x_star[i]^2 / (u_in - x_s_tup[i])^2 # Eq. (36) Byun(1990)
        else
            ϕ_in = recover_profile(
                param_set,
                z_in,
                x_star_tup[i],
                x_s_tup[i],
                (length(z0) > 1 ? z0_tup[i] : z0),
                L_MO,
                UF.HeatTransport(),
                scheme,
                universal_func,
            )
            C_i =
                x_star[1] * x_star[i] / (u_in - x_s_tup[1]) /
                (ϕ_in - x_s_tup[i]) # Eq. (37) Byun(1990)
        end
        C_i
    end
    return C
end

end # SurfaceFluxes module
