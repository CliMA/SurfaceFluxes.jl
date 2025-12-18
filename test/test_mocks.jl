# Shared mock parameter sets and utilities for SurfaceFluxes tests
#
# This module consolidates mock structs used across multiple test files
# to reduce code duplication and ensure consistency.

import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

"""
    SimpleMockThermoParams{FT}

Minimal thermodynamics parameters for testing.
"""
struct SimpleMockThermoParams{FT}
    cp_d::FT
    LH_v0::FT
    R_d::FT
end
Base.eltype(::SimpleMockThermoParams{FT}) where {FT} = FT

TDP.cp_d(tp::SimpleMockThermoParams) = tp.cp_d
TDP.LH_v0(tp::SimpleMockThermoParams) = tp.LH_v0
TD.cp_m(tp::SimpleMockThermoParams, q_tot, q_liq, q_ice) = tp.cp_d
TD.cp_m(tp::SimpleMockThermoParams, ts) = tp.cp_d
TD.cv_m(tp::SimpleMockThermoParams, q_tot, q_liq, q_ice) = tp.cp_d - tp.R_d
TD.gas_constant_air(tp::SimpleMockThermoParams, q_tot, q_liq, q_ice) = tp.R_d
TD.virtual_pottemp(tp::SimpleMockThermoParams, T, ρ, qt, ql, qi) = T
TD.dry_static_energy(tp::SimpleMockThermoParams, T, Φ) = tp.cp_d * T + Φ
TD.vapor_static_energy(tp::SimpleMockThermoParams, T, Φ) = tp.cp_d * T + Φ + tp.LH_v0

"""
    SimpleMockParamSet{FT, UFP, TP}

Mock SurfaceFluxes parameter set for unit testing.
"""
struct SimpleMockParamSet{FT, UFP, TP} <: SFP.AbstractSurfaceFluxesParameters{FT}
    von_karman_const::FT
    grav::FT
    gustiness_zi::FT
    gustiness_coeff::FT
    uf_params::UFP
    thermo_params::TP
end

SFP.von_karman_const(p::SimpleMockParamSet) = p.von_karman_const
SFP.grav(p::SimpleMockParamSet) = p.grav
SFP.gustiness_zi(p::SimpleMockParamSet) = p.gustiness_zi
SFP.gustiness_coeff(p::SimpleMockParamSet) = p.gustiness_coeff
SFP.uf_params(p::SimpleMockParamSet) = p.uf_params
SFP.thermodynamics_params(p::SimpleMockParamSet) = p.thermo_params
SFP.Rv_over_Rd(::SimpleMockParamSet{FT}) where {FT} = FT(1.61)
SFP.cp_d(p::SimpleMockParamSet) = p.thermo_params.cp_d

"""
    default_mock_param_set(FT)

Create a SimpleMockParamSet with Earth-like default values.
"""
function default_mock_param_set(::Type{FT}) where {FT}
    thermo = SimpleMockThermoParams(FT(1004.0), FT(2.5e6), FT(287.0))
    uf = UF.BusingerParams{FT}(
        Pr_0 = FT(0.74),
        a_m = FT(4.7),
        a_h = FT(4.7),
        b_m = FT(15.0),
        b_h = FT(9.0),
        ζ_a = FT(2.5),
        γ = FT(0.0),
    )
    return SimpleMockParamSet(
        FT(0.4),      # κ
        FT(9.81),     # g
        FT(1000.0),   # zi
        FT(1.2),      # gustiness_coeff
        uf,
        thermo,
    )
end

"""
    density_from_ideal_gas(thermo_params, T, p, q_tot)

Compute density from the ideal gas law.
"""
function density_from_ideal_gas(thermo_params, T, p::FT, q_tot) where {FT}
    R_m = TD.gas_constant_air(thermo_params, q_tot, FT(0), FT(0))
    return p / (R_m * T)
end
