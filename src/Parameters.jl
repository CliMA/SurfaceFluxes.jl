module Parameters

import Thermodynamics
import ..UniversalFunctions
import ..UniversalFunctions: universal_func_type

const TD = Thermodynamics
const TDPS = TD.Parameters.ThermodynamicsParameters

const UF = UniversalFunctions

abstract type AbstractSurfaceFluxesParameters end
const ASFP = AbstractSurfaceFluxesParameters

Base.broadcastable(ps::ASFP) = tuple(ps)

const AUFP{FT} = UF.AbstractUniversalFunctionParameters{FT} where FT

Base.@kwdef struct SurfaceFluxesParameters{FT, U<:AUFP{FT}, TP} <: ASFP
    von_karman_const::FT
    ufp::U
    thermo_params::TP
end

thermodynamics_params(ps::SurfaceFluxesParameters) = ps.thermo_params
uf_params(ps::SurfaceFluxesParameters) = ps.ufp
von_karman_const(ps::SurfaceFluxesParameters) = ps.von_karman_const

universal_func_type(::SurfaceFluxesParameters{FT, UFP}) where {FT, UFP} =
    universal_func_type(UFP)

for var in fieldnames(TDPS)
    @eval $var(ps::ASFP) = TD.Parameters.$var(thermodynamics_params(ps))
end

# Thermodynamics derived parameters
molmass_ratio(ps::ASFP) = TD.Parameters.molmass_ratio(thermodynamics_params(ps))
R_d(ps::ASFP) = TD.Parameters.R_d(thermodynamics_params(ps))
cp_d(ps::ASFP) = TD.Parameters.cp_d(thermodynamics_params(ps))

Pr_0(ps::ASFP) = Pr_0(uf_params(ps))
a_m(ps::ASFP) = a_m(uf_params(ps))
a_h(ps::ASFP) = a_h(uf_params(ps))
b_m(ps::ASFP) = b_m(uf_params(ps))
b_h(ps::ASFP) = b_h(uf_params(ps))
c_m(ps::ASFP) = c_m(uf_params(ps))
c_h(ps::ASFP) = c_h(uf_params(ps))
d_m(ps::ASFP) = d_m(uf_params(ps))
d_h(ps::ASFP) = d_h(uf_params(ps))
ζ_a(ps::ASFP) = ζ_a(uf_params(ps))
γ(ps::ASFP) = γ(uf_params(ps))

end
