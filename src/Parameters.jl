module Parameters

import Thermodynamics
const TD = Thermodynamics
const TDPS = TD.Parameters.ThermodynamicsParameters
import ..UniversalFunctions
const UF = UniversalFunctions

abstract type AbstractSurfaceFluxesParameters end
const ASFP = AbstractSurfaceFluxesParameters

Base.broadcastable(ps::ASFP) = tuple(ps)

Base.@kwdef struct SurfaceFluxesParameters{FT, AUFPS <: UF.AbstractUniversalFunctionParameters{FT}, TP} <:
                   AbstractSurfaceFluxesParameters
    von_karman_const::FT
    ufp::AUFPS
    thermo_params::TP
end

thermodynamics_params(ps::SurfaceFluxesParameters) = ps.thermo_params
uf_params(ps::SurfaceFluxesParameters) = ps.ufp
von_karman_const(ps::SurfaceFluxesParameters) = ps.von_karman_const

universal_func_type(::SurfaceFluxesParameters{FT, AUFPS}) where {FT, AUFPS} = UF.universal_func_type(AUFPS)

for var in fieldnames(TDPS)
    @eval $var(ps::ASFP) = TD.Parameters.$var(thermodynamics_params(ps))
end

# Thermodynamics derived parameters
molmass_ratio(ps::ASFP) = TD.Parameters.molmass_ratio(thermodynamics_params(ps))
R_d(ps::ASFP) = TD.Parameters.R_d(thermodynamics_params(ps))
cp_d(ps::ASFP) = TD.Parameters.cp_d(thermodynamics_params(ps))

Pr_0(ps::SurfaceFluxesParameters) = Pr_0(uf_params(ps))
a_m(ps::SurfaceFluxesParameters) = a_m(uf_params(ps))
a_h(ps::SurfaceFluxesParameters) = a_h(uf_params(ps))
b_m(ps::SurfaceFluxesParameters) = b_m(uf_params(ps))
b_h(ps::SurfaceFluxesParameters) = b_h(uf_params(ps))
ζ_a(ps::SurfaceFluxesParameters) = ζ_a(uf_params(ps))
γ(ps::SurfaceFluxesParameters) = γ(uf_params(ps))

end
