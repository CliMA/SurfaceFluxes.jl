module Parameters

import Thermodynamics
const TD = Thermodynamics
const TDPS = TD.Parameters.ThermodynamicsParameters
import ..UniversalFunctions
const UF = UniversalFunctions

"""
    AbstractSurfaceFluxesParameters{FT}

Abstract type for surface fluxes parameters.
"""
abstract type AbstractSurfaceFluxesParameters{FT} end
const ASFP = AbstractSurfaceFluxesParameters
Base.eltype(::ASFP{FT}) where {FT} = FT

Base.broadcastable(ps::ASFP) = tuple(ps)

"""
    SurfaceFluxesParameters{FT, AUFPS, TP}

Parameters for SurfaceFluxes.jl.

# Fields
- `von_karman_const`: Von Karman constant
- `ufp`: Universal function parameters
- `thermo_params`: Thermodynamics parameters
- `z0m_fixed`: [m] Fixed roughness length for momentum 
- `z0s_fixed`: [m] Fixed roughness length for scalars 
- `gustiness_coeff`: Scaling coefficient for Deardorff gustiness (often β)
- `gustiness_zi`: [m] Boundary layer height for gustiness (if fixed)
"""
Base.@kwdef struct SurfaceFluxesParameters{
    FT,
    AUFPS <: UF.AbstractUniversalFunctionParameters{FT},
    TP <: TDPS{FT},
} <: AbstractSurfaceFluxesParameters{FT}
    von_karman_const::FT
    ufp::AUFPS
    thermo_params::TP
    z0m_fixed::FT
    z0s_fixed::FT
    gustiness_coeff::FT
    gustiness_zi::FT
end

"Thermodynamics parameters"
thermodynamics_params(ps::SurfaceFluxesParameters) = ps.thermo_params
"Universal function parameters"
uf_params(ps::SurfaceFluxesParameters) = ps.ufp
"Von Karman constant"
von_karman_const(ps::SurfaceFluxesParameters) = ps.von_karman_const
"Fixed roughness length for momentum"
z0m_fixed(ps::SurfaceFluxesParameters) = ps.z0m_fixed
"Fixed roughness length for scalars"
z0s_fixed(ps::SurfaceFluxesParameters) = ps.z0s_fixed
"Gustiness coefficient"
gustiness_coeff(ps::SurfaceFluxesParameters) = ps.gustiness_coeff
"Fixed boundary layer height for Deardorff gustiness calculation"
gustiness_zi(ps::SurfaceFluxesParameters) = ps.gustiness_zi

for var in fieldnames(TDPS)
    @eval $var(ps::ASFP) = TD.Parameters.$var(thermodynamics_params(ps))
end

Pr_0(ps::SurfaceFluxesParameters) = UF.Pr_0(uf_params(ps))
a_m(ps::SurfaceFluxesParameters) = UF.a_m(uf_params(ps))
a_h(ps::SurfaceFluxesParameters) = UF.a_h(uf_params(ps))
b_m(ps::SurfaceFluxesParameters) = UF.b_m(uf_params(ps))
b_h(ps::SurfaceFluxesParameters) = UF.b_h(uf_params(ps))
c_m(ps::SurfaceFluxesParameters) = UF.c_m(uf_params(ps))
c_h(ps::SurfaceFluxesParameters) = UF.c_h(uf_params(ps))
d_m(ps::SurfaceFluxesParameters) = UF.d_m(uf_params(ps))
d_h(ps::SurfaceFluxesParameters) = UF.d_h(uf_params(ps))
ζ_a(ps::SurfaceFluxesParameters) = UF.ζ_a(uf_params(ps))
γ(ps::SurfaceFluxesParameters) = UF.γ(uf_params(ps))

end
