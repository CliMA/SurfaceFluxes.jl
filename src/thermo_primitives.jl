"""
    Thermo primitive helpers

Shared thermodynamic helper routines used by the functional surface-flux
solver. They operate directly on primitive scalars (`T`, `q`, `ρ`) and cached
thermodynamics parameters to avoid constructing intermediate state objects.
"""

@inline function dry_adiabatic_surface_density(
    cv_m::FT,
    R_m::FT,
    Tin::FT,
    ρin::FT,
    Ts::FT,
) where {FT}
    ratio = Ts / Tin
    return ρin * ratio^(cv_m / R_m)
end

@inline Δqt(q_in::FT, q_sfc::FT) where {FT} = q_in - q_sfc
@inline ΔT(T_in::FT, T_sfc::FT) where {FT} = T_in - T_sfc

@inline function θᵥ(param_set::APS, T::FT, ρ::FT, phase) where {FT}
    return TD.virtual_pottemp(param_set, T, ρ, phase)
end

@inline function virtual_temperature(param_set::APS, T::FT, phase) where {FT}
    return TD.virtual_temperature(param_set, T, phase)
end

@inline function Δθᵥ(
    param_set::APS,
    T_in::FT,
    ρ_in::FT,
    phase_in,
    T_sfc::FT,
    ρ_sfc::FT,
    phase_sfc,
) where {FT}
    return θᵥ(param_set, T_in, ρ_in, phase_in) -
           θᵥ(param_set, T_sfc, ρ_sfc, phase_sfc)
end

@inline function DSEᵥ(param_set::APS, T::FT, phase, Φ::FT) where {FT}
    cp_d = SFP.cp_d(param_set)
    return cp_d * virtual_temperature(param_set, T, phase) + Φ
end

@inline function ΔDSEᵥ(
    param_set::APS,
    T_in::FT,
    phase_in,
    Φ_in::FT,
    T_sfc::FT,
    phase_sfc,
    Φ_sfc::FT,
) where {FT}
    return DSEᵥ(param_set, T_in, phase_in, Φ_in) -
           DSEᵥ(param_set, T_sfc, phase_sfc, Φ_sfc)
end

