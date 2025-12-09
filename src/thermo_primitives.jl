"""
    Thermo primitive helpers

Shared thermodynamic helper routines used by the surface-flux solver. They operate directly 
on primitive scalars (`T`, `q`, `ρ`) and cached thermodynamics parameters.
"""

@inline function surface_density(
    cv_m,
    R_m,
    Tin,
    ρin,
    Ts,
)
    ratio = Ts / Tin
    return ρin * ratio^(cv_m / R_m)
end

@inline Δqt(q_in, q_sfc) = q_in - q_sfc
@inline ΔT(T_int, T_sfc) = T_int - T_sfc

@inline function θᵥ(param_set::APS, T, ρ, phase)
    return TD.virtual_pottemp(SFP.thermodynamics_params(param_set), T, ρ, phase)
end

@inline function virtual_temperature(param_set::APS, T, phase)
    return TD.virtual_temperature(SFP.thermodynamics_params(param_set), T, phase)
end

@inline function Δθᵥ(
    param_set::APS,
    T_int,
    ρ_int,
    phase_int,
    T_sfc,
    ρ_sfc,
    phase_sfc,
)
    return θᵥ(param_set, T_int, ρ_int, phase_int) -
           θᵥ(param_set, T_sfc, ρ_sfc, phase_sfc)
end

@inline function DSEᵥ(param_set::APS, T, phase, Φ)
    cp_d = SFP.cp_d(param_set)
    return cp_d * virtual_temperature(param_set, T, phase) + Φ
end

@inline function ΔDSEᵥ(
    param_set::APS,
    T_int,
    phase_int,
    Φ_int,
    T_sfc,
    phase_sfc,
    Φ_sfc,
)
    return DSEᵥ(param_set, T_int, phase_int, Φ_int) -
           DSEᵥ(param_set, T_sfc, phase_sfc, Φ_sfc)
end
