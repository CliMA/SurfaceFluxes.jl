# Roughness-length evaluation
#
# These helpers resolve the active momentum and scalar roughness specifications.
# They are intentionally lightweight so they can be inlined inside GPU kernels.
#

"""
    RoughnessLengths{FT}

Struct to hold momentum and scalar roughness lengths.
"""
struct RoughnessLengths{FT}
    momentum::FT
    scalar::FT
end

"""
    ConstantRoughnessParams{FT} <: AbstractRoughnessParams

Roughness lengths fixed to constant values.

# Fields
- `z0m`: Momentum roughness length [m]
- `z0s`: Scalar roughness length [m]

The default values specified here are used when constructing the struct manually. When loading
via `ClimaParams`, these values are overwritten by the parameters in the TOML file.
"""
Base.@kwdef struct ConstantRoughnessParams{FT} <: AbstractRoughnessParams
    z0m::FT = 2e-4
    z0s::FT = 2e-5
end

"""
    COARE3RoughnessParams{FT} <: AbstractRoughnessParams

COARE 3.0 roughness parameterization (Fairall et al. 2003).

The default values specified here are used when constructing the struct manually. When loading
via `ClimaParams`, these values are overwritten by the parameters in the TOML file.

"""
Base.@kwdef struct COARE3RoughnessParams{FT} <: AbstractRoughnessParams
    kinematic_visc::FT = 1.5e-5
    z0m_default::FT = 1e-4
    α_low::FT = 0.011
    α_high::FT = 0.018
    u_low::FT = 10.0
    u_high::FT = 18.0
end

"""
    RaupachRoughnessParams <: AbstractRoughnessParams

Raupach (1994) canopy roughness model.
"""
Base.@kwdef struct RaupachRoughnessParams{FT} <: AbstractRoughnessParams
    C_R::FT = 0.3
    C_S::FT = 0.003
    c_d1::FT = 7.5
end

@inline function default_surface_flux_config(::Type{FT}) where {FT}
    return SurfaceFluxConfig(ConstantRoughnessParams{FT}(), ConstantGustinessSpec(FT(1.0)))
end

"""
    roughness_lengths(z0m, z0s)

Helper to construct `ConstantRoughnessParams`.
"""
@inline function roughness_lengths(z0m, z0s)
    return ConstantRoughnessParams(z0m = z0m, z0s = z0s)
end

"""
    charnock_parameter(mag_u_10)

Compute the Charnock parameter `α` as a function of 10-m wind speed `mag_u_10` [m/s] 
computed from neutral profile. Piecewise linear interpolation based on COARE 3.0 
(Fairall et al. 2003).
"""
@inline function charnock_parameter(mag_u_10, α_low, α_high, u_low, u_high)
    return ifelse(
        mag_u_10 <= u_low,
        α_low,
        ifelse(
            mag_u_10 >= u_high,
            α_high,
            α_low + (α_high - α_low) * (mag_u_10 - u_low) / (u_high - u_low)
        )
    )
end

# Accessors
@inline function momentum_roughness(spec::ConstantRoughnessParams, u★, sfc_param_set, roughness_inputs)
    return spec.z0m
end

@inline function scalar_roughness(spec::ConstantRoughnessParams, u★, sfc_param_set, roughness_inputs)
    return spec.z0s
end

@inline function momentum_roughness(spec::COARE3RoughnessParams, u★, sfc_param_set, roughness_inputs)
    FT = float_type(sfc_param_set)
    grav = SFP.grav(sfc_param_set)
    kinematic_visc = spec.kinematic_visc

    # Recover 10m wind speed using neutral profile.
    z0_proxy = spec.z0m_default
    κ = SFP.von_karman_const(sfc_param_set)
    mag_u_10 = (u★ / κ) * log(FT(10) / z0_proxy)
    
    α = charnock_parameter(mag_u_10, spec.α_low, spec.α_high, spec.u_low, spec.u_high)

    # Smooth flow limit (Smith 1988)
    u★_safe = max(u★, eps(FT))
    z0_smooth = FT(0.11) * kinematic_visc / u★_safe

    # Rough flow limit (Charnock 1955)
    z0_rough = α * u★_safe^2 / grav

    return z0_smooth + z0_rough
end

@inline function scalar_roughness(spec::COARE3RoughnessParams, u★, sfc_param_set, roughness_inputs)
    FT = float_type(sfc_param_set)
    kinematic_visc = spec.kinematic_visc
    
    # Compute z0m for Reynolds number calculation.
    # Note: momentum_roughness is stateless, so recomputing it here is safe.
    z0m = momentum_roughness(spec, u★, sfc_param_set, roughness_inputs)

    # Roughness Reynolds number (Re_star)
    # Ensure positive u★ to prevent complex numbers in power law.
    u★_safe = max(u★, eps(FT))
    Re_star = z0m * u★_safe / kinematic_visc

    # Empirical fit to COARE and HEXOS data (Fairall et al. 2003)
    # Limits z0s to 1.1e-4 m (smooth) and decreases for rough flow.
    return min(FT(1.1e-4), FT(5.5e-5) * Re_star^FT(-0.6))
end

"""
    momentum_roughness(spec::RaupachRoughnessParams, u★, sfc_param_set, ctx, roughness_inputs)

Calculate momentum roughness length using the Raupach (1994) canopy roughness model.
Requires `roughness_inputs` to contain `LAI` (Leaf Area Index) and `h` (canopy height).
"""
@inline function momentum_roughness(spec::RaupachRoughnessParams, u★, sfc_param_set, roughness_inputs)
    FT = float_type(sfc_param_set)
    k = SFP.von_karman_const(sfc_param_set)
    
    LAI = roughness_inputs.LAI
    h = roughness_inputs.h
    
    # Roughness density lambda (frontal area index) approximated as LAI / 2
    λ = LAI / FT(2)
    
    C_R = spec.C_R
    C_S = spec.C_S
    c_d1 = spec.c_d1
    
    # Avoid division by zero if LAI is extremely small
    if λ <= eps(FT)
        return SFP.z0m_fixed(sfc_param_set)
    end
    
    # d/h (Eq 15 in Raupach 1994)
    sqrt_c_lambda = sqrt(c_d1 * λ)
    d_over_h = FT(1) - (FT(1) - exp(-sqrt_c_lambda)) / sqrt_c_lambda
    
    # u_star / U(h) (Eq 8 in Raupach 1994)
    ustar_over_Uh = sqrt(C_S + C_R * λ)
    Uh_over_ustar = FT(1) / ustar_over_Uh
    
    # Psi_h (Roughness sublayer influence function), approximated as 0.193
    Ψ_h = FT(0.193)
    
    z0m = h * (FT(1) - d_over_h) * exp(-k * Uh_over_ustar - Ψ_h)
    
    return z0m
end

@inline function scalar_roughness(::RaupachRoughnessParams, u★, sfc_param_set, roughness_inputs)
    return SFP.z0s_fixed(sfc_param_set) # Fallback to fixed scalar roughness
end

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.MomentumTransport,
)
    return momentum_roughness(inputs.roughness_model, u★, sfc_param_set, inputs.roughness_inputs)
end

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.HeatTransport,
)
    return scalar_roughness(inputs.roughness_model, u★, sfc_param_set, inputs.roughness_inputs)
end
