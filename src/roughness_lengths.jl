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
    ConstantRoughnessSpec{FT} <: AbstractRoughnessSpec

Roughness lengths fixed to constant values.

# Fields
- `z0m`: Momentum roughness length [m]
- `z0s`: Scalar roughness length [m]
"""
Base.@kwdef struct ConstantRoughnessSpec{FT} <: AbstractRoughnessSpec
    z0m::FT = 1e-3
    z0s::FT = 1e-4
end

"""
    COARE3RoughnessSpec{FT} <: AbstractRoughnessSpec

COARE 3.0 roughness parameterization (Fairall et al. 2003).

# Fields
- `kinematic_visc`: Kinematic viscosity of air [m^2/s]
"""
Base.@kwdef struct COARE3RoughnessSpec{FT} <: AbstractRoughnessSpec
    kinematic_visc::FT = 1.5e-5
    z0m_default::FT = 1e-4
    α_low::FT = 0.011
    α_high::FT = 0.018
end

"""
    RaupachRoughnessSpec <: AbstractRoughnessSpec

Raupach (1994) canopy roughness model.
"""
Base.@kwdef struct RaupachRoughnessSpec{FT} <: AbstractRoughnessSpec
    C_R::FT = 0.3
    C_S::FT = 0.003
    c_d1::FT = 7.5
end

SurfaceFluxConfig(; roughness = ConstantRoughnessSpec(), gustiness = ConstantGustinessSpec(1.0)) =
    SurfaceFluxConfig(roughness, gustiness)

"""
    roughness_lengths(z0m; scalar)

Helper to construct `ConstantRoughnessSpec`.
"""
function roughness_lengths(z0m; z0s)
    return ConstantRoughnessSpec(z0m = z0m, z0s = scalar)
end

"""
    charnock_parameter(mag_u_10)

Compute the Charnock parameter `α` as a function of 10m wind speed `mag_u_10` [m/s].
Piecewise linear interpolation based on COARE 3.0 (Fairall et al. 2003).
"""
@inline function charnock_parameter(mag_u_10, α_low, α_high)
    FT = eltype(mag_u_10)
    u_low = FT(10)
    u_high = FT(18)

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
@inline function momentum_roughness(spec::ConstantRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    return spec.z0m
end

@inline function scalar_roughness(spec::ConstantRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    return spec.z0s
end

@inline function momentum_roughness(spec::COARE3RoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    FT = float_type(sfc_param_set)
    grav = SFP.grav(sfc_param_set)
    kinematic_visc = spec.kinematic_visc

    # Recover 10m wind speed using profile recovery.
    # Use a constant proxy z0 for the restoration since z0m is not yet known.
    L_MO = ctx.L_MO
    z0_proxy = spec.z0m_default
    
    mag_u_10 = compute_profile_value(
        sfc_param_set,
        L_MO,
        z0_proxy,
        FT(10), # 10m height
        u★,
        FT(0), # Surface velocity
        UF.MomentumTransport(),
    )
    
    α = charnock_parameter(mag_u_10, spec.α_low, spec.α_high)

    # Smooth flow limit (Smith 1988)
    u★_safe = max(u★, FT(1e-9))
    z0_smooth = FT(0.11) * kinematic_visc / u★_safe

    # Rough flow limit (Charnock 1955)
    z0_rough = α * u★_safe^2 / grav

    return z0_smooth + z0_rough
end

@inline function scalar_roughness(spec::COARE3RoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    FT = float_type(sfc_param_set)
    kinematic_visc = spec.kinematic_visc
    
    # Compute z0m for Reynolds number calculation.
    # Note: momentum_roughness is stateless, so recomputing it here is safe.
    z0m = momentum_roughness(spec, u★, sfc_param_set, ctx, roughness_inputs)

    # Roughness Reynolds number (Re_star)
    # Ensure positive u★ to prevent complex numbers in power law.
    u★_safe = max(u★, FT(1e-9))
    Re_star = z0m * u★_safe / kinematic_visc

    # Empirical fit to COARE and HEXOS data (Fairall et al. 2003)
    # Limits z0s to 1.1e-4 m (smooth) and decreases for rough flow.
    return min(FT(1.1e-4), FT(5.5e-5) * Re_star^FT(-0.6))
end

"""
    momentum_roughness(spec::RaupachRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)

Calculate momentum roughness length using the Raupach (1994) canopy roughness model.
Requires `roughness_inputs` to contain `LAI` (Leaf Area Index) and `h` (canopy height).
"""
@inline function momentum_roughness(spec::RaupachRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
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

@inline function scalar_roughness(::RaupachRoughnessSpec, u★, sfc_param_set, ctx, roughness_inputs)
    return SFP.z0s_fixed(sfc_param_set) # Fallback to fixed scalar roughness
end

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.MomentumTransport,
    ctx,
)
    return momentum_roughness(inputs.roughness_model, u★, sfc_param_set, ctx, inputs.roughness_inputs)
end

@inline function compute_z0(
    u★,
    sfc_param_set,
    inputs::SurfaceFluxInputs,
    ::UF.HeatTransport,
    ctx,
)
    return scalar_roughness(inputs.roughness_model, u★, sfc_param_set, ctx, inputs.roughness_inputs)
end
