# Roughness-length evaluation
#
# These helpers resolve the momentum and scalar roughness length specifications.
# They are intentionally lightweight so they can be inlined inside GPU kernels.

"""
    ConstantRoughnessParams{FT} <: AbstractRoughnessParams

Roughness lengths fixed to constant values.

# Fields
- `z0m`: Momentum roughness length [m]
- `z0s`: Scalar roughness length [m]. Used for both heat (`z0h`) and humidity (`z0q`).

The default values specified here are used when constructing the struct manually. When loading
via `ClimaParams`, these values are overwritten by the parameters in the `ClimaParams` TOML file.
"""
Base.@kwdef struct ConstantRoughnessParams{FT} <: AbstractRoughnessParams
    z0m::FT = 2e-4
    z0s::FT = 2e-5
end

"""
    COARE3RoughnessParams{FT} <: AbstractRoughnessParams

COARE 3.0 roughness parameterization (Fairall et al. 2003).

# References
- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). 
    Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm.
    Journal of Climate, 16, 571–591.
    [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)

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

# References
- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane displacement 
    as functions of canopy height and area index.
    Boundary-Layer Meteorology, 71, 211–216.
    [DOI: 10.1007/BF00709229](https://doi.org/10.1007/BF00709229)

The default values specified here are used when constructing the struct manually. When loading
via `ClimaParams`, `stanton_number` is overwritten by the parameter in the TOML file, while
other parameters (`C_R`, `C_S`, `c_d1`) retain their default values.
"""
Base.@kwdef struct RaupachRoughnessParams{FT} <: AbstractRoughnessParams
    C_R::FT = 0.3
    C_S::FT = 0.003
    c_d1::FT = 7.5
    stanton_number::FT = 0.1
end

"""
    charnock_parameter(mag_u_10, α_low, α_high, u_low, u_high)

Compute the Charnock parameter `α` as a function of the 10-m wind speed `mag_u_10` [m/s] 
estimated from the neutral profile. Piecewise linear interpolation between the lower and 
upper bounds based on COARE 3.0 (Fairall et al. 2003).

# References
- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). 
    Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm.
    Journal of Climate, 16, 571–591.
    [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)
"""
@inline function charnock_parameter(mag_u_10, α_low, α_high, u_low, u_high)
    return ifelse(
        mag_u_10 <= u_low,
        α_low,
        ifelse(
            mag_u_10 >= u_high,
            α_high,
            α_low + (α_high - α_low) * (mag_u_10 - u_low) / (u_high - u_low),
        ),
    )
end

# Accessors
@inline function momentum_roughness(
    spec::ConstantRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    return spec.z0m
end

@inline function scalar_roughness(
    spec::ConstantRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    return spec.z0s
end

@inline function momentum_and_scalar_roughness(
    spec::ConstantRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    return (spec.z0m, spec.z0s)
end

"""
    momentum_roughness(spec::COARE3RoughnessParams, u★, sfc_param_set, roughness_inputs)

Calculate momentum roughness length using the COARE 3.0 algorithm (Fairall et al. 2003).

# Formulation
The momentum roughness length `z0m` is parameterized as the sum of a smooth flow limit
(Smith 1988) and a rough flow limit (Charnock 1955):
```math
z_{0m} = z_{0m,smooth} + z_{0m,rough}
```
- **Smooth flow**: Dominated by viscous limit, proportional to `ν / u★`.
- **Rough flow**: Dominated by wind stress, proportional to `α * u★^2 / g`.

The Charnock parameter `α` varies with the 10-m wind speed and is interpolated linearly between
lower and upper bounds defined in `spec`.

# Dependencies
- `u★`: Friction velocity [m/s]
- `kinematic_visc`: Kinematic viscosity of air [m^2/s], from `spec`
- `grav`: Gravitational acceleration [m/s^2], from `sfc_param_set`
- `mag_u_10`: 10m wind speed [m/s]. (internally recovered from u★ assuming neutral log profile)

# References
- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). 
    Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm.
    Journal of Climate, 16, 571–591.
    [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)
- Smith, S. D. (1988). Coefficients for sea surface wind stress, heat flux, and wind profiles 
    as a function of wind speed and temperature.
    Journal of Geophysical Research: Oceans, 93, 15467–15472.
    [DOI: 10.1029/JC093iC12p15467](https://doi.org/10.1029/JC093iC12p15467)
- Charnock, H. (1955). Wind stress on a water surface.
    Quarterly Journal of the Royal Meteorological Society, 81, 639–640.
    [DOI: 10.1002/qj.49708135027](https://doi.org/10.1002/qj.49708135027)
"""
@inline function momentum_roughness(
    spec::COARE3RoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    FT = eltype(sfc_param_set)
    grav = SFP.grav(sfc_param_set)
    kinematic_visc = spec.kinematic_visc

    # Recover 10-m wind speed using neutral profile with a proxy roughness length (to avoid 
    # circular dependency)
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

"""
    scalar_roughness(spec::COARE3RoughnessParams, u★, sfc_param_set, roughness_inputs)

Calculate scalar roughness length using the COARE 3.0 algorithm (Fairall et al. 2003).

# Formulation
The scalar roughness length `z0s` is parameterized as an empirical fit to COARE and HEXOS data.
It limits `z0s` to a smooth flow limit (`1.1e-4` m) and decreases for rough flow following a 
power law of the roughness Reynolds number (`Re_star`):
```math
z_{0s} = \\min(1.1 \\times 10^{-4}, 5.5 \\times 10^{-5} \\cdot R_{e*}^{-0.6})
```
where ``R_{e*} = z_{0m} u_* / \\nu``.

# Dependencies
- `u★`: Friction velocity [m/s]
- `kinematic_visc`: Kinematic viscosity of air [m^2/s]
- `z0m`: Momentum roughness length (calculated internally)

# References
- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). 
    Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm.
    Journal of Climate, 16, 571–591.
    [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)
"""
@inline function scalar_roughness(
    spec::COARE3RoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    # Forward to combined calculation to avoid code duplication
    _, z0s = momentum_and_scalar_roughness(spec, u★, sfc_param_set, roughness_inputs)
    return z0s
end

@inline function momentum_and_scalar_roughness(
    spec::COARE3RoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    FT = eltype(sfc_param_set)
    z0m = momentum_roughness(spec, u★, sfc_param_set, roughness_inputs)
    kinematic_visc = spec.kinematic_visc
    u★_safe = max(u★, eps(FT))
    Re_star = z0m * u★_safe / kinematic_visc
    z0s = min(FT(1.1e-4), FT(5.5e-5) * Re_star^FT(-0.6))
    return (z0m, z0s)
end


"""
    momentum_roughness(spec::RaupachRoughnessParams, u★, sfc_param_set, roughness_inputs)

Calculate momentum roughness length using the Raupach (1994) canopy roughness model.

# Formulation
This model estimates the aerodynamic roughness length `z0m` based on the geometry of the
roughness elements (canopy). It accounts for the division of surface drag between the
substrate (soil) and the roughness elements (plants).

Key features:
- **Displacement height (`d`)**: The height at which the mean drag appears to act.
- **Roughness density (`λ`)**: Characterized by the Frontal Area Index (FAI), approximated here as `LAI / 2`.
- **Wind attenuation**: Estimates `u★ / U(h)` (friction velocity ratio at canopy top).

# Dependencies
- `roughness_inputs.LAI`: Leaf Area Index [m^2/m^2]. Used to approximate frontal area index `λ`.
- `roughness_inputs.h`: Canopy height [m].
- `spec.C_R`, `spec.C_S`, `spec.c_d1`: Empirical model coefficients (defaults: 0.3, 0.003, 7.5).

# References
- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane 
    displacement as functions of canopy height and area index.
    Boundary-Layer Meteorology, 71, 211–216.
    [DOI: 10.1007/BF00709229](https://doi.org/10.1007/BF00709229)
"""
@inline function momentum_roughness(
    spec::RaupachRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    FT = eltype(sfc_param_set)
    k = SFP.von_karman_const(sfc_param_set)

    LAI = roughness_inputs.LAI
    h = roughness_inputs.h

    # Roughness density lambda (frontal area index) approximated as LAI / 2
    λ = LAI / FT(2)

    C_R = spec.C_R
    C_S = spec.C_S
    c_d1 = spec.c_d1

    # Avoid division by zero if LAI is extremely small
    # Small λ limit: Use fixed roughness
    z0m_small = SFP.z0m_fixed(sfc_param_set)

    # Large λ: Use Raupach model
    # Avoid div by zero by clamping λ
    λ_safe = max(λ, eps(FT))

    # d/h (Eq 15 in Raupach 1994)
    sqrt_c_lambda = sqrt(c_d1 * λ_safe)
    d_over_h = FT(1) - (FT(1) - exp(-sqrt_c_lambda)) / sqrt_c_lambda

    # u_star / U(h) (Eq 8 in Raupach 1994)
    ustar_over_Uh = sqrt(C_S + C_R * λ)
    Uh_over_ustar = FT(1) / ustar_over_Uh

    # Psi_h (roughness sublayer influence function), approximated by fixed value 0.193
    Ψ_h = FT(0.193)

    z0m_large = h * (FT(1) - d_over_h) * exp(-k * Uh_over_ustar - Ψ_h)

    return ifelse(λ <= eps(FT), z0m_small, z0m_large)
end

"""
    scalar_roughness(spec::RaupachRoughnessParams, u★, sfc_param_set, roughness_inputs)

Calculate scalar roughness length from the momentum roughness length using a fixed Stanton number.

# Formulation
Scaled from the momentum roughness length using a fixed Stanton number:
```math
z_{0s} = z_{0m} \\cdot St
```
where ``St`` is `spec.stanton_number`.

# Dependencies
- `spec.stanton_number`: Stanton number specific to the canopy/surface type.
"""
@inline function scalar_roughness(
    spec::RaupachRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    z0m = momentum_roughness(spec, u★, sfc_param_set, roughness_inputs)
    return z0m * spec.stanton_number
end

@inline function momentum_and_scalar_roughness(
    spec::RaupachRoughnessParams,
    u★,
    sfc_param_set,
    roughness_inputs,
)
    z0m = momentum_roughness(spec, u★, sfc_param_set, roughness_inputs)
    return (z0m, z0m * spec.stanton_number)
end

# =========================================================================================
# Iterative Solver for combined ustar and roughness
# =========================================================================================

struct UstarResidual{FT, PS, I, S} <: Function
    param_set::PS
    inputs::I
    scheme::S
    ζ::FT
end

function (ur::UstarResidual)(ustar)
    param_set = ur.param_set
    FT = eltype(param_set)
    inputs = ur.inputs
    scheme = ur.scheme
    ζ = ur.ζ

    # Ensure ustar is positive for physical consistency
    ustar_safe = max(ustar, FT(0))
    z0m, z0s = momentum_and_scalar_roughness(
        inputs.roughness_model,
        ustar_safe,
        param_set,
        inputs.roughness_inputs,
    )

    b_flux = buoyancy_flux(param_set, ζ, ustar_safe, inputs)
    gustiness_val = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    ustar_calc = compute_ustar(param_set, ζ, z0m, inputs, scheme, gustiness_val)

    return ustar - ustar_calc
end

"""
    compute_ustar_and_roughness(param_set, ζ, inputs, scheme)

Computes friction velocity `ustar` and roughness lengths `z0m`, `z0h` for a given stability `ζ`.

- If `inputs.ustar` is prescribed, it is returned directly.
- Uses a root solver to find `ustar` consistent with the roughness model.
- Converges in 1 iteration if roughness and gustiness are fixed (independent of `ustar`).
"""
function compute_ustar_and_roughness(
    param_set::APS,
    ζ,
    inputs,
    scheme,
)
    FT = eltype(param_set)

    if inputs.ustar !== nothing
        ustar = inputs.ustar
        z0m, z0s = momentum_and_scalar_roughness(
            inputs.roughness_model,
            ustar,
            param_set,
            inputs.roughness_inputs,
        )
        return ustar, z0m, z0s
    end

    rf = UstarResidual(param_set, inputs, scheme, ζ)

    ustar_min = FT(1e-4) # Sufficient for very calm conditions
    ustar_max = FT(4.0)  # Sufficient for hurricane-force winds
    maxiter = 4          # Relatively small maxiter for performance
    rtol = FT(0.02) # Relative tolerance on successive iterates of ustar

    sol = RS.find_zero(
        rf,
        RS.BrentsMethod(ustar_min, ustar_max),
        RS.CompactSolution(),
        RS.RelativeSolutionTolerance(rtol),
        maxiter,
    )
    ustar = sol.root

    z0m, z0s = momentum_and_scalar_roughness(
        inputs.roughness_model,
        ustar,
        param_set,
        inputs.roughness_inputs,
    )

    return ustar, z0m, z0s
end
