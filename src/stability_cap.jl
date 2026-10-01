# Stability caps for stable stratification.
#
# In very stable conditions, the Monin-Obukhov similarity theory (MOST) exchange
# coefficients decrease rapidly with the stability parameter ``ζ``. At fixed wind
# speed, the MOST sensible heat flux ``H ∝ ζ / F̂_m(ζ)^3`` (at fixed `ΔU`, from
# `u_* = κ ΔU / F̂_m` and `H ∝ u_*^3 ζ / Δz_eff`) reaches a maximum at a stability
# ``ζ_p`` that depends only on the geometry (`Δz_eff / z0m`) and the universal
# functions (the "maximum sustainable heat flux"; Derbyshire 1999; van de Wiel
# et al. 2012). On the branch beyond this maximum, the heat flux *decreases* as
# the surface–air temperature difference increases. In models, this positive
# feedback leads to runaway surface cooling and decoupling from the atmosphere
# at low wind speeds, more than is observed: turbulent exchange persists,
# intermittently and driven by submeso motions outside MOST (Mahrt 2014).
#
# A stability cap limits the stability parameter entering the flux-profile
# relations to `ζ ≤ ζ_cap`. Beyond the cap, the exchange coefficients are held
# at their values at `ζ_cap`, so the heat flux keeps increasing with the
# surface–air temperature difference, and the bulk Richardson number
# `Ri_b(ζ) = ζ F̂_h(ζ_cap) / F̂_m(ζ_cap)^2` increases linearly with `ζ`, so the
# MOST solve always has a root. Because the momentum and scalar scales are
# consistent with the capped coefficients, the returned `ζ` (and `L_MO`) are the
# Obukhov stability parameter (length) implied by the computed fluxes. Unstable
# conditions are unaffected.

"""
    NoStabilityCap()

No cap on the stability parameter: standard MOST (the default).
"""
struct NoStabilityCap <: AbstractStabilityCap end

"""
    ConstantStabilityCap(ζ_max)

Cap on the stability parameter entering the flux-profile relations at the constant
`ζ_max > 0`. Physick & Garratt (1995) limit `z/L` to 0.5 in their mesoscale model; caps
between 0.5 and 2 are common in land models. The cap must be positive, so that unstable
conditions are unaffected.

# Fields
- `ζ_max`: Maximum stability parameter entering the flux-profile relations [-].

# Examples
```julia
config = SurfaceFluxConfig(roughness, gustiness, MoistModel(), NoRoughnessSubLayer(),
    ConstantStabilityCap(0.5))
```

# References
- Physick, W. L., & Garratt, J. R. (1995). Incorporation of a high-roughness lower boundary
    into a mesoscale model for studies of dry deposition over complex terrain.
    Boundary-Layer Meteorology, 74, 55–71.
"""
struct ConstantStabilityCap{FT <: Real} <: AbstractStabilityCap
    ζ_max::FT
    function ConstantStabilityCap{FT}(ζ_max) where {FT}
        ζ_max > 0 ||
            throw(ArgumentError("ConstantStabilityCap: ζ_max must be > 0 (got $ζ_max)"))
        return new{FT}(ζ_max)
    end
end
ConstantStabilityCap(ζ_max::FT) where {FT <: Real} = ConstantStabilityCap{FT}(ζ_max)

"""
    MaxHeatFluxStabilityCap()

Cap on the stability parameter at the stability ``ζ_p`` at which the MOST sensible heat
flux at fixed wind speed, ``H ∝ ζ / F̂_m(ζ)^3``, is maximal.

``ζ_p`` depends only on the universal functions, the solver scheme, the ratio of
the effective height to the momentum roughness length (`Δz_eff / z0m`), and the
roughness sublayer correction (if any); it has no free parameters. For point
values, it satisfies ``F_m(ζ_p) = 3 [φ_m(ζ_p) - φ_m(ζ_p z_{0m}/Δz_{eff})]``; for
log-linear stable functions ``φ_m = 1 + b ζ``, ``ζ_p ≈ \\ln(Δz_{eff}/z_{0m}) / (2 b)``.
With the Gryanik et al. (2020) functions (point values, no RSL), ``ζ_p ≈ 0.15–0.3``
for `Δz_eff / z0m = 2–10` (tall canopies, forcing height a few roughness lengths above
the displacement height), ``≈ 0.4–0.8`` for `Δz_eff / z0m = 30–1000` (grass),
``≈ 1–1.2`` for `Δz_eff / z0m = 3000–10⁴` (bare soil), and ``≈ 1.4–1.6`` for
`Δz_eff / z0m = 3·10⁴–10⁵` (snow).

Up to ``ζ_p``, the fluxes are those of MOST. Beyond it, the exchange coefficients are
held at their values at ``ζ_p``, so the heat flux increases with the surface–air
temperature difference; MOST alone predicts a decreasing heat flux there, which is
dynamically unstable (runaway cooling). At the cap, the local flux Richardson number
``R_f = ζ_p / φ_m(ζ_p)`` is ≈ 0.09–0.22 for the Gryanik functions (increasing with
`Δz_eff / z0m`), i.e., at or below the critical value ``R_{f,cr} ≈ 0.20–0.25`` beyond
which local similarity theory ceases to apply (Grachev et al. 2013); the bulk ratio
``ζ_p / F_m(ζ_p)`` is ≈ 0.07–0.14 (≈ ``1/(3b)`` for log-linear functions).

The cap is computed once per surface flux solve from the momentum roughness
length at neutral stability, by golden-section maximization of ``ζ / F̂_m(ζ)^3``
in ``\\log ζ`` over ``ζ ∈ [10^{-2}, 20]`` with a fixed number of iterations, refined
by one parabolic step (see [`max_heat_flux_stability`](@ref)). The maximization evaluates
the dimensionless momentum profile 15 times per solve; when `z0m` depends on the friction
velocity (e.g., [`COARE3RoughnessParams`](@ref)), a neutral solve for the roughness length
precedes it. Since ``ζ_p`` depends only on the geometry `Δz_eff / z0m`, the scheme, and
the RSL parameters, it is constant in time for fixed roughness lengths; for such surfaces,
[`ConstantStabilityCap`](@ref)`(`[`max_heat_flux_stability`](@ref)`(...))`
computed once per column gives the same fluxes, with the maximization done once.

# Examples
```julia
config = SurfaceFluxConfig(roughness, gustiness, MoistModel(), NoRoughnessSubLayer(),
    MaxHeatFluxStabilityCap())
```

# References
- Derbyshire, S. H. (1999). Boundary-layer decoupling over cold surfaces as a physical
    boundary-instability. Boundary-Layer Meteorology, 90, 297–325.
- van de Wiel, B. J. H., et al. (2012). The minimum wind speed for sustainable turbulence
    in the nocturnal boundary layer. Journal of the Atmospheric Sciences, 69, 3116–3127.
- Grachev, A. A., Andreas, E. L, Fairall, C. W., Guest, P. S., & Persson, P. O. G. (2013).
    The critical Richardson number and limits of applicability of local similarity theory
    in the stable boundary layer. Boundary-Layer Meteorology, 147, 51–82.
- Gryanik, V. M., Lüpkes, C., Grachev, A., & Sidorenko, D. (2020). New modified and
    extended stability functions for the stable boundary layer based on SHEBA and
    parametrizations of bulk transfer coefficients for climate models. Journal of the
    Atmospheric Sciences, 77, 2687–2716.
"""
struct MaxHeatFluxStabilityCap <: AbstractStabilityCap end

Base.broadcastable(c::AbstractStabilityCap) = tuple(c)

"""
    capped_stability(ζ, ζ_cap)
    capped_stability(inputs, ζ)
    capped_stability(param_set, inputs, scheme, ζ)

Return the stability parameter entering the flux-profile relations:
`min(ζ, ζ_cap)`, or `ζ` if there is no cap (`ζ_cap === nothing`). The second
method reads the cap from `inputs.ζ_cap`, which is set by the MOST solver (see
[`with_stability_cap`](@ref)); inputs built without a solve have `ζ_cap = nothing`. The
third method computes the cap from `inputs.stability_cap` when `inputs.ζ_cap` is
`nothing` (see [`resolved_stability_cap`](@ref)), so that the exchange coefficients and
similarity scales computed from builder inputs agree with those of the solve.
"""
@inline capped_stability(ζ::Number, ::Nothing) = ζ
@inline capped_stability(ζ::Number, ζ_cap) = min(ζ, ζ_cap)
@inline capped_stability(inputs::NamedTuple, ζ) =
    capped_stability(ζ, get(inputs, :ζ_cap, nothing))
@inline capped_stability(param_set::APS, inputs::NamedTuple, scheme, ζ) =
    capped_stability(ζ, resolved_stability_cap(param_set, inputs, scheme))

"""
    resolved_stability_cap(param_set, inputs, scheme)

Return the numerical value of the stability cap for `inputs`, or `nothing`:
`inputs.ζ_cap` if set (inside the MOST solver and the prescribed-flux paths, see
[`with_stability_cap`](@ref)), otherwise computed from `inputs.stability_cap` with
[`stability_cap_value`](@ref). [`NoStabilityCap`](@ref) and
[`ConstantStabilityCap`](@ref) return their value directly;
[`MaxHeatFluxStabilityCap`](@ref) runs the maximization in
[`max_heat_flux_stability`](@ref) (and, when `z0m` depends on the friction velocity, a
neutral solve for the roughness length) on each call, so callers that evaluate several
quantities from the same inputs should set the cap once with [`with_stability_cap`](@ref).
"""
@inline resolved_stability_cap(param_set, inputs, scheme) =
    _resolved_stability_cap(get(inputs, :ζ_cap, nothing), param_set, inputs, scheme)
@inline _resolved_stability_cap(ζ_cap, param_set, inputs, scheme) = ζ_cap
@inline _resolved_stability_cap(::Nothing, param_set, inputs, scheme) = stability_cap_value(
    get(inputs, :stability_cap, NoStabilityCap()),
    param_set,
    inputs,
    scheme,
)

"""
    neutral_momentum_roughness(roughness_model, param_set, inputs, scheme)

Return the momentum roughness length `z0m` [m] at neutral stability. For a roughness
model independent of the friction velocity (see [`depends_on_ustar`](@ref)), this is the
model's roughness length; for the others (e.g., [`COARE3RoughnessParams`](@ref)), it comes
from a neutral solve with [`compute_ustar_and_roughness`](@ref).
"""
@inline function neutral_momentum_roughness(roughness_model, param_set, inputs, scheme)
    FT = eltype(param_set)
    if depends_on_ustar(roughness_model)
        # The cap is irrelevant at ζ = 0 and is disabled here, so that `compute_ustar`
        # (which resolves the cap from the inputs) does not recurse into
        # `stability_cap_value`.
        inputs_uncapped = (; inputs..., stability_cap = NoStabilityCap(), ζ_cap = nothing)
        _, z0m, _ =
            compute_ustar_and_roughness(param_set, zero(FT), inputs_uncapped, scheme)
        return z0m
    else
        return momentum_roughness(
            roughness_model,
            zero(FT),
            param_set,
            inputs.roughness_inputs,
        )
    end
end

"""
    stability_cap_value(stability_cap, param_set, inputs, scheme[, z0m])

Return the numerical value of the cap on the stability parameter (or
`nothing` for [`NoStabilityCap`](@ref)) for the given inputs. For
[`MaxHeatFluxStabilityCap`](@ref), the momentum roughness length `z0m` is evaluated at
neutral stability with [`neutral_momentum_roughness`](@ref) unless it is given (e.g., from
a prescribed friction velocity).
"""
@inline stability_cap_value(::NoStabilityCap, param_set, inputs, scheme) = nothing
@inline stability_cap_value(c::ConstantStabilityCap, param_set, inputs, scheme) = c.ζ_max
@inline stability_cap_value(
    c::MaxHeatFluxStabilityCap,
    param_set,
    inputs,
    scheme,
) = stability_cap_value(
    c,
    param_set,
    inputs,
    scheme,
    neutral_momentum_roughness(inputs.roughness_model, param_set, inputs, scheme),
)
@inline stability_cap_value(c::AbstractStabilityCap, param_set, inputs, scheme, z0m) =
    stability_cap_value(c, param_set, inputs, scheme)
@inline stability_cap_value(::MaxHeatFluxStabilityCap, param_set, inputs, scheme, z0m) =
    max_heat_flux_stability(
        param_set,
        effective_height(inputs),
        z0m,
        scheme,
        inputs.rsl_model,
    )

"""
    with_stability_cap(inputs, param_set, scheme[, z0m])

Return `inputs` with the field `ζ_cap` set to the numerical value of the stability cap
from [`stability_cap_value`](@ref) (or `nothing` for [`NoStabilityCap`](@ref)). The MOST
solver and the prescribed-flux paths call this once per solve, so that the cap is computed
once and read by [`capped_stability`](@ref) thereafter.
"""
@inline with_stability_cap(inputs, param_set, scheme, args...) = (;
    inputs...,
    ζ_cap = stability_cap_value(
        get(inputs, :stability_cap, NoStabilityCap()),
        param_set,
        inputs,
        scheme,
        args...,
    ),
)

"""
    effective_obukhov_length(inputs, ζ, Δz_eff, L_MO)

Return the effective Obukhov length `Δz_eff / min(ζ, ζ_cap)` for profile recovery: the
length at which the (capped) exchange coefficients and similarity scales were evaluated.
It equals `L_MO` when `ζ ≤ ζ_cap` or there is no cap.
"""
@inline function effective_obukhov_length(inputs, ζ, Δz_eff, L_MO)
    ζ_c = capped_stability(inputs, ζ)
    return ifelse(ζ_c < ζ, Δz_eff / ζ_c, L_MO)
end

"""
    _HeatFluxObjective(uf_params, rsl_model, Δz_eff, z0m, scheme, ϵ)

Logarithm of the heat-flux objective ``ζ / F̂_m(ζ)^3`` as a function of ``u = \\ln ζ``,
``f(u) = u - 3 \\ln \\max(F̂_m(e^u), ϵ)``, for [`max_heat_flux_stability`](@ref). ``F̂_m``
is positive whenever the MOST profile is; the floor `ϵ` only guards degenerate
geometries. The functor stores the value `ϵ` rather than the type `FT`: a stored type
is a `DataType` field, which defeats type inference and GPU compilation.
"""
struct _HeatFluxObjective{UFP, M, H, Z, S, E}
    uf_params::UFP
    rsl_model::M
    Δz_eff::H
    z0m::Z
    scheme::S
    ϵ::E
end
@inline function (f::_HeatFluxObjective)(u)
    F̂_m = rsl_corrected_profile(
        f.uf_params,
        f.rsl_model,
        f.Δz_eff,
        exp(u),
        f.z0m,
        UF.MomentumTransport(),
        f.scheme,
    )
    return u - 3 * log(max(F̂_m, f.ϵ))
end

"""
    max_heat_flux_stability(param_set, Δz_eff, z0m, scheme = PointValueScheme(), rsl_model = NoRoughnessSubLayer())

Return the stability parameter ``ζ_p`` at which ``ζ / F̂_m(ζ)^3`` (the MOST sensible heat
flux at fixed wind speed, up to a factor) is maximal, where `F̂_m = F_m + P_m` is the
roughness-sublayer-corrected dimensionless momentum profile
([`rsl_corrected_profile`](@ref)).

The maximum is located by golden-section search in ``\\ln ζ`` over ``[10^{-2}, 20]``
with a fixed number of iterations (no data-dependent control flow), followed by one
parabolic (Newton) step. The step refines the value to ``\\approx 10^{-4}`` relative
accuracy and makes ``ζ_p`` a smooth function of `Δz_eff`, `z0m`, and the RSL parameters
with the correct (implicit-function) derivatives under automatic differentiation. The
result lies within the search interval; if the maximum lies outside it, the nearest end
point is returned.

# Arguments
- `param_set`: Parameter set.
- `Δz_eff`: Effective height `Δz - d` above the displacement height [m].
- `z0m`: Momentum roughness length [m].
- `scheme`: [`PointValueScheme`](@ref) (default) or [`LayerAverageScheme`](@ref).
- `rsl_model`: Roughness sublayer model (default [`NoRoughnessSubLayer`](@ref)).

# Returns
The stability parameter ``ζ_p = Δz_eff / L`` of maximum heat flux [-].

See also [`MaxHeatFluxStabilityCap`](@ref).
"""
@inline function max_heat_flux_stability(
    param_set::APS,
    Δz_eff,
    z0m,
    scheme = UF.PointValueScheme(),
    rsl_model = NoRoughnessSubLayer(),
)
    FT = eltype(param_set)
    uf = SFP.uf_params(param_set)
    f = _HeatFluxObjective(uf, rsl_model, Δz_eff, z0m, scheme, eps(FT))
    g = FT(0.6180339887498949) # (√5 - 1)/2
    a = log(FT(1e-2))
    b = log(FT(20))
    a_0, b_0 = a, b
    c = b - g * (b - a)
    d = a + g * (b - a)
    fc = f(c)
    fd = f(d)
    # The width shrinks by g per iteration: 7.6 * g^10 ≈ 0.06 in ln ζ, so the midpoint
    # is within 0.03 of the maximum, which the parabolic step below then locates to
    # ~1e-4 relative accuracy (checked against brute-force maximization for the
    # Businger, Gryanik, and Grachev functions, both schemes, and RSL corrections).
    for _ in 1:10
        left = fc > fd # maximum in [a, d]
        a, b = ifelse(left, a, c), ifelse(left, d, b)
        c_new = ifelse(left, b - g * (b - a), d)
        d_new = ifelse(left, c, a + g * (b - a))
        fnew = f(ifelse(left, c_new, d_new))
        fc, fd = ifelse(left, fnew, fd), ifelse(left, fc, fnew)
        c, d = c_new, d_new
    end
    # The golden-section result is piecewise constant in the parameters (Δz_eff, z0m,
    # RSL parameters): its comparisons use values only, so its derivatives vanish
    # under AD and finite differences alike, although ζ_p varies smoothly with the
    # parameters. One parabolic (Newton) step on ∂f/∂u = 0, with central differences
    # of step h in u, makes the result smooth. Its derivatives then approximate the
    # implicit-function derivative -∂²f/∂u∂p / ∂²f/∂u² with O(h²) error, and the step
    # refines the value. The step is taken only if f is concave there and the vertex
    # lies within h, and the result is kept within the search interval, so a maximum
    # at an end point of the interval returns that end point.
    u_m = (a + b) / 2
    h = FT(0.05)
    f_0 = f(u_m)
    f_p = f(u_m + h)
    f_m = f(u_m - h)
    curv = f_p - 2 * f_0 + f_m  # ≈ h² ∂²f/∂u² (< 0 at a maximum)
    concave = curv < zero(curv)
    δ = -h * (f_p - f_m) / (2 * ifelse(concave, curv, -one(curv)))
    use_step = concave & (abs(δ) <= h)
    return exp(clamp(u_m + ifelse(use_step, δ, zero(δ)), a_0, b_0))
end
