import SurfaceFluxes
SurfaceFluxes.error_on_non_convergence() = false
using Random
const rseed = MersenneTwister(0)
const SF = SurfaceFluxes
import SurfaceFluxes.UniversalFunctions as UF
using Statistics
using StaticArrays
using Thermodynamics
using Thermodynamics.TemperatureProfiles
using Thermodynamics.TestedProfiles
Thermodynamics.print_warning() = false
using UnPack

using RootSolvers
const RS = RootSolvers

using CLIMAParameters
using CLIMAParameters: AbstractParameterSet
using CLIMAParameters.Planet
const CPP = CLIMAParameters.Planet
const CP = CLIMAParameters
const SFP = SF.Parameters
const APS = CP.AbstractParameterSet
const TP = Thermodynamics.TemperatureProfiles
import Thermodynamics.TestedProfiles: input_config, PhaseDryProfiles, PhaseEquilProfiles, shared_profiles

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
const FT = Float32;
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
param_set = create_parameters(toml_dict, UF.BusingerType())
thermo_params = SFP.thermodynamics_params(param_set)
uft = SFP.universal_func_type(param_set)

"""
    PhaseEquilProfiles(param_set, ::Type{ArrayType})

Returns a `ProfileSet` used to test moist states in thermodynamic equilibrium.
"""
function PhaseEquilProfiles(param_set::APS, ::Type{ArrayType}) where {ArrayType}
    phase_type = TD.PhaseEquil

    # Prescribe z_range, relative_sat, T_surface, T_min
    z_range, relative_sat, T_surface, T_min = input_config(ArrayType)

    # Compute T, p, from DecayingTemperatureProfile, (reshape RS)
    z, T_virt, p, RS = shared_profiles(param_set, z_range, relative_sat, T_surface, T_min)
    T = T_virt

    FT = eltype(T)
    R_d::FT = CPP.R_d(param_set)
    grav::FT = CPP.grav(param_set)
    # Compute total specific humidity from temperature, pressure
    # and relative saturation, and partition the saturation excess
    # according to temperature.
    ρ = p ./ (R_d .* T)
    q_tot = RS .* TD.q_vap_saturation.(Ref(param_set), T, ρ, Ref(phase_type))
    q_pt = TD.PhasePartition_equil.(Ref(param_set), T, ρ, q_tot, Ref(phase_type))

    # Extract phase partitioning and update pressure
    # to be thermodynamically consistent with T, ρ, q_pt
    q_liq = getproperty.(q_pt, :liq)
    q_ice = getproperty.(q_pt, :ice)
    p = TD.air_pressure.(Ref(param_set), T, ρ, q_pt)

    e_int = TD.internal_energy.(Ref(param_set), T, q_pt)
    θ_liq_ice = TD.liquid_ice_pottemp.(Ref(param_set), T, ρ, q_pt)
    RH = TD.relative_humidity.(Ref(param_set), T, p, Ref(phase_type), q_pt)
    e_pot = grav * z
    Random.seed!(15)
    u = rand(FT, size(T)) * 5
    v = rand(FT, size(T)) * 5
    w = rand(FT, size(T)) * 5
    e_kin = (u .^ 2 .+ v .^ 2 .+ w .^ 2) / 2

    return ProfileSet{typeof(T), typeof(q_pt), typeof(phase_type)}(
        z,
        T,
        p,
        RS,
        e_int,
        ρ,
        θ_liq_ice,
        q_tot,
        q_liq,
        q_ice,
        q_pt,
        RH,
        e_pot,
        u,
        v,
        w,
        e_kin,
        phase_type,
    )
end
function input_config(ArrayType; n = 10, n_RS1 = 20, n_RS2 = 20, T_surface = 290, T_min = 150)
    n_RS = n_RS1 + n_RS2
    z_range = ArrayType(range(0, stop = 5e2, length = n))
    relative_sat1 = ArrayType(range(0, stop = 1, length = n_RS1))
    relative_sat2 = ArrayType(range(1, stop = 1.02, length = n_RS2))
    relative_sat = vcat(relative_sat1, relative_sat2)
    return z_range, relative_sat, T_surface, T_min
end

"""
    shared_profiles(
        param_set::APS,
        z_range::AbstractArray,
        relative_sat::AbstractArray,
        T_surface,
        T_min,
    )

Compute profiles shared across `PhaseDry`,
`PhaseEquil` and `PhaseNonEquil` thermodynamic
states, including:
 - `z` altitude
 - `T` temperature
 - `p` pressure
 - `RS` relative saturation
"""
function shared_profiles(param_set::APS, z_range::AbstractArray, relative_sat::AbstractArray, T_surface, T_min)
    FT = eltype(z_range)
    n_RS = length(relative_sat)
    n = length(z_range)
    T = similar(z_range, n * n_RS)
    p = similar(z_range, n * n_RS)
    RS = similar(z_range, n * n_RS)
    z = similar(z_range, n * n_RS)
    linear_indices = LinearIndices((1:n, 1:n_RS))
    # We take the virtual temperature, returned here,
    # as the temperature, and then compute a thermodynamic
    # state consistent with that temperature. This profile
    # will not be in hydrostatic balance, but this does not
    # matter for the thermodynamic test profiles.
    profile = TP.DecayingTemperatureProfile{FT}(param_set, FT(T_surface), FT(T_min))
    for i in linear_indices.indices[1]
        for j in linear_indices.indices[2]
            k = linear_indices[i, j]
            z[k] = z_range[i]
            T[k], p[k] = profile(param_set, z[k])
            RS[k] = relative_sat[j]
        end
    end
    return z, T, p, RS
end

profiles_sfc = []
profiles_int = []
profiles = PhaseEquilProfiles(thermo_params, Array{FT});

## Properties contained in `profiles`
## :z, :T, :p, :RS, :e_int, :h, :ρ, :θ_liq_ice, :q_tot, :q_liq, :q_ice, :q_pt, :RH, :e_pot, :u, :v, :w, :e_kin, :phase_type

for prof in profiles
    if prof.z == FT(0)
        push!(profiles_sfc, prof)
    else
        push!(profiles_int, prof)
    end
end

z0_momentum = Array{FT}(range(1e-6, stop = 1e-2, length = 20))
z0_thermal = Array{FT}(range(1e-6, stop = 1e-2, length = 20))

function check_over_moist_states()
    @inbounds for (ii, pint) in enumerate(profiles_int)
        @inbounds for (jj, psfc) in enumerate(profiles_sfc)
            @inbounds for (kk, z0m) in enumerate(z0_momentum)
                @inbounds for (ll, z0b) in enumerate(z0_thermal)
                    Tsfc = psfc.T + rand()
                    Tint = pint.T + 2rand() * rand([-1, 1])
                    ts_int_test = Thermodynamics.PhaseEquil{FT}(pint.ρ, pint.p, pint.e_int, pint.q_tot, Tint)
                    ts_sfc_test = Thermodynamics.PhaseEquil{FT}(psfc.ρ, psfc.p, psfc.e_int, psfc.q_tot, Tsfc)
                    state_in = SF.InteriorValues(pint.z, (pint.u / 10, pint.v / 10), ts_int_test)
                    state_sfc = SF.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test)
                    sc = SF.ValuesOnly{FT}(; state_in, state_sfc, z0m, z0b)
                    sfc_output = SF.surface_conditions(param_set, sc; maxiter = 10, soltype = RS.VerboseSolution())
                    @test sfc_output.is_converged == true
                end
            end
        end
    end
end

function check_over_dry_states()
    nstates = length(profiles_int) * length(profiles_sfc)
    @inbounds for (ii, pint) in enumerate(profiles_int)
        @inbounds for (jj, psfc) in enumerate(profiles_sfc)
            @inbounds for (kk, z0m) in enumerate(z0_momentum)
                @inbounds for (ll, z0b) in enumerate(z0_thermal)
                    Tsfc = psfc.T + rand()
                    Tint = pint.T + 2rand() .* rand([-1, 1])
                    ts_int_test = Thermodynamics.PhaseDry{FT}(pint.e_int, pint.ρ)
                    ts_sfc_test = Thermodynamics.PhaseDry{FT}(psfc.e_int, psfc.ρ)
                    state_in = SF.InteriorValues(pint.z, (pint.u / 10, pint.v / 10), ts_int_test)
                    state_sfc = SF.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test)
                    sc = SF.ValuesOnly{FT}(; state_in, state_sfc, z0m, z0b)
                    sfc_output = SF.surface_conditions(param_set, sc; maxiter = 10)
                    @test sfc_output.is_converged == true
                end
            end
        end
    end
end

@testset "Check convergence (dry thermodynamic states)" begin
    check_over_dry_states()
end
@testset "Check convergence (moist thermodynamic states)" begin
    check_over_moist_states()
end
