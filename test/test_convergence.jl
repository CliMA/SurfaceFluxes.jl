import SurfaceFluxes
SurfaceFluxes.error_on_non_convergence() = true
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
import Thermodynamics.TestedProfiles: input_config, PhaseEquilProfiles
include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))

function input_config(ArrayType; n = 4, n_RS1 = 4, n_RS2 = 4, T_surface = 290, T_min = 150)
    n_RS = n_RS1 + n_RS2
    z_range = ArrayType(range(0, stop = 5e2, length = n))
    relative_sat1 = ArrayType(range(0, stop = 1, length = n_RS1))
    relative_sat2 = ArrayType(range(1, stop = 1.02, length = n_RS2))
    relative_sat = vcat(relative_sat1, relative_sat2)
    return z_range, relative_sat, T_surface, T_min
end

function generate_profiles(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    param_set = create_parameters(toml_dict, UF.BusingerType())
    thermo_params = SFP.thermodynamics_params(param_set)
    uft = SFP.universal_func_type(param_set)
    profiles = collect(Thermodynamics.TestedProfiles.PhaseEquilProfiles(thermo_params, Array{FT}))
    profiles_sfc = filter(p -> iszero(p.z), profiles)
    profiles_int = filter(p -> !iszero(p.z), profiles)
    ## Properties contained in `profiles_<sfc, int>`
    ## :z, :T, :p, :RS, :e_int, :h, :ρ, 
    ## :θ_liq_ice, :q_tot, :q_liq, :q_ice, :q_pt, :RH, 
    ## :e_pot, :u, :v, :w, :e_kin, :phase_type
    return profiles_sfc, profiles_int, param_set
end

function check_over_dry_states(
    param_set,
    ::Type{FT},
    profiles_int,
    profiles_sfc,
    scheme,
    z0_momentum,
    z0_thermal,
    maxiter,
) where {FT}
    counter = [0, 0, 0] # St, Unst, Neutral
    @inbounds for (ii, pint) in enumerate(profiles_int)
        @inbounds for (jj, psfc) in enumerate(profiles_sfc)
            @inbounds for (kk, z0m) in enumerate(z0_momentum)
                @inbounds for (ll, z0b) in enumerate(z0_thermal)
                    @inbounds for sch in scheme
                        if z0m / z0b >= FT(1e2) || z0b / z0m >= FT(1e2)
                            nothing
                        else
                            thermo_params = SFP.thermodynamics_params(param_set)
                            ts_sfc_test = Thermodynamics.PhaseDry{FT}(psfc.e_int, psfc.ρ)
                            ts_int_test = Thermodynamics.PhaseDry{FT}(pint.e_int, pint.ρ)
                            state_in = SF.InteriorValues(pint.z, (pint.u / 10, pint.v / 10), ts_int_test)
                            state_sfc = SF.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test)
                            sc = SF.ValuesOnly{FT}(; state_in, state_sfc, z0m, z0b)
                            DSEᵥ_int = TD.virtual_dry_static_energy(thermo_params, ts_int_test, pint.e_pot)
                            DSEᵥ_sfc = TD.virtual_dry_static_energy(thermo_params, ts_sfc_test, psfc.e_pot)
                            ΔDSEᵥ = DSEᵥ_int - DSEᵥ_sfc
                            if abs(ΔDSEᵥ) <= FT(SFP.cp_d(param_set) / 10)
                                counter[3] += 1
                            else
                                sign(ΔDSEᵥ) == 1 ? counter[1] += 1 : counter[2] += 1
                            end
                            @test try
                                sfcc =
                                    SF.surface_conditions(param_set, sc, sch; maxiter, soltype = RS.VerboseSolution())
                                true
                            catch
                                false
                            end
                        end
                    end
                end
            end
        end
    end
    @info "(Dry $FT) Stable BL: $(counter[1]), Unstable BL: $(counter[2]), Neutral BL: $(counter[3])"
end

function check_over_moist_states(
    param_set,
    FT::DataType,
    profiles_int,
    profiles_sfc,
    scheme,
    z0_momentum,
    z0_thermal,
    maxiter,
)
    counter = [0, 0, 0] # St, Unst, Neutral
    @inbounds for (ii, pint) in enumerate(profiles_int)
        @inbounds for (jj, psfc) in enumerate(profiles_sfc)
            @inbounds for (kk, z0m) in enumerate(z0_momentum)
                @inbounds for (ll, z0b) in enumerate(z0_thermal)
                    @inbounds for sch in scheme
                        if z0m / z0b >= FT(1e2) || z0b / z0m >= FT(1e2)
                            nothing
                        else
                            thermo_params = SFP.thermodynamics_params(param_set)
                            ts_sfc_test = Thermodynamics.PhaseEquil{FT}(psfc.ρ, psfc.p, psfc.e_int, psfc.q_tot, psfc.T)
                            ts_int_test = Thermodynamics.PhaseEquil{FT}(pint.ρ, pint.p, pint.e_int, pint.q_tot, pint.T)
                            state_in = SF.InteriorValues(pint.z, (pint.u / 10, pint.v / 10), ts_int_test)
                            state_sfc = SF.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test)
                            sc = SF.ValuesOnly{FT}(; state_in, state_sfc, z0m, z0b)
                            DSEᵥ_int = TD.virtual_dry_static_energy(thermo_params, ts_int_test, pint.e_pot)
                            DSEᵥ_sfc = TD.virtual_dry_static_energy(thermo_params, ts_sfc_test, psfc.e_pot)
                            ΔDSEᵥ = DSEᵥ_int - DSEᵥ_sfc
                            if abs(ΔDSEᵥ) <= FT(SFP.cp_d(param_set) / 10)
                                counter[3] += 1
                            else
                                sign(ΔDSEᵥ) == 1 ? counter[1] += 1 : counter[2] += 1
                            end
                            @test try
                                sfcc =
                                    SF.surface_conditions(param_set, sc, sch; maxiter, soltype = RS.VerboseSolution())
                                true
                            catch
                                false
                            end
                        end
                    end
                end
            end
        end
    end
    @info "(Moist $FT) Stable BL: $(counter[1]), Unstable BL: $(counter[2]), Neutral BL: $(counter[3])"
end

@testset "Check convergence (dry thermodynamic states): Stable/Unstable" begin
    for FT in [Float32, Float64]
        profiles_sfc, profiles_int = generate_profiles(FT)
        scheme = [SF.FVScheme(), SF.FDScheme()]
        z0_momentum = Array{FT}(range(1e-6, stop = 1e-1, length = 4))
        z0_thermal = Array{FT}(range(1e-6, stop = 1e-1, length = 4))
        maxiter = 40
        check_over_dry_states(param_set, FT, profiles_int, profiles_sfc, scheme, z0_momentum, z0_thermal, maxiter)
    end
end
#@testset "Check convergence (moist thermodynamic states): Stable/Unstable" begin
#    for FT in [Float32, Float64]
#        profiles_sfc, profiles_int = generate_profiles(FT)
#        scheme = [SF.FVScheme(), SF.FDScheme()]
#        z0_momentum = Array{FT}(range(1e-6, stop = 1e-1, length = 4))
#        z0_thermal = Array{FT}(range(1e-6, stop = 1e-1, length = 4))
#        maxiter = 10
#        check_over_moist_states(param_set, FT, profiles_int, profiles_sfc, scheme, z0_momentum, z0_thermal, maxiter)
#    end
#end
