include("../src/similarity_theory.jl")

using Test
import Thermodynamics as TD
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP
import ArtifactWrappers as AW
import NCDatasets as NC

# Assume constant values
# For method tests define functions for roughness that return constant values


FT = Float32 # TODO Check all floattypes
Σ₀ = SimilarityScales{FT, FT, FT}(1e-4,1e-4,1e-4)
Σₜ = SimilarityScales{FT, FT, FT}(1e-5,1e-5,1e-5)
ΔΣ = Σₜ- Σ₀

# Parameters (ClimaParams types)
param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
thermo_params = param_set.thermo_params

# Assign states

function z0test(surface_args, similarity_scales, atmos_state, param_set)
    u★ = similarity_scales.momentum
    FT =typeof(u★)
    return FT(0.015 * u★^2 / 9.81)
end

function generate_profiles(;FT=Float32, uf_params = UF.BusingerParams)
    param_set = SurfaceFluxesParameters(FT, uf_params)
    thermo_params = param_set.thermo_params
    profiles = collect(TD.TestedProfiles.PhaseEquilProfiles(thermo_params, Array{FT}))
    profiles_sfc = filter(p -> iszero(p.z), profiles)
    profiles_int = filter(p -> !iszero(p.z), profiles)
    ## Properties contained in `profiles_<sfc, int>`
    ## :z, :T, :p, :RS, :e_int, :h, :ρ, 
    ## :θ_liq_ice, :q_tot, :q_liq, :q_ice, :q_pt, :RH, 
    ## :e_pot, :u, :v, :w, :e_kin, :phase_type
    return profiles_sfc, profiles_int, param_set
end

function get_result(; FT=Float32)
    uf_params = UF.BusingerParams
    profiles_sfc, profiles_int = generate_profiles(; uf_params)
    function run_op(profiles_int, profiles_sfc)
        atmos_state = AtmosState(
                          profiles_int.u ./ 10,
                          profiles_int.v ./ 10,
                          profiles_int.q_tot,
                          profiles_int.θ_liq_ice,
                          profiles_int.z,
                          FT(1), # gustiness needs to be a function of u,v, ustar
                          FT(100),
                          (ρ=profiles_int.ρ, arg𝑏=FT(0.01), arg𝑐=z0test),
                        )
        surface_state = SurfaceState(
                          (𝑧0m=FT(0.01), 𝑧0θ=FT(0.01), 𝑧0q=z0test),
                          FT(0),
                          FT(0),
                          profiles_sfc.q_tot,
                          profiles_sfc.θ_liq_ice,
                          FT(0),
                          (arg𝑎=FT(0.01), arg𝑏=FT(0.01), arg𝑐=z0test),
                        )
        # Line by line debug and test for `refine` function
        Δstate = state_differences(surface_state, atmos_state, Σ₀, param_set); 
        (; 𝑧0m, 𝑧0θ, 𝑧0q) = surface_state.roughness_lengths
        # Generic info block 
        ζ₀ = FT(-10)
        L★ = Δstate.Δh ./ ζ₀
        sfc_params = SFP.uf_params(param_set)
        similarity_theory = SFP.universal_func_type(param_set)
        ufunc = UF.universal_func(similarity_theory, L★, sfc_params)
        Σ_est = (momentum=FT(0.1),temperature=FT(0.01),water_vapor=FT(0.001))
        ΔU_est = FT(10)
        # Initial guess
        similarity_profile = ufunc
        similarity_scales = refine_similarity_variables(Σ_est, ΔU_est, 
                                             similarity_profile,
                                             surface_state, 
                                             atmos_state, param_set)
        fluxes = compute_similarity_theory_fluxes(similarity_profile, 
                                         surface_state,
                                         atmos_state, 
                                         param_set)
        fluxes.scale_vars
    end
    [run_op(ii, jj) for ii in profiles_int, jj in profiles_sfc]
end

function run_profiles()
    PyCLES_output_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "PyCLES_output",
        AW.ArtifactFile[
        AW.ArtifactFile(url = "https://caltech.box.com/shared/static/johlutwhohvr66wn38cdo7a6rluvz708.nc", filename = "Rico.nc",),
        AW.ArtifactFile(url = "https://caltech.box.com/shared/static/zraeiftuzlgmykzhppqwrym2upqsiwyb.nc", filename = "Gabls.nc",),
        AW.ArtifactFile(url = "https://caltech.box.com/shared/static/toyvhbwmow3nz5bfa145m5fmcb2qbfuz.nc", filename = "DYCOMS_RF01.nc",),
        AW.ArtifactFile(url = "https://caltech.box.com/shared/static/jci8l11qetlioab4cxf5myr1r492prk6.nc", filename = "Bomex.nc",),
        ],
    )
    #! format: on
    PyCLES_output_dataset_path = AW.get_data_folder(PyCLES_output_dataset)
    files = ["DYCOMS_RF01.nc", "Bomex.nc", "Rico.nc", "Gabls.nc"];
    for f in files
        @info "Casename: $f"

        nt = NC.NCDataset(joinpath(PyCLES_output_dataset_path, f), "r") do data
            prof = data.group["profiles"]
            timeseries = data.group["timeseries"]
            z = Array(data["z_half"])
            u = Array(prof["u_mean"])
            v = Array(prof["v_mean"])
            qt = Array(prof["qt_mean"])
            ql = Array(prof["ql_mean"])
            ρ = Array(prof["rho"])
            p0 = Array(prof["p0"])
            b = Array(prof["buoyancy_mean"])
            θli = Array(prof["thetali_mean"]) # Care with variable names across cases ()?
            T = Array(prof["temperature_mean"]) # Care with variable names across cases ()?
            SHF = Array(timeseries["shf_surface_mean"]) # Care with variable names across cases ()?
            LHF = Array(timeseries["lhf_surface_mean"]) # Care with variable names across cases ()?
            (; z, u, v, qt, ql, ρ, p0, b, θli, T, SHF, LHF)
        end
        z, u, v, qt, ql, ρ, p0, b, θli, T, SHF, LHF = nt

        function getval(X) # Consider only the last timestep
            X[:, end]
        end

        # Data at first interior node (x_in)
        # TODO Make sure that the first node ii = 1 is at the "surface"
        # for the tests to be consistent
        ii = 2

        shf = mean(getval(SHF))
        lhf = mean(getval(LHF))

        z_in = z[ii]
        u_in = getval(u)[ii]
        v_in = getval(v)[ii]
        b_in = getval(b)[ii]
        ρ_in = getval(ρ)[ii]
        θ_in = getval(θli)[ii]
        T_in = getval(T)[ii]
        qt_in = getval(qt)[ii]

        # Surface values for variables
        jj = 1
        u_sfc = FT(0)
        v_sfc = FT(0)
        b_sfc = getval(b)[jj]
        ρ_sfc = getval(ρ)[jj]
        qt_sfc = getval(qt)[jj]
        θ_sfc = getval(θli)[jj]
        T_sfc = getval(T)[jj]
        z_sfc = FT(0)

        # Roughness
        ## Roughness lengths
        z0m = FT(0.001)
        z0b = FT(0.001)

        ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc, qt_sfc)
        ts_in = TD.PhaseEquil_ρθq(thermo_params, ρ_in, θ_in, qt_in)

        u_in = (FT(u_in), FT(v_in))
        u_sfc = (FT(u_sfc), FT(v_sfc))

        state_sfc = SF.StateValues(z_sfc, u_sfc, ts_sfc)
        state_in = SF.StateValues(z_in, u_in, ts_in)

        kwargs = (; state_in, state_sfc, z0m, z0b)

        # Initial guess
        sfc_params = SFP.uf_params(param_set)
        similarity_theory = SFP.universal_func_type(param_set)

        if f == "DYCOMS_RF01.nc"
            atmos_state = AtmosState(
                              u_in[1],
                              u_in[2],
                              qt_in,
                              θ_in,
                              z_in,
                              FT(1), # gustiness needs to be a function of u,v, ustar
                              FT(100),
                              (ρ=ρ_sfc, arg𝑏=FT(0.01), arg𝑐=z0test),
                            )
            surface_state = SurfaceState(
                              (𝑧0m=FT(0.001), 𝑧0θ=FT(0.001), 𝑧0q=z0test),
                              FT(0),
                              FT(0),
                              qt_sfc,
                              θ_sfc,
                              FT(0),
                              (arg𝑎=FT(0.01), arg𝑏=FT(0.01), arg𝑐=z0test),
                            )
            # Line by line debug and test for `refine` function
            Δstate = state_differences(surface_state, atmos_state, Σ₀, param_set); 
            ζ₀ = FT(-10)
            L★ = Δstate.Δh ./ ζ₀
            ufunc = UF.universal_func(similarity_theory, L★, sfc_params)
            Σ_est = (momentum=FT(0.1),temperature=FT(0.01),water_vapor=FT(0.001))
            ΔU_est = FT(10)
            similarity_profile = ufunc
            (; 𝑧0m, 𝑧0θ, 𝑧0q) = surface_state.roughness_lengths
            fluxes = compute_similarity_theory_fluxes(similarity_profile, 
                                             surface_state,
                                             atmos_state, 
                                             param_set)
            @show fluxes
        elseif f == "Bomex.nc"
            #sc = SF.FluxesAndFrictionVelocity(state_in, state_sfc, shf, lhf, FT(0.28), z0m, z0b)
        elseif f == "Rico.nc"
            #sc = SF.Coefficients(state_in, state_sfc, FT(0.001229), FT(0.001094), z0m, z0b)
        elseif f == "Gabls.nc"
            atmos_state = AtmosState(
                              u_in[1],
                              u_in[2],
                              qt_in,
                              θ_in,
                              z_in,
                              FT(1), # gustiness needs to be a function of u,v, ustar
                              FT(100),
                              (ρ=ρ_sfc, arg𝑏=FT(0.01), arg𝑐=z0test),
                            )
            surface_state = SurfaceState(
                              (𝑧0m=FT(0.001), 𝑧0θ=FT(0.001), 𝑧0q=z0test),
                              FT(0),
                              FT(0),
                              qt_sfc,
                              θ_sfc,
                              FT(0),
                              (arg𝑎=FT(0.01), arg𝑏=FT(0.01), arg𝑐=z0test),
                            )
            # Line by line debug and test for `refine` function
            Δstate = state_differences(surface_state, atmos_state, Σ₀, param_set); 
            (; 𝑧0m, 𝑧0θ, 𝑧0q) = surface_state.roughness_lengths
            # Initial guess
            ζ₀ = FT(-10)
            L★ = Δstate.Δh ./ ζ₀
            sfc_params = SFP.uf_params(param_set)
            similarity_theory = SFP.universal_func_type(param_set)
            ufunc = UF.universal_func(similarity_theory, L★, sfc_params)
            Σ_est = (momentum=FT(0.1),temperature=FT(0.01),water_vapor=FT(0.001))
            ΔU_est = FT(10)
            similarity_profile = ufunc
            fluxes = compute_similarity_theory_fluxes(similarity_profile, 
                                             surface_state,
                                             atmos_state, 
                                             param_set)
            @show fluxes
        end
    end
end
