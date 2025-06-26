include("../src/similarity_theory.jl")

using Test
import Thermodynamics as TD
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

FT = Float32
# Parameters (ClimaParams types)
param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
thermo_params = param_set.thermo_params

# Things that do not vary between test sets
h_a = FT(15)
u_a = FT(5)
v_a = FT(1)
gustiness = FT(1)
T_a = FT(298)
ρ_a = FT(1.2)
atmos_args = (;ρ = ρ_a,)
q_sat_a = TD.q_vap_saturation_generic(
    thermo_params,
    T_a,
    ρ_a,
    TD.Liquid(),
)
ζ₀ = FT(-10)
L★ = h_a ./ ζ₀
sfc_params = SFP.uf_params(param_set)
similarity_theory = SFP.universal_func_type(param_set)
similarity_profile = UF.universal_func(similarity_theory, L★, sfc_params)
similarity_profile = ufunc

@testset "q_atmos < q_land, β > 0, T_land = T_atmos" begin
    # Surface state
    function q_surface(surface_args, similarity_scales, atmos_state, param_set)
        q_a = atmos_state.q_a
        β = surface_args.β
        θ_s = surface_args.θ_s # instead pass surface_state as arg
        θ_a = atmos_state.θ_a
        ρ_a = atmos_state.args.ρ
        cv_m = TD.cv_m(thermo_params, TD.PhasePartition(q_a))
        R_m =  TD.gas_constant_air(thermo_params, TD.PhasePartition(q_a))
        ρ_s = ρ_a * (θ_s/θ_a)^(cv_m/R_m)
        q_sat = TD.q_vap_saturation_generic(
            thermo_params,
            θ_s,
            ρ_s,
            TD.Liquid(),
        )
        return q_sat*β + (1-β)*q_a
    end
    T_s = T_a
    surface_args = (; β = FT(0.5), θ_s = T_s)
    surface_state = SurfaceState(
        (𝑧0m=FT(0.01), 𝑧0θ=FT(0.01), 𝑧0q=FT(0.01)),
        FT(0), #u
        FT(0), # v
        q_surface, # q function
        T_a, # θ
        FT(0), # h
        surface_args,
    )

    # Atmos state (partial, we will vary q_a)
    q_a = [q_sat_a*0.01f0, q_sat_a*0.5f0, q_sat_a-eps(FT)]
    similarity_profile = ufunc
    lhf = []
    for atmos_q in q_a
        atmos_state = AtmosState(
            u_a, #u
            v_a, #v
            atmos_q, #q
            T_a, # θ
            h_a, # h
            gustiness, # gustiness
            FT(100), #h boundary layer
            atmos_args # extra args
        )
        fluxes = compute_similarity_theory_fluxes(similarity_profile, 
                                                  surface_state,
                                                  atmos_state, 
                                                  param_set)
        
        @test fluxes.sensible_heat ≈ 0
        push!(lhf, fluxes.latent_heat)
    end
    @test abs(lhf[2]/lhf[1] - FT(0.5)) < 0.05
    @test lhf[3] < 0.01
    @test all(lhf .> 0)
end


@testset "q_atmos > q_land, β <= 1, T_land < T_atmos" begin
    #Surface state (partial, surface_args specified below)
    T_s = FT(280)
    function q_surface(surface_args, similarity_scales, atmos_state, param_set)
        q_a = atmos_state.q_a
        β = surface_args.β
        θ_s = surface_args.θ_s # instead pass surface_state as arg
        θ_a = atmos_state.θ_a
        ρ_a = atmos_state.args.ρ
        cv_m = TD.cv_m(thermo_params, TD.PhasePartition(q_a))
        R_m =  TD.gas_constant_air(thermo_params, TD.PhasePartition(q_a))
        ρ_s = ρ_a * (θ_s/θ_a)^(cv_m/R_m)
        q_sat = TD.q_vap_saturation_generic(
            thermo_params,
            θ_s,
            ρ_s,
            TD.Liquid(),
        )
        return q_sat*β + (1-β)*q_a
    end

    # Atmos state (full)
    atmos_state = AtmosState(
        u_a, #u
        v_a, #v
        q_sat_a, #q
        T_a, # θ
        h_a, # h
        gustiness, # gustiness
        FT(100), #h boundary layer
        atmos_args # extra args
    )
    similarity_profile = ufunc
    lhf = []
    for β in [FT(0), FT(0.5), FT(1)]
        surface_args = (; β = β, θ_s = T_s)
        surface_state = SurfaceState(
            (𝑧0m=FT(0.01), 𝑧0θ=FT(0.01), 𝑧0q=FT(0.01)),
            FT(0), #u
            FT(0), # v
            q_surface, # q function
            T_s, # θ
            FT(0), # h
            surface_args,
        )
        fluxes = compute_similarity_theory_fluxes(similarity_profile, 
                                                  surface_state,
                                                  atmos_state, 
                                                  param_set)
        
        push!(lhf, fluxes.latent_heat)
    end
    @test abs(lhf[2]/lhf[3] - FT(0.5)) < 0.05
    @test all(lhf[2:3] .<0)
    @test lhf[1] ≈ 0
end

@testset "q_atmos = q_land, T_atmos = T_land, β = 0" begin
    #Surface state
    T_s = T_a
    # TODO: `surface_state` as arg in `surface_variable` function
    function q_surface(surface_args, similarity_scales, atmos_state, param_set)
        q_a = atmos_state.q_a
        β = surface_args.β
        θ_s = surface_args.θ_s # instead pass surface_state as arg
        θ_a = atmos_state.θ_a
        ρ_a = atmos_state.args.ρ
        cv_m = TD.cv_m(thermo_params, TD.PhasePartition(q_a))
        R_m =  TD.gas_constant_air(thermo_params, TD.PhasePartition(q_a))
        ρ_s = ρ_a * (θ_s/θ_a)^(cv_m/R_m)
        q_sat = TD.q_vap_saturation_generic(
            thermo_params,
            θ_s,
            ρ_s,
            TD.Liquid(),
        )
        return q_sat*β + (1-β)*q_a
    end
    surface_args = (; β = 0f0, θ_s = T_s)
    surface_state = SurfaceState(
        (𝑧0m=FT(0.01), 𝑧0θ=FT(0.01), 𝑧0q=FT(0.01)),
        FT(0), #u
        FT(0), # v
        q_surface, # q function
        T_s, # θ
        FT(0), # h
        surface_args,
    )
    
    # Atmos state
    atmos_state = AtmosState(
        u_a, #u
        v_a, #v
        q_a, #q
        T_a, # θ
        h_a, # h
        gustiness, # gustiness
        FT(100), #h boundary layer
        atmos_args # extra args
    )
    fluxes = compute_similarity_theory_fluxes(similarity_profile, 
                                              surface_state,
                                              atmos_state, 
                                              param_set)
    @test fluxes.latent_heat ≈ 0
    @test fluxes.sensible_heat ≈ 0
end
