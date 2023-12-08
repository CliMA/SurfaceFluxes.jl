# # Basic usage

# ```julia
# using SurfaceFluxes
# using Thermodynamics
# using StaticArrays
# import CLIMAParameters
# ```

# Next we define a helper function to estimate density of air at the surface
# given atmospheric conditions higher up, using the ideal gas
# law and an assumption of hydrostatic balance.
# This is useful when your surface model does not model the density of
# air.
function compute_ρ_sfc(thermo_params, ts_atmos, T_sfc)
    T_atmos = Thermodynamics.air_temperature(thermo_params, ts_atmos)
    Rm_atmos = Thermodynamics.gas_constant_air(thermo_params, ts_atmos)
    ρ_atmos = Thermodynamics.air_density(thermo_params, ts_atmos)
    ρ_sfc =
        ρ_atmos *
        (T_sfc / T_atmos)^(Thermodynamics.cv_m(thermo_params, ts_atmos) / Rm_atmos)
    return ρ_sfc
end

# # Set up thermodynamics and surface fluxes parameters

# ```julia
# const CP = CLIMAParameters
# include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
# FT = Float32
# toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
# uf_type = UF.BusingerType()
# param_set = create_parameters(toml_dict, uf_type)
# thermo_params = SurfaceFluxes.Parameters.thermodynamics_params(param_set)
# ```

# # Define some parameters that are intrinsic to your surface model
# Roughness lengths and displacement height
# ```julia
# z_0m = FT(0.01)
# z_0b = FT(0.001)
# d_sfc = FT(0)
# ```

# # Define conditions at the lowest level in the atmosphere
# ```julia
# h  = FT(10) # height at which measurements are made, in m
# u = FT(4) # horizontal windspeed
# P = FT(101350)
# T = FT(298.15)
# q = FT(0.006)
# ts_atmos = Thermodynamics.PhaseEquil_pTq(thermo_params, P, T, q) # you could use ρ instead
# ```


# # Repeat for the land
# ```julia
# u_sfc = FT(0)
# T_sfc = FT(290)
# ρ_sfc = compute_ρ_sfc(thermo_params, ts_atmos, T_sfc)
# q_sfc = FT(0.003)
# ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)
# h_sfc = FT(0) # height at which surface measurements are made
# ```

# # Compute fluxes

# SurfaceFluxes.jl expects a relative difference between where u = 0 
# and the atmosphere height. Here, we assume h and h_sfc are measured
# relative to a common reference. Then d_sfc + h_sfc + z_0m is the apparent
# source of momentum, and
# Δh ≈ h - d_sfc - h_sfc is the relative height difference between the
# apparent source of momentum and the atmosphere height.

# In this we have neglected z_0m and z_0b (i.e. assumed they are small
# compared to Δh).

# ```julia
# state_sfc = SurfaceFluxes.StateValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
# state_atmos = SurfaceFluxes.StateValues(
#    h - d_sfc - h_sfc,
#    SVector{2, FT}(u, 0),
#    ts_atmos,
#    )
# sc = SurfaceFluxes.ValuesOnly(
#     state_atmos,
#     state_sfc,
#     z_0m,
#     z_0b,
# )
# conditions = SurfaceFluxes.surface_conditions(
#     param_set,
#     sc;
#     )
# ```
