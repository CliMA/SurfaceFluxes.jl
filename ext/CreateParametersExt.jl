module CreateParametersExt

import Thermodynamics.Parameters.ThermodynamicsParameters
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions:
    AbstractUniversalFunctionParameters,
    BusingerParams,
    GryanikParams,
    GrachevParams,
    BeljaarsParams,
    ChengParams,
    HoltslagParams
import CLIMAParameters as CP

function SurfaceFluxesParameters(::Type{FT}, UFParams) where {FT <: AbstractFloat}
    SurfaceFluxesParameters(CP.create_toml_dict(FT), UFParams)
end

function SurfaceFluxesParameters(toml_dict::CP.AbstractTOMLDict, UFParams)
    name_map = Dict{String, String}("von_karman_constant" => "von_karman_const")
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    FT = CP.float_type(toml_dict)
    ufp = UFParams(toml_dict)
    thermo_params = ThermodynamicsParameters(toml_dict)
    return SurfaceFluxesParameters{FT, typeof(ufp), typeof(thermo_params)}(; parameters..., ufp, thermo_params)
end

BeljaarsParams(::Type{FT}) where {FT <: AbstractFloat} = BeljaarsParams(CP.create_toml_dict(FT))

function BeljaarsParams(toml_dict::CP.AbstractTOMLDict)
    name_map = Dict{String, String}(
        "prandtl_number_0_beljaars" => "Pr_0",
        "coefficient_a_m_beljaars" => "a_m",
        "coefficient_a_h_beljaars" => "a_h",
        "coefficient_b_m_beljaars" => "b_m",
        "coefficient_b_h_beljaars" => "b_h",
        "coefficient_c_m_beljaars" => "c_m",
        "coefficient_c_h_beljaars" => "c_h",
        "coefficient_d_m_beljaars" => "d_m",
        "coefficient_d_h_beljaars" => "d_h",
        "most_stability_parameter_beljaars" => "ζ_a",
        "most_stability_exponent_beljaars" => "γ",
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    FT = CP.float_type(toml_dict)
    return BeljaarsParams{FT}(; parameters...)
end

HoltslagParams(::Type{FT}) where {FT <: AbstractFloat} = HoltslagParams(CP.create_toml_dict(FT))

function HoltslagParams(toml_dict::CP.AbstractTOMLDict)
    name_map = Dict{String, String}(
        "prandtl_number_0_holtslag" => "Pr_0",
        "coefficient_a_m_holtslag" => "a_m",
        "coefficient_a_h_holtslag" => "a_h",
        "coefficient_b_m_holtslag" => "b_m",
        "coefficient_b_h_holtslag" => "b_h",
        "coefficient_c_m_holtslag" => "c_m",
        "coefficient_c_h_holtslag" => "c_h",
        "coefficient_d_m_holtslag" => "d_m",
        "coefficient_d_h_holtslag" => "d_h",
        "most_stability_parameter_holtslag" => "ζ_a",
        "most_stability_exponent_holtslag" => "γ",
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    FT = CP.float_type(toml_dict)
    return HoltslagParams{FT}(; parameters...)
end

BusingerParams(::Type{FT}) where {FT <: AbstractFloat} = BusingerParams(CP.create_toml_dict(FT))

function BusingerParams(toml_dict::CP.AbstractTOMLDict)
    name_map = Dict{String, String}(
        "prandtl_number_0_businger" => "Pr_0",
        "coefficient_a_m_businger" => "a_m",
        "coefficient_a_h_businger" => "a_h",
        "most_stability_parameter_businger" => "ζ_a",
        "most_stability_exponent_businger" => "γ",
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    FT = CP.float_type(toml_dict)
    return BusingerParams{FT}(; parameters...)
end

GryanikParams(::Type{FT}) where {FT <: AbstractFloat} = GryanikParams(CP.create_toml_dict(FT))

function GryanikParams(toml_dict::CP.AbstractTOMLDict)
    name_map = Dict{String, String}(
        "prandtl_number_0_gryanik" => "Pr_0",
        "coefficient_a_m_gryanik" => "a_m",
        "coefficient_a_h_gryanik" => "a_h",
        "coefficient_b_m_gryanik" => "b_m",
        "coefficient_b_h_gryanik" => "b_h",
        "most_stability_parameter_gryanik" => "ζ_a",
        "most_stability_exponent_gryanik" => "γ",
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    FT = CP.float_type(toml_dict)
    return GryanikParams{FT}(; parameters...)
end

GrachevParams(::Type{FT}) where {FT <: AbstractFloat} = GrachevParams(CP.create_toml_dict(FT))

function GrachevParams(toml_dict::CP.AbstractTOMLDict)
    name_map = Dict{String, String}(
        "prandtl_number_0_grachev" => "Pr_0",
        "coefficient_a_m_grachev" => "a_m",
        "coefficient_a_h_grachev" => "a_h",
        "coefficient_b_m_grachev" => "b_m",
        "coefficient_b_h_grachev" => "b_h",
        "coefficient_c_h_grachev" => "c_h",
        "most_stability_parameter_grachev" => "ζ_a",
        "most_stability_exponent_grachev" => "γ",
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    FT = CP.float_type(toml_dict)
    return GrachevParams{FT}(; parameters...)
end

ChengParams(::Type{FT}) where {FT <: AbstractFloat} = ChengParams(CP.create_toml_dict(FT))

function ChengParams(toml_dict::CP.AbstractTOMLDict)
    name_map = Dict{String, String}(
        "prandtl_number_0_cheng" => "Pr_0",
        "coefficient_a_m_cheng" => "a_m",
        "coefficient_a_h_cheng" => "a_h",
        "coefficient_b_m_cheng" => "b_m",
        "coefficient_b_h_cheng" => "b_h",
        "most_stability_parameter_cheng" => "ζ_a",
        "most_stability_exponent_cheng" => "γ",
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    FT = CP.float_type(toml_dict)
    return ChengParams{FT}(; parameters...)
end

end
