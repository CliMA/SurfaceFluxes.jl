import CLIMAParameters as CP
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD

function create_uf_parameters(toml_dict, ::UF.GryanikType)
    FT = CP.float_type(toml_dict)

    aliases = ["Pr_0_Gryanik", "a_m_Gryanik", "a_h_Gryanik", "b_m_Gryanik", "b_h_Gryanik"]

    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple

    pairs = (;
        Pr_0 = pairs.Pr_0_Gryanik,
        a_m = pairs.a_m_Gryanik,
        a_h = pairs.a_h_Gryanik,
        b_m = pairs.b_m_Gryanik,
        b_h = pairs.b_h_Gryanik,
    )
    return UF.GryanikParams{FT}(; pairs...)
end

function create_uf_parameters(toml_dict, ::UF.BusingerType)
    FT = CP.float_type(toml_dict)
    aliases = ["Pr_0_Businger", "a_m_Businger", "a_h_Businger"]

    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple

    pairs = (; Pr_0 = pairs.Pr_0_Businger, a_m = pairs.a_m_Businger, a_h = pairs.a_h_Businger)
    return UF.BusingerParams{FT}(; pairs...)
end

function create_uf_parameters(toml_dict, ::UF.GrachevType)
    FT = CP.float_type(toml_dict)
    aliases = ["Pr_0_Grachev", "a_m_Grachev", "a_h_Grachev", "b_m_Grachev", "b_h_Grachev", "c_h_Grachev"]

    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple

    pairs = (;
        Pr_0 = pairs.Pr_0_Grachev,
        a_m = pairs.a_m_Grachev,
        a_h = pairs.a_h_Grachev,
        b_m = pairs.b_m_Grachev,
        b_h = pairs.b_h_Grachev,
        c_h = pairs.c_h_Grachev,
    )
    return UF.GrachevParams{FT}(; pairs...)
end

function create_parameters(toml_dict, ufpt)
    FT = CP.float_type(toml_dict)

    ufp = create_uf_parameters(toml_dict, ufpt)
    AUFP = typeof(ufp)

    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TD.Parameters.ThermodynamicsParameters{FT}(; pairs...)
    TP = typeof(thermo_params)

    aliases = ["von_karman_const"]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "SurfaceFluxesParameters")
    return SFP.SurfaceFluxesParameters{FT, AUFP, TP}(; pairs..., ufp, thermo_params)
end
