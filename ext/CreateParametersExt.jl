module CreateParametersExt

import Thermodynamics.Parameters.ThermodynamicsParameters
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions:
    AbstractUniversalFunctionParameters,
    BusingerParams,
    GryanikParams,
    GrachevParams
import ClimaParams as CP

"""
    SurfaceFluxesParameters(::Type{FT}, UFParams)

Construct `SurfaceFluxesParameters` for the specified floating-point type `FT` and
universal function parameterization type `UFParams`.

Creates a TOML parameter dictionary using default values from ClimaParams and constructs
the parameter set.

# Arguments
- `FT`: Floating-point type (e.g., `Float32`, `Float64`).
- `UFParams`: Universal function parameterization type (e.g., `BusingerParams`, `GryanikParams`, `GrachevParams`).

# Returns
A `SurfaceFluxesParameters` instance containing surface flux parameters, universal function
parameters, and thermodynamics parameters.
"""
function SurfaceFluxesParameters(
    ::Type{FT},
    UFParams,
) where {FT <: AbstractFloat}
    SurfaceFluxesParameters(CP.create_toml_dict(FT), UFParams)
end

"""
    SurfaceFluxesParameters(toml_dict, UFParams)

Construct `SurfaceFluxesParameters` from a TOML parameter dictionary.

# Arguments
- `toml_dict`: A `ClimaParams.ParamDict` containing parameter values.
- `UFParams`: Universal function parameterization type (e.g., `BusingerParams`, `GryanikParams`, `GrachevParams`).

# Returns
A `SurfaceFluxesParameters` instance with parameters read from the TOML dictionary.
"""
function SurfaceFluxesParameters(toml_dict::CP.ParamDict{FT}, UFParams) where {FT}
    name_map = (; :von_karman_constant => :von_karman_const)
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    ufp = UFParams(toml_dict)
    thermo_params = ThermodynamicsParameters(toml_dict)
    return SurfaceFluxesParameters{FT, typeof(ufp), typeof(thermo_params)}(;
        parameters...,
        ufp,
        thermo_params,
    )
end

"""
    BusingerParams(::Type{FT})

Construct `BusingerParams` for the specified floating-point type using default parameter
values from ClimaParams.

# Arguments
- `FT`: Floating-point type (e.g., `Float32`, `Float64`).

# Returns
A `BusingerParams` instance with default parameter values from the TOML configuration.

See also: [`BusingerParams(toml_dict)`](@ref)
"""
BusingerParams(::Type{FT}) where {FT <: AbstractFloat} =
    BusingerParams(CP.create_toml_dict(FT))

"""
    BusingerParams(toml_dict)

Construct `BusingerParams` from a TOML parameter dictionary.

Reads Businger-specific parameters from the TOML dictionary and maps them to the
`BusingerParams` struct fields.

# Arguments
- `toml_dict`: A `ClimaParams.ParamDict` containing parameter values.

# Returns
A `BusingerParams` instance with parameters read from the TOML dictionary.

# Parameter Mappings
The following TOML keys are read from the "SurfaceFluxes" section:
- `prandtl_number_0_businger` → `Pr_0`
- `coefficient_a_m_businger` → `a_m`
- `coefficient_a_h_businger` → `a_h`
- `coefficient_b_m_businger` → `b_m`
- `coefficient_b_h_businger` → `b_h`
- `most_stability_parameter_businger` → `ζ_a`
- `most_stability_exponent_businger` → `γ`

See Businger et al. (1971) and Nishizawa & Kitamura (2018).
"""
function BusingerParams(toml_dict::CP.ParamDict{FT}) where {FT}
    name_map = (;
        :prandtl_number_0_businger => :Pr_0,
        :coefficient_a_m_businger => :a_m,
        :coefficient_a_h_businger => :a_h,
        :coefficient_b_m_businger => :b_m,
        :coefficient_b_h_businger => :b_h,
        :most_stability_parameter_businger => :ζ_a,
        :most_stability_exponent_businger => :γ,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    return BusingerParams{FT}(; parameters...)
end

"""
    _get_businger_unstable_params(toml_dict)

Helper function to get Businger b_m and b_h parameters for use in unstable branches
of Gryanik and Grachev parameterizations.
"""
function _get_businger_unstable_params(toml_dict::CP.ParamDict{FT}) where {FT}
    businger_name_map = (;
        :coefficient_b_m_businger => :b_m,
        :coefficient_b_h_businger => :b_h,
    )
    businger_params = CP.get_parameter_values(toml_dict, businger_name_map, "SurfaceFluxes")
    return (; b_m_unstable = businger_params.b_m, b_h_unstable = businger_params.b_h)
end

GryanikParams(::Type{FT}) where {FT <: AbstractFloat} =
    GryanikParams(CP.create_toml_dict(FT))

"""
    GryanikParams(toml_dict)

Construct `GryanikParams` from a TOML parameter dictionary.

Reads Gryanik-specific parameters from the TOML dictionary and maps them to the
`GryanikParams` struct fields. The unstable branch parameters (`b_m_unstable` and
`b_h_unstable`) are automatically set to the Businger parameter values.

# Arguments
- `toml_dict`: A `ClimaParams.ParamDict` containing parameter values.

# Returns
A `GryanikParams` instance with parameters read from the TOML dictionary.

# Parameter Mappings
The following TOML keys are read from the "SurfaceFluxes" section:
- `prandtl_number_0_gryanik` → `Pr_0`
- `coefficient_a_m_gryanik` → `a_m`
- `coefficient_a_h_gryanik` → `a_h`
- `coefficient_b_m_gryanik` → `b_m`
- `coefficient_b_h_gryanik` → `b_h`
- `most_stability_parameter_gryanik` → `ζ_a`
- `most_stability_exponent_gryanik` → `γ`

Additionally, `b_m_unstable` and `b_h_unstable` are set from the Businger parameters
(`coefficient_b_m_businger` and `coefficient_b_h_businger`) for use in unstable branches.

See Gryanik et al. (2020).
"""
function GryanikParams(toml_dict::CP.ParamDict{FT}) where {FT}
    name_map = (;
        :prandtl_number_0_gryanik => :Pr_0,
        :coefficient_a_m_gryanik => :a_m,
        :coefficient_a_h_gryanik => :a_h,
        :coefficient_b_m_gryanik => :b_m,
        :coefficient_b_h_gryanik => :b_h,
        :most_stability_parameter_gryanik => :ζ_a,
        :most_stability_exponent_gryanik => :γ,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    unstable_params = _get_businger_unstable_params(toml_dict)
    return GryanikParams{FT}(; parameters..., unstable_params...)
end

"""
    GrachevParams(::Type{FT})

Construct `GrachevParams` for the specified floating-point type using default parameter
values from ClimaParams.

# Arguments
- `FT`: Floating-point type (e.g., `Float32`, `Float64`).

# Returns
A `GrachevParams` instance with default parameter values from the TOML configuration.

See also: [`GrachevParams(toml_dict)`](@ref)
"""
GrachevParams(::Type{FT}) where {FT <: AbstractFloat} =
    GrachevParams(CP.create_toml_dict(FT))

"""
    GrachevParams(toml_dict)

Construct `GrachevParams` from a TOML parameter dictionary.

Reads Grachev-specific parameters from the TOML dictionary and maps them to the
`GrachevParams` struct fields. The unstable branch parameters (`b_m_unstable` and
`b_h_unstable`) are automatically set to the Businger parameter values.

# Arguments
- `toml_dict`: A `ClimaParams.ParamDict` containing parameter values.

# Returns
A `GrachevParams` instance with parameters read from the TOML dictionary.

# Parameter Mappings
The following TOML keys are read from the "SurfaceFluxes" section:
- `prandtl_number_0_grachev` → `Pr_0`
- `coefficient_a_m_grachev` → `a_m`
- `coefficient_a_h_grachev` → `a_h`
- `coefficient_b_m_grachev` → `b_m`
- `coefficient_b_h_grachev` → `b_h`
- `coefficient_c_h_grachev` → `c_h`
- `most_stability_parameter_grachev` → `ζ_a`
- `most_stability_exponent_grachev` → `γ`

Additionally, `b_m_unstable` and `b_h_unstable` are set from the Businger parameters
(`coefficient_b_m_businger` and `coefficient_b_h_businger`) for use in unstable branches.

See Grachev et al. (2007).
"""
function GrachevParams(toml_dict::CP.ParamDict{FT}) where {FT}
    name_map = (;
        :prandtl_number_0_grachev => :Pr_0,
        :coefficient_a_m_grachev => :a_m,
        :coefficient_a_h_grachev => :a_h,
        :coefficient_b_m_grachev => :b_m,
        :coefficient_b_h_grachev => :b_h,
        :coefficient_c_h_grachev => :c_h,
        :most_stability_parameter_grachev => :ζ_a,
        :most_stability_exponent_grachev => :γ,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "SurfaceFluxes")
    unstable_params = _get_businger_unstable_params(toml_dict)
    return GrachevParams{FT}(; parameters..., unstable_params...)
end

end
