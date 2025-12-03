# Tests for utility functions, including splattable methods

using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD

@testset "Splattable utility methods" begin
    FT = Float32
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)
    
    # Create a simple test case
    ρ_sfc = FT(1.15)
    ρ_in = FT(1.13)
    T_sfc = FT(280.0)
    T_in = FT(275.0)
    qt_sfc = FT(0.01)
    qt_in = FT(0.009)
    z0m = FT(1e-3)
    z0b = FT(1e-4)
    
    ts_sfc = TD.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, qt_sfc)
    ts_in = TD.PhaseEquil_ρTq(thermo_params, ρ_in, T_in, qt_in)
    state_sfc = SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc)
    state_in = SF.StateValues(FT(10), (FT(5), FT(0)), ts_in)
    sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
    
    @testset "surface_temperature splattable methods" begin
        # Base method with no args
        T_base = SF.surface_temperature(param_set, sc)
        @test T_base isa FT
        @test T_base ≈ T_sfc
        
        # Base method with args = nothing
        T_nothing = SF.surface_temperature(param_set, sc, nothing)
        @test T_nothing ≈ T_base
        
        # Verify that T_sfc and surface_temperature are the same when args = nothing
        T_sfc_result = SF.T_sfc(param_set, sc, nothing)
        @test T_sfc_result ≈ T_base
        @test T_sfc_result ≈ T_nothing
        
        # Splattable method with single argument (ignores args, returns base)
        u★ = FT(0.5)
        T_splat1 = SF.surface_temperature(param_set, sc, u★)
        @test T_splat1 isa FT
        @test T_splat1 ≈ T_base
        
        # Alternative formulation with multiple Real arguments (T1, T2) - returns T_sfc + T1 + T2
        T1 = FT(2.5)
        T2 = FT(1.5)
        T_splat2 = SF.surface_temperature(param_set, sc, T1, T2)
        @test T_splat2 isa FT
        # Verify alternative formulation: T_sfc + T1 + T2
        @test T_splat2 ≈ T_base + T1 + T2
        
        # Test with 3+ arguments (should match varargs method which ignores args)
        T3 = FT(0.5)
        T4 = FT(1.0)
        T_splat4 = SF.surface_temperature(param_set, sc, T1, T2, T3, T4)
        @test T_splat4 isa FT
        # Varargs method ignores arguments, so should return base
        @test T_splat4 ≈ T_base
        
        # Splattable method with tuple (ignores args, returns base)
        args_tuple = (; u★, q★ = FT(0.001))
        T_tuple = SF.surface_temperature(param_set, sc, args_tuple)
        @test T_tuple isa FT
        @test T_tuple ≈ T_base
    end
    
    @testset "surface_specific_humidity splattable methods" begin
        # Base method with no args
        q_base = SF.surface_specific_humidity(param_set, sc)
        @test q_base isa FT
        @test q_base ≈ qt_sfc
        
        # Base method with args = nothing
        q_nothing = SF.surface_specific_humidity(param_set, sc, nothing)
        @test q_nothing ≈ q_base
        
        # Verify that qt_sfc and surface_specific_humidity are the same when args = nothing
        qt_sfc_result = SF.qt_sfc(param_set, sc, nothing)
        @test qt_sfc_result ≈ q_base
        @test qt_sfc_result ≈ q_nothing
        
        # Splattable method with single argument (ignores args, returns base)
        u★ = FT(0.5)
        q_splat1 = SF.surface_specific_humidity(param_set, sc, u★)
        @test q_splat1 isa FT
        @test q_splat1 ≈ q_base
        
        # Alternative formulation with multiple Real arguments (Q1, Q2) - returns qt_sfc + Q1 + Q2
        Q1 = FT(0.002)
        Q2 = FT(0.001)
        q_splat2 = SF.surface_specific_humidity(param_set, sc, Q1, Q2)
        @test q_splat2 isa FT
        @test q_splat2 ≈ q_base + Q1 + Q2
        
        # Test with 3+ arguments (should match varargs method which ignores args)
        Q3 = FT(0.0005)
        Q4 = FT(0.0003)
        q_splat4 = SF.surface_specific_humidity(param_set, sc, Q1, Q2, Q3, Q4)
        @test q_splat4 isa FT
        # Varargs method ignores arguments, so should return base
        @test q_splat4 ≈ q_base
        
        # Splattable method with tuple (ignores args, returns base)
        args_tuple = (; u★, q★ = FT(0.001))
        q_tuple = SF.surface_specific_humidity(param_set, sc, args_tuple)
        @test q_tuple isa FT
        @test q_tuple ≈ q_base
    end
    
    @testset "Consistency between base and splattable methods" begin
        # Test that varargs method ignores args and returns base result
        # Use 3+ arguments to ensure we call the varargs method (not the 2-arg alternative formulation)
        u★ = FT(0.5)
        q★ = FT(0.001)
        extra_arg = FT(1.0)
        
        T_base = SF.surface_temperature(param_set, sc, nothing)
        # Call with 3 args to match varargs method, which should ignore them
        T_splat = SF.surface_temperature(param_set, sc, u★, q★, extra_arg)
        @test T_base ≈ T_splat
        
        q_base = SF.surface_specific_humidity(param_set, sc, nothing)
        # Call with 3 args to match varargs method, which should ignore them
        q_splat = SF.surface_specific_humidity(param_set, sc, u★, q★, extra_arg)
        @test q_base ≈ q_splat
        
        # Test alternative formulation with Real arguments (exactly 2)
        T1 = FT(2.5)
        T2 = FT(1.5)
        T_alt = SF.surface_temperature(param_set, sc, T1, T2)
        @test T_alt ≈ T_base + T1 + T2
        
        Q1 = FT(0.002)
        Q2 = FT(0.001)
        q_alt = SF.surface_specific_humidity(param_set, sc, Q1, Q2)
        @test q_alt ≈ q_base + Q1 + Q2
    end
    
    @testset "Type stability for splattable methods" begin
        u★ = FT(0.5)
        q★ = FT(0.001)
        
        # Test that return types are consistent
        T1 = SF.surface_temperature(param_set, sc)
        T2 = SF.surface_temperature(param_set, sc, u★)
        T3 = SF.surface_temperature(param_set, sc, u★, q★)
        @test typeof(T1) == typeof(T2) == typeof(T3) == FT
        
        q1 = SF.surface_specific_humidity(param_set, sc)
        q2 = SF.surface_specific_humidity(param_set, sc, u★)
        q3 = SF.surface_specific_humidity(param_set, sc, u★, q★)
        @test typeof(q1) == typeof(q2) == typeof(q3) == FT
    end
    
    @testset "surface_virtual_pottemp with args = nothing" begin
        # Base method with no args
        θ_base = SF.surface_virtual_pottemp(param_set, sc)
        @test θ_base isa FT
        
        # Base method with args = nothing
        θ_nothing = SF.surface_virtual_pottemp(param_set, sc, nothing)
        @test θ_nothing ≈ θ_base
        
        # Verify that θᵥ_sfc and surface_virtual_pottemp are the same when args = nothing
        θᵥ_sfc_result = SF.θᵥ_sfc(param_set, sc)
        @test θᵥ_sfc_result ≈ θ_base
        @test θᵥ_sfc_result ≈ θ_nothing
    end
end

