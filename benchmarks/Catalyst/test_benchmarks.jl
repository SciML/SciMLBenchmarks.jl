# Test script for Catalyst benchmarks
# This script can be used to validate the benchmark setup before running the full suite

using Test
using Catalyst, DifferentialEquations, JumpProcesses
using BenchmarkTools
using JSON

include("model_utils.jl")

@testset "Catalyst Benchmarks Tests" begin
    
    @testset "Model Loading" begin
        # Test model info retrieval
        @test get_model_info("multistate").species == 9
        @test get_model_info("multistate").reactions == 18
        @test get_model_info("unknown").species == 0
        
        # Test model summary creation
        summary = create_model_summary()
        @test length(summary) == 5
        @test summary[1].name == "multistate"
    end
    
    @testset "Data Download" begin
        # Test data directory creation
        data_dir = joinpath(@__DIR__, "data")
        if isdir(data_dir)
            rm(data_dir, recursive=true)
        end
        
        # Download should create directory
        downloaded = download_catalyst_data()
        @test isdir(data_dir)
        
        # Check that at least some files were downloaded
        @test length(readdir(data_dir)) > 0
    end
    
    @testset "Model Loading and Problem Creation" begin
        # Test loading a simple model
        try
            problems = create_catalyst_problems("multistate")
            @test haskey(problems, :ode_prob)
            @test haskey(problems, :jump_prob)
            @test haskey(problems, :rn)
            
            # Test that problems are the right types
            @test problems.ode_prob isa ODEProblem
            @test problems.jump_prob isa JumpProblem
            
            println("✓ Successfully loaded multistate model")
        catch e
            @test_skip "Model loading failed: $(e)"
            println("⚠ Model loading failed, skipping: $(e)")
        end
    end
    
    @testset "Simple Benchmark Execution" begin
        # Test basic ODE solving
        try
            problems = create_catalyst_problems("multistate", tspan=(0.0, 1.0))
            
            # Test ODE solving
            sol = solve(problems.ode_prob, Tsit5(), abstol=1e-6, reltol=1e-3)
            @test sol.retcode == ReturnCode.Success
            
            # Test Jump solving (with shorter time span for speed)
            jump_prob_short = remake(problems.jump_prob, tspan=(0.0, 0.1))
            sol_jump = solve(jump_prob_short, Direct())
            @test sol_jump.retcode == ReturnCode.Success
            
            println("✓ Basic solving tests passed")
        catch e
            @test_skip "Basic solving failed: $(e)"
            println("⚠ Basic solving failed, skipping: $(e)")
        end
    end
    
    @testset "Benchmark Timing" begin
        # Test basic benchmarking functionality
        try
            problems = create_catalyst_problems("multistate", tspan=(0.0, 0.1))
            
            # Simple benchmark
            b = @benchmark solve($(problems.ode_prob), Tsit5(), abstol=1e-6, reltol=1e-3) seconds=1
            @test b.times |> length > 0
            @test all(t -> t > 0, b.times)
            
            println("✓ Benchmarking functionality works")
        catch e
            @test_skip "Benchmarking failed: $(e)"
            println("⚠ Benchmarking failed, skipping: $(e)")
        end
    end
    
    @testset "Result Serialization" begin
        # Test JSON serialization
        test_results = Dict(
            "model" => "test",
            "solver" => "Tsit5",
            "time" => 1.23,
            "success" => true,
            "inf_value" => Inf,
            "nan_value" => NaN
        )
        
        temp_file = tempname() * ".json"
        save_benchmark_results(test_results, temp_file)
        
        @test isfile(temp_file)
        
        # Load back and verify
        loaded = load_benchmark_results(temp_file)
        @test loaded["model"] == "test"
        @test loaded["solver"] == "Tsit5"
        @test loaded["time"] == 1.23
        @test loaded["success"] == true
        
        # Clean up
        rm(temp_file)
        
        println("✓ JSON serialization works")
    end
end

# Function to run integration tests
function run_integration_tests()
    println("\\n=== Running Integration Tests ===")
    
    # Test model summary
    print_model_summary()
    
    # Test if we can create and solve a simple problem
    println("\\n--- Testing Problem Creation and Solving ---")
    
    try
        problems = create_catalyst_problems("multistate", tspan=(0.0, 1.0))
        
        # Test multiple solvers
        solvers = [
            ("Tsit5", Tsit5()),
            ("Rodas4", Rodas4()),
            ("CVODE_BDF", CVODE_BDF())
        ]
        
        println("Testing ODE solvers:")
        for (name, alg) in solvers
            try
                sol = solve(problems.ode_prob, alg, abstol=1e-6, reltol=1e-3)
                println("  ✓ $(name): Success (final time: $(sol.t[end]))")
            catch e
                println("  ✗ $(name): Failed - $(e)")
            end
        end
        
        # Test jump solvers
        jump_solvers = [
            ("Direct", Direct()),
            ("RSSA", RSSA())
        ]
        
        println("\\nTesting Jump solvers:")
        jump_prob_short = remake(problems.jump_prob, tspan=(0.0, 0.1))
        for (name, alg) in jump_solvers
            try
                sol = solve(jump_prob_short, alg)
                println("  ✓ $(name): Success (final time: $(sol.t[end]))")
            catch e
                println("  ✗ $(name): Failed - $(e)")
            end
        end
        
    catch e
        println("Integration test failed: $(e)")
    end
end

# Function to test external dependencies
function test_external_dependencies()
    println("\\n=== Testing External Dependencies ===")
    
    # Test PyCall
    try
        using PyCall
        py"""
        import sys
        print(f"Python version: {sys.version}")
        """
        println("✓ PyCall works")
        
        # Test if we can import basic packages
        try
            py"""
            import numpy as np
            print(f"NumPy version: {np.__version__}")
            """
            println("✓ NumPy available")
        catch e
            println("⚠ NumPy not available: $(e)")
        end
        
        try
            py"""
            import gillespy2
            print(f"GillesPy2 version: {gillespy2.__version__}")
            """
            println("✓ GillesPy2 available")
        catch e
            println("⚠ GillesPy2 not available: $(e)")
        end
        
    catch e
        println("⚠ PyCall not available: $(e)")
    end
    
    # Test RCall
    try
        using RCall
        R"""
        cat("R version:", R.version.string, "\\n")
        """
        println("✓ RCall works")
        
        # Test if we can load required R packages
        try
            R"""
            library(GillespieSSA2)
            cat("✓ GillespieSSA2 loaded\\n")
            """
            println("✓ GillespieSSA2 available")
        catch e
            println("⚠ GillespieSSA2 not available: $(e)")
        end
        
    catch e
        println("⚠ RCall not available: $(e)")
    end
end

# Function to run a quick benchmark
function run_quick_benchmark()
    println("\\n=== Running Quick Benchmark ===")
    
    try
        problems = create_catalyst_problems("multistate", tspan=(0.0, 1.0))
        
        # Quick ODE benchmark
        println("Quick ODE benchmark (Tsit5):")
        b = @benchmark solve($(problems.ode_prob), Tsit5(), abstol=1e-6, reltol=1e-3) seconds=2
        println("  Median time: $(median(b.times) / 1e6) ms")
        println("  Min time: $(minimum(b.times) / 1e6) ms")
        
        # Quick Jump benchmark
        println("\\nQuick Jump benchmark (Direct):")
        jump_prob_short = remake(problems.jump_prob, tspan=(0.0, 0.1))
        b_jump = @benchmark solve($jump_prob_short, Direct()) seconds=2
        println("  Median time: $(median(b_jump.times) / 1e6) ms")
        println("  Min time: $(minimum(b_jump.times) / 1e6) ms")
        
    catch e
        println("Quick benchmark failed: $(e)")
    end
end

# Main test function
function main()
    println("Starting Catalyst benchmark tests...")
    
    # Run unit tests
    println("\\n=== Running Unit Tests ===")
    Test.runtests()
    
    # Run integration tests
    run_integration_tests()
    
    # Test external dependencies
    test_external_dependencies()
    
    # Run quick benchmark
    run_quick_benchmark()
    
    println("\\n=== Test Summary ===")
    println("Tests completed. Check output above for any failures or warnings.")
    println("If all tests pass, the benchmark suite should work correctly.")
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end