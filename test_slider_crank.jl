#!/usr/bin/env julia
# Test script for slider_crank.jmd

using OrdinaryDiffEq, Sundials, LinearAlgebra

try
    println("="^70)
    println("TESTING SLIDER CRANK DAE SOLVER")
    println("="^70)
    
    # Load benchmark code
    include("benchmarks/DAE/slider_crank.jmd")
    
    println("\n" * "="^70)
    println("✓ CODE EXECUTED SUCCESSFULLY!")
    println("="^70)
    println("\nVerifying solution...")
    println("Final time step has $(length(sol.t)) time points")
    println("Final state vector has $(length(sol.u[end])) components")
    
    # Extract reference component
    x3_final = sol.u[end][3]
    println("\nFinal state (t=$(sol.t[end])):")
    println("  x₃(0.1) = $x3_final m")
    println("  Reference (Table II.19.2) ≈ 0.169737 m")
    println("  Error ≈ $(abs(x3_final - 0.169737)) m")
    
    lambda2_final = sol.u[end][23]
    println("  λ₂(0.1) = $lambda2_final")
    
    # Check constraint residuals
    println("\nConstraint verification (should be near zero):")
    residual = similar(sol.u[end])
    slider_crank_resmbs!(residual, sol.u[end], similar(sol.u[end]), sol.t[end])
    println("  Constraint residuals (rows 22-24):")
    println("    |residual[22]| = $(abs(residual[22]))")
    println("    |residual[23]| = $(abs(residual[23]))")
    println("    |residual[24]| = $(abs(residual[24]))")
    
catch e
    println("\n" * "="^70)
    println("✗ ERROR OCCURRED")
    println("="^70)
    println("\n$e\n")
    println(stacktrace(catch_backtrace()))
    exit(1)
end
