# Regression test for Brusselator autodiff stability
# This test ensures that ForwardDiff doesn't fail with NaN errors
# and that the performance is reasonable compared to FiniteDiff

using Test
using OrdinaryDiffEq, ForwardDiff, FiniteDiff, SciMLSensitivity
using LinearAlgebra

# Include the makebrusselator function from BrussScaling.jmd
include("BrussScaling.jmd")

function auto_sen_l2_test(f, u0, tspan, p, t, alg=Tsit5(); diffalg=ForwardDiff.gradient, kwargs...)
    test_f(p) = begin
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(f,convert.(eltype(p),u0),tspan,p)
        sol = solve(prob,alg,saveat=t; kwargs...)
        sum(sol.u) do x
            sum(z->(1-z)^2/2, x)
        end
    end
    diffalg(test_f, p)
end

@testset "Brusselator Autodiff Regression Tests" begin
    
    # Test parameters - use smaller sizes for faster testing
    bt = 0:0.5:1.0  # Reduced time points
    tspan = (0.0, 1.0)
    tols = (abstol=1e-5, reltol=1e-7)
    
    # Test small size first (N=4)
    @testset "Small Brusselator (N=4)" begin
        N = 4
        bfun, b_u0, b_p, brusselator_jac = makebrusselator(N)
        
        # Test that ForwardDiff doesn't throw errors
        @testset "ForwardDiff Stability" begin
            local result_fd
            @test_nowarn begin
                result_fd = auto_sen_l2_test(bfun, b_u0, tspan, b_p, bt, Rodas5(); 
                                           diffalg=ForwardDiff.gradient, tols...)
            end
            
            # Ensure result is finite (not NaN or Inf)
            @test all(isfinite, result_fd)
            @test !any(isnan, result_fd)
        end
        
        # Test that FiniteDiff works as baseline
        @testset "FiniteDiff Baseline" begin
            local result_finite
            @test_nowarn begin
                result_finite = auto_sen_l2_test(bfun, b_u0, tspan, b_p, bt, Rodas5(); 
                                                diffalg=FiniteDiff.finite_difference_gradient, tols...)
            end
            
            @test all(isfinite, result_finite)
            @test !any(isnan, result_finite)
        end
    end
    
    # Test medium size (N=6) 
    @testset "Medium Brusselator (N=6)" begin
        N = 6
        bfun, b_u0, b_p, brusselator_jac = makebrusselator(N)
        
        # Test ForwardDiff doesn't have SingularException or NaN issues
        @testset "ForwardDiff Medium Scale" begin
            local result_fd
            @test_nowarn begin
                result_fd = auto_sen_l2_test(bfun, b_u0, tspan, b_p, bt, Rodas5(); 
                                           diffalg=ForwardDiff.gradient, tols...)
            end
            
            @test all(isfinite, result_fd)
            @test !any(isnan, result_fd)
        end
    end
    
    # Performance regression test
    @testset "Performance Sanity Check" begin
        N = 4
        bfun, b_u0, b_p, brusselator_jac = makebrusselator(N)
        
        # Time both methods
        t_fd = @elapsed auto_sen_l2_test(bfun, b_u0, tspan, b_p, bt, Rodas5(); 
                                       diffalg=ForwardDiff.gradient, tols...)
        
        t_finite = @elapsed auto_sen_l2_test(bfun, b_u0, tspan, b_p, bt, Rodas5(); 
                                           diffalg=FiniteDiff.finite_difference_gradient, tols...)
        
        # ForwardDiff should be competitive (not 10x+ slower)
        # This is a sanity check to catch severe performance regressions
        ratio = t_finite / t_fd
        @test ratio > 0.1  # ForwardDiff shouldn't be more than 10x slower than FiniteDiff
        
        @info "Performance ratio (FiniteDiff/ForwardDiff): $ratio"
        @info "ForwardDiff time: $t_fd seconds"
        @info "FiniteDiff time: $t_finite seconds"
    end
    
    # Dependency compatibility test
    @testset "Package Loading" begin
        @test_nowarn using OrdinaryDiffEq
        @test_nowarn using ForwardDiff
        @test_nowarn using FiniteDiff
        @test_nowarn using SciMLSensitivity
        @test_nowarn using LinearAlgebra
        
        # Test that problematic methods are not accessible
        @test_throws UndefVarError MooncakeVJP()
    end
end

println("All regression tests completed successfully!")