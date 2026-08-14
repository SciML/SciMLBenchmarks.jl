#!/usr/bin/env julia
# Quick feedback for Allen-Cahn - now single final plot with baseline+best by default, tuning optional via PLOT_TUNING
ENV["GKSwstype"] = "100"
OPENBLAS_NUM_THREADS = get(ENV,"OPENBLAS_NUM_THREADS","1")

using OrdinaryDiffEq, OrdinaryDiffEqBDF, OrdinaryDiffEqExponentialRK
using DiffEqDevTools, SciMLOperators, LinearAlgebra, SparseArrays, Plots
using SummationByPartsOperators
const SBP=SummationByPartsOperators
using ADTypes: AutoFiniteDiff
gr(show=false)

function forcing_term(du,u,p,t); du .= p .* @. (u - u^3); du[1]=0.0; du[end]=0.0; end
function allen_cahn(N,L)
    eps=1e-3; D2=derivative_operator(MattssonSvärdNordström2004(); derivative_order=2, accuracy_order=2, xmin=-L, xmax=L, N=N)
    x=LinRange(-L,L,N); u0=@. cos(2π*x); p=3.0
    SplitODEProblem(MatrixOperator(eps*sparse(D2)), forcing_term, u0, (0.0,1.0), p) |> (prob->(x,prob))
end
N=256; L=2.0; xs, prob = allen_cahn(N,L)
struct CountingMatrix <: AbstractMatrix{Float64}; A::SparseMatrixCSC{Float64,Int}; count::Base.RefValue{Int}; end
CountingMatrix(A)=CountingMatrix(A,Ref(0)); Base.size(C::CountingMatrix)=size(C.A); Base.getindex(C::CountingMatrix,i,j)=getindex(C.A,i,j)
LinearAlgebra.mul!(y,C::CountingMatrix,x)=(C.count[]+=1; mul!(y,C.A,x)); LinearAlgebra.mul!(y,C::CountingMatrix,x,a,b)=(C.count[]+=1; mul!(y,C.A,x,a,b))
Base.:*(C::CountingMatrix,x::AbstractVector)=(C.count[]+=1; C.A*x)
LinearAlgebra.mul!(y::AbstractVector, C::CountingMatrix, x::AbstractVector, a::Number, b::Number)=(C.count[]+=1; mul!(y,C.A,x,a,b))

function full_branched_f!(du,u,p,t); prob.f.f1(du,u,p,t); tmp=similar(du); prob.f.f2(tmp,u,p,t); @. du+=tmp; du[1]=0.0; du[end]=0.0; end
prob_full = ODEProblem(full_branched_f!, prob.u0, prob.tspan, prob.p)

println("Reference solve N=$N L=$L FBDF 1e-12...")
@time sol = solve(prob_full, FBDF(autodiff=AutoFiniteDiff()), abstol=1e-12, reltol=1e-12, maxiters=Int(1e6), save_everystep=false)
test_sol = TestSolution(sol); println("Ref done steps=$(sol.stats.naccept)")

numruns = parse(Int, get(ENV,"NUMRUNS","1"))
plot_tuning = get(ENV,"PLOT_TUNING","false")=="true"
best_m=15; best_iop=2
println("Best from prior tuning data: m=$best_m iop=$best_iop (PLOT_TUNING=$plot_tuning) - m=15 31% matvec save (10.4 vs 30), iop=2 15x orth save")

if plot_tuning
    println("=== Stage 1 Tuning: m & IOP sweep (Exp-only, selects best) ===")
    abstols=0.1 .^ (5:8); reltols=0.1 .^ (1:4); multipliers=[0.05,0.025,0.0125,0.005]
    setups=[Dict(:alg=>ETDRK2(krylov=true,m=15),:dts=>multipliers), Dict(:alg=>ETDRK2(krylov=true,m=30),:dts=>multipliers), Dict(:alg=>ETDRK2(krylov=true,m=60),:dts=>multipliers),
            Dict(:alg=>ETDRK4(krylov=true,m=30,iop=0),:dts=>multipliers), Dict(:alg=>ETDRK4(krylov=true,m=30,iop=2),:dts=>multipliers), Dict(:alg=>ETDRK4(krylov=true,m=30,iop=10),:dts=>multipliers)]
    labels=hcat("ETDRK2 m=15","ETDRK2 m=30","ETDRK2 m=60","ETDRK4 m=30 iop=0 full","ETDRK4 m=30 iop=2 best","ETDRK4 m=30 iop=10")
    @time wp_tune = WorkPrecisionSet(prob, abstols, reltols, setups; print_names=true, names=labels, numruns=numruns, error_estimate=:l2, save_everystep=false, appxsol=test_sol, maxiters=Int(1e5))
    savefig(plot(wp_tune, label=labels, markershape=:auto, title="Stage1 Tuning m&IOP N=$N selects m=$best_m iop=$best_iop"), "benchmarks/SimpleHandwrittenPDE/m_iop_sweep.png")
    println("Saved m_iop_sweep.png")
else
    println("Skipping tuning plots (PLOT_TUNING=false) - using best m=$best_m iop=$best_iop")
end

println("=== Stage 2 Final: single plot with baseline + best Exp (addresses reviewer) ===")
abstols_ad=0.1 .^ (5:8); reltols_ad=0.1 .^ (2:5)
setups_ad=[Dict(:alg=>FBDF()), Dict(:alg=>Exprb32(autodiff=AutoFiniteDiff(), m=best_m, iop=best_iop)), Dict(:alg=>Exprb43(autodiff=AutoFiniteDiff(), m=best_m, iop=best_iop))]
labels_ad=hcat("FBDF (Newton baseline)","Exprb32 m=$best_m iop=$best_iop best","Exprb43 m=$best_m iop=$best_iop best")
@time wp_ad = WorkPrecisionSet(prob_full, abstols_ad, reltols_ad, setups_ad; print_names=true, names=labels_ad, numruns=numruns, error_estimate=:l2, save_everystep=false, appxsol=test_sol, maxiters=Int(1e6))

abstols_dt=0.1 .^ (5:8); reltols_dt=0.1 .^ (1:4); multipliers=[0.05,0.025,0.0125,0.005]
setups_dt=[Dict(:alg=>ETDRK2(krylov=true,m=best_m,iop=best_iop),:dts=>multipliers), Dict(:alg=>ETDRK4(krylov=true,m=best_m,iop=best_iop),:dts=>multipliers), Dict(:alg=>EPIRK4s3A(adaptive_krylov=true,m=best_m),:dts=>multipliers)]
labels_dt=hcat("ETDRK2 m=$best_m iop=$best_iop best","ETDRK4 m=$best_m iop=$best_iop best final","EPIRK4s3A m=$best_m block")
@time wp_dt = WorkPrecisionSet(prob, abstols_dt, reltols_dt, setups_dt; print_names=true, names=labels_dt, numruns=numruns, error_estimate=:l2, save_everystep=false, appxsol=test_sol, maxiters=Int(1e6))

p_final = plot(wp_ad, label=labels_ad, markershape=:auto, title="Allen-Cahn N=$N Final - Baseline vs Best Exp m=$best_m iop=$best_iop\nDT vs TOL fair, spatial floor 1e-3")
plot!(p_final, wp_dt, label=labels_dt, markershape=:auto)
savefig(p_final, "benchmarks/SimpleHandwrittenPDE/dt_vs_tol_fair.png")
savefig(p_final, "dt_vs_tol_fair.png")
println("Saved single final plot dt_vs_tol_fair.png with baseline FBDF + best m=$best_m iop=$best_iop - addresses reviewer feedback")
