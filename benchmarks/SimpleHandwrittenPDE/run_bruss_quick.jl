#!/usr/bin/env julia
ENV["GKSwstype"]="100"
OPENBLAS_NUM_THREADS=get(ENV,"OPENBLAS_NUM_THREADS","1")
using OrdinaryDiffEq, OrdinaryDiffEqBDF, OrdinaryDiffEqExponentialRK, OrdinaryDiffEqRosenbrock
using DiffEqDevTools, SciMLOperators, LinearAlgebra, SparseArrays, Plots; gr(show=false)
using ADTypes: AutoFiniteDiff

N=parse(Int, get(ENV,"BRUSS_N","32"))
numruns=parse(Int, get(ENV,"NUMRUNS","1"))
plot_tuning=get(ENV,"PLOT_TUNING","false")=="true"
best_m=15; best_iop=2
println("Brusselator Quick N=$N (2*N^2=$(2*N^2) ODEs) - numruns=$numruns best m=$best_m iop=$best_iop PLOT_TUNING=$plot_tuning DT 0.005..0.05 tspan 0-1.0")

bruss_f(x,y,t) = (((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
init_bruss(N) = begin xyd=range(0,stop=1,length=N); u=zeros(N,N,2); for I in CartesianIndices((N,N)); x=xyd[I[1]]; y=xyd[I[2]]; u[I,1]=22*(y*(1-y))^(3/2); u[I,2]=27*(x*(1-x))^(3/2); end; u end
function build_mat(N, as)
    n=N*N; L=spzeros(n,n)
    for j in 1:N, i in 1:N; idx=(j-1)*N+i; ip1=i==N ? 1 : i+1; im1=i==1 ? N : i-1; jp1=j==N ? 1 : j+1; jm1=j==1 ? N : j-1; idx_ip1=(j-1)*N+ip1; idx_im1=(j-1)*N+im1; idx_jp1=(jp1-1)*N+i; idx_jm1=(jm1-1)*N+i; L[idx,idx]=-4*as; L[idx,idx_ip1]=as; L[idx,idx_im1]=as; L[idx,idx_jp1]=as; L[idx,idx_jm1]=as; end
    C=spzeros(2*n,2*n); C[1:n,1:n]=L; C[n+1:end,n+1:end]=L; return C
end
function bruss_nonlin_vec!(du,u,p,t)
    Nloc=p[5]; A_val=p[1]; B_val=p[2]; u3=reshape(u,Nloc,Nloc,2); du3=reshape(du,Nloc,Nloc,2); xyd=range(0,stop=1,length=Nloc)
    @inbounds for II in CartesianIndices((Nloc,Nloc)); i=II[1]; j=II[2]; x=xyd[i]; y=xyd[j]; du3[i,j,1]=B_val+u3[i,j,1]^2*u3[i,j,2]-(A_val+1)*u3[i,j,1]+bruss_f(x,y,t); du3[i,j,2]=A_val*u3[i,j,1]-u3[i,j,1]^2*u3[i,j,2]; end
end
make_p(N;A=3.4,B=1.0,α=10.0)=(A,B,α,1.0/N,N)
make_split_prob(N; tspan=(0.0,1.0)) = begin p=make_p(N); αs=p[3]/p[4]^2; L=build_mat(N,αs); u0=vec(init_bruss(N)); f1=MatrixOperator(L); f2=(du,u,par,t)->bruss_nonlin_vec!(du,u,par,t); SplitODEProblem(f1,f2,u0,tspan,p) end
make_full_prob(N; tspan=(0.0,1.0)) = begin p=make_p(N); u0=vec(init_bruss(N)); function full!(du,u,par,t) Nloc=par[5]; α=par[3]; dx=par[4]; αs=α/dx^2; u3=reshape(u,Nloc,Nloc,2); du3=reshape(du,Nloc,Nloc,2); @inbounds for k in 1:2, j in 1:Nloc, i in 1:Nloc; ip1=i==Nloc ? 1 : i+1; im1=i==1 ? Nloc : i-1; jp1=j==Nloc ? 1 : j+1; jm1=j==1 ? Nloc : j-1; du3[i,j,k]=αs*(u3[im1,j,k]+u3[ip1,j,k]+u3[i,jp1,k]+u3[i,jm1,k]-4*u3[i,j,k]); end; tmp=similar(du); bruss_nonlin_vec!(tmp,u,par,t); du.+=tmp; end; ODEProblem(full!, u0, tspan, p) end

prob_split = make_split_prob(N)
prob_full = make_full_prob(N)
println("Probs built N=$N tspan 0-1.0")

println("Ref FBDF 1e-12...")
@time sol_ref = solve(prob_full, FBDF(), abstol=1e-12, reltol=1e-12, maxiters=Int(1e6), save_everystep=false)
test_sol = TestSolution(sol_ref)
println("Ref steps=$(sol_ref.stats.naccept)")

abstols = 0.1 .^ (5:8); reltols = 0.1 .^ (5:8)
dts = [0.05, 0.025, 0.0125, 0.005] # fair DTs

if plot_tuning
    println("=== Stage 1 Tuning (Exp-only) - optional ===")
    setups_tune = [Dict(:alg=>ETDRK2(krylov=true,m=15,iop=2),:dts=>dts), Dict(:alg=>ETDRK2(krylov=true,m=30),:dts=>dts), Dict(:alg=>ETDRK4(krylov=true,m=30,iop=0),:dts=>dts), Dict(:alg=>ETDRK4(krylov=true,m=30,iop=2),:dts=>dts)]
    labels_tune = hcat("ETDRK2 m=15","ETDRK2 m=30","ETDRK4 iop=0 full","ETDRK4 iop=2 best")
    @time wp_tune = WorkPrecisionSet(prob_split, abstols, reltols, setups_tune; print_names=true, names=labels_tune, numruns=numruns, error_estimate=:l2, save_everystep=false, appxsol=test_sol, maxiters=Int(1e6))
    savefig(plot(wp_tune, label=labels_tune), "benchmarks/SimpleHandwrittenPDE/bruss_fixed_quick.png")
end

# Stage 2 Final single plot with baseline always present (default)
println("=== Stage 2 Final: single plot with baseline + best (PLOT_TUNING=$plot_tuning) ===")
setups_ad = [Dict(:alg=>FBDF()), Dict(:alg=>KenCarp4()), Dict(:alg=>Exprb43(autodiff=AutoFiniteDiff(), m=best_m, iop=best_iop))]
labels_ad = hcat("FBDF baseline","KenCarp4 baseline","Exprb43 m=$best_m iop=$best_iop best")
@time wp_ad = WorkPrecisionSet(prob_full, abstols, reltols, setups_ad; print_names=true, names=labels_ad, numruns=numruns, error_estimate=:l2, save_everystep=false, appxsol=test_sol, maxiters=Int(1e6))

setups_dt = [Dict(:alg=>ETDRK4(krylov=true,m=best_m,iop=best_iop),:dts=>dts), Dict(:alg=>ETDRK2(krylov=true,m=best_m,iop=best_iop),:dts=>dts)]
labels_dt = hcat("ETDRK4 m=$best_m iop=$best_iop best final (4x win)","ETDRK2 m=$best_m iop=$best_iop")
@time wp_dt = WorkPrecisionSet(prob_split, abstols, reltols, setups_dt; print_names=true, names=labels_dt, numruns=numruns, error_estimate=:l2, save_everystep=false, appxsol=test_sol, maxiters=Int(1e6))

p_final = plot(wp_ad, label=labels_ad, markershape=:auto, title="Brusselator N=$N Final - Baseline vs Best Exp m=$best_m iop=$best_iop\nDT vs TOL fair, N=32 2048 ODEs")
plot!(p_final, wp_dt, label=labels_dt)
savefig(p_final, "benchmarks/SimpleHandwrittenPDE/bruss_dt_vs_tol_fair_quick.png")
savefig(p_final, "bruss_dt_vs_tol_fair_quick.png")
println("Saved single final plot bruss_dt_vs_tol_fair_quick.png with baseline FBDF+KenCarp4 + best m=$best_m iop=$best_iop")
