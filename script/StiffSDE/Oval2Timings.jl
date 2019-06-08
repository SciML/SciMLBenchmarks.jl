
using Distributed
addprocs()

@everywhere begin
  using StochasticDiffEq, DiffEqProblemLibrary, ParallelDataTransfer, Random, Distributed
  using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
  prob = DiffEqProblemLibrary.SDEProblemLibrary.oval2ModelExample(largeFluctuations=true,useBigs=false)
  Random.seed!(99 + myid())
  prob = remake(prob,tspan=(0.0,1.0))
  println("Solve once to compile.")
  sol = solve(prob,EM(),dt=1/2^(18),adaptive=false,save_everystep=false)
  sol = solve(prob,RKMil(),dt=1/2^(18),adaptive=false,save_everystep=false)
  sol = solve(prob,SRIW1(),dt=1/2^(18),adaptive=false,save_everystep=false)
  sol = solve(prob,SRI(),dt=1/2^(18),adaptive=false,save_everystep=false)
  sol = solve(prob,SOSRI(),dt=1/2^(18),adaptive=false,save_everystep=false)
  sol = solve(prob,SOSRI2(),dt=1/2^(18),adaptive=false,save_everystep=false)
  Int(sol.u[1]!=NaN)
  println("Compilation complete.")
  js = 16:21
  dts = 1.0 ./ 2.0 .^ (js)
  fails = Array{Int}(undef,length(dts),3)
  times = Array{Float64}(undef,length(dts),3)
  numRuns = 10000
end
println("Setup Complete")


## Timing Runs

@everywhere function runAdaptiveSRIW1(i)
  sol = solve(prob,SRIW1(),abstol=2.0^(-13),reltol=2.0^(-7),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIW1,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRI(i)
    sol = solve(prob,SRI(error_terms=2),abstol=2.0^(-13),reltol=2.0^(-7),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRI,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRI(i)
  sol = solve(prob,SRI(),abstol=2.0^(-14),reltol=2.0^(-18),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRI,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt1(i)
  sol = solve(prob,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),abstol=2.0^(-7),reltol=2.0^(-4),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt1,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt1(i)
  sol = solve(prob,SOSRI(),abstol=2.0^(-7),reltol=2.0^(-4),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt1,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt1(i)
  sol = solve(prob,SOSRI(),abstol=2.0^(-7),reltol=2.0^(-6),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt1,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt1(i)
  sol = solve(prob,SOSRI(),abstol=2.0^(-12),reltol=2.0^(-15),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt1,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt1(i)
  sol = solve(prob,SOSRI(),abstol=2.0^(-13),reltol=2.0^(-7),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt1,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt1(i)
  sol = solve(prob,SOSRI(),abstol=2.0^(-12),reltol=2.0^(-15),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt1,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt2(i)
    sol = solve(prob,SOSRI2(),abstol=2.0^(-12),reltol=2.0^(-15),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
    Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt2,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt2(i)
    sol = solve(prob,SOSRI2(),abstol=2.0^(-13),reltol=2.0^(-11),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
    Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt2,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


## Timing Runs

@everywhere function runAdaptiveSRIOpt2(i)
    sol = solve(prob,SOSRI2(),abstol=2.0^(-16),reltol=2.0^(-9),maxIters=Int(1e11),qmax=1.125,save_everystep=false)
    Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
@everywhere Random.seed!(99 + myid())
adaptiveTime = @elapsed numFails = sum(pmap(runAdaptiveSRIOpt2,1:numRuns))
println("The number of Adaptive Fails is $numFails. Elapsed time was $adaptiveTime")


@everywhere function runEM(i,j)
  sol =solve(prob,EM(),dt=dts[j],maxIters=Int(1e11),save_everystep=false,verbose=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
for j in eachindex(js)
  println("j = $j")
  sendto(workers(), j=j)
  @everywhere Random.seed!(99 + myid())
  t1 = @elapsed numFails = sum(pmap((i)->runEM(i,j),1:numRuns))
  println("The number of Euler-Maruyama Fails is $numFails. Elapsed time was $t1")
  fails[j,1] = numFails
  times[j,1] = t1
end


@everywhere function runSRI(i,j)
  sol =solve(prob,SRIW1(),dt=dts[j],maxIters=Int(1e11),adaptive=false,save_everystep=false,verbose=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
for j in 1:4
  println("j = $j")
  sendto(workers(), j=j)
  @everywhere Random.seed!(99 + myid())
  t2 = @elapsed numFails = sum(pmap((i)->runSRI(i,j),1:numRuns))
  println("The number of Rossler-SRI Fails is $numFails. Elapsed time was $t2")
  fails[j,2] = numFails
  times[j,2] = t2
end


@everywhere js = 17:21
@everywhere dts = 1.0 ./2.0 .^ (js)
@everywhere function runIEM(i,j)
  sol =solve(prob,ImplicitEM(),dt=dts[j],maxIters=Int(1e11),save_everystep=false,verbose=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
for j in 1:6
  println("j = $j")
  sendto(workers(), j=j)
  @everywhere Random.seed!(99 + myid())
  t2 = @elapsed numFails = sum(pmap((i)->runIEM(i,j),1:numRuns))
    println("The number of Implicit-EM Fails is $numFails. Elapsed time was $t2")
  fails[j,2] = numFails
  times[j,2] = t2
end


@everywhere js = 17:21
@everywhere dts = 1.0 ./ 2.0 .^(js)
@everywhere function runIRM(i,j)
    sol =solve(prob,ImplicitRKMil(),dt=dts[j],maxIters=Int(1e11),save_everystep=false,verbose=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
for j in 1:4
  println("j = $j")
  sendto(workers(), j=j)
  @everywhere Random.seed!(99 + myid())
  t2 = @elapsed numFails = sum(pmap((i)->runIRM(i,j),1:numRuns))
    println("The number of Implicit-RKMil Fails is $numFails. Elapsed time was $t2")
  fails[j,2] = numFails
  times[j,2] = t2
end


@everywhere function runMil(i,j)
  sol =solve(prob,RKMil(),dt=dts[j],maxIters=Int(1e11),save_everystep=false,verbose=false)
  Int(any(isnan,sol[end]) || sol.t[end] != 1)
end
for j in eachindex(js)
  println("j = $j")
  sendto(workers(), j=j)
  @everywhere Random.seed!(99 + myid())
  t3 = @elapsed numFails = sum(pmap((i)->runMil(i,j),1:numRuns))
  println("The number of RK-Milstein Fails is $numFails. Elapsed time was $t3")
  fails[j,3] = numFails
  times[j,3] = t3
end


using Plots
lw = 3
p2 = plot(dts,times,xscale=:log2,yscale=:log2,guidefont=font(16),tickfont=font(14),yguide="Elapsed Time (s)",xguide=L"Chosen $\Delta t$",top_margin=50px,linewidth=lw,lab=["Euler-Maruyama" "RK-Mil" "RosslerSRI"],legendfont=font(14))
plot!(dts,repmat([adaptiveTime],11),linewidth=lw,line=:dash,lab="ESRK+RSwM3",left_margin=75px)
scatter!([2.0^(-20);2.0^(-20);2.0^(-18)],[times[5,1];times[5,2];times[3,3]],markersize=20,c=:red,lab="")
plot(p2,size=(800,800))


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

