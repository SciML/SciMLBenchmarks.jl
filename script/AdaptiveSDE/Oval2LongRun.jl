
using StochasticDiffEq, DiffEqProblemLibrary, Random
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
prob = DiffEqProblemLibrary.SDEProblemLibrary.oval2ModelExample(largeFluctuations=true,useBigs=false)


Random.seed!(250)
prob = remake(prob,tspan = (0.0,500.0))
sol = solve(prob,SRIW1(),dt=(1/2)^(18),progress=true,qmax=1.125,
          saveat=0.1,abstol=1e-5,reltol=1e-3,maxiters=1e7);
Random.seed!(250)
prob = remake(prob,tspan = (0.0,500.0))
@time sol = solve(prob,SRIW1(),dt=(1/2)^(18),progress=true,qmax=1.125,
    saveat=0.1,abstol=1e-5,reltol=1e-3,maxiters=1e7);


println(maximum(sol[:,2]))
using Plots; gr()
lw = 2
lw2 = 3
p1 = plot(sol,vars=(0,16),
          title="(A) Timeseries of Ecad Concentration",xguide="Time (s)",
          yguide="Concentration",guidefont=font(16),tickfont=font(16),
          linewidth=lw,leg=false)


p2 = plot(sol,vars=(0,17),
          title="(B) Timeseries of Vim Concentration",xguide="Time (s)",
          yguide="Concentration",guidefont=font(16),
          tickfont=font(16),linewidth=lw,leg=false)


prob = remake(prob,tspan = (0.0,1.0))
## Little Run
sol = solve(prob,EM(),dt=(1/2)^(20),
          progressbar=true,saveat=0.1)
println("EM")
@time sol = solve(prob,EM(),dt=(1/2)^(20),
          progressbar=true,saveat=0.1)

sol = solve(prob,SRI(),dt=(1/2)^(18),adaptive=false,
            progressbar=true,save_everystep=false)
println("SRI")
@time sol = solve(prob,SRI(),dt=(1/2)^(18),adaptive=false,
          progressbar=true,save_everystep=false)

sol = solve(prob,SRIW1(),dt=(1/2)^(18),adaptive=false,
          adaptivealg=:RSwM3,progressbar=false,qmax=4,saveat=0.1)
println("SRIW1")
@time sol = solve(prob,SRIW1(),dt=(1/2)^(18),adaptive=false,
          adaptivealg=:RSwM3,progressbar=false,qmax=4,saveat=0.1)

sol = solve(prob,SRI(),dt=(1/2)^(18),
          adaptivealg=:RSwM3,progressbar=false,qmax=1.125,
          saveat=0.1,abstol=1e-6,reltol=1e-4)
println("SRI Adaptive")
@time sol = solve(prob,SRI(),dt=(1/2)^(18),
          adaptivealg=:RSwM3,progressbar=false,qmax=1.125,
          saveat=0.1,abstol=1e-6,reltol=1e-4)
@show length(sol.t)

sol = solve(prob,SRIW1(),dt=(1/2)^(18),
          adaptivealg=:RSwM3,progressbar=false,qmax=1.125,
          saveat=0.1,abstol=1e-6,reltol=1e-4)
println("SRIW1 Adaptive")
@time sol = solve(prob,SRIW1(),dt=(1/2)^(18),
          adaptivealg=:RSwM3,progressbar=false,qmax=1.125,
          saveat=0.1,abstol=1e-6,reltol=1e-4)
@show length(sol.t)


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

