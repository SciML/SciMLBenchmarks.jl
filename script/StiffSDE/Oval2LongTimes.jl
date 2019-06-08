
using StochasticDiffEq, DiffEqProblemLibrary, Random
Random.seed!(200)
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
prob = DiffEqProblemLibrary.SDEProblemLibrary.oval2ModelExample(largeFluctuations=true,useBigs=false)


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SRIW1(),dt=(1/2)^(18),qmax=1.125,
        saveat=0.1,maxiters=1e7,abstol=1e-5,reltol=1e-3)
end


Random.seed!(200)
@time for i in 1:10
    @show i
    sol = solve(prob,ImplicitEM(),dt=1/60000)
end


Random.seed!(200)
@time for i in 1:10
    @show i
    sol = solve(prob,ImplicitRKMil(),dt=1/50000)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-4,reltol=1e-2)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI2(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-4,reltol=1e-4)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI2(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-5,reltol=1e-3)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-3,reltol=1e-2)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-4,reltol=1e-4)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-2,reltol=1e-2)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-5,reltol=1e-3)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-2,reltol=1e-1)
end


Random.seed!(200)
@time for i in 1:10
    sol = solve(prob,SOSRI2(),dt=(1/2)^(18),qmax=1.125,
          saveat=0.1,maxiters=1e7,abstol=1e-4,reltol=1e-1)
end


Random.seed!(200)
@time for i in 1:10
    @show i
    sol = solve(prob,ImplicitEM(),dt=1/50000)
end


Random.seed!(200)
@time for i in 1:10
    @show i
    sol = solve(prob,ImplicitRKMil(),dt=1/40000)
end


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

