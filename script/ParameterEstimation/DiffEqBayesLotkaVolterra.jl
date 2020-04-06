
using DiffEqBayes, CmdStan, DynamicHMC


using Distributions, BenchmarkTools
using OrdinaryDiffEq, RecursiveArrayTools, ParameterizedFunctions
using Plots


gr(fmt=:png)


f = @ode_def LotkaVolterraTest begin
    dx = a*x - b*x*y
    dy = -c*y + d*x*y
end a b c d


u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1,0]


prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob,Tsit5())


t = collect(range(1,stop=10,length=10))
sig = 0.49
data = convert(Array, VectorOfArray([(sol(t[i]) + sig*randn(2)) for i in 1:length(t)]))


scatter(t, data[1,:], lab="#prey (data)")
scatter!(t, data[2,:], lab="#predator (data)")
plot!(sol)


priors = [Truncated(Normal(1.5,0.5),0.5,2.5),Truncated(Normal(1.2,0.5),0,2),Truncated(Normal(3.0,0.5),1,4),Truncated(Normal(1.0,0.5),0,2)]


@btime bayesian_result_stan = stan_inference(prob,t,data,priors,num_samples=10_000,printsummary=false)


@btime bayesian_result_turing = turing_inference(prob,Tsit5(),t,data,priors,num_samples=10_000)


@btime bayesian_result_dynamichmc = dynamichmc_inference(prob,Tsit5(),t,data,priors,num_samples=10_000)


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

