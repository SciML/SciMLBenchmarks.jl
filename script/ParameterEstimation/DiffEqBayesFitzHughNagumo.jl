
using DiffEqBayes, BenchmarkTools


using OrdinaryDiffEq, RecursiveArrayTools, Distributions, ParameterizedFunctions, Mamba, BenchmarkTools
using Plots


gr(fmt=:png)


fitz = @ode_def FitzhughNagumo begin
  dv = v - v^3/3 -w + l
  dw = τinv*(v +  a - b*w)
end a b τinv l


prob_ode_fitzhughnagumo = ODEProblem(fitz,[1.0,1.0],(0.0,10.0),[0.7,0.8,1/12.5,0.5])
sol = solve(prob_ode_fitzhughnagumo, Tsit5())


t = collect(range(1,stop=10,length=10))
sig = 0.20
data = convert(Array, VectorOfArray([(sol(t[i]) + sig*randn(2)) for i in 1:length(t)]))


scatter(t, data[1,:])
scatter!(t, data[2,:])
plot!(sol)


priors = [Truncated(Normal(1.0,0.5),0,1.5),Truncated(Normal(1.0,0.5),0,1.5),Truncated(Normal(0.0,0.5),-0.5,0.5),Truncated(Normal(0.5,0.5),0,1)]


@time bayesian_result_stan = stan_inference(prob_ode_fitzhughnagumo,t,data,priors;reltol=1e-5,abstol=1e-5,vars =(StanODEData(),InverseGamma(3,2)))


plot_chain(bayesian_result_stan)


@time bayesian_result_turing = turing_inference(prob_ode_fitzhughnagumo,Tsit5(),t,data,priors)


plot_chain(bayesian_result_turing)


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

