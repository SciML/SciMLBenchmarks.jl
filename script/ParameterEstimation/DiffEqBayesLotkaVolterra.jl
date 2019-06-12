
using DiffEqBayes


using Distributions
using OrdinaryDiffEq, RecursiveArrayTools, ParameterizedFunctions, DiffEqUncertainity
using Plots
using DiffEqMonteCarlo


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


cb = AdaptiveProbIntsUncertainty(5)
monte_prob = MonteCarloProblem(prob)
sim = solve(monte_prob,Tsit5(),num_monte=100,callback=cb,reltol=1e-5,abstol=1e-5)
plot(sim,vars=(0,1),linealpha=0.4)


@time bayesian_result_stan = stan_inference(prob,t,data,priors;reltol=1e-5,abstol=1e-5,vars =(StanODEData(),InverseGamma(3,2)))


plot_chain(bayesian_result_stan)


@time bayesian_result_turing = turing_inference(prob,Tsit5(),t,data,priors)


plot_chain(bayesian_result_turing)

