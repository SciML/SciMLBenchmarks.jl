
using DiffEqBayes
using Distributions
using OrdinaryDiffEq, RecursiveArrayTools, ParameterizedFunctions
using Plots


gr(fmt=:png)


g1 = @ode_def LorenzExample begin
  dx = σ*(y-x)
  dy = x*(ρ-z) - y
  dz = x*y - β*z
end σ ρ β


r0 = [1.0; 0.0; 0.0]
tspan = (0.0, 30.0)
p = [10.0,28.0,2.66]


prob = ODEProblem(g1,r0,tspan,p)
sol = solve(prob,Tsit5())


t = collect(range(1,stop=30,length=30))
sig = 0.49
data = convert(Array, VectorOfArray([(sol(t[i]) + sig*randn(3)) for i in 1:length(t)]))


Plots.scatter(t, data[1,:],markersize=4,color=:purple)
Plots.scatter!(t, data[2,:],markersize=4,color=:yellow)
Plots.scatter!(t, data[3,:],markersize=4,color=:black)
plot!(sol)


cb = AdaptiveProbIntsUncertainty(5)
monte_prob = MonteCarloProblem(prob)
sim = solve(monte_prob,Tsit5(),num_monte=100,callback=cb,reltol=1e-5,abstol=1e-5)
plot(sim,vars=(0,1),linealpha=0.4)


cb = AdaptiveProbIntsUncertainty(5)
monte_prob = MonteCarloProblem(prob)
sim = solve(monte_prob,Tsit5(),num_monte=100,callback=cb,reltol=1e-6,abstol=1e-6)
plot(sim,vars=(0,1),linealpha=0.4)


cb = AdaptiveProbIntsUncertainty(5)
monte_prob = MonteCarloProblem(prob)
sim = solve(monte_prob,Tsit5(),num_monte=100,callback=cb,reltol=1e-8,abstol=1e-8)
plot(sim,vars=(0,1),linealpha=0.4)


priors = [Truncated(Normal(10,2),1,15),Truncated(Normal(30,5),1,45),Truncated(Normal(2.5,0.5),1,4)]


@time bayesian_result = stan_inference(prob,t,data,priors;reltol=1e-8,abstol=1e-8,vars=(StanODEData(),InverseGamma(3,2)))


plot_chain(bayesian_result)


@time bayesian_result_turing = turing_inference(prob,Tsit5(),t,data,priors)


plot_chain(bayesian_result_turing)

