
using OrdinaryDiffEq, ODE, ODEInterfaceDiffEq, DiffEqDevTools
tspan = (0.0,1.0)
# Linear ODE
f = (u,p,t) -> (1.01*u)
(::typeof(f))(::Type{Val{:analytic}},u₀,p,t) = u₀*exp(1.01*t)
probnum = ODEProblem(f,1/2,tspan)

# 2D Linear ODE
f = (du,u,p,t) -> begin
  @inbounds for i in eachindex(u)
    du[i] = 1.01*u[i]
  end
end
(::typeof(f))(::Type{Val{:analytic}},u₀,p,t) = u₀*exp(1.01*t)
prob = ODEProblem(f,rand(100,100),tspan)
using Plots; gr()


setups = [Dict(:alg=>RK4());Dict(:alg=>RK4(),:dense=>false);Dict(:alg=>RK4(),:dense=>false,:save_everystep=>false);Dict(:alg=>ode4())]
names = ["RK4 Continuous";"RK4 Non-Dense";"RK4 Save End";"ODE.jl rk4"]
shoot = Shootout(probnum,setups;dt=1/2^(6),names=names);


@show shoot.times # Times
@show shoot.errors # Errors
@show shoot.effs # Efficiencies
@show shoot.effratios[shoot.bestidx,:] # Efficiencies


plot(shoot)


setups = [Dict(:alg=>RK4());Dict(:alg=>RK4(),:dense=>false);Dict(:alg=>RK4(),:dense=>false,:save_everystep=>false);Dict(:alg=>ode4())]
names = ["RK4 Continuous";"RK4 Non-Dense";"RK4 Save End";"ode4"]
shoot = Shootout(prob,setups;dt=1/2^(6),names=names)
println(shoot.effratios[shoot.bestidx,:])
plot(shoot)


setups = [Dict(:alg=>DP5())
          Dict(:abstol=>1e-3,:reltol=>1e-6,:alg=>ode45()); # Fix ODE to be normal
          Dict(:alg=>dopri5())]
shoot = Shootout(probnum,setups)
println(shoot.times)
println(shoot.errors)
println(shoot.effratios[shoot.bestidx,:])
plot(shoot)


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

