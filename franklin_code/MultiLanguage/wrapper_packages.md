+++
author = "Chris Rackauckas"
title = "ODE Solver Multi-Language Wrapper Package Work-Precision Benchmarks (MATLAB, SciPy, Julia, deSolve (R))"
+++


The following benchmarks demonstrate the performance differences due to using
similar algorithms from wrapper packages in the main scripting languages across
a range of stiff and non-stiff ODEs. It takes into account solver time and
error in order to ensure correctness of interpretations. These were ran with
Julia 1.7, MATLAB 2019B, deSolve 1.3.0, and SciPy 1.6.1.

These benchmarks are generated using the following bindings:

- [MATLABDiffEq.jl](https://github.com/JuliaDiffEq/MATLABDiffEq.jl) (MATLAB)
- [SciPyDiffEq.jl](https://github.com/JuliaDiffEq/SciPyDiffEq.jl) (SciPy)
- [deSolveDiffEq.jl](https://github.com/JuliaDiffEq/deSolveDiffEq.jl) (deSolve)
- [OrdinaryDiffEq.jl](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl) (OrdinaryDiffEq.jl)
- [Sundials.jl](https://github.com/JuliaDiffEq/Sundials.jl) (Sundials)
- [ODEInterfaceDiffEq.jl](https://github.com/JuliaDiffEq/ODEInterfaceDiffEq.jl) (Hairer and Netlib)

The respective repos verify negligible overhead on interop (MATLAB, ODEInterface,
and Sundials overhead are negligable, SciPy is accelerated 3x over SciPy+Numba
setups due to the Julia JIT on the ODE function, deSolve sees a 3x overhead
over the pure-R version). Error and timing is compared together to ensure
the methods are solving to the same accuracy when compared.

More wrappers will continue to be added as necessary.

## Setup

```julia
using ParameterizedFunctions, MATLABDiffEq, OrdinaryDiffEq,
      ODEInterfaceDiffEq, Plots, Sundials, SciPyDiffEq, deSolveDiffEq
using DiffEqDevTools
using LinearAlgebra, StaticArrays
```




#### Non-Stiff Problem 1: Lotka-Volterra

```julia
f = @ode_def_bare LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
p = [1.5,1,3,1]
tspan = (0.0,10.0)
u0 = [1.0,1.0]
prob = ODEProblem(f,u0,tspan,p)
staticprob = ODEProblem{false}(f,SVector{2}(u0),tspan,SVector{4}(p))

sol = solve(prob,Vern7(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

setups = [
          Dict(:alg=>DP5())
          Dict(:alg=>Tsit5())
          Dict(:alg=>Vern7())
          Dict(:prob_choice => 2, :alg=>DP5())
          Dict(:prob_choice => 2, :alg=>Tsit5())
          Dict(:prob_choice => 2, :alg=>Vern7())
          Dict(:alg=>dopri5())
          Dict(:alg=>MATLABDiffEq.ode45())
          Dict(:alg=>MATLABDiffEq.ode113())
          Dict(:alg=>SciPyDiffEq.RK45())
          Dict(:alg=>SciPyDiffEq.LSODA())
          Dict(:alg=>SciPyDiffEq.odeint())
          Dict(:alg=>deSolveDiffEq.lsoda())
          Dict(:alg=>deSolveDiffEq.ode45())
          Dict(:alg=>CVODE_Adams())
  ]

labels = [
  "Julia: DP5"
  "Julia: Tsit5"
  "Julia: Vern7"
  "Julia: DP5 Static"
  "Julia: Tsit5 Static"
  "Julia: Vern7 Static"
  "Hairer: dopri5"
  "MATLAB: ode45"
  "MATLAB: ode113"
  "SciPy: RK45"
  "SciPy: LSODA"
  "SciPy: odeint"
  "deSolve: lsoda"
  "deSolve: ode45"
  "Sundials: Adams"
  ]

abstols = 1.0 ./ 10.0 .^ (6:13)
reltols = 1.0 ./ 10.0 .^ (3:10)
wp = WorkPrecisionSet([prob,staticprob],abstols,reltols,setups;
                      names = labels,print_names = true,
                      appxsol=[test_sol,test_sol],dense=false,
                      save_everystep=false,numruns=100,maxiters=10000000,
                      timeseries_errors=false,verbose=false)
plot(wp,title="Non-stiff 1: Lotka-Volterra",legend=:outertopleft,
     color=permutedims([repeat([:LightGreen],3)...,repeat([:DarkGreen],3)...,
     :Red,repeat([:Orange],2)...,repeat([:Yellow],3)...,
     repeat([:Blue],2)...,:Purple]),size = (800,350),
     xticks = 10.0 .^ (-12:1:5),
     yticks = 10.0 .^ (-6:0.5:5),
     bottom_margin=5Plots.mm)
```

```
Julia: DP5
Julia: Tsit5
Julia: Vern7
Julia: DP5 Static
Julia: Tsit5 Static
Julia: Vern7 Static
Hairer: dopri5
MATLAB: ode45
MATLAB: ode113
SciPy: RK45
SciPy: LSODA
SciPy: odeint
deSolve: lsoda
deSolve: ode45
Sundials: Adams
```


![](figures/wrapper_packages_2_1.png)



#### Non-Stiff Problem 2: Rigid Body

```julia
f = @ode_def_bare RigidBodyBench begin
  dy1  = -2*y2*y3
  dy2  = 1.25*y1*y3
  dy3  = -0.5*y1*y2 + 0.25*sin(t)^2
end
u0 = [1.0;0.0;0.9]
prob = ODEProblem(f,u0,(0.0,100.0))
staticprob = ODEProblem{false}(f,SVector{3}(u0),(0.0,100.0))
sol = solve(prob,Vern7(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

setups = [Dict(:alg=>DP5())
          Dict(:alg=>Tsit5())
          Dict(:alg=>Vern7())
          Dict(:prob_choice => 2, :alg=>DP5())
          Dict(:prob_choice => 2, :alg=>Tsit5())
          Dict(:prob_choice => 2, :alg=>Vern7())
          Dict(:alg=>dopri5())
          Dict(:alg=>MATLABDiffEq.ode45())
          Dict(:alg=>MATLABDiffEq.ode113())
          Dict(:alg=>SciPyDiffEq.RK45())
          Dict(:alg=>SciPyDiffEq.LSODA())
          Dict(:alg=>SciPyDiffEq.odeint())
          Dict(:alg=>deSolveDiffEq.lsoda())
          Dict(:alg=>deSolveDiffEq.ode45())
          Dict(:alg=>CVODE_Adams())
  ]

labels = [
  "Julia: DP5"
  "Julia: Tsit5"
  "Julia: Vern7"
  "Julia: DP5 Static"
  "Julia: Tsit5 Static"
  "Julia: Vern7 Static"
  "Hairer: dopri5"
  "MATLAB: ode45"
  "MATLAB: ode113"
  "SciPy: RK45"
  "SciPy: LSODA"
  "SciPy: odeint"
  "deSolve: lsoda"
  "deSolve: ode45"
  "Sundials: Adams"
  ]

abstols = 1.0 ./ 10.0 .^ (6:13)
reltols = 1.0 ./ 10.0 .^ (3:10)

wp = WorkPrecisionSet([prob,staticprob],abstols,reltols,setups;
                      names = labels,print_names = true,
                      appxsol=[test_sol,test_sol],dense=false,
                      save_everystep=false,numruns=100,maxiters=10000000,
                      timeseries_errors=false,verbose=false)
plot(wp,title="Non-stiff 2: Rigid-Body",legend=:outertopleft,
    color=permutedims([repeat([:LightGreen],3)...,repeat([:DarkGreen],3)...,
    :Red,repeat([:Orange],2)...,repeat([:Yellow],3)...,
    repeat([:Blue],2)...,:Purple]),size = (800,350),
    xticks = 10.0 .^ (-12:1:5),
    yticks = 10.0 .^ (-6:0.5:5),
    bottom_margin=5Plots.mm)
```

```
Julia: DP5
Julia: Tsit5
Julia: Vern7
Julia: DP5 Static
Julia: Tsit5 Static
Julia: Vern7 Static
Hairer: dopri5
MATLAB: ode45
MATLAB: ode113
SciPy: RK45
SciPy: LSODA
SciPy: odeint
deSolve: lsoda
deSolve: ode45
Sundials: Adams
```


![](figures/wrapper_packages_3_1.png)



#### Stiff Problem 1: ROBER

```julia
rober = @ode_def begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃
u0 = [1.0,0.0,0.0]
p = [0.04,3e7,1e4]
prob = ODEProblem(rober,u0,(0.0,1e5),p)
staticprob = ODEProblem{false}(rober,SVector{3}(u0),(0.0,1e5),SVector{3}(p))
sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (7:12)
reltols = 1.0 ./ 10.0 .^ (3:8);

setups = [Dict(:alg=>Rosenbrock23())
          Dict(:alg=>Rodas4())
          Dict(:alg=>Rodas5())
          Dict(:prob_choice => 2, :alg=>Rosenbrock23())
          Dict(:prob_choice => 2, :alg=>Rodas4())
          Dict(:prob_choice => 2, :alg=>Rodas5())
          Dict(:alg=>rodas())
          Dict(:alg=>radau())
          Dict(:alg=>MATLABDiffEq.ode23s())
          Dict(:alg=>MATLABDiffEq.ode15s())
          Dict(:alg=>SciPyDiffEq.LSODA())
          Dict(:alg=>SciPyDiffEq.BDF())
          Dict(:alg=>SciPyDiffEq.odeint())
          Dict(:alg=>deSolveDiffEq.lsoda())
          Dict(:alg=>CVODE_BDF())
          ]

labels = [
  "Julia: Rosenbrock23"
  "Julia: Rodas4"
  "Julia: Rodas5"
  "Julia: Rosenbrock23 Static"
  "Julia: Rodas4 Static"
  "Julia: Rodas5 Static"
  "Hairer: rodas"
  "Hairer: radau"
  "MATLAB: ode23s"
  "MATLAB: ode15s"
  "SciPy: LSODA"
  "SciPy: BDF"
  "SciPy: odeint"
  "deSolve: lsoda"
  "Sundials: CVODE"
  ]

wp = WorkPrecisionSet([prob,staticprob],abstols,reltols,setups;
                      names = labels,print_names = true,
                      dense=false,verbose = false,
                      save_everystep=false,appxsol=[test_sol,test_sol],
                      maxiters=Int(1e5))
plot(wp,title="Stiff 1: ROBER", legend=:outertopleft,
     color=permutedims([repeat([:LightGreen],3)...,repeat([:DarkGreen],3)...,
     :Red,:Red,repeat([:Orange],2)...,repeat([:Yellow],3)...,
     repeat([:Blue],1)...,:Purple]),size = (800,350),
     xticks = 10.0 .^ (-12:1:5),
     yticks = 10.0 .^ (-6:0.5:5),
     bottom_margin=5Plots.mm)
```

```
Julia: Rosenbrock23
Julia: Rodas4
Julia: Rodas5
Julia: Rosenbrock23 Static
Julia: Rodas4 Static
Julia: Rodas5 Static
Hairer: rodas
Hairer: radau
MATLAB: ode23s
MATLAB: ode15s
SciPy: LSODA
SciPy: BDF
SciPy: odeint
deSolve: lsoda
Sundials: CVODE
```


![](figures/wrapper_packages_4_1.png)



#### Stiff Problem 2: HIRES

```julia
f = @ode_def Hires begin
  dy1 = -1.71*y1 + 0.43*y2 + 8.32*y3 + 0.0007
  dy2 = 1.71*y1 - 8.75*y2
  dy3 = -10.03*y3 + 0.43*y4 + 0.035*y5
  dy4 = 8.32*y2 + 1.71*y3 - 1.12*y4
  dy5 = -1.745*y5 + 0.43*y6 + 0.43*y7
  dy6 = -280.0*y6*y8 + 0.69*y4 + 1.71*y5 -
           0.43*y6 + 0.69*y7
  dy7 = 280.0*y6*y8 - 1.81*y7
  dy8 = -280.0*y6*y8 + 1.81*y7
end

u0 = zeros(8)
u0[1] = 1
u0[8] = 0.0057
prob = ODEProblem(f,u0,(0.0,321.8122))
staticprob = ODEProblem{false}(f,SVector{8}(u0),(0.0,321.8122))

sol = solve(prob,Rodas5(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (5:10)
reltols = 1.0 ./ 10.0 .^ (1:6);

setups = [Dict(:alg=>Rosenbrock23())
          Dict(:alg=>Rodas4())
          Dict(:alg=>RadauIIA5())
          Dict(:prob_choice => 2, :alg=>Rosenbrock23())
          Dict(:prob_choice => 2, :alg=>Rodas4())
          Dict(:prob_choice => 2, :alg=>RadauIIA5())
          Dict(:alg=>rodas())
          Dict(:alg=>radau())
          Dict(:alg=>MATLABDiffEq.ode23s())
          Dict(:alg=>MATLABDiffEq.ode15s())
          Dict(:alg=>SciPyDiffEq.LSODA())
          Dict(:alg=>SciPyDiffEq.BDF())
          Dict(:alg=>SciPyDiffEq.odeint())
          Dict(:alg=>deSolveDiffEq.lsoda())
          Dict(:alg=>CVODE_BDF())
          ]

labels = [
  "Julia: Rosenbrock23"
  "Julia: Rodas4"
  "Julia: radau"
  "Julia: Rosenbrock23 Static"
  "Julia: Rodas4 Static"
  "Julia: radau Static"
  "Hairer: rodas"
  "Hairer: radau"
  "MATLAB: ode23s"
  "MATLAB: ode15s"
  "SciPy: LSODA"
  "SciPy: BDF"
  "SciPy: odeint"
  "deSolve: lsoda"
  "Sundials: CVODE"
  ]

wp = WorkPrecisionSet([prob,staticprob],abstols,reltols,setups;
                      names = labels,print_names = true,
                      dense=false,verbose = false,
                      save_everystep=false,appxsol=[test_sol,test_sol],
                      maxiters=Int(1e5),numruns=100)
plot(wp,title="Stiff 2: Hires",legend=:outertopleft,
     color=permutedims([repeat([:LightGreen],3)...,repeat([:DarkGreen],3)...,
     :Red,:Red,repeat([:Orange],2)...,repeat([:Yellow],3)...,
     repeat([:Blue],1)...,:Purple]),size = (800,350),
     xticks = 10.0 .^ (-12:1:5),
     yticks = 10.0 .^ (-6:0.5:5),
     bottom_margin=5Plots.mm)
```

```
Julia: Rosenbrock23
Julia: Rodas4
Julia: radau
Julia: Rosenbrock23 Static
Julia: Rodas4 Static
Julia: radau Static
Hairer: rodas
Hairer: radau
MATLAB: ode23s
MATLAB: ode15s
SciPy: LSODA
SciPy: BDF
SciPy: odeint
deSolve: lsoda
Sundials: CVODE
```


![](figures/wrapper_packages_5_1.png)
