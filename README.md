# DiffEqBenchmarks.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

DiffEqBenchmarks.jl holds Jupyter notebooks showing the benchmarks for the
JuliaDiffEq ecosystem.

## Viewing the Notebooks Locally

To view the notebooks locally and interact with the contents, use the following
commands (requires [IJulia](https://github.com/JuliaLang/IJulia.jl)):

```julia
Pkg.clone("https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl")
using IJulia
notebook(dir=Pkg.dir("DiffEqBenchmarks"))
```

## Table of Contents

The notebooks can be viewed remotely on Github or via [nbviewer](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/tree/master/)

- Non-stiff ODEs
  - [Linear Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffODE/Linear%20Work-Precision%20Diagrams.ipynb)
  - [Runge-Kutta Benchmarks on Linear ODEs](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffODE/Runge-Kutta%20Benchmarks%20on%20Linear%20ODEs.ipynb)
  - [Three-Body Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffODE/ThreeBody%20Work-Precision%20Diagrams.ipynb)
  - [Pleides Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffODE/Pleiades%20Work-Precision%20Diagrams.ipynb)
  - [Rigid Body Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffODE/RigidBody%20Work-Precision%20Diagrams.ipynb)
  - [Fizhugh-Nagumo Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffODE/FitzhughNagumo%20Work-Precision%20Diagrams.ipynb)
  - [Lotka-Volterra Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffODE/LotkaVolterra%20Work-Precision%20Diagrams.ipynb)
  - [Direct vs MATLAB Benchmark](https://github.com/JuliaDiffEq/MATLABDiffEq.jl#benchmark)
- Stiff ODEs
  - [Van der Pol Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/VanDerPol.ipynb)
  - [ROBER Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/ROBER.ipynb)
  - [Orego Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/Orego.ipynb)
  - [Hires Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/Hires.ipynb)
  - [Pollution Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/Pollution.ipynb)
  - [Filament PDE Discretization Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/Filament.ipynb)
- Dynamical ODEs
  - [Single Pendulum Comparison Benchmark](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/DynamicalODE/single_pendulums.ipynb)
  - [Henon-Heiles Energy Conservation Benchmark](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/DynamicalODE/Henon-Heiles_energy_conservation_benchmark.ipynb)
  - [Quadrupole Boson Hamiltonian Energy Conservation Benchmark](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/DynamicalODE/Quadrupole_boson_Hamiltonian_energy_conservation_benchmark.ipynb)
- Parallelism
  - [Multithreaded Runge-Kutta Methods](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/Parallelism/Multithreaded%20Runge-Kutta%20Methods.ipynb)
- Nonstiff SDEs
  - [Simple Nonstiff SDE Strong Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffSDE/BasicSDEWorkPrecision.ipynb)
  - [Simple Nonstiff SDE Weak Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffSDE/BasicSDEWeakWorkPrecision.ipynb)
  - [Lotka-Volterra SDE Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonStiffSDE/LotkaVolterraSDE.ipynb)
- Stiff SDEs
  - [Stochastic Heat Equation Investigation](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffSDE/StochasticHeat.ipynb)
  - [Quadratic Diffusion Noise Investigation](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffSDE/QuadraticStiffness.ipynb)
  - [Oval2 Long Run](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/Oval2LongRun.ipynb)
  - [Oval2 Long Timings](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/Oval2LongTimes.ipynb)
  - [Oval2 Timings](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/Oval2Timings.ipynb)
- Nonstiff DDEs
  - Constant Delay DDEs
    - [Mackey and Glass Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonstiffDDE/Mackey%20and%20Glass.ipynb)
    - [Wheldon, Kirk, and Finlay Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/NonstiffDDE/Wheldon%2C%20Kirk%2C%20and%20Finlay.ipynb)
- Stiff DDEs
  - [Quorum Sensing Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffDDE/Quorum%20Sensing.ipynb)
- Parameter Estimation
  - [Chaotic Lorenz Equation Parameter Estimation](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/ParameterEstimation/LorenzParameterEstimation.ipynb)
  - [Bayesian Lotka-Volterra Parameter Estimation](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/ParameterEstimation/DiffEqBayesLotkaVolterra.ipynb)
  - [Bayesian Lorenz Equation Estimation](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/ParameterEstimation/DiffEqBayesLorenz.ipynb)

The following tests were developed for the paper *Adaptive Methods for Stochastic Differential Equations via Natural Embeddings and Rejection Sampling with Memory*. These notebooks track their latest developments.

- SDE Adaptivity

  - [qmax Determination Tests](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/qmaxDetermination.ipynb)
  - [Adaptive Efficiency Tests](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/AdaptiveEfficiencyTests.ipynb)

## Current Summary

The following is a quick summary of the benchmarks. These paint broad strokes
over the set of tested equations and some specific examples may differ.

### Non-Stiff ODEs

- OrdinaryDiffEq.jl's methods are the most efficient by a good amount
- The `Vern` methods tend to do the best in every benchmark of this category
- At lower tolerances, `Tsit5` does well consistently.
- ARKODE and Hairer's `dopri5`/`dop853` perform very similarly, but are both
  far less efficient than the `Vern` methods.
- The multistep methods, `CVODE_Adams` and `lsoda`, tend to not do very well.
- The ODEInterface multistep method `ddeabm` does not do as well as the other
  multistep methods.
- ODE.jl's methods are not able to consistently solve the problems.
- Fixed time step methods are less efficient than the adaptive methods.

### Stiff ODEs

- In this category, the best methods are much more problem dependent.
- For smaller problems:
  - `Rosenbrock23`, `lsoda`, and `TRBDF2` tend to be the most efficient at high  
    tolerances.
  - `Rodas4` and `Rodas5` tend to be the most efficient at low tolerances.
- For larger problems (Filament PDE):
  - `CVODE_BDF` does the best at all tolerances.
  - The ESDIRK methods like `TRBDF2` and `KenCarp4` can come close.
- `radau` is always the most efficient when tolerances go to the low extreme
  (`1e-13`)
- Fixed time step methods tend to diverge on every tested problem because the
  high stiffness results in divergence of the Newton solvers.
- ARKODE is very inconsistent and requires a lot of tweaking in order to not
  diverge on many of the tested problems. When it doesn't diverge, the similar
  algorithms in OrdinaryDiffEq.jl (`KenCarp4`) are much more efficient in most
  cases.
- ODE.jl and GeometricIntegrators.jl fail to converge on any of the tested
  problems.

### Dynamical ODEs

- Higher order (generally order >=6) symplectic integrators are much more
  efficient than the lower order counterparts.
- For high accuracy, using a symplectic integrator is not preferred. Their extra
  cost is not necessary since the other integrators are able to not drift simply
  due to having low enough error.
- In this class, the `DPRKN` methods are by far the most efficient. The `Vern`
  methods do well for not being specific to the domain.

### Non-Stiff SDEs

- For simple 1-dimensional SDEs at low accuracy, the `EM` and `RKMil` methods
  can do well. Beyond that, they are simply outclassed.
- The `SRA` and `SRI` methods both are very similar within-class on the simple
  SDEs.
- `SRA3` is the most efficient when applicable and the tolerances are low.
- Generally, only low accuracy is necessary to get to sampling error of the mean.
- The adaptive method is very conservative with error estimates.

### Stiff SDEs

- The high order adaptive methods (`SRIW1`) generally do well on stiff problems.
- The "standard" low-order implicit methods, `ImplicitEM` and `ImplicitRK`, do
  not do well on all stiff problems. Some exceptions apply to well-behaved
  problems like the Stochastic Heat Equation.

### Non-Stiff DDEs

- The efficiency ranking tends to match the ODE Tests, but the cutoff from
  low to high tolerance is lower.
- `Tsit5` does well in a large class of problems here.
- The `Vern` methods do well in low tolerance cases.

### Stiff DDEs

- The Rosenbrock methods, specifically `Rodas5`, perform well.
