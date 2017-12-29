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
- Constant Delay DDEs
  - [Mackey and Glass Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/DDE/Mackey%20and%20Glass.ipynb)
  - [Stiff Quorum Sensing Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/DDE/Quorum%20Sensing.ipynb)
  - [Wheldon, Kirk, and Finlay Work-Precision Diagrams](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/DDE/Wheldon%2C%20Kirk%2C%20and%20Finlay.ipynb)

The following tests were developed for the paper *Adaptive Methods for Stochastic Differential Equations via Natural Embeddings and Rejection Sampling with Memory*. These notebooks track their latest developments.

- SDE Adaptivity
  - [Oval2 Long Run](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/Oval2LongRun.ipynb)
  - [Oval2 Long Timings](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/Oval2LongTimes.ipynb))
  - [Oval2 Timings](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/Oval2Timings.ipynb)
  - [qmax Determination Tests](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/qmaxDetermination.ipynb)
  - [Adaptive Efficiency Tests](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/AdaptiveEfficiencyTests.ipynb)
