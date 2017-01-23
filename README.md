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
  - [Van der Pol Work-Precision Diagrams](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/VanDerPol.ipynb)
  - [ROBER Work-Precision Diagrams](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/ROBER.ipynb)
- Parallelism
  - [Multithreaded Runge-Kutta Methods](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/Parallelism/Multithreaded%20Runge-Kutta%20Methods.ipynb)
  
The following tests were developed for the paper *Adaptive Methods for Stochastic Differential Equations via Natural Embeddings and Rejection Sampling with Memory*. These notebooks track their latest developments.

- SDE Adaptivity
  - [Oval2 Long Run](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/Oval2LongRun.ipynb)
  - [Oval2 Timings](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/Oval2Timings.ipynb)
  - [qmax Determination Tests](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/qmaxDetermination.ipynb)
  - [Adaptive Efficiency Tests](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/AdaptiveSDE/AdaptiveEfficiencyTests.ipynb)
