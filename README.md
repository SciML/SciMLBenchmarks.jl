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

The notebooks can be viewed remotely on Github or via [nbviewer]()

- Non-stiff ODEs
  - Linear Work-Precision Diagrams
  - Runge-Kutta Benchmarks on Linear ODEs
  - Three-Body Work Benchmarks
  - Three-Body Work-Precision Diagrams
  - Pleides Work-Precision Diagrams
  - Rigid Body Work-Precision Diagrams
  - Fizhugh-Nagumo Work-Precision Diagrams
  - Lotka-Volterra Work-Precision Diagrams
- Parallelism
  - Multithreaded Runge-Kutta Methods
