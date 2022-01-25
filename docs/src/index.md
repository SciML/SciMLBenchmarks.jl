# LinearSolve.jl: High-Performance Unified Linear Solvers

LinearSolve.jl is a unified interface for the linear solving packages of
Julia. It interfaces with other packages of the Julia ecosystem
to make it easy to test alternative solver packages and pass small types to
control algorithm swapping. It also interfaces with the
[ModelingToolkit.jl](https://mtk.sciml.ai/dev/) world of symbolic modeling to
allow for automatically generating high-performance code.

Performance is key: the current methods are made to be highly performant on
scalar and statically sized small problems, with options for large-scale systems.
If you run into any performance issues, please file an issue.

## Installation

To install LinearSolve.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("LinearSolve")
```

## Contributing

- Please refer to the
  [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
  for guidance on PRs, issues, and other matters relating to contributing to ModelingToolkit.
- There are a few community forums:
    - the #diffeq-bridged channel in the [Julia Slack](https://julialang.org/slack/)
    - [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
    - on the [Julia Discourse forums](https://discourse.julialang.org)
    - see also [SciML Community page](https://sciml.ai/community/)

## Roadmap

Wrappers for every linear solver in the Julia language is on the roadmap. If
there are any important ones that are missing that you would like to see added,
please open an issue. The current algorithms should support automatic differentiation.
Pre-defined preconditioners would be a welcome addition.
