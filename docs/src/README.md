# SciMLBenchmarks.jl: Benchmarks for Scientific Machine Learning (SciML) and Differential Equation Solver Software

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build status](https://badge.buildkite.com/2f4b5708bf098c75ce193f04b3f3c4047f993f0e363e314c61.svg)](https://buildkite.com/julialang/scimlbenchmarks-dot-jl)

SciMLBenchmarks.jl holds webpages, pdfs, and notebooks showing the benchmarks
for the SciML Scientific Machine Learning Software ecosystem, including cross-language
benchmarks of differential equation solvers and methods for parameter estimation,
training universal differential equations (and subsets like neural ODEs), tests
of physics-informed neural networks (PINNs), and more.

## Results
Static outputs in pdf, markdown, and html reside in [SciMLBenchmarksOutput](https://github.com/SciML/SciMLBenchmarksOutput).

## Table of Contents

- Multi-Language Wrapper Benchmarks
  - [ODE Solver Multi-Language Wrapper Package Work-Precision Benchmarks (MATLAB, SciPy, Julia, deSolve (R))](https://benchmarks.sciml.ai/html/MultiLanguage/wrapper_packages.html)
  - [Torchdiffeq vs DifferentialEquations.jl (/ DiffEqFlux.jl) Benchmarks](https://gist.github.com/ChrisRackauckas/cc6ac746e2dfd285c28e0584a2bfd320)
  - [torchdiffeq vs Julia DiffEqFlux Neural ODE Training Benchmark](https://gist.github.com/ChrisRackauckas/4a4d526c15cc4170ce37da837bfc32c4)
  - [torchsde vs DifferentialEquations.jl / DiffEqFlux.jl](https://gist.github.com/ChrisRackauckas/6a03e7b151c86b32d74b41af54d495c6)
  - [JITCODE vs SciPy vs DifferentialEquations.jl on large network dynamics](https://github.com/PIK-ICoN/NetworkDynamicsBenchmarks)
  - [DifferentialEquations.jl vs Mujuco and DiffTaichi](https://arxiv.org/abs/2012.06684)
  - [DiffEqFlux.jl / DifferentialEquations.jl vs Jax on an epidemic model](https://gist.github.com/ChrisRackauckas/62a063f23cccf3a55a4ac9f6e497739a)
  - [DifferentialEquations.jl vs SciPy vs NumbaLSODA on a stiff ODE](https://gist.github.com/ChrisRackauckas/fd62e005c4c86520306338b6bdae6b79)
  - [DifferentialEquations.jl vs SciPy vs NumbaLSODA](https://github.com/Nicholaswogan/NumbaLSODA/tree/main/benchmark)
- Non-stiff Ordinary Differential Equations (ODEs)
  - [Linear Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffODE/linear_wpd.html)
  - [Three-Body Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffODE/ThreeBody_wpd.html)
  - [Pleides Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffODE/Pleiades_wpd.html)
  - [Rigid Body Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffODE/RigidBody_wpd.html)
  - [Fizhugh-Nagumo Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffODE/FitzhughNagumo_wpd.html)
  - [Lotka-Volterra Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffODE/LotkaVolterra_wpd.html)
  - [Direct vs MATLAB Benchmark](https://github.com/JuliaDiffEq/MATLABDiffEq.jl#benchmark)
  - [Runge-Kutta vs Taylor Integration on Pleides](https://gist.github.com/ChrisRackauckas/1301b23aa12ad83de7138d8e41d64dd6)
- Stiff Ordinary Differential Equations (ODEs)
  - [Van der Pol Work-Precision Diagrams](https://benchmarks.sciml.ai/html/StiffODE/VanDerPol.html)
  - [ROBER Work-Precision Diagrams](https://benchmarks.sciml.ai/html/StiffODE/ROBER.html)
  - [Orego Work-Precision Diagrams](https://benchmarks.sciml.ai/html/StiffODE/Orego.html)
  - [Hires Work-Precision Diagrams](https://benchmarks.sciml.ai/html/StiffODE/Hires.html)
  - [Pollution Work-Precision Diagrams](https://benchmarks.sciml.ai/html/StiffODE/Pollution.html)
  - [BCR (1122 ODE) Work-Precision Diagrams](https://benchmarks.sciml.ai/html/Bio/BCR.html)
- Differential-Algebraic Equations (DAEs)
  - [ROBER DAE Work-Precision Diagrams](https://benchmarks.sciml.ai/html/DAE/ROBERDAE.html)
  - [OREGO DAE Work-Precision Diagrams](https://benchmarks.sciml.ai/html/DAE/OregoDAE.html)
  - [Chemical Akzo Nobel Differential-Algebraic Equation (DAE) Work-Precision Diagrams](https://benchmarks.sciml.ai/html/DAE/ChemicalAkzoNobel.html)
- Method of Lines PDEs
  - [Filament PDE Discretization Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/Filament.html)
  - [Allen-Cahn Finite Difference Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/allen_cahn_fdm_wpd.html)
  - [Allen-Cahn Pseudospectral Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/allen_cahn_spectral_wpd.html)
  - [Burger's Finite Difference Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/burgers_fdm_wpd.html)
  - [Burger's Pseudospectral Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/burgers_spectral_wpd.html)
  - [KdV Finite Difference Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/kdv_fdm_wpd.html)
  - [KdV Pseudospectral Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/kdv_spectral_wpd.html)
  - [Kuramoto–Sivashinsky Finite Difference Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/ks_fdm_wpd.html)
  - [Kuramoto–Sivashinsky Pseudospectral Work-Precision Diagrams](https://benchmarks.sciml.ai/html/MOLPDE/ks_spectral_wpd.html)
- Dynamical ODEs
  - [Single Pendulum Comparison Benchmark](https://benchmarks.sciml.ai/html/DynamicalODE/single_pendulums.html)
  - [Henon-Heiles Energy Conservation Benchmark](https://benchmarks.sciml.ai/html/DynamicalODE/Henon-Heiles_energy_conservation_benchmark.html)
  - [Quadrupole Boson Hamiltonian Energy Conservation Benchmark](https://benchmarks.sciml.ai/html/DynamicalODE/Quadrupole_boson_Hamiltonian_energy_conservation_benchmark.html)
- N-body problems
  - [Acceleration function benchmarks](https://benchmarks.sciml.ai/html/NBodySimulator/acceleration_functions.html)
  - [Liquid argon benchmarks](https://benchmarks.sciml.ai/html/NBodySimulator/liquid_argon.html)
- Nonstiff SDEs
  - [Simple Nonstiff SDE Strong Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffSDE/BasicSDEWorkPrecision.html)
  - [Simple Nonstiff SDE Weak Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffSDE/BasicSDEWeakWorkPrecision.html)
  - [Lotka-Volterra SDE Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffSDE/LotkaVolterraSDE.html	)
- Stiff SDEs
  - [Stochastic Heat Equation Investigation](https://benchmarks.sciml.ai/html/StiffSDE/StochasticHeat.html)
  - [Quadratic Diffusion Noise Investigation](https://benchmarks.sciml.ai/html/StiffSDE/QuadraticStiffness.html)
  - [Oval2 Long Run](https://benchmarks.sciml.ai/html/StiffSDE/Oval2LongRun.html)
  - [Oval2 Long Timings](https://benchmarks.sciml.ai/html/StiffSDE/Oval2LongTimes.html)
  - [Oval2 Timings](https://benchmarks.sciml.ai/html/StiffSDE/Oval2Timings.html)
- Nonstiff DDEs
  - Constant Delay DDEs
    - [Mackey and Glass Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffDDE/Mackey_Glass_wpd.html)
    - [Wheldon, Kirk, and Finlay Work-Precision Diagrams](https://benchmarks.sciml.ai/html/NonStiffDDE/Wheldon_Kirk_Finlay_wpd.html)
- Stiff DDEs
  - [Quorum Sensing Work-Precision Diagrams](https://benchmarks.sciml.ai/html/StiffDDE/QuorumSensing.html)
- Jump Equations
  - [Diffusion Model](https://benchmarks.sciml.ai/html/Jumps/Diffusion_CTRW.html)
  - [Mendes Multistate Model](https://benchmarks.sciml.ai/html/Jumps/Mendes_multistate_example.html)
  - [Negative Feedback Gene Expression Model](https://benchmarks.sciml.ai/html/Jumps/NegFeedback_GeneExpr.html)
  - [Negative Feedback Marchetti Model](https://benchmarks.sciml.ai/html/Jumps/NegFeedback_GeneExpr_Marchetti.html)
- Parameter Estimation
  - [Lorenz Equation Parameter Estimation by Optimization Methods](https://benchmarks.sciml.ai/html/ParameterEstimation/LorenzParameterEstimation.html)
  - [Bayesian Lotka-Volterra Parameter Estimation](https://benchmarks.sciml.ai/html/ParameterEstimation/DiffEqBayesLotkaVolterra.html)
  - [Bayesian Lorenz Equation Estimation](https://benchmarks.sciml.ai/html/ParameterEstimation/DiffEqBayesLorenz.html)
  - [Bayesian FitzHugh-Nagumo Equation Estimation](https://benchmarks.sciml.ai/html/ParameterEstimation/DiffEqBayesFitzHughNagumo.html)
  - [Lotka Volterra Equation Parameter Estimation by Optimization Methods](https://benchmarks.sciml.ai/html/ParameterEstimation/LotkaVolterraParameterEstimation.html)
  - [FitzHugh-Nagumo Equation Parameter Estimation by Optimization Methods](https://benchmarks.sciml.ai/html/ParameterEstimation/FitzHughNagumoParameterEstimation.html)
- Physics-Informed Neural Network (Neural Network PDE Solver) Cost Function Benchmarks
  - [Allen-Cahn PDE Physics-Informed Neural Network (PINN) Loss Function Error vs Time Benchmarks](https://benchmarks.sciml.ai/html/PINNErrorsVsTime/allen_cahn_et.html)
  - [Diffusion Equation Physics-Informed Neural Network (PINN) Loss Function Error vs Time Benchmarks](https://benchmarks.sciml.ai/html/PINNErrorsVsTime/diffusion_et.html)
  - [Hamilton-Jacobi PDE Physics-Informed Neural Network (PINN) Loss Function Error vs Time Benchmarks](https://benchmarks.sciml.ai/html/PINNErrorsVsTime/hamilton_jacobi_et.html)
  - [Level Set PDE Physics-Informed Neural Network (PINN) Loss Function Error vs Time Benchmarks](https://benchmarks.sciml.ai/html/PINNErrorsVsTime/level_set_et.html)
  - [Nernst-Planck PDE Physics-Informed Neural Network (PINN) Loss Function Error vs Time Benchmarks](https://benchmarks.sciml.ai/html/PINNErrorsVsTime/nernst_planck_et.html)
- Physics-Informed Neural Network (Neural Network PDE Solver) Optimizer Benchmarks
  - [Diffusion Equation Physics-Informed Neural Network (PINN) Optimizer Benchmarks](https://benchmarks.sciml.ai/html/PINNOptimizers/1d_diffusion.html)
  - [1D Nernst-Planck Equation Physics-Informed Neural Network (PINN) Optimizer Benchmarks](https://benchmarks.sciml.ai/html/PINNOptimizers/1d_poisson_nernst_planck.html)
  - [Allen-Cahn Equation Physics-Informed Neural Network (PINN) Optimizer Benchmarks](https://benchmarks.sciml.ai/html/PINNOptimizers/allen_cahn.html)
  - [Berger's Equation Physics-Informed Neural Network (PINN) Optimizer Benchmarks](https://benchmarks.sciml.ai/html/PINNOptimizers/burgers_equation.html)
  - [Hamilton-Jacobi Equation Physics-Informed Neural Network (PINN) Optimizer Benchmarks](https://benchmarks.sciml.ai/html/PINNOptimizers/hamilton_jacobi.html)
  - [Poisson Equation Physics-Informed Neural Network (PINN) Optimizer Benchmarks](https://benchmarks.sciml.ai/html/PINNOptimizers/poisson.html)

The following tests were developed for the paper *Adaptive Methods for Stochastic Differential Equations via Natural Embeddings and Rejection Sampling with Memory*. These notebooks track their latest developments.

- SDE Adaptivity

  - [qmax Determination Tests](https://benchmarks.sciml.ai/html/AdaptiveSDE/qmaxDetermination.html)
  - [Adaptive Efficiency Tests](https://benchmarks.sciml.ai/html/AdaptiveSDE/AdaptiveEfficiencyTests.html)

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
  - `QNDF` and `FBDF` does the best at all normal tolerances.
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

### Parameter Estimation

- Broadly two different approaches have been used, Bayesian Inference and Optimisation
  algorithms.
- In general it seems that the optimisation algorithms perform more accurately but that can be
  attributed to the larger number of data points being used in the optimisation cases, Bayesian
  approach tends to be slower of the two and hence lesser data points are used, accuracy can
  increase if proper data is used.
- Within the different available optimisation algorithms, BBO from the BlackBoxOptim package and GN_CRS2_LM
  for the global case while LD_SLSQP,LN_BOBYQA and LN_NELDERMEAD for the local case from the NLopt package
  perform the best.
- Another algorithm being used is the [QuadDIRECT](https://github.com/timholy/QuadDIRECT.jl) algorithm, it gives very good results in the shorter problem case
  but doesn't do very well in the case of the longer problems.
- The choice of global versus local optimization make a huge difference in the timings. BBO tends to find
  the correct solution for a global optimization setup. For local optimization, most methods in NLopt,
  like :LN_BOBYQA, solve the problem very fast but require a good initial condition.
- The different backends options available for Bayesian method offer some tradeoffs beteween
  time, accuracy and control. It is observed that sufficiently high accuracy can be observed with
  any of the backends with the fine tuning of stepsize, constraints on the parameters, tightness of the
  priors and number of iterations being passed.

## Interactive Notebooks

To run the tutorials interactively via Jupyter notebooks and benchmark on your
own machine
1. Run Weave for the file (or folder) you are interested in
2. Activate the appropriate environment
3. Open and run the notebook.

Note: Since notebooks default to looking for a Project.toml file at the same level or parent folder, you might need to move the notebook to the folder with the appropriate Project.toml.

### Example (starting from the project root folder)
```julia
]activate .
]instantiate
using SciMLBenchmarks
SciMLBenchmarks.weave_file("benchmarks/Jumps", "Diffusion_CTRW.jmd", [:notebook])
]activate benchmarks/Jumps
```

Then move `Diffusion_CTRW.ipynb` to "benchmarks/Jumps" and open the notebook.

## Contributing

All of the files are generated from the Weave.jl files in the `benchmarks` folder. To run the generation process, do for example:

```julia
]activate SciMLBenchmarks # Get all of the packages
using SciMLBenchmarks
SciMLBenchmarks.weave_file("NonStiffODE","linear_wpd.jmd")
```

To generate all of the files in a folder, for example, run:

```julia
SciMLBenchmarks.weave_folder("NonStiffODE")
```

To generate all of the notebooks, do:

```julia
SciMLBenchmarks.weave_all()
```

Each of the benchmarks displays the computer characteristics at the bottom of
the benchmark. Since performance-necessary computations are normally performed on
compute clusters, the official benchmarks use a workstation with an
Intel Xeon CPU E5-2680 v4 @ 2.40GHz to match the performance characteristics of
a standard node in a high performance computing (HPC) cluster or cloud computing
setup.

### Inspecting Benchmark Results

To see benchmark results before merging, click into the BuildKite, click onto
Artifacts, and then investigate the trained results.

![](https://user-images.githubusercontent.com/1814174/118359358-02ddc980-b551-11eb-8a9b-24de947cefee.PNG)
