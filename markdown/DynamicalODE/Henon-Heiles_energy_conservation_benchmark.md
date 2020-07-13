---
author: "Sebastian Micluța-Câmpeanu, Chris Rackauckas"
title: "Hénon–Heiles Energy Conservation"
---


In this notebook we will study the energy conservation properties of several high-order methods
for the Hénon–Heiles system.
We will se how the energy error behaves at very thight tolerances and how different techniques,
such as using symplectic solvers or manifold projections, benchmark against each other.
The Hamiltonian for this system is given by:
$$
\mathcal{H}=\frac{1}{2}(p_1^2 + p_2^2) + \frac{1}{2}\left(q_1^2 + q_2^2 + 2q_1^2 q_2 - \frac{2}{3}q_2^3\right)
$$

We will also compare the in place apporach with the out of place approach by using `Array`s
(for the in place version) and `StaticArrays` (for out of place versions).
In order to separate these two, we will define the relevant functions and initial conditions
in the `:inplace` and `OutofPlace` modules. In this way the rest of the code will work for both.

````julia
using OrdinaryDiffEq, Plots, DiffEqCallbacks
using TaylorIntegration, LinearAlgebra

T(p) = 1//2 * norm(p)^2
V(q) = 1//2 * (q[1]^2 + q[2]^2 + 2q[1]^2 * q[2]- 2//3 * q[2]^3)
H(p,q, params) = T(p) + V(q)

function iip_q̇(dq,p,q,params,t)
    dq[1] = p[1]
    dq[2] = p[2]
end

function iip_ṗ(dp,p,q,params,t)
    dp[1] = -q[1] * (1 + 2q[2])
    dp[2] = -q[2] - (q[1]^2 - q[2]^2)
end

iip_q0 = [0.1, 0.]
iip_p0 = [0., 0.5]


function oop_q̇(p, q, params, t)
    p
end

function oop_ṗ(p, q, params, t)
    dp1 = -q[1] * (1 + 2q[2])
    dp2 = -q[2] - (q[1]^2 - q[2]^2)
    @SVector [dp1, dp2]
end
````


````
Error: LoadError: UndefVarError: @SVector not defined
in expression starting at none:5
````



````julia

oop_q0 = @SVector [0.1, 0.]
````


````
Error: LoadError: UndefVarError: @SVector not defined
in expression starting at none:2
````



````julia
oop_p0 = @SVector [0., 0.5]
````


````
Error: LoadError: UndefVarError: @SVector not defined
in expression starting at none:1
````



````julia


function g(resid, u)
    resid[1] = H([u[1],u[2]], [u[3],u[4]], nothing) - E
    resid[2:4] .= 0
end

const cb = ManifoldProjection(g, nlopts=Dict(:ftol=>1e-13))

const E = H(iip_p0, iip_q0, nothing)
````


````
0.13
````





For the comparison we will use the following function

````julia
energy_err(sol) = map(i->H([sol[1,i], sol[2,i]], [sol[3,i], sol[4,i]], nothing)-E, 1:length(sol.u))
abs_energy_err(sol) = [abs.(H([sol[1,j], sol[2,j]], [sol[3,j], sol[4,j]], nothing) - E) for j=1:length(sol.u)]

function compare(mode=:inplace, all=true, plt=nothing; tmax=1e2)
    if mode == :inplace
        prob = DynamicalODEProblem(iip_ṗ, iip_q̇, iip_p0, iip_q0, (0., tmax))
    else
        prob = DynamicalODEProblem(oop_ṗ, oop_q̇, oop_p0, oop_q0, (0., tmax))
    end

    GC.gc()
    (mode == :inplace && all) && @time sol1 = solve(prob, Vern9(), callback=cb, abstol=1e-14, reltol=1e-14)
    GC.gc()
    @time sol2 = solve(prob, KahanLi8(), dt=1e-2, maxiters=1e10)
    GC.gc()
    @time sol3 = solve(prob, SofSpa10(), dt=1e-2, maxiters=1e8)
    GC.gc()
    @time sol4 = solve(prob, Vern9(), abstol=1e-14, reltol=1e-14)
    GC.gc()
    @time sol5 = solve(prob, DPRKN12(), abstol=1e-14, reltol=1e-14)
    GC.gc()
    (mode == :inplace && all) && @time sol6 = solve(prob, TaylorMethod(50), abstol=1e-20)

    (mode == :inplace && all) && println("Vern9 + ManifoldProjection max energy error:\t"*
        "$(maximum(abs_energy_err(sol1)))\tin\t$(length(sol1.u))\tsteps.")
    println("KahanLi8 max energy error:\t\t\t$(maximum(abs_energy_err(sol2)))\tin\t$(length(sol2.u))\tsteps.")
    println("SofSpa10 max energy error:\t\t\t$(maximum(abs_energy_err(sol3)))\tin\t$(length(sol3.u))\tsteps.")
    println("Vern9 max energy error:\t\t\t\t$(maximum(abs_energy_err(sol4)))\tin\t$(length(sol4.u))\tsteps.")
    println("DPRKN12 max energy error:\t\t\t$(maximum(abs_energy_err(sol5)))\tin\t$(length(sol5.u))\tsteps.")
    (mode == :inplace && all) && println("TaylorMethod max energy error:\t\t\t$(maximum(abs_energy_err(sol6)))"*
        "\tin\t$(length(sol6.u))\tsteps.")

    if plt === nothing
        plt = plot(xlabel="t", ylabel="Energy error")
    end
    (mode == :inplace && all) && plot!(sol1.t, energy_err(sol1), label="Vern9 + ManifoldProjection")
    plot!(sol2.t, energy_err(sol2), label="KahanLi8", ls=mode==:inplace ? :solid : :dash)
    plot!(sol3.t, energy_err(sol3), label="SofSpa10", ls=mode==:inplace ? :solid : :dash)
    plot!(sol4.t, energy_err(sol4), label="Vern9", ls=mode==:inplace ? :solid : :dash)
    plot!(sol5.t, energy_err(sol5), label="DPRKN12", ls=mode==:inplace ? :solid : :dash)
    (mode == :inplace && all) && plot!(sol6.t, energy_err(sol6), label="TaylorMethod")

    return plt
end
````


````
compare (generic function with 4 methods)
````





The `mode` argument choses between the in place approach
and the out of place one. The `all` parameter is used to compare only the integrators that support both
the in place and the out of place versions (we reffer here only to the 6 high order methods chosen bellow).
The `plt` argument can be used to overlay the results over a previous plot and the `tmax` keyword determines
the simulation time.

Note:
1. The `Vern9` method is used with `ODEProblem` because of performance issues with `ArrayPartition` indexing which manifest for `DynamicalODEProblem`.
2. The `NLsolve` call used by `ManifoldProjection` was modified to use `ftol=1e-13` in order to obtain a very low energy error.

Here are the results of the comparisons between the in place methods:

````julia
compare(tmax=1e2)
````


````
28.477321 seconds (42.07 M allocations: 1.597 GiB, 1.31% gc time)
  2.993858 seconds (9.36 M allocations: 451.505 MiB, 1.53% gc time)
  4.171649 seconds (11.93 M allocations: 551.610 MiB, 1.79% gc time)
  0.196334 seconds (1.07 M allocations: 43.356 MiB)
  2.035021 seconds (5.77 M allocations: 315.643 MiB, 1.19% gc time)
Error: type DynamicalODEFunction has no field f
````



````julia
compare(tmax=1e3)
````


````
0.141850 seconds (6.11 M allocations: 193.736 MiB, 25.09% gc time)
  0.063475 seconds (1.50 M allocations: 87.744 MiB)
  0.098909 seconds (1.50 M allocations: 87.744 MiB)
  0.064437 seconds (4.98 M allocations: 133.042 MiB)
  0.003351 seconds (48.86 k allocations: 2.924 MiB)
Error: type DynamicalODEFunction has no field f
````



````julia
compare(tmax=1e4)
````


````
1.898340 seconds (61.25 M allocations: 1.899 GiB, 35.43% gc time)
  1.336070 seconds (15.00 M allocations: 831.609 MiB, 40.89% gc time)
  1.736581 seconds (15.00 M allocations: 831.609 MiB, 39.65% gc time)
  6.773495 seconds (49.77 M allocations: 1.297 GiB, 89.40% gc time)
  0.035806 seconds (486.99 k allocations: 28.546 MiB)
Error: type DynamicalODEFunction has no field f
````



````julia
compare(tmax=5e4)
````


````
12.047769 seconds (305.15 M allocations: 9.435 GiB, 49.59% gc time)
 20.483846 seconds (75.00 M allocations: 4.061 GiB, 80.40% gc time)
 32.531808 seconds (75.00 M allocations: 4.061 GiB, 83.93% gc time)
132.424261 seconds (248.84 M allocations: 6.482 GiB, 97.20% gc time)
  1.258019 seconds (2.43 M allocations: 140.316 MiB, 86.65% gc time)
Error: type DynamicalODEFunction has no field f
````





We can see that as the simulation time increases, the energy error increases. For this particular example
the energy error for all the methods is comparable. For relatively short simulation times,
if a highly accurate solution is required, the symplectic method is not recommended as
its energy error fluctuations are larger than for other methods.
An other thing to notice is the fact that the two versions of `Vern9` behave identically, as expected,
untill the energy error set by `ftol` is reached.

We will now compare the in place with the out of place versions. In the plots bellow we will use
a dashed line for the out of place versions.

````julia
function in_vs_out(;all=false, tmax=1e2)
    println("In place versions:")
    plt = compare(:inplace, all, tmax=tmax)
    println("\nOut of place versions:")
    plt = compare(:oop, false, plt; tmax=tmax)
end
````


````
in_vs_out (generic function with 1 method)
````





First, here is a summary of all the available methods for `tmax = 1e3`:

````julia
in_vs_out(all=true, tmax=1e3)
````


````
In place versions:
  0.121460 seconds (6.11 M allocations: 193.736 MiB, 10.67% gc time)
  0.063232 seconds (1.50 M allocations: 87.744 MiB)
  0.097984 seconds (1.50 M allocations: 87.744 MiB)
  0.065835 seconds (4.98 M allocations: 133.042 MiB)
  0.003351 seconds (48.86 k allocations: 2.924 MiB)
Error: type DynamicalODEFunction has no field f
````





Now we will compare the in place and the out of place versions, but only for the integrators
that are compatible with `StaticArrays`

````julia
in_vs_out(tmax=1e2)
````


````
In place versions:
  0.006880 seconds (150.07 k allocations: 8.321 MiB)
  0.008175 seconds (150.07 k allocations: 8.321 MiB)
  0.004847 seconds (498.57 k allocations: 13.350 MiB)
  0.000302 seconds (5.05 k allocations: 299.688 KiB)
KahanLi8 max energy error:			4.9404924595819466e-15	in	10001	steps.
SofSpa10 max energy error:			5.440092820663267e-15	in	10001	steps.
Vern9 max energy error:				1.4432899320127035e-15	in	853	steps.
DPRKN12 max energy error:			2.220446049250313e-16	in	329	steps.

Out of place versions:
Error: UndefVarError: oop_ṗ not defined
````



````julia
in_vs_out(tmax=1e3)
````


````
In place versions:
  0.098971 seconds (1.50 M allocations: 87.744 MiB, 38.78% gc time)
  0.088299 seconds (1.50 M allocations: 87.744 MiB)
  0.067235 seconds (4.98 M allocations: 133.042 MiB)
  0.003286 seconds (48.86 k allocations: 2.924 MiB)
KahanLi8 max energy error:			1.815214645262131e-14	in	100002	steps.
SofSpa10 max energy error:			2.8033131371785203e-14	in	100002	steps.
Vern9 max energy error:				1.1796119636642288e-14	in	8472	steps.
DPRKN12 max energy error:			3.0531133177191805e-16	in	3240	steps.

Out of place versions:
Error: UndefVarError: oop_ṗ not defined
````



````julia
in_vs_out(tmax=1e4)
````


````
In place versions:
  1.456041 seconds (15.00 M allocations: 831.609 MiB, 46.57% gc time)
  1.431695 seconds (15.00 M allocations: 831.609 MiB, 27.18% gc time)
  1.060642 seconds (49.77 M allocations: 1.297 GiB, 28.32% gc time)
  0.034983 seconds (486.99 k allocations: 28.546 MiB)
KahanLi8 max energy error:			3.161360062620133e-14	in	1000001	steps.
SofSpa10 max energy error:			1.136590821460004e-13	in	1000001	steps.
Vern9 max energy error:				1.2015388684005757e-13	in	84654	steps.
DPRKN12 max energy error:			1.1879386363489175e-14	in	32360	steps.

Out of place versions:
Error: UndefVarError: oop_ṗ not defined
````



````julia
in_vs_out(tmax=5e4)
````


````
In place versions:
  8.836682 seconds (75.00 M allocations: 4.061 GiB, 52.88% gc time)
 18.178625 seconds (75.00 M allocations: 4.061 GiB, 71.20% gc time)
101.562255 seconds (248.84 M allocations: 6.482 GiB, 96.14% gc time)
  1.106361 seconds (2.43 M allocations: 140.316 MiB, 84.79% gc time)
KahanLi8 max energy error:			1.2331802246023926e-13	in	5000001	steps.
SofSpa10 max energy error:			1.5035195310986182e-13	in	5000001	steps.
Vern9 max energy error:				6.044331701815508e-13	in	423231	steps.
DPRKN12 max energy error:			5.789813073420191e-14	in	161763	steps.

Out of place versions:
Error: UndefVarError: oop_ṗ not defined
````





As we see from the above comparisons, the `StaticArray` versions are significantly faster and use less memory.
The speedup provided for the out of place version is more proeminent at larger values for `tmax`.
We can see again that if the simulation time is increased, the energy error of the symplectic methods
is less noticeable compared to the rest of the methods.

The benchmarks were performed on a machine with

````julia
using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])
````



## Appendix

These benchmarks are a part of the DiffEqBenchmarks.jl repository, found at: [https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl)

To locally run this tutorial, do the following commands:

```
using DiffEqBenchmarks
DiffEqBenchmarks.weave_file("DynamicalODE","Henon-Heiles_energy_conservation_benchmark.jmd")
```

Computer Information:

```
Julia Version 1.4.2
Commit 44fa15b150* (2020-05-23 18:35 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
Environment:
  JULIA_DEPOT_PATH = /builds/JuliaGPU/DiffEqBenchmarks.jl/.julia
  JULIA_CUDA_MEMORY_LIMIT = 2147483648
  JULIA_PROJECT = @.
  JULIA_NUM_THREADS = 8

```

Package Information:

```
Status: `/builds/JuliaGPU/DiffEqBenchmarks.jl/benchmarks/DynamicalODE/Project.toml`
[459566f4-90b8-5000-8ac3-15dfb0a30def] DiffEqCallbacks 2.13.3
[055956cb-9e8b-5191-98cc-73ae4a59e68a] DiffEqPhysics 3.2.0
[b305315f-e792-5b7a-8f41-49f472929428] Elliptic 1.0.1
[1dea7af3-3e70-54e6-95c3-0bf5283fa5ed] OrdinaryDiffEq 5.41.0
[65888b18-ceab-5e60-b2b9-181511a3b968] ParameterizedFunctions 5.4.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.5.4
[d330b81b-6aea-500a-939a-2ce795aea3ee] PyPlot 2.9.0
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.4
[92b13dbe-c966-51a2-8445-caca9f8a7d42] TaylorIntegration 0.8.3
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra 
[de0858da-6303-5e67-8744-51eddeeeb8d7] Printf 
[10745b16-79ce-11e8-11f9-7d13ad32a3b2] Statistics 
```

