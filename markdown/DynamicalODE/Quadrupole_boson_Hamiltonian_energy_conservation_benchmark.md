---
author: "Sebastian Micluța-Câmpeanu, Chris Rackauckas"
title: "Quadruple Boson Energy Conservation"
---


In this notebook we will study the energy conservation properties of several high-order methods for a system with the following Hamiltonian:
$$
\mathcal{H}\left(q_0,q_2,p_0,p_2\right) = \frac{A}{2} \left(p_0^2 + p_2^2 + q_0^2 + q_2^2\right) + \frac{B}{\sqrt{2}} q_0 \left(3q_2^2 - q_0^2\right) + \frac{D}{4} \left(q_0^2+q_2^2\right)^2
$$
This Hamiltonian resembles the Hénon-Heiles one, but it has an additional fourth order term.
The aim of this benchmark is to see what happens with the energy error when highly accurate solutions are needed and how the results compare with the Hénon-Heiles case.

````julia
using OrdinaryDiffEq, Plots, DiffEqCallbacks, LinearAlgebra
using TaylorIntegration

T(p) = A / 2 * norm(p)^2
V(q) = A / 2 * (q[1]^2 + q[2]^2) + B / √2 * q[1] * (3 * q[2]^2 - q[1]^2) + D / 4 * (q[1]^2 + q[2]^2)^2
H(p, q, params) = T(p) + V(q)

const A, B, D = 1., 0.55, 0.4

module InPlace
    using ParameterizedFunctions

    function q̇(dq, p, q, params, t)
        dq[1] = A * p[1]
        dq[2] = A * p[2]
    end

    function ṗ(dp, p, q, params, t)
        dp[1] = -A * q[1] - 3 * B / √2 * (q[2]^2 - q[1]^2) - D * q[1] * (q[1]^2 + q[2]^2)
        dp[2] = -q[2] * (A + 3 * √2 * B * q[1] + D * (q[1]^2 + q[2]^2))
    end

    const A, B, D = 1., 0.55, 0.4

    const q0 = [4.919080920016389, 2.836942666663649]
    const p0 = [0., 0.]
    const u0 = vcat(p0,q0)

    h_eqs = @ode_def HamiltEqs begin
      dp₀ = -A * q₀ - 3 * B / √2 * (q₂^2 - q₀^2) - D * q₀ * (q₀^2 + q₂^2)
      dp₂ = -q₂ * (A + 3 * √2 * B * q₀ + D * (q₀^2 + q₂^2))
      dq₀ = A * p₀
      dq₂ = A * p₂
    end A B D

    p = [1.0,0.55,0.4]

end
````


````
Error: LoadError: UndefVarError: ModelingToolkit not defined
in expression starting at <macrocall>:0
````



````julia

module OutOfPlace
    using StaticArrays

    function q̇(p, q, params, t)
        p
    end

    function ṗ(p, q, params, t)
        dp1 = -A * q[1] - 3 * B / √2 * (q[2]^2 - q[1]^2) - D * q[1] * (q[1]^2 + q[2]^2)
        dp2 = -q[2] * (A + 3 * √2 * B * q[1] + D * (q[1]^2 + q[2]^2))
        @SVector [dp1, dp2]
    end

    const A, B, D = 1., 0.55, 0.4

    const q0 = @SVector [4.919080920016389, 2.836942666663649]
    const p0 = @SVector [0., 0.]
    const u0 = vcat(p0,q0)

    h_eqs(z, params, t) = SVector(
        -A * z[3] - 3 * B / √2 * (z[4]^2 - z[3]^2) - D * z[3] * (z[3]^2 + z[4]^2),
        -z[4] * (A + 3 * √2 * B * z[3] + D * (z[3]^2 + z[4]^2)),
        z[1],
        z[2]
    )

    p = nothing

end

function g(resid,u)
    resid[1] = H([u[1],u[2]],[u[3],u[4]],nothing) - E
    resid[2:4] .= 0
end

const E = H(InPlace.p0, InPlace.q0, InPlace.p)
````


````
Error: UndefVarError: p not defined
````



````julia
const cb = ManifoldProjection(g, nlopts=Dict(:ftol=>1e-13))
````


````
DiffEqBase.DiscreteCallback{DiffEqCallbacks.var"#7#8",DiffEqCallbacks.Manif
oldProjection{true,DiffEqCallbacks.NonAutonomousFunction{typeof(Main.##Weav
eSandBox#311.g),true},DiffEqCallbacks.NLSOLVEJL_SETUP{0,true},Dict{Symbol,F
loat64}},typeof(DiffEqCallbacks.Manifold_initialize)}(DiffEqCallbacks.var"#
7#8"(), DiffEqCallbacks.ManifoldProjection{true,DiffEqCallbacks.NonAutonomo
usFunction{typeof(Main.##WeaveSandBox#311.g),true},DiffEqCallbacks.NLSOLVEJ
L_SETUP{0,true},Dict{Symbol,Float64}}(DiffEqCallbacks.NonAutonomousFunction
{typeof(Main.##WeaveSandBox#311.g),true}(Main.##WeaveSandBox#311.g, 0, 0), 
DiffEqCallbacks.NonAutonomousFunction{typeof(Main.##WeaveSandBox#311.g),tru
e}(Main.##WeaveSandBox#311.g, 0, 0), DiffEqCallbacks.NLSOLVEJL_SETUP{0,true
}(), Dict(:ftol => 1.0e-13)), DiffEqCallbacks.Manifold_initialize, Bool[0, 
1])
````





For the comparison we will use the following function

````julia
energy_err(sol) = map(i->H([sol[1,i], sol[2,i]], [sol[3,i], sol[4,i]],nothing)-E, 1:length(sol.u))
abs_energy_err(sol) = [abs.(H([sol[1,j], sol[2,j]], [sol[3,j], sol[4,j]],nothing) - E) for j=1:length(sol.u)]

function compare(mode=InPlace, all=true, plt=nothing; tmax=1e2)
    prob1 = DynamicalODEProblem(mode.ṗ,  mode.q̇, mode.p0, mode.q0, (0., tmax))
    prob2 = ODEProblem(mode.h_eqs, mode.u0, (0., tmax), mode.p)

    GC.gc()
    (mode == InPlace  && all) && @time sol1 = solve(prob2, Vern9(), callback=cb, abstol=1e-14, reltol=1e-14)
    GC.gc()
    @time sol2 = solve(prob1, KahanLi8(), dt=1e-2, maxiters=1e10)
    GC.gc()
    @time sol3 = solve(prob1, SofSpa10(), dt=1e-2, maxiters=1e8)
    GC.gc()
    @time sol4 = solve(prob2, Vern9(), abstol=1e-14, reltol=1e-14)
    GC.gc()
    @time sol5 = solve(prob1, DPRKN12(), abstol=1e-14, reltol=1e-14)
    GC.gc()
    (mode == InPlace && all) && @time sol6 = solve(prob2, TaylorMethod(50), abstol=1e-20)

    (mode == InPlace && all) && println("Vern9 + ManifoldProjection max energy error:\t"*
        "$(maximum(abs_energy_err(sol1)))\tin\t$(length(sol1.u))\tsteps.")
    println("KahanLi8 max energy error:\t\t\t$(maximum(abs_energy_err(sol2)))\tin\t$(length(sol2.u))\tsteps.")
    println("SofSpa10 max energy error:\t\t\t$(maximum(abs_energy_err(sol3)))\tin\t$(length(sol3.u))\tsteps.")
    println("Vern9 max energy error:\t\t\t\t$(maximum(abs_energy_err(sol4)))\tin\t$(length(sol4.u))\tsteps.")
    println("DPRKN12 max energy error:\t\t\t$(maximum(abs_energy_err(sol5)))\tin\t$(length(sol5.u))\tsteps.")
    (mode == InPlace && all) && println("TaylorMethod max energy error:\t\t\t$(maximum(abs_energy_err(sol6)))"*
        "\tin\t$(length(sol6.u))\tsteps.")

    if plt == nothing
        plt = plot(xlabel="t", ylabel="Energy error")
    end

    (mode == InPlace && all) && plot!(sol1.t, energy_err(sol1), label="Vern9 + ManifoldProjection")
    plot!(sol2.t, energy_err(sol2), label="KahanLi8", ls=mode==InPlace ? :solid : :dash)
    plot!(sol3.t, energy_err(sol3), label="SofSpa10", ls=mode==InPlace ? :solid : :dash)
    plot!(sol4.t, energy_err(sol4), label="Vern9", ls=mode==InPlace ? :solid : :dash)
    plot!(sol5.t, energy_err(sol5), label="DPRKN12", ls=mode==InPlace ? :solid : :dash)
    (mode == InPlace && all) && plot!(sol6.t, energy_err(sol6), label="TaylorMethod")

    return plt
end
````


````
compare (generic function with 4 methods)
````





The `mode` argument choses between the in place approach
and the out of place one. The `all` parameter is used to compare only the integrators that support both the in place and the out of place versions (we reffer here only to the 6 high order methods chosen bellow).
The `plt` argument can be used to overlay the results over a previous plot and the `tmax` keyword determines the simulation time.

Note:
1. The `Vern9` method is used with `ODEProblem` because of performance issues with `ArrayPartition` indexing which manifest for `DynamicalODEProblem`.
2. The `NLsolve` call used by `ManifoldProjection` was modified to use `ftol=1e-13` in order to obtain a very low energy error.

Here are the results of the comparisons between the in place methods:

````julia
compare(tmax=1e2)
````


````
Error: UndefVarError: h_eqs not defined
````



````julia
compare(tmax=1e3)
````


````
Error: UndefVarError: h_eqs not defined
````



````julia
compare(tmax=1e4)
````


````
Error: UndefVarError: h_eqs not defined
````



````julia
compare(tmax=2e4)
````


````
Error: UndefVarError: h_eqs not defined
````





As we can see from the above plots, we can achieve a very low energy error for long time simulation by manifold projection and with very high order Taylor methods. In comparison with the Hénon-Heiles system we see that as the Hamiltonian got more complex, the energy error for the other integration methods increased significantly.

We will now compare the in place with the out of place versions. In the plots bellow we will use a dashed line for the out of place versions.

````julia
function in_vs_out(;all=false, tmax=1e2)
    println("In place versions:")
    plt = compare(InPlace, all, tmax=tmax)
    println("\nOut of place versions:")
    plt = compare(OutOfPlace, false, plt; tmax=tmax)
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
Error: UndefVarError: h_eqs not defined
````





Now we will compare the in place and the out of place versions, but only for the integrators that are compatible with `StaticArrays`

````julia
in_vs_out(tmax=1e2)
````


````
In place versions:
Error: UndefVarError: h_eqs not defined
````



````julia
in_vs_out(tmax=1e3)
````


````
In place versions:
Error: UndefVarError: h_eqs not defined
````



````julia
in_vs_out(tmax=1e4)
````


````
In place versions:
Error: UndefVarError: h_eqs not defined
````



````julia
in_vs_out(tmax=2e4)
````


````
In place versions:
Error: UndefVarError: h_eqs not defined
````





As we see from the above comparisons, the `StaticArray` versions are significantly faster and use less memory. The speedup provided for the out of place version is more proeminent at larger values for `tmax`.
We can see again that if the simulation time is increased, the energy error of the symplectic methods is less noticeable compared to the rest of the methods.
In comparison with the Henon-Heiles case, we see that the symplectic methods are more competitive with `DPRKN12`.

````julia
using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])
````



## Appendix

These benchmarks are a part of the DiffEqBenchmarks.jl repository, found at: [https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl)

To locally run this tutorial, do the following commands:

```
using DiffEqBenchmarks
DiffEqBenchmarks.weave_file("DynamicalODE","Quadrupole_boson_Hamiltonian_energy_conservation_benchmark.jmd")
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
[65888b18-ceab-5e60-b2b9-181511a3b968] ParameterizedFunctions 5.3.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.5.3
[d330b81b-6aea-500a-939a-2ce795aea3ee] PyPlot 2.9.0
[90137ffa-7385-5640-81b9-e52037218182] StaticArrays 0.12.3
[92b13dbe-c966-51a2-8445-caca9f8a7d42] TaylorIntegration 0.8.3
[37e2e46d-f89d-539d-b4ee-838fcccc9c8e] LinearAlgebra 
[de0858da-6303-5e67-8744-51eddeeeb8d7] Printf 
[10745b16-79ce-11e8-11f9-7d13ad32a3b2] Statistics 
```

