+++
author = "Kirill Zubov, Zoe McCarthy, Yingbo Ma, Francesco Calisto, Valerio Pagliarino, Simone Azeglio, Luca Bottero, Emmanuel Luján, Valentin Sulzer, Ashutosh Bharambe, Nand Vinchhi, Kaushik Balakrishnan, Devesh Upadhyay, Chris Rackauckas"
title = "Poisson PDE Physics-Informed Neural Network (PINN) Optimizer Benchmarks"
+++


Adapted from [NeuralPDE: Automating Physics-Informed Neural Networks (PINNs) with Error Approximations](https://arxiv.org/abs/2107.09443).
Uses the [NeuralPDE.jl](https://neuralpde.sciml.ai/dev/) library from the
[SciML Scientific Machine Learning Open Source Organization](https://sciml.ai/)
for the implementation of physics-informed neural networks (PINNs) and other
science-guided AI techniques.

## Setup Code

```julia
using NeuralPDE, Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Quadrature,Cubature,Cuba
using Plots
```


```julia
function solve(opt)
    strategy = QuadratureTraining()

    @parameters x y
    @variables u(..)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2

    # 2D PDE
    eq  = Dxx(u(x,y)) + Dyy(u(x,y)) ~ -sin(pi*x)*sin(pi*y)

    # Boundary conditions
    bcs = [u(0,y) ~ 0.f0, u(1,y) ~ -sin(pi*1)*sin(pi*y),
           u(x,0) ~ 0.f0, u(x,1) ~ -sin(pi*x)*sin(pi*1)]
    # Space and time domains
    domains = [x ∈ IntervalDomain(0.0,1.0),
               y ∈ IntervalDomain(0.0,1.0)]

    # Neural network
    dim = 2 # number of dimensions
    chain = FastChain(FastDense(dim,16,Flux.σ),FastDense(16,16,Flux.σ),FastDense(16,1))

    discretization = PhysicsInformedNN(chain,strategy)

    indvars = [x, y]   #phisically independent variables
    depvars = [u]       #dependent (target) variable

    loss = []
    initial_time = 0

    times = []

    cb = function (p,l)
        if initial_time == 0
            initial_time = time()
        end
        push!(times, time() - initial_time)
        #println("Current loss for $opt is: $l")
        push!(loss, l)
        return false
    end

    @named pde_system = PDESystem(eq, bcs, domains, indvars, depvars)
    prob = discretize(pde_system, discretization)

    if opt == "both"
        res = GalacticOptim.solve(prob, ADAM(); cb = cb, maxiters=50)
        prob = remake(prob,u0=res.minimizer)
        res = GalacticOptim.solve(prob, BFGS(); cb = cb, maxiters=150)
    else
        res = GalacticOptim.solve(prob, opt; cb = cb, maxiters=200)
    end

    times[1] = 0.001

    return loss, times #add numeric solution
end
```

```
solve (generic function with 1 method)
```



```julia
opt1 = ADAM()
opt2 = ADAM(0.005)
opt3 = ADAM(0.05)
opt4 = RMSProp()
opt5 = RMSProp(0.005)
opt6 = RMSProp(0.05)
opt7 = Optim.BFGS()
opt8 = Optim.LBFGS()
```

```
Optim.LBFGS{Nothing, LineSearches.InitialStatic{Float64}, LineSearches.Hage
rZhang{Float64, Base.RefValue{Bool}}, Optim.var"#17#19"}(10, LineSearches.I
nitialStatic{Float64}
  alpha: Float64 1.0
  scaled: Bool false
, LineSearches.HagerZhang{Float64, Base.RefValue{Bool}}
  delta: Float64 0.1
  sigma: Float64 0.9
  alphamax: Float64 Inf
  rho: Float64 5.0
  epsilon: Float64 1.0e-6
  gamma: Float64 0.66
  linesearchmax: Int64 50
  psi3: Float64 0.1
  display: Int64 0
  mayterminate: Base.RefValue{Bool}
, nothing, Optim.var"#17#19"(), Optim.Flat(), true)
```





## Solve

```julia
loss_1, times_1 = solve(opt1)
loss_2, times_2 = solve(opt2)
loss_3, times_3 = solve(opt3)
loss_4, times_4 = solve(opt4)
loss_5, times_5 = solve(opt5)
loss_6, times_6 = solve(opt6)
loss_7, times_7 = solve(opt7)
loss_8, times_8 = solve(opt8)
loss_9, times_9 = solve("both")
```

```
(Any[104.2684735341356, 101.07773885593438, 97.94397489905161, 94.867988074
16685, 91.85047362340062, 88.89199402025474, 85.9930383485208, 83.154064923
7131, 80.37550967081319, 77.65774472404973  …  0.0003188404442962145, 0.000
2753817434812261, 0.00025846908100736715, 0.00022924711159700477, 0.0002119
7987118970702, 0.0001984253703549455, 0.00019508461948719214, 0.00019226966
042200556, 0.00018558273727754988, 0.00017566174663987923], Any[0.001, 0.01
857614517211914, 0.03719305992126465, 0.055773019790649414, 0.0743460655212
4023, 0.09290003776550293, 0.1114511489868164, 0.1299901008605957, 0.202987
1940612793, 0.22084498405456543  …  8.621284008026123, 8.6758131980896, 8.7
12090015411377, 8.74848198890686, 8.78489899635315, 8.821274042129517, 8.85
7802152633667, 8.894743204116821, 8.950502157211304, 9.006330013275146])
```





## Results

```julia
p = plot([times_1, times_2, times_3, times_4, times_5, times_6, times_7, times_8, times_9], [loss_1, loss_2, loss_3, loss_4, loss_5, loss_6, loss_7, loss_8, loss_9],xlabel="time (s)", ylabel="loss", xscale=:log10, yscale=:log10, labels=["ADAM(0.001)" "ADAM(0.005)" "ADAM(0.05)" "RMSProp(0.001)" "RMSProp(0.005)" "RMSProp(0.05)" "BFGS()" "LBFGS()" "ADAM + BFGS"], legend=:bottomleft, linecolor=["#2660A4" "#4CD0F4" "#FEC32F" "#F763CD" "#44BD79" "#831894" "#A6ED18" "#980000" "#FF912B"])
```

![](figures/poisson_5_1.png)

```julia
p = plot(1:201, [loss_1, loss_2, loss_3, loss_4, loss_5, loss_6, loss_7, loss_8, loss_9[2:end]], xlabel="iterations", ylabel="loss", yscale=:log10, labels=["ADAM(0.001)" "ADAM(0.005)" "ADAM(0.05)" "RMSProp(0.001)" "RMSProp(0.005)" "RMSProp(0.05)" "BFGS()" "LBFGS()" "ADAM + BFGS"], legend=:bottomleft, linecolor=["#2660A4" "#4CD0F4" "#FEC32F" "#F763CD" "#44BD79" "#831894" "#A6ED18" "#980000" "#FF912B"])
```

![](figures/poisson_6_1.png)

```julia
@show loss_1[end], loss_2[end], loss_3[end], loss_4[end], loss_5[end], loss_6[end], loss_7[end], loss_8[end], loss_9[end]
```

```
(loss_1[end], loss_2[end], loss_3[end], loss_4[end], loss_5[end], loss_6[en
d], loss_7[end], loss_8[end], loss_9[end]) = (3.732579140083826, 3.83871068
57725737, 0.0734829146939961, 3.8387786739370595, 2.3763009329942824, 0.090
72786483935653, 0.00010257991924182217, 0.02300612987128965, 0.000175661746
63987923)
(3.732579140083826, 3.8387106857725737, 0.0734829146939961, 3.8387786739370
595, 2.3763009329942824, 0.09072786483935653, 0.00010257991924182217, 0.023
00612987128965, 0.00017566174663987923)
```




## Appendix

These benchmarks are a part of the SciMLBenchmarks.jl repository, found at: [https://github.com/SciML/SciMLBenchmarks.jl](https://github.com/SciML/SciMLBenchmarks.jl). For more information on high-performance scientific machine learning, check out the SciML Open Source Software Organization [https://sciml.ai](https://sciml.ai).

To locally run this benchmark, do the following commands:

```
using SciMLBenchmarks
SciMLBenchmarks.weave_file("benchmarks/PINNOptimizers","poisson.jmd")
```

Computer Information:

```
Julia Version 1.6.5
Commit 9058264a69 (2021-12-19 12:30 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: AMD EPYC 7502 32-Core Processor
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, znver2)
Environment:
  BUILDKITE_PLUGIN_JULIA_CACHE_DIR = /cache/julia-buildkite-plugin
  JULIA_DEPOT_PATH = /cache/julia-buildkite-plugin/depots/5b300254-1738-4989-ae0a-f4d2d937f953

```

Package Information:

```
      Status `/cache/build/exclusive-amdci3-0/julialang/scimlbenchmarks-dot-jl/benchmarks/PINNOptimizers/Project.toml`
  [8a292aeb] Cuba v2.2.0
  [667455a9] Cubature v1.5.1
  [aae7a2af] DiffEqFlux v1.44.0
  [587475ba] Flux v0.12.8
  [a75be94c] GalacticOptim v2.2.0
  [961ee093] ModelingToolkit v6.7.1
  [315f7962] NeuralPDE v4.0.1
  [429524aa] Optim v1.5.0
  [91a5bcdd] Plots v1.24.2
  [67601950] Quadrature v1.12.0
  [31c91b34] SciMLBenchmarks v0.1.0
```

And the full manifest:

```
      Status `/cache/build/exclusive-amdci3-0/julialang/scimlbenchmarks-dot-jl/benchmarks/PINNOptimizers/Manifest.toml`
  [621f4979] AbstractFFTs v1.0.1
  [1520ce14] AbstractTrees v0.3.4
  [79e6a3ab] Adapt v3.3.1
  [ec485272] ArnoldiMethod v0.1.0
  [4fba245c] ArrayInterface v3.1.40
  [13072b0f] AxisAlgorithms v1.0.1
  [ab4f0b2a] BFloat16s v0.2.0
  [e2ed5e7c] Bijections v0.1.3
  [62783981] BitTwiddlingConvenienceFunctions v0.1.1
  [fa961155] CEnum v0.4.1
  [2a0fbf3d] CPUSummary v0.1.6
  [00ebfdb7] CSTParser v3.3.0
  [052768ef] CUDA v3.5.0
  [7057c7e9] Cassette v0.3.9
  [082447d4] ChainRules v1.14.0
  [d360d2e6] ChainRulesCore v1.11.1
  [9e997f8a] ChangesOfVariables v0.1.1
  [fb6a15b2] CloseOpenIntervals v0.1.4
  [944b1d66] CodecZlib v0.7.0
  [35d6a980] ColorSchemes v3.15.0
  [3da002f7] ColorTypes v0.11.0
  [5ae59095] Colors v0.12.8
  [861a8166] Combinatorics v1.0.2
  [a80b9123] CommonMark v0.8.3
  [38540f10] CommonSolve v0.2.0
  [bbf7d656] CommonSubexpressions v0.3.0
  [34da2185] Compat v3.40.0
  [b152e2b5] CompositeTypes v0.1.2
  [8f4d0f93] Conda v1.5.2
  [88cd18e8] ConsoleProgressMonitor v0.1.2
  [187b0558] ConstructionBase v1.3.0
  [d38c429a] Contour v0.5.7
  [a8cc5b0e] Crayons v4.0.4
  [8a292aeb] Cuba v2.2.0
  [667455a9] Cubature v1.5.1
  [754358af] DEDataArrays v0.2.0
  [9a962f9c] DataAPI v1.9.0
  [82cc6244] DataInterpolations v3.6.1
  [864edb3b] DataStructures v0.18.10
  [e2d170a0] DataValueInterfaces v1.0.0
  [b429d917] DensityInterface v0.4.0
  [2b5f629d] DiffEqBase v6.76.0
  [459566f4] DiffEqCallbacks v2.17.0
  [aae7a2af] DiffEqFlux v1.44.0
  [c894b116] DiffEqJump v7.3.1
  [77a26b50] DiffEqNoiseProcess v5.9.0
  [41bf760c] DiffEqSensitivity v6.60.3
  [163ba53b] DiffResults v1.0.3
  [b552c78f] DiffRules v1.5.0
  [b4f34e82] Distances v0.10.6
  [31c24e10] Distributions v0.25.34
  [ced4e74d] DistributionsAD v0.6.34
  [ffbed154] DocStringExtensions v0.8.6
  [5b8099bc] DomainSets v0.5.9
  [7c1d4256] DynamicPolynomials v0.3.21
  [da5c29d0] EllipsisNotation v1.1.1
  [7da242da] Enzyme v0.7.2
  [d4d017d3] ExponentialUtilities v1.10.2
  [e2ba6199] ExprTools v0.1.6
  [c87230d0] FFMPEG v0.4.1
  [7a1cc6ca] FFTW v1.4.5
  [7034ab61] FastBroadcast v0.1.11
  [9aa1b823] FastClosures v0.3.2
  [1a297f60] FillArrays v0.12.7
  [6a86dc24] FiniteDiff v2.8.1
  [53c48c17] FixedPointNumbers v0.8.4
  [587475ba] Flux v0.12.8
  [59287772] Formatting v0.4.2
  [f6369f11] ForwardDiff v0.10.23
  [069b7b12] FunctionWrappers v1.1.2
  [d9f16b24] Functors v0.2.7
  [0c68f7d7] GPUArrays v8.1.2
  [61eb1bfa] GPUCompiler v0.13.8
  [28b8d3ca] GR v0.62.1
  [a75be94c] GalacticOptim v2.2.0
  [5c1252a2] GeometryBasics v0.4.1
  [d7ba0133] Git v1.2.1
  [af5da776] GlobalSensitivity v1.2.2
  [86223c79] Graphs v1.4.1
  [42e2da0e] Grisu v1.0.2
  [19dc6840] HCubature v1.5.0
  [cd3eb016] HTTP v0.9.17
  [eafb193a] Highlights v0.4.5
  [3e5b6fbb] HostCPUFeatures v0.1.5
  [0e44f5e4] Hwloc v2.0.0
  [7073ff75] IJulia v1.23.2
  [7869d1d1] IRTools v0.4.4
  [615f187c] IfElse v0.1.1
  [d25df0c9] Inflate v0.1.2
  [83e8ac13] IniFile v0.5.0
  [a98d9a8b] Interpolations v0.13.4
  [8197267c] IntervalSets v0.5.3
  [3587e190] InverseFunctions v0.1.2
  [92d709cd] IrrationalConstants v0.1.1
  [c8e1da08] IterTools v1.3.0
  [42fd0dbc] IterativeSolvers v0.9.2
  [82899510] IteratorInterfaceExtensions v1.0.0
  [692b3bcd] JLLWrappers v1.3.0
  [682c06a0] JSON v0.21.2
  [98e50ef6] JuliaFormatter v0.18.1
  [e5e0dc1b] Juno v0.8.4
  [5ab0869b] KernelDensity v0.6.3
  [929cbde3] LLVM v4.7.0
  [b964fa9f] LaTeXStrings v1.3.0
  [2ee39098] LabelledArrays v1.6.7
  [23fbe1c1] Latexify v0.15.9
  [a5e1c1ea] LatinHypercubeSampling v1.8.0
  [73f95e8e] LatticeRules v0.0.1
  [10f19ff3] LayoutPointers v0.1.4
  [1d6d02ad] LeftChildRightSiblingTrees v0.1.2
  [093fc24a] LightGraphs v1.3.5
  [d3d80556] LineSearches v7.1.1
  [2ab3a3ac] LogExpFunctions v0.3.5
  [e6f89c97] LoggingExtras v0.4.7
  [bdcacae8] LoopVectorization v0.12.98
  [1914dd2f] MacroTools v0.5.9
  [d125e4d3] ManualMemory v0.1.6
  [739be429] MbedTLS v1.0.3
  [442fdcdd] Measures v0.3.1
  [e89f7d12] Media v0.5.0
  [e1d29d7a] Missings v1.0.2
  [961ee093] ModelingToolkit v6.7.1
  [4886b29c] MonteCarloIntegration v0.0.3
  [46d2c3a1] MuladdMacro v0.2.2
  [102ac46a] MultivariatePolynomials v0.3.18
  [ffc61752] Mustache v1.0.12
  [d8a4904e] MutableArithmetics v0.2.22
  [d41bc354] NLSolversBase v7.8.2
  [2774e3e8] NLsolve v4.5.1
  [872c559c] NNlib v0.7.31
  [a00861dc] NNlibCUDA v0.1.10
  [77ba4419] NaNMath v0.3.5
  [315f7962] NeuralPDE v4.0.1
  [8913a72c] NonlinearSolve v0.3.11
  [d8793406] ObjectFile v0.3.7
  [6fe1bfb0] OffsetArrays v1.10.8
  [429524aa] Optim v1.5.0
  [bac558e1] OrderedCollections v1.4.1
  [1dea7af3] OrdinaryDiffEq v5.68.0
  [90014a1f] PDMats v0.11.5
  [d96e819e] Parameters v0.12.3
  [69de0a69] Parsers v2.1.2
  [ccf2f8ad] PlotThemes v2.0.1
  [995b91a9] PlotUtils v1.0.15
  [91a5bcdd] Plots v1.24.2
  [e409e4f3] PoissonRandom v0.4.0
  [f517fe37] Polyester v0.5.4
  [1d0040c9] PolyesterWeave v0.1.2
  [85a6dd25] PositiveFactorizations v0.2.4
  [d236fae5] PreallocationTools v0.2.0
  [21216c6a] Preferences v1.2.2
  [33c8b6b6] ProgressLogging v0.1.4
  [92933f4c] ProgressMeter v1.7.1
  [1fd47b50] QuadGK v2.4.2
  [67601950] Quadrature v1.12.0
  [8a4e6c94] QuasiMonteCarlo v0.2.3
  [74087812] Random123 v1.4.2
  [e6cf234a] RandomNumbers v1.5.3
  [c84ed2f1] Ratios v0.4.2
  [c1ae055f] RealDot v0.1.0
  [3cdcf5f2] RecipesBase v1.2.1
  [01d81517] RecipesPipeline v0.4.1
  [731186ca] RecursiveArrayTools v2.20.0
  [f2c3362d] RecursiveFactorization v0.2.5
  [189a3867] Reexport v1.2.2
  [05181044] RelocatableFolders v0.1.3
  [ae029012] Requires v1.1.3
  [ae5879a3] ResettableStacks v1.1.1
  [37e2e3b7] ReverseDiff v1.10.0
  [79098fc4] Rmath v0.7.0
  [7e49a35a] RuntimeGeneratedFunctions v0.5.3
  [3cdde19b] SIMDDualNumbers v0.1.0
  [94e857df] SIMDTypes v0.1.0
  [476501e8] SLEEFPirates v0.6.28
  [1bc83da4] SafeTestsets v0.0.1
  [0bca4576] SciMLBase v1.19.5
  [31c91b34] SciMLBenchmarks v0.1.0
  [6c6a2e73] Scratch v1.1.0
  [efcf1570] Setfield v0.8.0
  [992d4aef] Showoff v1.0.3
  [699a6c99] SimpleTraits v0.9.4
  [ed01d8cd] Sobol v1.5.0
  [b85f4697] SoftGlobalScope v1.1.0
  [a2af1166] SortingAlgorithms v1.0.1
  [47a9eef4] SparseDiffTools v1.18.1
  [276daf66] SpecialFunctions v1.8.1
  [860ef19b] StableRNGs v1.0.0
  [aedffcd0] Static v0.4.0
  [90137ffa] StaticArrays v1.2.13
  [82ae8749] StatsAPI v1.1.0
  [2913bbd2] StatsBase v0.33.13
  [4c63d2b9] StatsFuns v0.9.14
  [789caeaf] StochasticDiffEq v6.41.0
  [7792a7ef] StrideArraysCore v0.2.9
  [69024149] StringEncodings v0.3.5
  [09ab397b] StructArrays v0.6.3
  [53d494c1] StructIO v0.3.0
  [d1185830] SymbolicUtils v0.16.0
  [0c5d862f] Symbolics v3.5.1
  [3783bdb8] TableTraits v1.0.1
  [bd369af6] Tables v1.6.0
  [8ea1fca8] TermInterface v0.1.8
  [5d786b92] TerminalLoggers v0.1.5
  [8290d209] ThreadingUtilities v0.4.6
  [a759f4b9] TimerOutputs v0.5.13
  [0796e94c] Tokenize v0.5.21
  [9f7883ad] Tracker v0.2.16
  [3bb67fe8] TranscodingStreams v0.9.6
  [592b5752] Trapz v2.0.3
  [a2a6695c] TreeViews v0.3.0
  [d5829a12] TriangularSolve v0.1.8
  [5c2747f8] URIs v1.3.0
  [3a884ed6] UnPack v1.0.2
  [1cfade01] UnicodeFun v0.4.1
  [1986cc42] Unitful v1.9.2
  [3d5dd08c] VectorizationBase v0.21.21
  [81def892] VersionParsing v1.2.1
  [19fa3120] VertexSafeGraphs v0.2.0
  [44d3d7a6] Weave v0.10.10
  [efce3f68] WoodburyMatrices v0.5.5
  [ddb6d928] YAML v0.4.7
  [c2297ded] ZMQ v1.2.1
  [a5390f91] ZipFile v0.9.4
  [e88e6eb3] Zygote v0.6.32
  [700de1a5] ZygoteRules v0.2.2
  [6e34b625] Bzip2_jll v1.0.6+5
  [83423d85] Cairo_jll v1.16.0+6
  [3bed1096] Cuba_jll v4.2.2+0
  [7bc98958] Cubature_jll v1.0.5+0
  [5ae413db] EarCut_jll v2.2.3+0
  [7cc45869] Enzyme_jll v0.0.22+0
  [2e619515] Expat_jll v2.2.10+0
  [b22a6f82] FFMPEG_jll v4.3.1+4
  [f5851436] FFTW_jll v3.3.10+0
  [a3f928ae] Fontconfig_jll v2.13.1+14
  [d7e528f0] FreeType2_jll v2.10.1+5
  [559328eb] FriBidi_jll v1.0.10+0
  [0656b61e] GLFW_jll v3.3.5+1
  [d2c73de3] GR_jll v0.58.1+0
  [78b55507] Gettext_jll v0.20.1+7
  [f8c6e375] Git_jll v2.31.0+0
  [7746bdde] Glib_jll v2.59.0+4
  [e33a78d0] Hwloc_jll v2.5.0+0
  [1d5cc7b8] IntelOpenMP_jll v2018.0.3+2
  [aacddb02] JpegTurbo_jll v2.1.0+0
  [c1c5ebd0] LAME_jll v3.100.1+0
  [dad2f222] LLVMExtra_jll v0.0.13+0
  [dd4b983a] LZO_jll v2.10.1+0
  [dd192d2f] LibVPX_jll v1.10.0+0
  [e9f186c6] Libffi_jll v3.2.2+1
  [d4300ac3] Libgcrypt_jll v1.8.7+0
  [7e76a0d4] Libglvnd_jll v1.3.0+3
  [7add5ba3] Libgpg_error_jll v1.42.0+0
  [94ce4f54] Libiconv_jll v1.16.1+1
  [4b2f31a3] Libmount_jll v2.35.0+0
  [89763e89] Libtiff_jll v4.3.0+0
  [38a345b3] Libuuid_jll v2.36.0+0
  [856f044c] MKL_jll v2021.1.1+2
  [e7412a2a] Ogg_jll v1.3.5+0
  [458c3c95] OpenSSL_jll v1.1.10+0
  [efe28fd5] OpenSpecFun_jll v0.5.5+0
  [91d4177d] Opus_jll v1.3.2+0
  [2f80f16e] PCRE_jll v8.44.0+0
  [30392449] Pixman_jll v0.40.1+0
  [ea2cea3b] Qt5Base_jll v5.15.2+0
  [f50d1b31] Rmath_jll v0.3.0+0
  [a2964d1f] Wayland_jll v1.19.0+0
  [2381bf8a] Wayland_protocols_jll v1.23.0+0
  [02c8fc9c] XML2_jll v2.9.12+0
  [aed1982a] XSLT_jll v1.1.34+0
  [4f6342f7] Xorg_libX11_jll v1.6.9+4
  [0c0b7dd1] Xorg_libXau_jll v1.0.9+4
  [935fb764] Xorg_libXcursor_jll v1.2.0+4
  [a3789734] Xorg_libXdmcp_jll v1.1.3+4
  [1082639a] Xorg_libXext_jll v1.3.4+4
  [d091e8ba] Xorg_libXfixes_jll v5.0.3+4
  [a51aa0fd] Xorg_libXi_jll v1.7.10+4
  [d1454406] Xorg_libXinerama_jll v1.1.4+4
  [ec84b674] Xorg_libXrandr_jll v1.5.2+4
  [ea2f1a96] Xorg_libXrender_jll v0.9.10+4
  [14d82f49] Xorg_libpthread_stubs_jll v0.1.0+3
  [c7cfdc94] Xorg_libxcb_jll v1.13.0+3
  [cc61e674] Xorg_libxkbfile_jll v1.1.0+4
  [12413925] Xorg_xcb_util_image_jll v0.4.0+1
  [2def613f] Xorg_xcb_util_jll v0.4.0+1
  [975044d2] Xorg_xcb_util_keysyms_jll v0.4.0+1
  [0d47668e] Xorg_xcb_util_renderutil_jll v0.3.9+1
  [c22f9ab0] Xorg_xcb_util_wm_jll v0.4.1+1
  [35661453] Xorg_xkbcomp_jll v1.4.2+4
  [33bec58e] Xorg_xkeyboard_config_jll v2.27.0+4
  [c5fb5394] Xorg_xtrans_jll v1.4.0+3
  [8f1865be] ZeroMQ_jll v4.3.4+0
  [3161d3a3] Zstd_jll v1.5.0+0
  [0ac62f75] libass_jll v0.14.0+4
  [f638f0a6] libfdk_aac_jll v0.1.6+4
  [b53b4c65] libpng_jll v1.6.38+0
  [a9144af2] libsodium_jll v1.0.20+0
  [f27f6e37] libvorbis_jll v1.3.7+0
  [1270edf5] x264_jll v2020.7.14+2
  [dfaa095f] x265_jll v3.0.0+3
  [d8fb68d0] xkbcommon_jll v0.9.1+5
  [0dad84c5] ArgTools
  [56f22d72] Artifacts
  [2a0f44e3] Base64
  [ade2ca70] Dates
  [8bb1440f] DelimitedFiles
  [8ba89e20] Distributed
  [f43a241f] Downloads
  [7b1f6079] FileWatching
  [9fa8497b] Future
  [b77e0a4c] InteractiveUtils
  [4af54fe1] LazyArtifacts
  [b27032c2] LibCURL
  [76f85450] LibGit2
  [8f399da3] Libdl
  [37e2e46d] LinearAlgebra
  [56ddb016] Logging
  [d6f4376e] Markdown
  [a63ad114] Mmap
  [ca575930] NetworkOptions
  [44cfe95a] Pkg
  [de0858da] Printf
  [9abbd945] Profile
  [3fa0cd96] REPL
  [9a3f8284] Random
  [ea8e919c] SHA
  [9e88b42a] Serialization
  [1a1011a3] SharedArrays
  [6462fe0b] Sockets
  [2f01184e] SparseArrays
  [10745b16] Statistics
  [4607b0f0] SuiteSparse
  [fa267f1f] TOML
  [a4e569a6] Tar
  [8dfed614] Test
  [cf7118a7] UUIDs
  [4ec0a83e] Unicode
  [e66e0078] CompilerSupportLibraries_jll
  [deac9b47] LibCURL_jll
  [29816b5a] LibSSH2_jll
  [c8ffd9c3] MbedTLS_jll
  [14a3606d] MozillaCACerts_jll
  [05823500] OpenLibm_jll
  [efcefdf7] PCRE2_jll
  [83775a58] Zlib_jll
  [8e850ede] nghttp2_jll
  [3f19e933] p7zip_jll
```

