+++
author = "Chris Rackauckas"
title = "OREGO Differential-Algebraic Equation (DAE) Work-Precision Diagrams"
+++
```julia
using OrdinaryDiffEq, DiffEqDevTools, Sundials, ModelingToolkit, ODEInterfaceDiffEq,
      Plots, DASSL, DASKR
using LinearAlgebra

@variables t y1(t)=1.0 y2(t)=2.0 y3(t)=3.0
@parameters p1=77.27 p2=8.375e-6 p3=0.161
D = Differential(t)
eqs = [
  D(y1) ~ p1*(y2+y1*(1-p2*y1-y2))
  D(y2) ~ (y3-(1+y1)*y2)/p1
  D(y3) ~ p3*(y1-y3)
]
@named sys = ODESystem(eqs)
simpsys = structural_simplify(sys)
mmprob = ODEProblem(sys,[],(0.0,30.0))
daeprob = DAEProblem(sys,[D(y1)=>77.26935286375,
                              D(y2)=>-0.012941633234114146,
                              D(y3)=>-0.322],[],(0.0,30.0))
odaeprob = ODAEProblem(simpsys,[],(0.0,30.0))

ref_sol = solve(daeprob,IDA(),abstol=1/10^14,reltol=1/10^14);
ode_ref_sol = solve(odaeprob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14);

probs = [mmprob,daeprob,odaeprob]
refs = [ref_sol,ref_sol,ode_ref_sol];
```


```julia
plot(ref_sol)
```

![](figures/OregoDAE_2_1.png)



## High Tolerances

```julia
abstols = 1.0 ./ 10.0 .^ (6:9)
reltols = 1.0 ./ 10.0 .^ (2:5);
setups = [Dict(:prob_choice => 1, :alg=>Rosenbrock23()),
          Dict(:prob_choice => 1, :alg=>Rodas4()),
          #Dict(:prob_choice => 1, :alg=>FBDF()),
          Dict(:prob_choice => 1, :alg=>QNDF()),
          Dict(:prob_choice => 1, :alg=>rodas()),
          Dict(:prob_choice => 1, :alg=>radau()),
          Dict(:prob_choice => 1, :alg=>RadauIIA5()),
          Dict(:prob_choice => 2, :alg=>DFBDF()),
          Dict(:prob_choice => 2, :alg=>IDA()),
]

wp = WorkPrecisionSet(probs,abstols,reltols,setups;print_names=true,
                      save_everystep=false,appxsol=refs,maxiters=Int(1e5),numruns=10)
plot(wp)
```

```
Rosenbrock23
Rodas4
QNDF
rodas
radau
RadauIIA5
DFBDF
IDA
```


![](figures/OregoDAE_3_1.png)

```julia
setups = [Dict(:prob_choice => 1, :alg=>Rosenbrock23()),
          Dict(:prob_choice => 1, :alg=>Rodas4()),
          Dict(:prob_choice => 2, :alg=>IDA()),
          Dict(:prob_choice => 3, :alg=>Rosenbrock23()),
          Dict(:prob_choice => 3, :alg=>Rodas4()),
          Dict(:prob_choice => 3, :alg=>CVODE_BDF()),
          Dict(:prob_choice => 3, :alg=>TRBDF2()),
          Dict(:prob_choice => 3, :alg=>KenCarp4()),
          ]
wp = WorkPrecisionSet(probs,abstols,reltols,setups;
                      save_everystep=false,appxsol=refs,maxiters=Int(1e5),numruns=10)
plot(wp)
```

![](figures/OregoDAE_4_1.png)

```julia
abstols = 1.0 ./ 10.0 .^ (6:8)
reltols = 1.0 ./ 10.0 .^ (2:4);
setups = [Dict(:prob_choice => 3, :alg=>Rosenbrock23()),
          Dict(:prob_choice => 3, :alg=>Rodas4()),
          Dict(:prob_choice => 2, :alg=>IDA()),
          Dict(:prob_choice => 2, :alg=>DASSL.dassl()),
          Dict(:prob_choice => 2, :alg=>DASKR.daskr()),
          ]
wp = WorkPrecisionSet(probs,abstols,reltols,setups;
                      save_everystep=false,appxsol=refs,maxiters=Int(1e5),numruns=10)
plot(wp)
```

![](figures/OregoDAE_5_1.png)



### Timeseries Errors

```julia
abstols = 1.0 ./ 10.0 .^ (6:9)
reltols = 1.0 ./ 10.0 .^ (2:5);
setups = [Dict(:prob_choice => 1, :alg=>Rosenbrock23()),
          Dict(:prob_choice => 1, :alg=>Rodas4()),
          #Dict(:prob_choice => 1, :alg=>FBDF()), # too slow
          Dict(:prob_choice => 1, :alg=>QNDF()),
          Dict(:prob_choice => 1, :alg=>rodas()),
          Dict(:prob_choice => 1, :alg=>radau()),
          Dict(:prob_choice => 1, :alg=>RadauIIA5()),
          Dict(:prob_choice => 2, :alg=>DFBDF()),
          Dict(:prob_choice => 2, :alg=>IDA()),
          ]
gr()
wp = WorkPrecisionSet(probs,abstols,reltols,setups;error_estimate = :l2,
                      save_everystep=false,appxsol=refs,maxiters=Int(1e5),numruns=10)
plot(wp)
```

![](figures/OregoDAE_6_1.png)

```julia
abstols = 1.0 ./ 10.0 .^ (6:9)
reltols = 1.0 ./ 10.0 .^ (2:5);
setups = [Dict(:prob_choice => 1, :alg=>Rosenbrock23()),
          Dict(:prob_choice => 1, :alg=>Rodas4()),
          Dict(:prob_choice => 2, :alg=>IDA()),
          Dict(:prob_choice => 3, :alg=>Rosenbrock23()),
          Dict(:prob_choice => 3, :alg=>Rodas4()),
          Dict(:prob_choice => 3, :alg=>CVODE_BDF()),
          Dict(:prob_choice => 3, :alg=>TRBDF2()),
          Dict(:prob_choice => 3, :alg=>KenCarp4()),
          ]
wp = WorkPrecisionSet(probs,abstols,reltols,setups;error_estimate = :l2,
                      save_everystep=false,appxsol=refs,maxiters=Int(1e5),numruns=10)
plot(wp)
```

![](figures/OregoDAE_7_1.png)



### Low Tolerances

This is the speed at lower tolerances, measuring what's good when accuracy is needed.

```julia
abstols = 1.0 ./ 10.0 .^ (7:12)
reltols = 1.0 ./ 10.0 .^ (4:9)

setups = [Dict(:prob_choice => 1, :alg=>Rodas5()),
          Dict(:prob_choice => 3, :alg=>Rodas5()),
          Dict(:prob_choice => 1, :alg=>Rodas4()),
          Dict(:prob_choice => 3, :alg=>Rodas4()),
          #Dict(:prob_choice => 1, :alg=>FBDF()),
          Dict(:prob_choice => 1, :alg=>QNDF()),
          Dict(:prob_choice => 1, :alg=>rodas()),
          Dict(:prob_choice => 1, :alg=>radau()),
          Dict(:prob_choice => 1, :alg=>RadauIIA5()),
          Dict(:prob_choice => 2, :alg=>DFBDF()),
          Dict(:prob_choice => 2, :alg=>IDA()),
          ]
gr()
wp = WorkPrecisionSet(probs,abstols,reltols,setups;
                      save_everystep=false,appxsol=refs,maxiters=Int(1e5),numruns=10)
plot(wp)
```

![](figures/OregoDAE_8_1.png)

```julia
wp = WorkPrecisionSet(probs,abstols,reltols,setups;error_estimate = :l2,
                      save_everystep=false,appxsol=refs,maxiters=Int(1e5),numruns=10)
plot(wp)
```

![](figures/OregoDAE_9_1.png)



### Conclusion


## Appendix

These benchmarks are a part of the SciMLBenchmarks.jl repository, found at: [https://github.com/SciML/SciMLBenchmarks.jl](https://github.com/SciML/SciMLBenchmarks.jl). For more information on high-performance scientific machine learning, check out the SciML Open Source Software Organization [https://sciml.ai](https://sciml.ai).

To locally run this benchmark, do the following commands:

```
using SciMLBenchmarks
SciMLBenchmarks.weave_file("benchmarks/DAE","OregoDAE.jmd")
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
      Status `/cache/build/exclusive-amdci3-0/julialang/scimlbenchmarks-dot-jl/benchmarks/DAE/Project.toml`
  [165a45c3] DASKR v2.8.0
  [e993076c] DASSL v2.6.0
  [f3b72e0c] DiffEqDevTools v2.28.0
  [961ee093] ModelingToolkit v7.1.3
  [09606e27] ODEInterfaceDiffEq v3.10.1
  [1dea7af3] OrdinaryDiffEq v5.70.0
  [91a5bcdd] Plots v1.25.2
  [31c91b34] SciMLBenchmarks v0.1.0
  [c3572dad] Sundials v4.6.0
```

And the full manifest:

```
      Status `/cache/build/exclusive-amdci3-0/julialang/scimlbenchmarks-dot-jl/benchmarks/DAE/Manifest.toml`
  [1520ce14] AbstractTrees v0.3.4
  [79e6a3ab] Adapt v3.3.1
  [dce04be8] ArgCheck v2.1.0
  [ec485272] ArnoldiMethod v0.2.0
  [4fba245c] ArrayInterface v3.2.1
  [15f4f7f2] AutoHashEquals v0.2.0
  [198e06fe] BangBang v0.3.32
  [9718e550] Baselet v0.1.1
  [e2ed5e7c] Bijections v0.1.3
  [62783981] BitTwiddlingConvenienceFunctions v0.1.1
  [fa961155] CEnum v0.4.1
  [2a0fbf3d] CPUSummary v0.1.6
  [00ebfdb7] CSTParser v3.3.0
  [d360d2e6] ChainRulesCore v1.11.2
  [9e997f8a] ChangesOfVariables v0.1.2
  [fb6a15b2] CloseOpenIntervals v0.1.4
  [35d6a980] ColorSchemes v3.15.0
  [3da002f7] ColorTypes v0.11.0
  [5ae59095] Colors v0.12.8
  [861a8166] Combinatorics v1.0.2
  [a80b9123] CommonMark v0.8.4
  [38540f10] CommonSolve v0.2.0
  [bbf7d656] CommonSubexpressions v0.3.0
  [34da2185] Compat v3.40.0
  [b152e2b5] CompositeTypes v0.1.2
  [a33af91c] CompositionsBase v0.1.1
  [8f4d0f93] Conda v1.6.0
  [187b0558] ConstructionBase v1.3.0
  [d38c429a] Contour v0.5.7
  [a8cc5b0e] Crayons v4.0.4
  [165a45c3] DASKR v2.8.0
  [e993076c] DASSL v2.6.0
  [754358af] DEDataArrays v0.2.0
  [9a962f9c] DataAPI v1.9.0
  [864edb3b] DataStructures v0.18.11
  [e2d170a0] DataValueInterfaces v1.0.0
  [244e2a9f] DefineSingletons v0.1.1
  [b429d917] DensityInterface v0.4.0
  [2b5f629d] DiffEqBase v6.77.0
  [459566f4] DiffEqCallbacks v2.18.0
  [f3b72e0c] DiffEqDevTools v2.28.0
  [c894b116] DiffEqJump v8.0.0
  [77a26b50] DiffEqNoiseProcess v5.9.0
  [163ba53b] DiffResults v1.0.3
  [b552c78f] DiffRules v1.6.0
  [b4f34e82] Distances v0.10.7
  [31c24e10] Distributions v0.25.35
  [ffbed154] DocStringExtensions v0.8.6
  [5b8099bc] DomainSets v0.5.9
  [7c1d4256] DynamicPolynomials v0.3.21
  [da5c29d0] EllipsisNotation v1.1.3
  [d4d017d3] ExponentialUtilities v1.10.2
  [e2ba6199] ExprTools v0.1.6
  [c87230d0] FFMPEG v0.4.1
  [7034ab61] FastBroadcast v0.1.11
  [9aa1b823] FastClosures v0.3.2
  [1a297f60] FillArrays v0.12.7
  [6a86dc24] FiniteDiff v2.8.1
  [53c48c17] FixedPointNumbers v0.8.4
  [59287772] Formatting v0.4.2
  [f6369f11] ForwardDiff v0.10.23
  [069b7b12] FunctionWrappers v1.1.2
  [28b8d3ca] GR v0.62.1
  [5c1252a2] GeometryBasics v0.4.1
  [d7ba0133] Git v1.2.1
  [86223c79] Graphs v1.4.1
  [42e2da0e] Grisu v1.0.2
  [cd3eb016] HTTP v0.9.17
  [eafb193a] Highlights v0.4.5
  [3e5b6fbb] HostCPUFeatures v0.1.5
  [0e44f5e4] Hwloc v2.0.0
  [7073ff75] IJulia v1.23.2
  [615f187c] IfElse v0.1.1
  [d25df0c9] Inflate v0.1.2
  [83e8ac13] IniFile v0.5.0
  [22cec73e] InitialValues v0.2.11
  [8197267c] IntervalSets v0.5.3
  [3587e190] InverseFunctions v0.1.2
  [92d709cd] IrrationalConstants v0.1.1
  [c8e1da08] IterTools v1.4.0
  [42fd0dbc] IterativeSolvers v0.9.2
  [82899510] IteratorInterfaceExtensions v1.0.0
  [692b3bcd] JLLWrappers v1.3.0
  [682c06a0] JSON v0.21.2
  [98e50ef6] JuliaFormatter v0.19.2
  [b964fa9f] LaTeXStrings v1.3.0
  [2ee39098] LabelledArrays v1.6.7
  [23fbe1c1] Latexify v0.15.9
  [10f19ff3] LayoutPointers v0.1.4
  [d3d80556] LineSearches v7.1.1
  [2ab3a3ac] LogExpFunctions v0.3.5
  [bdcacae8] LoopVectorization v0.12.99
  [1914dd2f] MacroTools v0.5.9
  [d125e4d3] ManualMemory v0.1.6
  [739be429] MbedTLS v1.0.3
  [442fdcdd] Measures v0.3.1
  [e9d8d322] Metatheory v1.3.2
  [128add7d] MicroCollections v0.1.1
  [e1d29d7a] Missings v1.0.2
  [961ee093] ModelingToolkit v7.1.3
  [46d2c3a1] MuladdMacro v0.2.2
  [102ac46a] MultivariatePolynomials v0.3.18
  [ffc61752] Mustache v1.0.12
  [d8a4904e] MutableArithmetics v0.2.22
  [d41bc354] NLSolversBase v7.8.2
  [2774e3e8] NLsolve v4.5.1
  [77ba4419] NaNMath v0.3.5
  [8913a72c] NonlinearSolve v0.3.12
  [54ca160b] ODEInterface v0.5.0
  [09606e27] ODEInterfaceDiffEq v3.10.1
  [6fe1bfb0] OffsetArrays v1.10.8
  [429524aa] Optim v1.5.0
  [bac558e1] OrderedCollections v1.4.1
  [1dea7af3] OrdinaryDiffEq v5.70.0
  [90014a1f] PDMats v0.11.5
  [d96e819e] Parameters v0.12.3
  [69de0a69] Parsers v2.1.2
  [ccf2f8ad] PlotThemes v2.0.1
  [995b91a9] PlotUtils v1.0.15
  [91a5bcdd] Plots v1.25.2
  [e409e4f3] PoissonRandom v0.4.0
  [f517fe37] Polyester v0.5.4
  [1d0040c9] PolyesterWeave v0.1.2
  [85a6dd25] PositiveFactorizations v0.2.4
  [d236fae5] PreallocationTools v0.2.1
  [21216c6a] Preferences v1.2.2
  [1fd47b50] QuadGK v2.4.2
  [74087812] Random123 v1.4.2
  [e6cf234a] RandomNumbers v1.5.3
  [3cdcf5f2] RecipesBase v1.2.1
  [01d81517] RecipesPipeline v0.4.1
  [731186ca] RecursiveArrayTools v2.20.0
  [f2c3362d] RecursiveFactorization v0.2.5
  [189a3867] Reexport v1.2.2
  [42d2dcc6] Referenceables v0.1.2
  [05181044] RelocatableFolders v0.1.3
  [ae029012] Requires v1.1.3
  [ae5879a3] ResettableStacks v1.1.1
  [79098fc4] Rmath v0.7.0
  [47965b36] RootedTrees v1.2.1
  [7e49a35a] RuntimeGeneratedFunctions v0.5.3
  [3cdde19b] SIMDDualNumbers v0.1.0
  [94e857df] SIMDTypes v0.1.0
  [476501e8] SLEEFPirates v0.6.28
  [1bc83da4] SafeTestsets v0.0.1
  [0bca4576] SciMLBase v1.22.0
  [31c91b34] SciMLBenchmarks v0.1.0
  [6c6a2e73] Scratch v1.1.0
  [efcf1570] Setfield v0.8.1
  [992d4aef] Showoff v1.0.3
  [699a6c99] SimpleTraits v0.9.4
  [b85f4697] SoftGlobalScope v1.1.0
  [a2af1166] SortingAlgorithms v1.0.1
  [47a9eef4] SparseDiffTools v1.18.3
  [276daf66] SpecialFunctions v1.8.1
  [171d559e] SplittablesBase v0.1.14
  [aedffcd0] Static v0.4.0
  [90137ffa] StaticArrays v1.2.13
  [82ae8749] StatsAPI v1.1.0
  [2913bbd2] StatsBase v0.33.13
  [4c63d2b9] StatsFuns v0.9.14
  [7792a7ef] StrideArraysCore v0.2.9
  [69024149] StringEncodings v0.3.5
  [09ab397b] StructArrays v0.6.3
  [c3572dad] Sundials v4.6.0
  [d1185830] SymbolicUtils v0.18.1
  [0c5d862f] Symbolics v4.2.0
  [3783bdb8] TableTraits v1.0.1
  [bd369af6] Tables v1.6.0
  [8ea1fca8] TermInterface v0.2.3
  [8290d209] ThreadingUtilities v0.4.6
  [ac1d9e8a] ThreadsX v0.1.8
  [a759f4b9] TimerOutputs v0.5.13
  [0796e94c] Tokenize v0.5.21
  [28d57a85] Transducers v0.4.67
  [a2a6695c] TreeViews v0.3.0
  [d5829a12] TriangularSolve v0.1.8
  [5c2747f8] URIs v1.3.0
  [3a884ed6] UnPack v1.0.2
  [1cfade01] UnicodeFun v0.4.1
  [1986cc42] Unitful v1.9.2
  [3d5dd08c] VectorizationBase v0.21.22
  [81def892] VersionParsing v1.2.1
  [19fa3120] VertexSafeGraphs v0.2.0
  [44d3d7a6] Weave v0.10.10
  [ddb6d928] YAML v0.4.7
  [c2297ded] ZMQ v1.2.1
  [700de1a5] ZygoteRules v0.2.2
  [6e34b625] Bzip2_jll v1.0.6+5
  [83423d85] Cairo_jll v1.16.0+6
  [655fdf9c] DASKR_jll v1.0.0+0
  [5ae413db] EarCut_jll v2.2.3+0
  [2e619515] Expat_jll v2.2.10+0
  [b22a6f82] FFMPEG_jll v4.3.1+4
  [a3f928ae] Fontconfig_jll v2.13.1+14
  [d7e528f0] FreeType2_jll v2.10.1+5
  [559328eb] FriBidi_jll v1.0.10+0
  [0656b61e] GLFW_jll v3.3.5+1
  [d2c73de3] GR_jll v0.58.1+0
  [78b55507] Gettext_jll v0.20.1+7
  [f8c6e375] Git_jll v2.31.0+0
  [7746bdde] Glib_jll v2.59.0+4
  [e33a78d0] Hwloc_jll v2.5.0+0
  [aacddb02] JpegTurbo_jll v2.1.0+0
  [c1c5ebd0] LAME_jll v3.100.1+0
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
  [c771fb93] ODEInterface_jll v0.0.1+0
  [e7412a2a] Ogg_jll v1.3.5+0
  [458c3c95] OpenSSL_jll v1.1.10+0
  [efe28fd5] OpenSpecFun_jll v0.5.5+0
  [91d4177d] Opus_jll v1.3.2+0
  [2f80f16e] PCRE_jll v8.44.0+0
  [30392449] Pixman_jll v0.40.1+0
  [ea2cea3b] Qt5Base_jll v5.15.2+0
  [f50d1b31] Rmath_jll v0.3.0+0
  [fb77eaff] Sundials_jll v5.2.0+1
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
  [4536629a] OpenBLAS_jll
  [05823500] OpenLibm_jll
  [efcefdf7] PCRE2_jll
  [bea87d4a] SuiteSparse_jll
  [83775a58] Zlib_jll
  [8e850ede] nghttp2_jll
  [3f19e933] p7zip_jll
```

