This is a test of the builder system.

````julia
using Plots
plot(rand(10,10))
````


{{< figure src="../figures/test_1_1.png"  >}}


```math
\begin{equation}
 y'(t) = \frac{0.2y(t-14)}{1 + y(t-14)^{10}} - 0.1y(t)
\end{equation}
```

$\alpha$

``u_0``

## Appendix

````julia
using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])
````



<div class="markdown"><h2>Appendix</h2>
<p>These benchmarks are a part of the DiffEqBenchmarks.jl repository, found at: <a href="https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl">https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl</a></p>
</div>
<div class="markdown"><p>To locally run this tutorial, do the following commands:</p>
<pre><code>using DiffEqBenchmarks
DiffEqBenchmarks.weave_file&#40;&quot;.&quot;,&quot;test.jmd&quot;&#41;</code></pre>
</div>
<div class="markdown"><p>Computer Information:</p>
</div>
<div class="markdown"><pre><code>Julia Version 1.3.0
Commit 46ce4d7933 &#40;2019-11-26 06:09 UTC&#41;
Platform Info:
  OS: Windows &#40;x86_64-w64-mingw32&#41;
  CPU: Intel&#40;R&#41; Core&#40;TM&#41; i7-8550U CPU @ 1.80GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 &#40;ORCJIT, skylake&#41;
Environment:
  JULIA_EDITOR &#61; &quot;C:\Users\accou\AppData\Local\atom\app-1.41.0\atom.exe&quot;  -a
  JULIA_NUM_THREADS &#61; 4
</code></pre>
</div>
<div class="markdown"><p>Package Information:</p>
</div>
<div class="markdown"><pre><code>Status: &#96;C:\Users\accou\.julia\environments\v1.3\Project.toml&#96;
&#91;0bf59076-c3b1-5ca4-86bd-e02cd72cde3d&#93; AdvancedHMC 0.2.16
&#91;c52e3926-4ff0-5f6e-af25-54175e0327b1&#93; Atom 0.11.3
&#91;a134a8b2-14d6-55f6-9291-3336d3ab0209&#93; BlackBoxOptim 0.5.0
&#91;49dc2e85-a5d0-5ad3-a950-438e2897f1b9&#93; Calculus 0.5.1
&#91;3a865a2d-5b23-5a0f-bc46-62713ec82fae&#93; CuArrays 1.5.0
&#91;2445eb08-9709-466a-b3fc-47e12bd697a2&#93; DataDrivenDiffEq 0.1.0
&#91;bcd4f6db-9728-5f36-b5f7-82caef46ccdb&#93; DelayDiffEq 5.18.0
&#91;2b5f629d-d688-5b77-993f-72d75c75574e&#93; DiffEqBase 6.9.0
&#91;31c91b34-3c75-11e9-0341-95557aab0344&#93; DiffEqBenchmarks 0.1.0
&#91;aae7a2af-3d4f-5e19-a356-7da93b79d9d0&#93; DiffEqFlux 0.10.0
&#91;071ae1c0-96b5-11e9-1965-c90190d839ea&#93; DiffEqGPU 1.1.0
&#91;41bf760c-e81c-5289-8e54-58b1f1f8abe2&#93; DiffEqSensitivity 5.0.0
&#91;ef61062a-5684-51dc-bb67-a0fcdec5c97d&#93; DiffEqUncertainty 1.4.0
&#91;587475ba-b771-5e3f-ad9e-33799f191a9c&#93; Flux 0.10.0
&#91;f6369f11-7733-5829-9624-2563aa707210&#93; ForwardDiff 0.10.7
&#91;d1acc4aa-44c8-5952-acd4-ba5d80a2a253&#93; IntervalArithmetic 0.16.2
&#91;e5e0dc1b-0480-54bc-9374-aad01c23163d&#93; Juno 0.7.2
&#91;7f56f5a3-f504-529b-bc02-0b1fe5e64312&#93; LSODA 0.6.1
&#91;10e44e05-a98a-55b3-a45b-ba969058deb6&#93; MATLAB 0.7.0
&#91;1914dd2f-81c6-5fcd-8719-6d5c9610ff09&#93; MacroTools 0.5.3
&#91;eff96d63-e80a-5855-80a2-b1b0885c5ab7&#93; Measurements 2.1.1
&#91;961ee093-0014-501f-94e3-6117800e7a78&#93; ModelingToolkit 1.0.3
&#91;429524aa-4258-5aef-a3af-852621145aeb&#93; Optim 0.19.7
&#91;1dea7af3-3e70-54e6-95c3-0bf5283fa5ed&#93; OrdinaryDiffEq 5.26.4
&#91;d96e819e-fc66-5662-9728-84c9c7592b0a&#93; Parameters 0.12.0
&#91;14b8a8f1-9102-5b29-a752-f990bacb7fe1&#93; PkgTemplates 0.6.3
&#91;91a5bcdd-55d7-5caf-9e0b-520d859cae80&#93; Plots 0.28.3
&#91;1fd47b50-473d-5c70-9696-f719f8f3bcdc&#93; QuadGK 2.3.1
&#91;8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b&#93; QuasiMonteCarlo 0.1.0
&#91;731186ca-8d62-57ce-b412-fbd966d074cd&#93; RecursiveArrayTools 1.2.0
&#91;37e2e3b7-166d-5795-8a7a-e32c996b4267&#93; ReverseDiff 0.3.1
&#91;789caeaf-c7a9-5a7d-9973-96adeb23e2a0&#93; StochasticDiffEq 6.15.0
&#91;c3572dad-4567-51f8-b174-8c6c989267f4&#93; Sundials 3.8.1
&#91;6fc51010-71bc-11e9-0e15-a3fcc6593c49&#93; Surrogates 0.5.0
&#91;9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c&#93; Tracker 0.2.5
&#91;44d3d7a6-8a23-5bf8-98c5-b353f8df5ec9&#93; Weave 0.9.1
&#91;e88e6eb3-aa80-5325-afca-941959d7151f&#93; Zygote 0.4.1</code></pre>
</div>
