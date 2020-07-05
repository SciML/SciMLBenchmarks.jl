---
author: "Vaibhav Dixit, Chris Rackauckas"
title: "Lotka-Volterra Bayesian Parameter Estimation Benchmarks"
---


## Parameter Estimation of Lotka-Volterra Equation using DiffEqBayes.jl

````julia
using DiffEqBayes, CmdStan, DynamicHMC
````



````julia
using Distributions, BenchmarkTools
using OrdinaryDiffEq, RecursiveArrayTools, ParameterizedFunctions
using Plots
````



````julia
gr(fmt=:png)
````


````
Plots.GRBackend()
````





#### Initializing the problem

````julia
f = @ode_def LotkaVolterraTest begin
    dx = a*x - b*x*y
    dy = -c*y + d*x*y
end a b c d
````


````
(::Main.##WeaveSandBox#371.LotkaVolterraTest{Main.##WeaveSandBox#371.var"##
#ParameterizedDiffEqFunction#391",Main.##WeaveSandBox#371.var"###Parameteri
zedTGradFunction#392",Main.##WeaveSandBox#371.var"###ParameterizedJacobianF
unction#393",Nothing,Nothing,ModelingToolkit.ODESystem}) (generic function 
with 1 method)
````



````julia
u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1,0]
````


````
5-element Array{Float64,1}:
 1.5
 1.0
 3.0
 1.0
 0.0
````



````julia
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob,Tsit5())
````


````
retcode: Success
Interpolation: specialized 4th order "free" interpolation
t: 34-element Array{Float64,1}:
  0.0
  0.0776084743154256
  0.23264513699277584
  0.4291185174543143
  0.6790821776882875
  0.9444045910389707
  1.2674601253261835
  1.6192913723304114
  1.9869755337814992
  2.264090367186479
  ⋮
  7.584862904164952
  7.978068388305894
  8.483164907244102
  8.719247868929038
  8.949206527971544
  9.200184813643565
  9.438028630962807
  9.711807852444823
 10.0
u: 34-element Array{Array{Float64,1},1}:
 [1.0, 1.0]
 [1.0454942346944578, 0.8576684823217127]
 [1.1758715885138267, 0.639459570317544]
 [1.4196809607170826, 0.4569962601282084]
 [1.8767193485546056, 0.32473343696185236]
 [2.5882499852859384, 0.26336255804531]
 [3.860708771268753, 0.2794458027885767]
 [5.750812903389158, 0.5220073140479389]
 [6.814978737433837, 1.917783300239219]
 [4.3929977807914105, 4.194671536988031]
 ⋮
 [2.614252575185928, 0.26416950055716665]
 [4.2410731685818694, 0.30512345857554246]
 [6.791122470590543, 1.1345265418479897]
 [6.26537352594436, 2.741690196017545]
 [3.78076791065078, 4.431164786168439]
 [1.8164212283793362, 4.0640577258289365]
 [1.1465027088171469, 2.791172606389902]
 [0.9557986534742364, 1.6235632025270912]
 [1.03375813933372, 0.9063703701433561]
````





#### We take the solution data obtained and add noise to it to obtain data for using in the Bayesian Inference of the parameters

````julia
t = collect(range(1,stop=10,length=10))
sig = 0.49
data = convert(Array, VectorOfArray([(sol(t[i]) + sig*randn(2)) for i in 1:length(t)]))
````


````
2×10 Array{Float64,2}:
 2.12986   6.86243  1.54252  1.23979   …  4.05179   3.50023  0.897005
 0.231067  1.4358   1.76601  0.205241     0.702947  4.85362  0.554
````





#### Plots of the actual data and generated data

````julia
scatter(t, data[1,:], lab="#prey (data)")
scatter!(t, data[2,:], lab="#predator (data)")
plot!(sol)
````


![](figures/DiffEqBayesLotkaVolterra_8_1.png)

````julia
priors = [truncated(Normal(1.5,0.5),0.5,2.5),truncated(Normal(1.2,0.5),0,2),truncated(Normal(3.0,0.5),1,4),truncated(Normal(1.0,0.5),0,2)]
````


````
4-element Array{Distributions.Truncated{Distributions.Normal{Float64},Distr
ibutions.Continuous,Float64},1}:
 Truncated(Distributions.Normal{Float64}(μ=1.5, σ=0.5), range=(0.5, 2.5))
 Truncated(Distributions.Normal{Float64}(μ=1.2, σ=0.5), range=(0.0, 2.0))
 Truncated(Distributions.Normal{Float64}(μ=3.0, σ=0.5), range=(1.0, 4.0))
 Truncated(Distributions.Normal{Float64}(μ=1.0, σ=0.5), range=(0.0, 2.0))
````





### Stan.jl backend

The solution converges for tolerance values lower than 1e-3, lower tolerance leads to better accuracy in result but is accompanied by longer warmup and sampling time, truncated normal priors are used for preventing Stan from stepping into negative values.

````julia
@btime bayesian_result_stan = stan_inference(prob,t,data,priors,num_samples=10_000,printsummary=false)
````


````
File /builds/JuliaGPU/DiffEqBenchmarks.jl/tmp/parameter_estimation_model.st
an will be updated.

Error: IOError: chdir : no such file or directory (ENOENT)
````





### Turing.jl backend

````julia
@btime bayesian_result_turing = turing_inference(prob,Tsit5(),t,data,priors,num_samples=10_000)
````


````
14.363 s (107025365 allocations: 7.75 GiB)
Object of type Chains, with data of type 9000×17×1 Array{Float64,3}

Iterations        = 1:9000
Thinning interval = 1
Chains            = 1
Samples per chain = 9000
internals         = acceptance_rate, hamiltonian_energy, hamiltonian_energy
_error, is_accept, log_density, lp, max_hamiltonian_energy_error, n_steps, 
nom_step_size, numerical_error, step_size, tree_depth
parameters        = theta[1], theta[2], theta[3], theta[4], σ[1]

2-element Array{MCMCChains.ChainDataFrame,1}

Summary Statistics
  parameters    mean     std  naive_se    mcse        ess   r_hat
  ──────────  ──────  ──────  ────────  ──────  ─────────  ──────
    theta[1]  1.4155  0.0772    0.0008  0.0018  1816.3494  1.0004
    theta[2]  1.0540  0.1074    0.0011  0.0020  3071.1881  0.9999
    theta[3]  3.3235  0.2742    0.0029  0.0065  1757.6657  1.0003
    theta[4]  1.0912  0.0951    0.0010  0.0022  1816.5819  1.0001
        σ[1]  0.5279  0.1035    0.0011  0.0018  3210.6111  1.0000

Quantiles
  parameters    2.5%   25.0%   50.0%   75.0%   97.5%
  ──────────  ──────  ──────  ──────  ──────  ──────
    theta[1]  1.2854  1.3606  1.4074  1.4626  1.5855
    theta[2]  0.8774  0.9805  1.0416  1.1150  1.3050
    theta[3]  2.7812  3.1349  3.3268  3.5114  3.8519
    theta[4]  0.9047  1.0273  1.0920  1.1557  1.2785
        σ[1]  0.3678  0.4544  0.5117  0.5857  0.7625
````





### DynamicHMC.jl backend

````julia
@btime bayesian_result_dynamichmc = dynamichmc_inference(prob,Tsit5(),t,data,priors,num_samples=10_000)
````


````
17.077 s (89167748 allocations: 9.05 GiB)
(posterior = NamedTuple{(:parameters, :σ),Tuple{Array{Float64,1},Array{Floa
t64,1}}}[(parameters = [1.4202887326680922, 1.0113300793606446, 3.284902175
4156613, 1.077262513738021], σ = [0.5497845366150373, 0.20521750485193588])
, (parameters = [1.4308002868915999, 1.005958878622269, 3.3612526472844246,
 1.0553425083055064], σ = [0.7041400932957184, 0.16134502282402802]), (para
meters = [1.4322983995451768, 0.9783664704139745, 3.176347229907761, 1.0806
268922386297], σ = [0.5471395780512464, 0.3406634577526517]), (parameters =
 [1.4688314539174536, 1.0058925104621452, 3.218995032887886, 1.002413324759
798], σ = [0.7703862127542543, 0.24044124906653097]), (parameters = [1.4498
665610550476, 0.9782528137580607, 3.2365361492321907, 1.029044971458564], σ
 = [0.8676359657062644, 0.2927168421758142]), (parameters = [1.415366910998
7734, 1.05628208894169, 3.3385457045281033, 1.1014931867201445], σ = [0.438
8172599013083, 0.2626132394246403]), (parameters = [1.4823942499278913, 1.0
12573169943293, 3.2176462793681915, 1.0021115443656021], σ = [0.83738968413
04317, 0.29191299482976885]), (parameters = [1.4404231640116802, 1.04764896
4993562, 3.2210794712626196, 1.0528415702518492], σ = [0.40171569997655693,
 0.27805158197077606]), (parameters = [1.335140877058643, 0.969385152432012
4, 3.6046237258027714, 1.1963446007435405], σ = [0.43458012472142127, 0.247
28060529472637]), (parameters = [1.3470365636412858, 0.9130970289254319, 3.
505078706670574, 1.2002480031109874], σ = [0.4211954602201436, 0.1617894830
5826252])  …  (parameters = [1.4428912681953772, 1.0560597529101257, 3.3163
87360145428, 1.0492223492779102], σ = [0.832586663097917, 0.225460537259522
45]), (parameters = [1.4245863882813397, 1.042304589282836, 3.2160740819606
635, 1.075668052119757], σ = [0.5640854222315207, 0.2842984976317355]), (pa
rameters = [1.4459792989083908, 1.004285777581193, 3.130879483936254, 1.086
4667714625802], σ = [0.6397960258835007, 0.3070716525103144]), (parameters 
= [1.424997334481176, 1.0175717884180464, 3.2790001073895985, 1.09082655220
75777], σ = [0.6314977610077022, 0.2879940185548185]), (parameters = [1.411
9576997682337, 1.0433944764315561, 3.307208549411129, 1.077528783739846], σ
 = [0.5894245572723399, 0.3210575092012224]), (parameters = [1.348864501921
7864, 0.905503887598705, 3.507541804094798, 1.210411269510027], σ = [0.5397
838386042715, 0.23982305976332438]), (parameters = [1.3454971361074275, 0.9
07999967453259, 3.5829955133923406, 1.2233878696659166], σ = [0.57542576059
49225, 0.1923065810840884]), (parameters = [1.3224206262868137, 0.928556563
1447762, 3.5131986427355946, 1.2483530430110097], σ = [0.6295166833081953, 
0.20990132290236987]), (parameters = [1.3940924210249286, 0.970337498600050
8, 3.544934030645019, 1.1137471044402818], σ = [0.6681589578612176, 0.18393
28088206675]), (parameters = [1.3765185627615781, 0.9900869970544621, 3.480
9713142348393, 1.1440445808024102], σ = [0.564104036128537, 0.1958265426895
4577])], chain = [[0.3508601838086978, 0.011266374745451967, 1.189336872319
8507, 0.07442311381953264, -0.5982288291193515, -1.583684863027914], [0.358
2339289140279, 0.005941194721091153, 1.2123137161898945, 0.0538653666621162
76, -0.35077794646251304, -1.8242102080375546], [0.35928042640110097, -0.02
1870965010405977, 1.1557318667249208, 0.07754132851334554, -0.6030513389985
44, -1.0768602164437249], [0.38446715535894616, 0.005875217521673809, 1.169
0692092548673, 0.0024104173683135183, -0.2608633148755797, -1.4252795059048
1], [0.37147152533916405, -0.02198714157777506, 1.1745036681770142, 0.02863
1159938283267, -0.14198304663981492, -1.228549546195547], [0.34738879854169
56, 0.05475527929092019, 1.2055352943607933, 0.09666670178506578, -0.823672
2170307794, -1.3370729016872893], [0.3936585170877308, 0.012494783998111162
, 1.1686501231278008, 0.002109318189027722, -0.1774657444000903, -1.2312994
840253741], [0.3649369343153716, 0.04654857272099013, 1.1697165428986565, 0
.051492766230976234, -0.9120106545455605, -1.2799486358456842], [0.28903681
21713721, -0.031093271933243596, 1.282217389645227, 0.1792707417396126, -0.
8333749445538808, -1.397231533350182], [0.2979070415596627, -0.090913129197
75939, 1.2542129757752516, 0.18252820469996142, -0.8646582769851892, -1.821
4592761102242]  …  [0.3666489257367335, 0.05454476787843919, 1.198876045934
7971, 0.04804927003982272, -0.18321796272823052, -1.4896101364202967], [0.3
538815177511596, 0.04143421278792391, 1.168161386486938, 0.0729419123957016
1, -0.5725495810953553, -1.2577305447860838], [0.368786807526483, 0.0042766
19792628258, 1.1413139503662386, 0.0829309371473469, -0.4466058629841342, -
1.1806741628284325], [0.3541699431793388, 0.017419189573625645, 1.187538530
613227, 0.08693571365084453, -0.45966088271220407, -1.244815567968748], [0.
34497718095108776, 0.04247931777534795, 1.1961044950324733, 0.0746702560866
2084, -0.5286085446799571, -1.1361350155086098], [0.299263128831542, -0.099
2637083507068, 1.2549154514109235, 0.19096019401544315, -0.6165865184457804
, -1.4278538785294512], [0.296763562685663, -0.09651093622527282, 1.2761991
861289936, 0.2016239525130845, -0.5526450589764549, -1.6486644039668168], [
0.2794638650185937, -0.07412398120287625, 1.2565269167279471, 0.22182511697
029766, -0.4628029233736292, -1.5611177496391082], [0.33224360929769076, -0
.03011133129956243, 1.265519550684005, 0.10773009997035021, -0.403229172883
60216, -1.693184757516646], [0.3195575309884039, -0.009962463901766474, 1.2
473113681273897, 0.1345698614312015, -0.5725165832741471, -1.63052599786099
76]], tree_statistics = DynamicHMC.TreeStatisticsNUTS[DynamicHMC.TreeStatis
ticsNUTS(-15.916747375810925, 5, turning at positions -22:9, 0.999227010022
3319, 31, DynamicHMC.Directions(0x9053eb69)), DynamicHMC.TreeStatisticsNUTS
(-16.802628003560184, 5, turning at positions 20:35, 0.9619622978644781, 47
, DynamicHMC.Directions(0xab4ac633)), DynamicHMC.TreeStatisticsNUTS(-17.175
794074699787, 5, turning at positions -3:28, 0.6793465443139691, 31, Dynami
cHMC.Directions(0xf72820bc)), DynamicHMC.TreeStatisticsNUTS(-17.31400258754
8114, 4, turning at positions 15:30, 0.9991959832162177, 31, DynamicHMC.Dir
ections(0xc869b27e)), DynamicHMC.TreeStatisticsNUTS(-16.930157636210488, 4,
 turning at positions 1:16, 0.9695392525802309, 31, DynamicHMC.Directions(0
x3509e130)), DynamicHMC.TreeStatisticsNUTS(-18.870480715176335, 5, turning 
at positions -28:3, 0.9220568312299521, 31, DynamicHMC.Directions(0x6aedcb4
3)), DynamicHMC.TreeStatisticsNUTS(-19.768456323638, 5, turning at position
s -30:1, 0.9916901133681929, 31, DynamicHMC.Directions(0x8addff81)), Dynami
cHMC.TreeStatisticsNUTS(-17.683361997572995, 5, turning at positions -26:5,
 0.9867581932952088, 31, DynamicHMC.Directions(0xae6bc365)), DynamicHMC.Tre
eStatisticsNUTS(-15.637982449015434, 6, turning at positions -2:61, 0.99087
92759508681, 63, DynamicHMC.Directions(0x9a72e3bd)), DynamicHMC.TreeStatist
icsNUTS(-18.856131784920304, 5, turning at positions -22:9, 0.9852049789438
235, 31, DynamicHMC.Directions(0xc5dd4949))  …  DynamicHMC.TreeStatisticsNU
TS(-19.131640143610767, 5, turning at positions 36:39, 0.7246814945456048, 
39, DynamicHMC.Directions(0x4b6ae6ff)), DynamicHMC.TreeStatisticsNUTS(-16.5
71416886436115, 5, turning at positions -25:-56, 0.9860182976208569, 63, Dy
namicHMC.Directions(0xe2aa6ac7)), DynamicHMC.TreeStatisticsNUTS(-18.6687508
394307, 4, turning at positions 10:25, 0.9630465022662542, 31, DynamicHMC.D
irections(0x1f373d39)), DynamicHMC.TreeStatisticsNUTS(-17.678366940789758, 
4, turning at positions -12:3, 0.9941799719562847, 15, DynamicHMC.Direction
s(0x8d611213)), DynamicHMC.TreeStatisticsNUTS(-17.194198525338347, 5, turni
ng at positions -17:14, 0.9940337991859914, 31, DynamicHMC.Directions(0x469
9fc4e)), DynamicHMC.TreeStatisticsNUTS(-18.419567029597733, 5, turning at p
ositions -28:-59, 0.9099737275745559, 63, DynamicHMC.Directions(0x88ab2784)
), DynamicHMC.TreeStatisticsNUTS(-16.983507572348476, 5, turning at positio
ns 12:27, 0.9685844542640334, 47, DynamicHMC.Directions(0x23b6c82b)), Dynam
icHMC.TreeStatisticsNUTS(-16.612136145790323, 4, turning at positions -7:8,
 0.9841065368964172, 15, DynamicHMC.Directions(0x8ffce8b8)), DynamicHMC.Tre
eStatisticsNUTS(-16.91516767903357, 5, turning at positions -13:-28, 0.8987
158893316309, 47, DynamicHMC.Directions(0x62b21413)), DynamicHMC.TreeStatis
ticsNUTS(-15.8153287132636, 4, turning at positions -7:-22, 0.9914429318616
386, 31, DynamicHMC.Directions(0xceea4389))], κ = Gaussian kinetic energy (
LinearAlgebra.Diagonal), √diag(M⁻¹): [0.041825173276473, 0.057837169951327,
 0.055618463151112074, 0.08254112097343654, 0.2657686653610789, 0.304300715
24605873], ϵ = 0.09114840376662581)
````






## Conclusion

Lotka-Volterra Equation is a "predator-prey" model, it models population of two species in which one is the predator (wolf) and the other is the prey (rabbit). 
It depicts a cyclic behaviour, which is also seen in its Uncertainity Quantification Plots. This behaviour makes it easy to estimate even at very high tolerance values (1e-3).



````julia
using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])
````



## Appendix

These benchmarks are a part of the DiffEqBenchmarks.jl repository, found at: [https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl](https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl)

To locally run this tutorial, do the following commands:

```
using DiffEqBenchmarks
DiffEqBenchmarks.weave_file("ParameterEstimation","DiffEqBayesLotkaVolterra.jmd")
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
Status: `/builds/JuliaGPU/DiffEqBenchmarks.jl/benchmarks/ParameterEstimation/Project.toml`
[6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf] BenchmarkTools 0.5.0
[a134a8b2-14d6-55f6-9291-3336d3ab0209] BlackBoxOptim 0.5.0
[593b3428-ca2f-500c-ae53-031589ec8ddd] CmdStan 6.0.6
[ebbdde9d-f333-5424-9be2-dbf1e9acfb5e] DiffEqBayes 2.16.0
[1130ab10-4a5a-5621-a13d-e4788d82bd4c] DiffEqParamEstim 1.15.0
[ef61062a-5684-51dc-bb67-a0fcdec5c97d] DiffEqUncertainty 1.4.1
[31c24e10-a181-5473-b8eb-7969acd0382f] Distributions 0.23.4
[bbc10e6e-7c05-544b-b16e-64fede858acb] DynamicHMC 2.1.5
[76087f3c-5699-56af-9a33-bf431cd00edd] NLopt 0.6.0
[1dea7af3-3e70-54e6-95c3-0bf5283fa5ed] OrdinaryDiffEq 5.41.0
[65888b18-ceab-5e60-b2b9-181511a3b968] ParameterizedFunctions 5.3.0
[91a5bcdd-55d7-5caf-9e0b-520d859cae80] Plots 1.5.3
[731186ca-8d62-57ce-b412-fbd966d074cd] RecursiveArrayTools 2.5.0
```

