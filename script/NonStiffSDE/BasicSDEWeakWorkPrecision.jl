
using StochasticDiffEq, DiffEqDevTools, ParameterizedFunctions, DiffEqProblemLibrary
using Plots; gr()
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems; importsdeproblems()
import DiffEqProblemLibrary.SDEProblemLibrary: prob_sde_additive,
            prob_sde_linear, prob_sde_wave
const N = 1000


prob = prob_sde_additive

reltols = 1.0 ./ 10.0 .^ (1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [
          Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
          Dict(:alg=>SRIW1(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
          Dict(:alg=>SRA1(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
          Dict(:alg=>SRA1())
          Dict(:alg=>SRIW1())
          ]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns_error=N,
                      save_everystep = false,
                      parallel_type = :none,
                      error_estimate=:weak_final)#
plot(wp)


sample_size = Int[10;1e2;1e3;1e4]
se = get_sample_errors(prob,setups[6],numruns=sample_size,
                                      sample_error_runs = 100_000,solution_runs=100)


times = [wp[i].times for i in 1:length(wp)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!([se[end];se[end]],times,color=:red,linestyle=:dash,label="Sample Error: 10000",lw=3)


prob = prob_sde_additive

reltols = 1.0 ./ 10.0 .^ (1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [
          Dict(:alg=>SRA1())
          Dict(:alg=>SRA2())
          Dict(:alg=>SRA3())
          Dict(:alg=>SOSRA())
          Dict(:alg=>SOSRA2())
          ]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns_error=N,
                      save_everystep = false,
                      maxiters = 1e7,
                      parallel_type = :none,
                      error_estimate=:weak_final)
plot(wp)


sample_size = Int[10;1e2;1e3;1e4]
se = get_sample_errors(prob,setups[4],numruns=sample_size,
                                      sample_error_runs = 100_000,solution_runs=100)


times = [wp[i].times for i in 1:length(wp)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!([se[end];se[end]],times,color=:red,linestyle=:dash,label="Sample Error: 10000",lw=3)


prob = prob_sde_linear

reltols = 1.0 ./ 10.0 .^ (1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [Dict(:alg=>SRIW1())
          Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
          Dict(:alg=>SRIW1(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
          ]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns_error=N,
                      save_everystep = false,
                      maxiters = 1e7,
                      parallel_type = :none,
                      error_estimate=:weak_final)
plot(wp)


sample_size = Int[10;1e2;1e3;1e4]
se = get_sample_errors(prob,setups[1],numruns=sample_size,
                                      sample_error_runs = 100_000,solution_runs=100)


times = [wp[i].times for i in 1:length(wp)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!([se[end];se[end]],times,color=:red,linestyle=:dash,label="Sample Error: 10000",lw=3)


prob = prob_sde_linear

reltols = 1.0 ./ 10.0 .^ (1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 2))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 2),:adaptive=>false)
          Dict(:alg=>SRI())
          Dict(:alg=>SRIW1())
          Dict(:alg=>SRIW2())
          Dict(:alg=>SOSRI())
          Dict(:alg=>SOSRI2())
          ]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns_error=N,
                      save_everystep = false,
                      maxiters = 1e7,
                      parallel_type = :none,
                      error_estimate=:weak_final)
plot(wp)


sample_size = Int[10;1e2;1e3;1e4]
se = get_sample_errors(prob,setups[6],numruns=sample_size,
                                      sample_error_runs = 100_000,solution_runs=100)


times = [wp[i].times for i in 1:length(wp)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!([se[end];se[end]],times,color=:red,linestyle=:dash,label="Sample Error: 10000",lw=3)


prob = prob_sde_wave

reltols = 1.0 ./ 10.0 .^ (1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [
          Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
          Dict(:alg=>SRIW1(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
          Dict(:alg=>SRIW1())
          ]

wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns_error=N,
                      save_everystep = false,
                      maxiters = 1e7,
                      parallel_type = :none,
                      error_estimate=:weak_final)
plot(wp)


sample_size = Int[10;1e2;1e3;1e4]
se = get_sample_errors(prob,setups[4],numruns=sample_size,
                                      sample_error_runs = 100_000,solution_runs=100)


times = [wp[i].times for i in 1:length(wp)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!([se[end];se[end]],times,color=:red,linestyle=:dash,label="Sample Error: 10000",lw=3)


prob = prob_sde_wave

reltols = 1.0 ./ 10.0 .^ (1:5)
abstols = reltols#[0.0 for i in eachindex(reltols)]

setups = [Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 2))
          Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 2),:adaptive=>false)
          Dict(:alg=>SRI())
          Dict(:alg=>SRIW1())
          Dict(:alg=>SRIW2())
          Dict(:alg=>SOSRI())
          Dict(:alg=>SOSRI2())
          ]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns_error=N,
                      save_everystep = false,
                      maxiters = 1e7,
                      parallel_type = :none,
                      error_estimate=:weak_final)
plot(wp)


sample_size = Int[10;1e2;1e3;1e4]
se = get_sample_errors(prob,setups[6],numruns=sample_size,
                                      sample_error_runs = 100_000,solution_runs=100)


times = [wp[i].times for i in 1:length(wp)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!([se[end];se[end]],times,color=:red,linestyle=:dash,label="Sample Error: 10000",lw=3)


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

