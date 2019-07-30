
using StochasticDiffEq, DiffEqDevTools, ParameterizedFunctions
using Plots; gr()
const N = 100

f = @ode_def LotkaVolterraTest begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d

p = [1.5,1.0,3.0,1.0]

function g(du,u,p,t)
  du .= 0.1u
end
u0 = [1.0;1.0]
tspan = (0.0,10.0)
prob = SDEProblem(f,g,u0,tspan,p);


sol = solve(prob,SRIW1(),abstol=1e-4,reltol=1e-4)
plot(sol)


reltols = 1.0 ./ 4.0 .^ (2:4)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [Dict(:alg=>SRIW1())
          Dict(:alg=>EM(),:dts=>1.0./12.0.^((1:length(reltols)) .+ 1.5))
          Dict(:alg=>RKMil(),:dts=>1.0./12.0.^((1:length(reltols)) .+ 1.5),:adaptive=>false)
          Dict(:alg=>SRIW1(),:dts=>1.0./4.0.^((1:length(reltols)) .+ 5),:adaptive=>false)
          Dict(:alg=>SRIW2())
          Dict(:alg=>SOSRI())
          Dict(:alg=>SOSRI2())
          ]
test_dt = 1/10^2
appxsol_setup = Dict(:alg=>SRIW1(),:abstol=>1e-4,:reltol=>1e-4)
wp = WorkPrecisionSet(prob,abstols,reltols,setups,test_dt;
                                     maxiters = 1e7,
                                     verbose=false,save_everystep=false,
                                     parallel_type = :threads,
                                     appxsol_setup = appxsol_setup,
                                     numruns_error=N,error_estimate=:final)
plot(wp)


reltols = 1.0 ./ 4.0 .^ (2:4)
abstols = reltols#[0.0 for i in eachindex(reltols)]
setups = [Dict(:alg=>SRIW1())
          Dict(:alg=>EM(),:dts=>1.0./12.0.^((1:length(reltols)) .+ 1.5))
          Dict(:alg=>RKMil(),:dts=>1.0./12.0.^((1:length(reltols)) .+ 1.5),:adaptive=>false)
          Dict(:alg=>SRIW1(),:dts=>1.0./4.0.^((1:length(reltols)) .+ 5),:adaptive=>false)
          Dict(:alg=>SRIW2())
          Dict(:alg=>SOSRI())
          Dict(:alg=>SOSRI2())
          ]
test_dt = 1e-2
appxsol_setup = Dict(:alg=>SRIW1(),:abstol=>1e-4,:reltol=>1e-4)
wp = WorkPrecisionSet(prob,abstols,reltols,setups,test_dt;
                                     maxiters = 1e7,
                                     verbose=false,save_everystep=false,
                                     parallel_type = :none,
                                     appxsol_setup = appxsol_setup,
                                     numruns_error=N,error_estimate=:weak_final)
plot(wp;legend=:topleft)


sample_size = Int[10;1e2;1e3]
se = get_sample_errors(prob,setups[6],test_dt,numruns=sample_size,
                                      appxsol_setup = appxsol_setup,
                                      sample_error_runs = 100_000,solution_runs=20)


plot(wp;legend=:topleft)
times = [wp[i].times for i in 1:length(wp)]
times = [minimum(minimum(t) for t in times),maximum(maximum(t) for t in times)]
plot!([se[end];se[end]],times,color=:orange,linestyle=:dash,label="Sample Error: 1000",lw=3)


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

