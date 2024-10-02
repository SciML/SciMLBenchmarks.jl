```julia
using JumpProcesses, Plots, StableRNGs, Random, BenchmarkTools

u0 = [1.]
p = (λ = 1.0,)
tspan = (0.0, 100.0)
rng = StableRNG(53124)

function generateJumpProblem(method, N::Int64) 
    depgraph = Vector{Vector{Int64}}()
    jumpvec = Vector{ConstantRateJump}()
    for i in 1:N
        rate(u, p, t) = p.λ * u[1]
        # We make each reaction affect the rates of 10 others, so there are 10 updates per jump. 
        push!(depgraph, rand(1:N, 10))
        # Simplest way to get random updates on the times, since we don't resample times
        # when the propensities are updated. 
        affect!(integrator) = begin 
            integrator.u[1] = randexp(rng)
        end
        push!(jumpvec, ConstantRateJump(rate, affect!))
    end

    jset = JumpSet((), jumpvec, nothing, nothing)
    dprob = DiscreteProblem(u0, tspan, p)
    jump_prob = JumpProblem(dprob, method, jset; dep_graph = depgraph, save_positions = (false, false), rng)
    jump_prob
end

@kwdef mutable struct EventCallback
    n::Int = 0
end

function (ecb::EventCallback)(u, t, integ)
    ecb.n += 1
    ecb.n == 10^7
end

function (ecb::EventCallback)(integ) 
    terminate!(integ)
    nothing
end

algs = [NRM(), CCNRM(), DirectCR()]

function benchmark_and_save!(bench_dict, end_times, Nv, algs)
    @assert length(end_times) == length(Nv)

    # callback for terminating simulations
    ecb = EventCallback()
    cb = DiscreteCallback(ecb, ecb)

    for (end_time, N) in zip(end_times, Nv)
        names = ["$s"[1:end-2] for s in algs]
        # benchmarking and saving
        benchmarks = Vector{BenchmarkTools.Trial}(undef, length(algs))

        # callback for terminating simulations
        for (i, alg) in enumerate(algs)
            name = names[i]
            println("benchmarking $name")
            jp = generateJumpProblem(alg, N)
            b = @benchmarkable solve($jp, SSAStepper(); saveat = tspan, callback) setup = (callback = deepcopy($cb)) samples = 1 seconds = 360
            bench_dict[name, N] = run(b)
        end
    end
end

function fetch_and_plot(bench_dict)
    names = unique([key[1] for key in keys(bench_dict)])
    Nv = sort(unique([key[2] for key in keys(bench_dict)]))

    plt1 = plot()
    plt2 = plot()

    medtimes = [Float64[] for i in 1:length(names)]
    for (i,name) in enumerate(names)
        for N in Nv
            try
                push!(medtimes[i], median(bench_dict[name, N]).time/1e9)
            catch
                break
            end
        end
        len = length(medtimes[i])
        plot!(plt1, Nv[1:len], medtimes[i], marker = :hex, label = name, lw = 2)
    end

    plot!(plt1, xlabel = "Number of reaction channels", ylabel = "Median time in seconds",
                legend = :bottomright)
end

N = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000]
bench_dict = Dict{Tuple{String, Int}, BenchmarkTools.Trial}()
end_times = 2000. * ones(length(N))
benchmark_and_save!(bench_dict, end_times, N, algs)
plt = fetch_and_plot(bench_dict)
```
