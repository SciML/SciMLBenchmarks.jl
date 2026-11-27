function _weave_subprocess(folder, file, build_list, use_benchmark_project)
    caller_project = something(Base.active_project(), dirname(@__DIR__))
    project = use_benchmark_project ? abspath(folder) : caller_project
    worker = joinpath(@__DIR__, "weave_worker.jl")
    threads = string(Threads.nthreads())
    if VERSION >= v"1.9"
        interactive = Threads.nthreads(:interactive)
        # Julia 1.12 accepts an explicit zero for the interactive pool.
        if interactive > 0 || VERSION >= v"1.12"
            threads = "$threads,$interactive"
        end
    end
    gc_flags = if VERSION >= v"1.10"
        # Threads.ngcthreads() reports the total. Julia exposes the sweep count through its startup options.
        sweep = Int(Base.JLOptions().nsweepthreads)
        mark = Threads.ngcthreads() - sweep
        `--gcthreads=$mark,$sweep`
    else
        ``
    end
    outputs = string.(collect(build_list))
    command = `$(Base.julia_cmd()) --startup-file=no --threads=$threads $gc_flags --project=$project $worker $folder $file $caller_project $use_benchmark_project $(length(outputs)) $outputs $LOAD_PATH`
    inherited_threads = get(ENV, "JULIA_NUM_THREADS", string(Threads.nthreads()))
    return run(addenv(command, "JULIA_NUM_THREADS" => inherited_threads))
end
