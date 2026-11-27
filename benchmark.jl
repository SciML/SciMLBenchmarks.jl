ENV["JULIA_NUM_THREADS"] = get(ENV, "JULIA_NUM_THREADS", string(Threads.nthreads()))
target = ARGS[1]
if isdir(target)
    if !isfile(joinpath(target, "Project.toml"))
        error("Cannot benchmark folder $(target) without Project.toml!")
    end
    println("Benchmarking the $(target) folder")
    using SciMLBenchmarks
    SciMLBenchmarks.weave_folder(target)
elseif isfile(target)
    folder = dirname(target)
    file = basename(target)
    println("Benchmarking $(folder)/$(file)")
    if isfile(joinpath(folder, "Project.toml"))
        include("src/weave_subprocess.jl")
        _weave_subprocess(folder, file, (:script, :github), true)
    else
        using SciMLBenchmarks
        SciMLBenchmarks.weave_file(folder, file)
    end
else
    error("Unable to find benchmarking target $(target)!")
end
