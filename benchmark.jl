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
        caller_project = something(Base.active_project(), @__DIR__)
        include("src/weave_bootstrap.jl")
        Pkg.activate(abspath(folder))
        _bootstrap_weave(folder, file, (:script, :github), caller_project, true)
    else
        using SciMLBenchmarks
        SciMLBenchmarks.weave_file(folder, file)
    end
else
    error("Unable to find benchmarking target $(target)!")
end
