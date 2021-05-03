using SciMLBenchmarks
target = ARGS[1]
if isdir(target)
    if !isfile(joinpath(target, "Project.toml"))
        error("Cannot benchmark folder $(target) without Project.toml!")
    end
    println("Benchmarking the $(target) folder")
    SciMLBenchmarks.weave_folder(target)
elseif isfile(target)
    folder = dirname(target)
    file = basename(target)
    println("Benchmarking $(folder)/$(file)")
    SciMLBenchmarks.weave_file(folder, file)
else
    error("Unable to find benchmarking target $(target)!")
end
