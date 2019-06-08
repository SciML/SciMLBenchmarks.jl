using DiffEqBenchmarks
if length(ARGS) == 1
    println("Benchmarking the $(ARGS[1]) folder")
    DiffEqBenchmarks.weave_folder(ARGS[1])
elseif length(ARGS) == 2
    println("Benchmarking $(ARGS[1])/$(ARGS[2])")
    DiffEqBenchmarks.weave_file(ARGS[1],ARGS[2])
else
    error("Only 1 or 2 arguments allowed!")
end
