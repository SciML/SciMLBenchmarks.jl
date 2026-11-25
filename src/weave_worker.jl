using Pkg: Pkg

folder, file, caller_project, use_benchmark_project = ARGS[1:4]
build_list = Symbol.(ARGS[5:end])
pushfirst!(LOAD_PATH, "@")

@info("Weaving $(joinpath(folder, file))")
if use_benchmark_project == "true"
    @info("Instantiating", folder)
    withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do
        Pkg.instantiate()
    end
    Pkg.precompile()
end

# The harness project exposes Weave by name. The benchmark manifest takes
# precedence when resolving Weave and its dependencies.
insert!(LOAD_PATH, 2, caller_project)
insert!(LOAD_PATH, 3, joinpath(@__DIR__, ".."))
using Weave: Weave, tangle, weave
deleteat!(LOAD_PATH, 2:3)

include("weave.jl")
_render_file(folder, file, build_list)
