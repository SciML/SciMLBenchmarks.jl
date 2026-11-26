push!(LOAD_PATH, "@stdlib")
using Pkg: Pkg
pop!(LOAD_PATH)

include("weave.jl")

function _bootstrap_weave(folder, file, build_list, caller_project, use_benchmark_project)
    @info("Weaving $(joinpath(folder, file))")
    if use_benchmark_project
        if isempty(LOAD_PATH) || first(LOAD_PATH) != "@"
            pushfirst!(LOAD_PATH, "@")
        end
        @info("Instantiating", folder)
        withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do
            Pkg.instantiate()
        end
        Pkg.precompile()
    end

    # The benchmark manifest takes precedence over the projects that expose Weave by name.
    bootstrap_index = use_benchmark_project ? 2 : 1
    insert!(LOAD_PATH, bootstrap_index, caller_project)
    insert!(LOAD_PATH, bootstrap_index + 1, repo_directory)
    try
        @eval using Weave: Weave, tangle, weave
    finally
        deleteat!(LOAD_PATH, bootstrap_index:(bootstrap_index + 1))
    end
    return Base.invokelatest(_render_file, folder, file, build_list)
end
