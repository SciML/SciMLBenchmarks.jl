module SciMLBenchmarks

using InteractiveUtils: InteractiveUtils
using Markdown: Markdown, @md_str
using Pkg: Pkg
using Weave: Weave, tangle, weave
using PrecompileTools: @compile_workload, @setup_workload

include("weave.jl")

macro subprocess(ex, wait = true)
    return quote
        local project = Base.active_project()
        local ex_str = $(string(ex))
        run(`$(Base.julia_cmd()) --project=$(project) -e "$(ex_str)"`; wait = $(wait))
    end
end

function _benchmark_priority(path)
    return open(path) do io
        eof(io) && return 0
        strip(readline(io)) == "---" || return 0
        for line in eachline(io)
            strip(line) == "---" && return 0
            key_value = split(line, ':'; limit = 2)
            if length(key_value) == 2 && strip(first(key_value)) == "priority"
                return parse(Int, strip(last(key_value)))
            end
        end
        return 0
    end
end

function _uses_benchmark_project(folder, build_list)
    return isfile(joinpath(folder, "Project.toml")) && build_list != (:notebook,)
end

function weave_file(folder, file, build_list = (:script, :github))
    if _uses_benchmark_project(folder, build_list)
        return _weave_subprocess(folder, file, build_list, true)
    end
    @info("Weaving $(joinpath(folder, file))")
    return _render_file(folder, file, build_list)
end

function _weave_subprocess(folder, file, build_list, use_benchmark_project = _uses_benchmark_project(folder, build_list))
    caller_project = something(Base.active_project(), repo_directory)
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
    gc_flags = VERSION >= v"1.10" ? `--gcthreads=$(Threads.ngcthreads())` : ``
    outputs = string.(collect(build_list))
    run(
        `$(Base.julia_cmd()) --startup-file=no --threads=$threads $gc_flags --project=$project $worker $folder $file $caller_project $use_benchmark_project $(length(outputs)) $outputs $LOAD_PATH`
    )
    return :notebook ∈ build_list ? _notebook_path(folder, file) : nothing
end

function weave_all(build_list = (:script, :github))
    for folder in readdir(joinpath(repo_directory, "benchmarks"))
        folder == "test.jmd" && continue
        weave_folder(joinpath(repo_directory, "benchmarks", folder), build_list)
    end
    return
end

function weave_folder(folder, build_list = (:script, :github))
    weave_files = String[]
    priorities = Int[]
    for file in readdir(folder)
        # Skip non-`.jmd` files
        endswith(file, ".jmd") || continue
        push!(weave_files, file)
        push!(priorities, _benchmark_priority(joinpath(folder, file)))
    end

    weave_files = weave_files[sortperm(priorities; rev = true)]

    failures = Tuple{String, String}[]
    for file in weave_files
        try
            _weave_subprocess(folder, file, build_list)
        catch e
            @show folder, file
            @error(e)
            push!(failures, (folder, file))
        end
    end
    isempty(failures) || error("weave_folder: $(length(failures)) file(s) failed: $failures")
    return
end

function bench_footer(folder = nothing, file = nothing)
    display(
        md"""
        ## Appendix

        These benchmarks are a part of the SciMLBenchmarks.jl repository, found at: <https://github.com/SciML/SciMLBenchmarks.jl>.
        For more information on high-performance scientific machine learning, check out the SciML Open Source Software Organization <https://sciml.ai>.

        """
    )
    if folder !== nothing && file !== nothing
        display(
            "text/markdown",
            """
            To locally run this benchmark, do the following commands:
            ```
            using SciMLBenchmarks
            SciMLBenchmarks.weave_file("$folder","$file")
            ```
            """,
        )
    end
    display(md"Computer Information:")
    vinfo = sprint(InteractiveUtils.versioninfo)
    display(
        "text/markdown",
        """
        ```
        $(vinfo)
        ```
        """,
    )

    display(
        md"""
        Package Information:
        """
    )

    proj = sprint(io -> Pkg.status(io = io))
    mani = sprint(io -> Pkg.status(io = io, mode = Pkg.PKGMODE_MANIFEST))

    md = """
    ```
    $(chomp(proj))
    ```

    And the full manifest:

    ```
    $(chomp(mani))
    ```
    """
    return display("text/markdown", md)
end

"""
    open_notebooks()

Open the SciMLBenchmarks notebooks in IJulia/Jupyter.

Requires IJulia to be loaded. If IJulia is not loaded, this function will error
with instructions to load it first.
"""
function open_notebooks(args...)
    error("IJulia is required for open_notebooks(). Please run `using IJulia` first.")
end

# Precompilation workload to improve TTFX for common operations
@setup_workload begin
    # Setup code - this runs at precompile time but is not precompiled
    @compile_workload begin
        # Precompile Markdown operations commonly used in bench_footer
        md_text = md"""
        ## Test Header
        This is a test markdown string.
        """

        markdown_text = """
        Test markdown content
        ```
        code block
        ```
        """
        repr("text/markdown", markdown_text)

        # Precompile sprint with versioninfo (used in bench_footer)
        vinfo = sprint(InteractiveUtils.versioninfo)

        # Precompile Pkg.status operations (used in bench_footer)
        proj_status = sprint(io -> Pkg.status(io = io))
        manifest_status = sprint(io -> Pkg.status(io = io, mode = Pkg.PKGMODE_MANIFEST))

        # Precompile string operations commonly used
        chomp(proj_status)
        chomp(manifest_status)
    end
end

end # module SciMLBenchmarks
