@testset "folder subprocess inherits configured environment" begin
    repository_root = dirname(@__DIR__)
    mktempdir() do folder
        benchmark = joinpath(folder, "thread_probe.jmd")
        probe = joinpath(folder, "environment.txt")
        probe_code = """
        values = (
            Threads.nthreads(),
            ENV["SCIMLBENCHMARKS_REVIEW_SENTINEL"],
            ENV["SCIMLBENCHMARKS_OUTER_SENTINEL"],
            ENV["TMPDIR"],
            ENV["OPENBLAS_NUM_THREADS"],
            ENV["PATH"],
        )
        write(ENV["SCIMLBENCHMARKS_REVIEW_PROBE"], join(string.(values), "\\n"))
        """
        write(benchmark, "```julia\n", probe_code, "\n```\n")
        output_dir = joinpath(repository_root, "markdown", basename(folder))
        code = "using SciMLBenchmarks; SciMLBenchmarks.weave_folder($(repr(folder)), (:github,))"
        command = `$(Base.julia_cmd()) --threads=2 --project=$(Base.active_project()) -e $code`
        inherited_tmpdir = get(ENV, "TMPDIR", tempdir())
        inherited_blas_threads = get(ENV, "OPENBLAS_NUM_THREADS", "unset")

        inherited_environment = (
            "JULIA_NUM_THREADS" => "auto",
            "SCIMLBENCHMARKS_REVIEW_SENTINEL" => "kept",
            "SCIMLBENCHMARKS_REVIEW_PROBE" => probe,
            "TMPDIR" => inherited_tmpdir,
            "OPENBLAS_NUM_THREADS" => inherited_blas_threads,
            "PATH" => ENV["PATH"],
        )

        try
            withenv("SCIMLBENCHMARKS_OUTER_SENTINEL" => "outer") do
                run(addenv(command, inherited_environment...))
            end
            @test read(probe, String) == join(
                ("2", "kept", "outer", inherited_tmpdir, inherited_blas_threads, ENV["PATH"]),
                "\n",
            )
        finally
            rm(output_dir; recursive = true, force = true)
        end
    end
end
