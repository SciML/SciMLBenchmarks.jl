@testset "benchmark entrypoints inherit configured environment" begin
    repository_root = dirname(@__DIR__)
    mktempdir() do folder
        write(joinpath(folder, "Project.toml"), "[deps]\n")
        benchmark = joinpath(folder, "thread_probe.jmd")
        probe = joinpath(folder, "environment.txt")
        probe_code = """
        values = (
            Threads.nthreads(),
            ENV["JULIA_NUM_THREADS"],
            read(`\$(Base.julia_cmd()) -e "print(Threads.nthreads())"`, String),
            read(`\$(Base.julia_cmd()) --threads=1 -e "print(Threads.nthreads())"`, String),
            ENV["SCIMLBENCHMARKS_REVIEW_SENTINEL"],
            ENV["SCIMLBENCHMARKS_OUTER_SENTINEL"],
            ENV["TMPDIR"],
            ENV["OPENBLAS_NUM_THREADS"],
            ENV["PATH"],
        )
        write(ENV["SCIMLBENCHMARKS_REVIEW_PROBE"], join(string.(values), "\\n"))
        """
        write(benchmark, "```julia\n", probe_code, "\n```\n")
        entrypoint = joinpath(repository_root, "benchmark.jl")
        inherited_tmpdir = get(ENV, "TMPDIR", tempdir())
        inherited_blas_threads = get(ENV, "OPENBLAS_NUM_THREADS", "unset")

        inherited_environment = (
            "JULIA_NUM_THREADS" => "3",
            "SCIMLBENCHMARKS_REVIEW_SENTINEL" => "kept",
            "SCIMLBENCHMARKS_REVIEW_PROBE" => probe,
            "TMPDIR" => inherited_tmpdir,
            "OPENBLAS_NUM_THREADS" => inherited_blas_threads,
            "PATH" => ENV["PATH"],
        )

        try
            for (label, target) in (("folder API", nothing), ("page", benchmark), ("folder", folder))
                @testset "$label" begin
                    command = if isnothing(target)
                        code = "using SciMLBenchmarks; SciMLBenchmarks.weave_folder($(repr(folder)), (:github,))"
                        `$(Base.julia_cmd()) --threads=2 --project=$(Base.active_project()) -e $code`
                    else
                        `$(Base.julia_cmd()) --threads=2 --project=$(Base.active_project()) $entrypoint $target`
                    end
                    withenv("SCIMLBENCHMARKS_OUTER_SENTINEL" => "outer") do
                        run(addenv(command, inherited_environment...))
                    end
                    @test read(probe, String) == join(
                        ("2", "2", "2", "1", "kept", "outer", inherited_tmpdir, inherited_blas_threads, ENV["PATH"]),
                        "\n",
                    )
                end
            end
        finally
            for format in ("markdown", "script")
                rm(joinpath(repository_root, format, basename(folder)); recursive = true, force = true)
            end
        end
    end
end
