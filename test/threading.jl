@testset "folder subprocess inherits configured thread count" begin
    repository_root = dirname(@__DIR__)
    mktempdir() do folder
        benchmark = joinpath(folder, "thread_probe.jmd")
        write(benchmark, "```julia\nThreads.nthreads()\n```\n")
        output_dir = joinpath(repository_root, "markdown", basename(folder))
        output = joinpath(output_dir, "thread_probe.md")
        code = "using SciMLBenchmarks; SciMLBenchmarks.weave_folder($(repr(folder)), (:github,))"
        command = `$(Base.julia_cmd()) --threads=2 --project=$(Base.active_project()) -e $code`

        try
            run(setenv(command, "JULIA_NUM_THREADS" => "auto"))
            @test occursin("```\n2\n```", read(output, String))
        finally
            rm(output_dir; recursive = true, force = true)
        end
    end
end
