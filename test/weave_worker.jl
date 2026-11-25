using Pkg: Pkg
using SciMLBenchmarks: SciMLBenchmarks
using Test: @test, @testset

@testset "weave worker environment" begin
    mktempdir() do directory
        directory = realpath(directory)
        caller = joinpath(directory, "caller")
        benchmark = joinpath(directory, basename(directory))
        package_directory = joinpath(directory, "extra packages")
        support = mkpath(joinpath(package_directory, "WeaveLoadPathFixture", "src"))
        write(
            joinpath(support, "WeaveLoadPathFixture.jl"), """
            __precompile__(false)
            module WeaveLoadPathFixture
            const MESSAGE = "Loaded from the caller's package directory."
            end
            """
        )

        uuid = "0479c93c-fce4-45d0-a4d9-9aaa8049f391"
        for (project, version) in ((caller, "1.0.0"), (benchmark, "2.0.0"))
            package = joinpath(project, "fixture")
            mkpath(joinpath(package, "src"))
            write(joinpath(project, "Project.toml"), "[deps]\nWeaveProjectFixture = \"$uuid\"\n")
            write(
                joinpath(project, "Manifest.toml"), """
                julia_version = "$(VERSION)"
                manifest_format = "2.0"

                [[deps.WeaveProjectFixture]]
                path = "fixture"
                uuid = "$uuid"
                version = "$version"
                """
            )
            write(
                joinpath(package, "Project.toml"), """
                name = "WeaveProjectFixture"
                uuid = "$uuid"
                version = "$version"
                """
            )
            write(
                joinpath(package, "src", "WeaveProjectFixture.jl"),
                "__precompile__(false)\nmodule WeaveProjectFixture\nconst VERSION_LABEL = \"$version\"\nend\n"
            )
        end

        write(
            joinpath(benchmark, "environment.jmd"), """
            # Worker environment

            ```julia
            import WeaveProjectFixture
            open(joinpath(WEAVE_ARGS[:folder], "environment.txt"), "w") do io
                println(io, Base.active_project())
                println(io, WeaveProjectFixture.VERSION_LABEL)
                println(io, repr(LOAD_PATH))
                println(io, Threads.nthreads())
                println(io, VERSION >= v"1.9" ? Threads.nthreads(:interactive) : 0)
                println(io, VERSION >= v"1.10" ? Threads.ngcthreads() : 0)
                println(io, VERSION >= v"1.10" ? Base.JLOptions().nsweepthreads : 0)
            end
            println("Benchmark resolved ", WeaveProjectFixture.VERSION_LABEL)
            ```

            ```julia, echo = false
            import WeaveLoadPathFixture
            println(WeaveLoadPathFixture.MESSAGE)
            ```
            """
        )
        output = joinpath(dirname(@__DIR__), "markdown", basename(benchmark), "environment.md")
        original_project = Base.active_project()
        original_load_path = copy(LOAD_PATH)
        try
            Pkg.activate(caller)
            pushfirst!(LOAD_PATH, package_directory)
            push!(LOAD_PATH, original_project)
            caller_version = @eval begin
                import WeaveProjectFixture
                WeaveProjectFixture.VERSION_LABEL
            end
            @test caller_version == "1.0.0"

            state = joinpath(benchmark, "environment.txt")
            for render in (
                    () -> SciMLBenchmarks.weave_file(benchmark, "environment.jmd", (:github,)),
                    () -> SciMLBenchmarks.weave_folder(benchmark, (:github,)),
                )
                render()
                @test Base.active_project() == joinpath(caller, "Project.toml")
                project, version, load_path, threads, interactive, gc_threads, sweep_threads = readlines(state)
                @test project == joinpath(benchmark, "Project.toml")
                @test version == "2.0.0"
                @test load_path == repr(["@"; LOAD_PATH])
                @test threads == string(Threads.nthreads())
                @test interactive == string(VERSION >= v"1.9" ? Threads.nthreads(:interactive) : 0)
                sweep = VERSION >= v"1.10" ? Base.JLOptions().nsweepthreads : 0
                mark = VERSION >= v"1.10" ? Threads.ngcthreads() - sweep : 0
                @test parse(Int, gc_threads) - parse(Int, sweep_threads) == mark
                @test parse(Int, sweep_threads) == sweep
                @test occursin("Loaded from the caller's package directory.", read(output, String))
                rm(state)
                rm(output)
            end

            @testset "CLI isolates preloaded packages" begin
                preload = joinpath(directory, "preload.jl")
                write(
                    preload, """
                    import WeaveProjectFixture
                    println("CLI preload resolved ", WeaveProjectFixture.VERSION_LABEL)
                    """
                )
                entrypoint = joinpath(dirname(@__DIR__), "benchmark.jl")
                page = joinpath(benchmark, "environment.jmd")
                command = `$(Base.julia_cmd()) --startup-file=no --project=$caller -L $preload $entrypoint $page`
                command = addenv(command, "JULIA_LOAD_PATH" => join(LOAD_PATH, Sys.iswindows() ? ';' : ':'))
                stdout = IOBuffer()
                @test success(pipeline(command; stdout))
                @test occursin("CLI preload resolved 1.0.0", String(take!(stdout)))
                @test occursin("Benchmark resolved 2.0.0", read(output, String))
            end
        finally
            append!(empty!(LOAD_PATH), original_load_path)
            Pkg.activate(original_project)
            rm(dirname(output); recursive = true, force = true)
            rm(joinpath(dirname(@__DIR__), "script", basename(benchmark)); recursive = true, force = true)
        end
    end
end

if VERSION >= v"1.10"
    @testset "Weave thread settings" begin
        mktempdir() do directory
            write(joinpath(directory, "Project.toml"), "[deps]\n")
            page = joinpath(directory, "threads.jmd")
            write(
                page, """
                # Weave thread settings

                ```julia, echo = false
                sweep = Int(Base.JLOptions().nsweepthreads)
                mark = Threads.ngcthreads() - sweep
                print("Weave thread settings: ", (Threads.nthreads(), Threads.nthreads(:interactive), mark, sweep))
                ```
                """
            )
            root = dirname(@__DIR__)
            entrypoint = joinpath(root, "benchmark.jl")
            output = joinpath(root, "markdown", basename(directory), "threads.md")
            script = joinpath(root, "script", basename(directory))
            try
                for (gc_threads, sweep) in (("1", 0), ("1,1", 1))
                    julia = `$(Base.julia_cmd()) --startup-file=no --threads=3,1 --gcthreads=$gc_threads --project=$(Base.active_project())`
                    for (entry, command) in (
                            ("single-file CLI", `$julia $entrypoint $page`),
                            ("folder CLI", `$julia $entrypoint $directory`),
                            ("weave_file", `$julia -e 'using SciMLBenchmarks; SciMLBenchmarks.weave_file(ARGS[1], "threads.jmd", (:github,))' $directory`),
                            ("weave_folder", `$julia -e 'using SciMLBenchmarks; SciMLBenchmarks.weave_folder(ARGS[1], (:github,))' $directory`),
                        )
                        @testset "$entry with --gcthreads=$gc_threads" begin
                            run(command)
                            @test occursin("Weave thread settings: (3, 1, 1, $sweep)", read(output, String))
                            rm(output)
                        end
                    end
                end
            finally
                rm(dirname(output); recursive = true, force = true)
                rm(script; recursive = true, force = true)
            end
        end
    end
end
