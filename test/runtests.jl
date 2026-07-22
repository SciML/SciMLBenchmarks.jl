using SciMLBenchmarks, Test

@testset "Explicit Imports" begin
    include("explicit_imports.jl")
end

@testset "weave_file resolves stale manifests" begin
    mktempdir() do folder
        write(
            joinpath(folder, "Project.toml"),
            """
            [deps]
            Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
            """,
        )
        write(
            joinpath(folder, "Manifest.toml"),
            """
            julia_version = "$(VERSION)"
            manifest_format = "2.0"
            project_hash = "0000000000000000000000000000000000000000"
            """,
        )

        @test SciMLBenchmarks.weave_file(folder, "unused.jmd", ()) === nothing
        @test occursin("[[deps.Dates]]", read(joinpath(folder, "Manifest.toml"), String))
    end
end

@testset "weave_file" begin
    benchmarks_dir = joinpath(dirname(@__DIR__), "benchmarks")
    SciMLBenchmarks.weave_file(joinpath(benchmarks_dir, "Testing"), "test.jmd")

    #@test isfile(joinpath(dirname(@__DIR__), "script", "Testing", "test.jl"))
    #@test isfile(joinpath(dirname(@__DIR__), "html", "Testing", "test.html"))
    @test isfile(joinpath(dirname(@__DIR__), "markdown", "Testing", "test.md"))
end
