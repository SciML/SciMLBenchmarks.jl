using SciMLBenchmarks
using Test
using TOML

@testset "benchmark dependencies" begin
    project = TOML.parsefile(
        joinpath(dirname(@__DIR__), "benchmarks", "MultiLanguage", "Project.toml")
    )
    for package in (
            "ODEInterfaceDiffEq", "Plots", "SciPyDiffEq", "Sundials", "deSolveDiffEq",
        )
        @test haskey(project["deps"], package)
    end
end

@testset "IJulia extension" begin
    fallback_method = which(SciMLBenchmarks.open_notebooks, Tuple{})
    @test fallback_method.module === SciMLBenchmarks

    @eval using IJulia
    extension = Base.get_extension(SciMLBenchmarks, :SciMLBenchmarksIJuliaExt)
    @test extension !== nothing
    @test fallback_method in methods(SciMLBenchmarks.open_notebooks)
    @test which(SciMLBenchmarks.open_notebooks, Tuple{}).module === extension
end

@testset "benchmark priority" begin
    mktempdir() do directory
        prioritized = joinpath(directory, "prioritized.jmd")
        write(prioritized, "---\ntitle: Prioritized\npriority: 42\n---\n")
        @test SciMLBenchmarks._benchmark_priority(prioritized) == 42

        unprioritized = joinpath(directory, "unprioritized.jmd")
        write(unprioritized, "---\ntitle: Unprioritized\n---\n")
        @test SciMLBenchmarks._benchmark_priority(unprioritized) == 0

        no_header = joinpath(directory, "no_header.jmd")
        write(no_header, "# No front matter\n")
        @test SciMLBenchmarks._benchmark_priority(no_header) == 0
    end
end

@testset "subprocess" begin
    process = SciMLBenchmarks.@subprocess exit()
    @test success(process)
end

@testset "weave_file" begin
    benchmarks_dir = joinpath(dirname(@__DIR__), "benchmarks")
    SciMLBenchmarks.weave_file(joinpath(benchmarks_dir, "Testing"), "test.jmd")

    output = joinpath(dirname(@__DIR__), "markdown", "Testing", "test.md")
    @test isfile(output)
    @test occursin("To locally run this benchmark", read(output, String))
end
