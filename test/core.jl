using SciMLBenchmarks
using Test

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
