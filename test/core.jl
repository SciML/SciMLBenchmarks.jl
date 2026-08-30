using SciMLBenchmarks
using Test

@testset "documentation workflow consistency" begin
    repository_root = dirname(@__DIR__)
    @test !isfile(joinpath(repository_root, ".github", "workflows", "DocPreviewCleanup.yml"))
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

@testset "benchmark publication" begin
    workflow = read(joinpath(dirname(@__DIR__), ".github", "workflows", "benchmarks.yml"), String)
    benchmark_job = match(
        r"(?ms)^  benchmark:\n(?<body>.*?)(?=^  [a-zA-Z][a-zA-Z0-9_-]*:\n|\z)", workflow
    )
    @test !isnothing(benchmark_job)
    if !isnothing(benchmark_job)
        @test occursin("- name: Publish to SciMLBenchmarksOutput", benchmark_job[:body])
        @test occursin("if: success() && github.ref == 'refs/heads/master'", benchmark_job[:body])
    end
end

@testset "superseded pull request cancellation" begin
    workflow = read(joinpath(dirname(@__DIR__), ".github", "workflows", "benchmarks.yml"), String)
    @test occursin("types: [opened, synchronize, reopened, closed]", workflow)
    @test occursin("format('pr-{0}', github.event.pull_request.number)", workflow)
    @test occursin("github.event_name == 'pull_request' || github.ref != 'refs/heads/master'", workflow)
    @test occursin("github.event.action != 'closed'", workflow)

    cancellation_path = joinpath(
        dirname(@__DIR__), ".github", "workflows", "cancel-superseded-benchmarks.yml"
    )
    @test isfile(cancellation_path)
    if isfile(cancellation_path)
        cancellation = read(cancellation_path, String)
        @test occursin("workflow_run:", cancellation)
        @test occursin("workflows: [Benchmarks]", cancellation)
        @test occursin("types: [requested]", cancellation)
        @test occursin("actions: write", cancellation)
        @test occursin("run.run_number >= source.run_number", cancellation)
        @test occursin("run.pull_requests", cancellation)
        @test occursin("/force-cancel", cancellation)
    end
end

@testset "Jumps spatial benchmark runtime" begin
    benchmark = read(
        joinpath(dirname(@__DIR__), "benchmarks", "Jumps", "Spatial_Signaling_Sanft.jmd"),
        String
    )
    @test occursin("samples=3 seconds=300", benchmark)
    @test !occursin("seconds=3600", benchmark)
end

@testset "StiffODE work-precision names" begin
    benchmark = read(
        joinpath(dirname(@__DIR__), "benchmarks", "StiffODE", "Bruss.jmd"), String
    )
    checked = 0
    for chunk in eachmatch(r"(?ms)^```julia[^\n]*\n(?<code>.*?)^```", benchmark)
        code = join(
            (first(split(line, "#"; limit = 2)) for line in eachline(IOBuffer(chunk[:code]))),
            "\n"
        )
        setups = match(r"(?ms)setups = \[(?<values>.*?)\]", code)
        names = match(r"(?ms)names = \[(?<values>.*?)\]", code)
        if !isnothing(setups) && !isnothing(names) && occursin("WorkPrecisionSet", code)
            @test count("Dict(", setups[:values]) == count("\"", names[:values]) ÷ 2
            checked += 1
        end
    end
    @test checked == count("names = names", benchmark)
end

@testset "weave_file" begin
    benchmarks_dir = joinpath(dirname(@__DIR__), "benchmarks")
    SciMLBenchmarks.weave_file(joinpath(benchmarks_dir, "Testing"), "test.jmd")

    output = joinpath(dirname(@__DIR__), "markdown", "Testing", "test.md")
    @test isfile(output)
    @test occursin("To locally run this benchmark", read(output, String))
end
