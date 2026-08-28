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

function benchmark_assignment(path, variable)
    prefix = string(variable, " = ")
    line = only(line for line in eachline(path) if startswith(strip(line), prefix))
    return Meta.parse(strip(split(line, "="; limit = 2)[2]))
end

@testset "DiffEqGPU singleton parameter grids" begin
    benchmarks_dir = joinpath(dirname(@__DIR__), "benchmarks", "DiffEqGPU")

    crn_path = joinpath(benchmarks_dir, "crn_sde.jmd")
    crn_source = read(crn_path, String)
    @test occursin("ps = crn_parameters(2)", crn_source)
    @test !occursin("ps = crn_parameters(1)", crn_source)
    @test occursin("GPUEM(), KERNEL; trajectories = 2", crn_source)
    @test !occursin("GPUEM(), KERNEL; trajectories = 1", crn_source)
    for variable in ("S_grid", "D_grid")
        crn_expression = benchmark_assignment(crn_path, variable)
        crn_grid = Core.eval(@__MODULE__, :((N, T) -> $crn_expression))
        @test crn_grid(1, Float32) == Float32[0.1]
        @test crn_grid(2, Float32) == Float32[0.1, 100]
    end

    lorenz_path = joinpath(benchmarks_dir, "lorenz_ensemble.jmd")
    lorenz_source = read(lorenz_path, String)
    @test occursin("make_ensemble(2)", lorenz_source)
    @test !occursin("make_ensemble(1)", lorenz_source)
    lorenz_expression = benchmark_assignment(lorenz_path, "plist")
    lorenz_grid = Core.eval(@__MODULE__, :((n, T) -> $lorenz_expression))
    @test collect(lorenz_grid(1, Float32)) == Float32[0]
    @test collect(lorenz_grid(4, Float32)) == Float32[0, 7, 14, 21]
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

@testset "V100 CUDA runtime preferences" begin
    benchmarks_dir = joinpath(dirname(@__DIR__), "benchmarks")
    cuda_runtime = "CUDA_Runtime_jll = \"76a88914-d11a-5bdc-97e0-2f5a05c973a2\""
    for folder in ("DiffEqGPU", "PINNErrorsVsTime", "PINNOptimizers")
        preferences = read(joinpath(benchmarks_dir, folder, "LocalPreferences.toml"), String)
        project = read(joinpath(benchmarks_dir, folder, "Project.toml"), String)
        @test occursin("version = \"12.9\"", preferences)
        @test occursin(cuda_runtime, project)
    end
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
