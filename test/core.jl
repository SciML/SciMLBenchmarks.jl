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

function workflow_job(workflow, name)
    return match(Regex("(?ms)^  $(name):\\n(?<body>.*?)(?=^  [a-zA-Z][a-zA-Z0-9_-]*:\\n|\\z)"), workflow)
end

@testset "AdaptiveSDE plot labels" begin
    benchmark = joinpath(
        dirname(@__DIR__), "benchmarks", "AdaptiveSDE", "AdaptiveEfficiencyTests.jmd"
    )
    assignment = only(
        line for line in eachline(benchmark) if startswith(strip(line), "leg=") ||
            startswith(strip(line), "leg =")
    )
    labels = Core.eval(@__MODULE__, Meta.parse(split(assignment, "="; limit = 2)[2]))
    @test size(labels) == (1, 3)
    @test vec(labels) == ["RSwM1", "RSwM2", "RSwM3"]
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
        @test Base.invokelatest(crn_grid, 1, Float32) == Float32[0.1]
        @test Base.invokelatest(crn_grid, 2, Float32) == Float32[0.1, 100]
    end

    lorenz_path = joinpath(benchmarks_dir, "lorenz_ensemble.jmd")
    lorenz_source = read(lorenz_path, String)
    @test occursin("make_ensemble(2)", lorenz_source)
    @test !occursin("make_ensemble(1)", lorenz_source)
    lorenz_expression = benchmark_assignment(lorenz_path, "plist")
    lorenz_grid = Core.eval(@__MODULE__, :((n, T) -> $lorenz_expression))
    @test collect(Base.invokelatest(lorenz_grid, 1, Float32)) == Float32[0]
    @test collect(Base.invokelatest(lorenz_grid, 4, Float32)) == Float32[0, 7, 14, 21]
end

@testset "subprocess" begin
    process = SciMLBenchmarks.@subprocess exit()
    @test success(process)
end

@testset "benchmark publication" begin
    workflow = read(joinpath(dirname(@__DIR__), ".github", "workflows", "benchmarks.yml"), String)
    benchmark_job = workflow_job(workflow, "benchmark")
    @test !isnothing(benchmark_job)
    if !isnothing(benchmark_job)
        @test occursin("- name: Publish to SciMLBenchmarksOutput", benchmark_job[:body])
        @test occursin("if: success() && github.ref == 'refs/heads/master'", benchmark_job[:body])
    end
end

@testset "NeuralNetworks V100 environment" begin
    folder = joinpath(dirname(@__DIR__), "benchmarks", "NeuralNetworks")
    preferences_path = joinpath(folder, "LocalPreferences.toml")
    project = read(joinpath(folder, "Project.toml"), String)
    manifest = read(joinpath(folder, "Manifest.toml"), String)
    python_dependencies = read(joinpath(folder, "CondaPkg.toml"), String)
    benchmark = read(joinpath(folder, "simple_networks.jmd"), String)

    @test isfile(preferences_path)
    if isfile(preferences_path)
        preferences = read(preferences_path, String)
        @test occursin(
            "[CUDA_Runtime_jll]\n__clear__ = [\"local\"]\nversion = \"12.8\"", preferences
        )
        @test occursin(
            "[Reactant_jll]\ncuda_version = \"12.8\"\ngpu = \"cuda\"", preferences
        )
    end
    @test occursin(
        "CUDA_Runtime_jll = \"76a88914-d11a-5bdc-97e0-2f5a05c973a2\"", project
    )
    @test occursin(
        "Reactant_jll = \"0192cb87-2b54-54ad-80e0-3be72ad8a3c0\"", project
    )
    @test occursin("CUDNN_jll = \"62b44479-cb7b-5706-934f-f13b2eb2e645\"", project)
    @test occursin("CUDNN_jll = \"=9.10.0\"", project)
    @test occursin("CUDA = \"5\"", project)
    @test occursin("cuDNN = \"=1.4.4\"", project)
    @test occursin("Reactant = \"=0.2.171\"", project)
    @test occursin(
        r"(?s)\[\[deps\.CUDNN_jll\]\].*?version = \"9\.10\.0\+0\"", manifest
    )
    @test occursin(
        r"(?s)\[\[deps\.Reactant\]\].*?version = \"0\.2\.171\"", manifest
    )
    @test occursin(
        r"(?s)\[\[deps\.Reactant_jll\]\].*?version = \"0\.0\.251\+0\"", manifest
    )
    @test occursin(
        r"(?s)\[\[deps\.GPUCompiler\]\].*?pinned = true.*?version = \"1\.9\.1\"", manifest
    )
    @test occursin("torch = \">=2.0,<2.11\"", python_dependencies)
    @test occursin(
        "ENV[\"XLA_REACTANT_GPU_MEM_FRACTION\"] = \"0.25\"\n" *
            "ENV[\"XLA_REACTANT_GPU_PREALLOCATE\"] = \"false\"\n" *
            "ENV[\"XLA_PYTHON_CLIENT_PREALLOCATE\"] = \"false\"",
        benchmark
    )
    worker_call = findfirst("python_output = read(", benchmark)
    cuda_import = findfirst("using CUDA, LuxCUDA", benchmark)
    @test !isnothing(worker_call)
    @test !isnothing(cuda_import)
    if !isnothing(worker_call) && !isnothing(cuda_import)
        @test first(worker_call) < first(cuda_import)
    end
    @test occursin("delete!(python_env, \"LD_LIBRARY_PATH\")", benchmark)
    @test occursin("PythonCall.python_executable_path()", benchmark)
    @test !occursin("pyimport(", benchmark)
    @test !occursin("nn_utils.", benchmark)

    python_helper = read(joinpath(folder, "nn_benchmark_utils.py"), String)
    @test occursin("if __name__ == \"__main__\":", python_helper)
    @test occursin("TIMING\\t{framework}\\t{model}\\t{operation}", python_helper)
    @test occursin("require_jax_gpu()", python_helper)
    @test occursin("require_torch_cuda()", python_helper)
    reactant_steps = (
        ("reactant_step", "ts"),
        ("reactant_gelu_step", "ts_gelu"),
        ("reactant_bn_step", "ts_bn"),
        ("reactant_lenet_step", "ts_lenet"),
        ("reactant_resnet_step", "ts_resnet"),
    )
    for (step, state) in reactant_steps
        @test occursin("$state, _ = $step($state, x_ra, y_ra)", benchmark)
    end
    @test occursin(
        r"Training\.TrainState\(\s*lux_mlp_bn,\s*ps_bn_ra,\s*st_bn_train_ra", benchmark
    )
    @test occursin(
        r"Training\.TrainState\(\s*lux_resnet,\s*ps_resnet_ra,\s*st_resnet_train_ra",
        benchmark,
    )
    @test !occursin("REACTANT_TRAINING_COMPILE_OPTIONS", benchmark)
    @test !occursin("compile_options=", benchmark)
    @test count("sync=true", benchmark) == 5
    @test count(r"@benchmark \$reactant_", benchmark) == 5
    @test !occursin("compiled_reactant_step", benchmark)
    @test !occursin(r"@compile reactant_", benchmark)
end

@testset "superseded pull request cancellation" begin
    workflow = read(joinpath(dirname(@__DIR__), ".github", "workflows", "benchmarks.yml"), String)
    @test occursin("types: [opened, synchronize, reopened, closed]", workflow)
    @test occursin("format('pr-{0}', github.event.pull_request.number)", workflow)
    @test occursin("github.event_name == 'pull_request' || github.ref != 'refs/heads/master'", workflow)
    @test occursin("github.event.action != 'closed'", workflow)
    pull_request_trigger = match(
        r"(?ms)^  pull_request:\n(?<body>.*?)^  workflow_dispatch:", workflow
    )
    @test !isnothing(pull_request_trigger)
    if !isnothing(pull_request_trigger)
        @test !occursin("paths:", pull_request_trigger[:body])
    end

    benchmark_job = workflow_job(workflow, "benchmark")
    @test !isnothing(benchmark_job)
    if !isnothing(benchmark_job)
        @test occursin(
            "format('pr-{0}', github.event.pull_request.number) || github.ref",
            benchmark_job[:body],
        )
        @test occursin("\${{ matrix.target }}", benchmark_job[:body])
        @test occursin("cancel-in-progress: true", benchmark_job[:body])
    end

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
