#!/usr/bin/env julia

using YAML, Pkg

"""
Generate dynamic Buildkite pipeline based on queue configuration and changed files.
"""

# Read queue configuration
config_path = joinpath(@__DIR__, "queue_config.yml")
if !isfile(config_path)
    error("Queue configuration file not found: $config_path")
end

config = YAML.load_file(config_path)
default_queue = get(config, "default_queue", "juliaecosystem")
queues = get(config, "queues", Dict())
benchmark_queues = get(config, "benchmark_queues", Dict())

# Get changed files from environment or arguments
changed_files = if haskey(ENV, "BUILDKITE_CHANGED_FILES")
    split(ENV["BUILDKITE_CHANGED_FILES"], "\n")
elseif length(ARGS) > 0
    ARGS
else
    # Default to testing folder for demonstration
    ["benchmarks/Testing/test.jmd"]
end

# Process changed files using the existing project-coalescing logic
coalescing_script = joinpath(@__DIR__, "path_processors", "project-coalescing")
if isfile(coalescing_script) && haskey(ENV, "BUILDKITE_AGENT_ACCESS_TOKEN")
    # Only use project-coalescing if we're in a real Buildkite environment
    try
        build_targets = strip(read(`$coalescing_script $changed_files`, String))
        targets = split(build_targets)
    catch e
        @warn "project-coalescing failed, using fallback" exception=e
        targets = changed_files
    end
else
    # Fallback: use the provided targets directly or extract from file paths
    targets = if any(startswith(f, "benchmarks/") for f in changed_files)
        unique([startswith(f, "benchmarks/") ? f : "benchmarks/$f" for f in changed_files])
    else
        changed_files
    end
end

# Generate pipeline steps
steps = []

for target in targets
    if isfile(target)
        # Single file target
        folder = dirname(target)
        benchmark_name = basename(folder)
    else
        # Folder target  
        benchmark_name = basename(target)
    end
    
    # Determine queue for this benchmark
    # Priority: 1. Project.toml metadata, 2. queue_config.yml, 3. default
    queue_name = default_queue
    
    # Check if benchmark has its own Project.toml with queue specification
    project_path = if isfile(target)
        joinpath(dirname(target), "Project.toml")
    else
        joinpath(target, "Project.toml")
    end
    
    if isfile(project_path)
        try
            project_content = read(project_path, String)
            # Simple regex to extract queue from [buildkite] section
            queue_match = match(r"\[buildkite\].*?queue\s*=\s*\"([^\"]+)\""s, project_content)
            if queue_match !== nothing
                queue_name = queue_match.captures[1]
            end
        catch e
            @warn "Could not parse Project.toml for queue info: $(project_path)" exception=e
        end
    end
    
    # Fallback to configuration file
    if queue_name == default_queue
        queue_name = get(benchmark_queues, benchmark_name, default_queue)
    end
    
    queue_info = get(queues, queue_name, Dict("arch" => "x86_64", "os" => "linux"))
    
    # Create step configuration
    step = Dict(
        "label" => ":julia: $(benchmark_name) on $(queue_name)",
        "command" => "julia benchmark.jl $(target)",
        "plugins" => [
            Dict("JuliaCI/julia#v1" => Dict("version" => "1.10")),
            Dict("JuliaCI/julia-test#v1" => nothing)
        ],
        "timeout_in_minutes" => 60,  # Increased timeout for potential GPU jobs
        "artifact_paths" => [
            "html/$(benchmark_name)/*.html",
            "markdown/$(benchmark_name)/*.md", 
            "notebook/$(benchmark_name)/*.ipynb",
            "pdf/$(benchmark_name)/*.pdf",
            "script/$(benchmark_name)/*.jl"
        ],
        "agents" => Dict(
            "queue" => queue_name,
            "arch" => get(queue_info, "arch", "x86_64"),
            "os" => get(queue_info, "os", "linux")
        )
    )
    
    # Add GPU-specific environment variables if needed
    if queue_name == "gpu"
        step["env"] = Dict(
            "JULIA_CUDA_USE_BINARYBUILDER" => "false",
            "JULIA_GPU_ALLOW_DEFAULT" => "true"
        )
    end
    
    push!(steps, step)
end

# Generate final pipeline
pipeline = Dict("steps" => steps)

# Output as YAML
println(YAML.write(pipeline))