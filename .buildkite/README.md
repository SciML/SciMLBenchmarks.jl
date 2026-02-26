# Multi-Queue Buildkite Configuration

This directory contains the configuration for running SciMLBenchmarks on multiple Buildkite queues, enabling GPU-intensive benchmarks to run on GPU-enabled compute resources while maintaining CPU-only benchmarks on the standard queue.

## Files

- `queue_config.yml`: Central configuration mapping benchmarks to queues
- `generate_pipeline.jl`: Dynamic pipeline generator that creates appropriate Buildkite steps
- `test_sciml.yml`: Main Buildkite pipeline configuration
- `path_processors/project-coalescing`: Script to determine which benchmarks need rebuilding

## How It Works

1. **Queue Detection**: When a benchmark runs, the system determines which queue to use based on:
   - Benchmark-specific `[buildkite]` section in `Project.toml` (highest priority)
   - Central mapping in `queue_config.yml`
   - Default queue (`juliaecosystem`) as fallback

2. **Dynamic Pipeline Generation**: The `generate_pipeline.jl` script:
   - Reads the queue configuration
   - Processes changed files to determine which benchmarks to run
   - Generates appropriate Buildkite steps with correct queue assignments
   - Uploads the dynamic pipeline to Buildkite

3. **Queue-Specific Configuration**: Each queue can have different:
   - Architecture requirements (`arch`)
   - Operating system (`os`)
   - Environment variables (e.g., GPU-specific settings)
   - Timeout values

## Available Queues

### `juliaecosystem` (Default)
- **Purpose**: CPU-only compute queue
- **Architecture**: x86_64 Linux
- **Usage**: Standard Julia benchmarks, ODE solvers, linear algebra

### `gpu`
- **Purpose**: GPU-enabled compute queue  
- **Architecture**: x86_64 Linux with GPU support
- **Usage**: Neural PDE benchmarks, GPU-accelerated simulations
- **Environment**: Includes CUDA-related environment variables

## Configuring Benchmarks for Specific Queues

### Method 1: Project.toml (Recommended)

Add a `[buildkite]` section to your benchmark's `Project.toml`:

```toml
[deps]
# ... your dependencies

[buildkite]
# Specify the queue this benchmark should use
queue = "gpu"
```

### Method 2: Central Configuration

Edit `.buildkite/queue_config.yml` and add your benchmark folder:

```yaml
benchmark_queues:
  YourBenchmarkFolder: "gpu"
```

## Adding New Queues

1. **Define the queue** in `queue_config.yml`:
```yaml
queues:
  your_new_queue:
    arch: "x86_64"
    os: "linux" 
    description: "Description of queue capabilities"
```

2. **Configure benchmarks** to use the new queue using either method above

3. **Set up the Buildkite infrastructure** to have agents with the `your_new_queue` tag

## Testing

To test the pipeline generation locally:

```bash
# Test with specific files
julia .buildkite/generate_pipeline.jl benchmarks/PINNOptimizers/poisson.jmd

# Test with benchmark folders
julia .buildkite/generate_pipeline.jl benchmarks/PINNOptimizers
```

## Backward Compatibility

- Benchmarks without queue specifications continue to use `juliaecosystem`
- Existing CI behavior is preserved for benchmarks that don't opt into GPU queues
- The system gracefully handles missing configuration files or invalid queue names

## Environment Variables

For GPU benchmarks, the following environment variables are automatically set:

- `JULIA_CUDA_USE_BINARYBUILDER=false`: Use system CUDA instead of artifacts
- `JULIA_GPU_ALLOW_DEFAULT=true`: Allow default GPU selection

## Troubleshooting

### Pipeline Not Generated
- Check that `YAML.jl` is available in the Julia environment
- Verify the `queue_config.yml` syntax is valid YAML
- Ensure changed files are properly detected

### Wrong Queue Assignment
- Check the `[buildkite]` section in your benchmark's `Project.toml`
- Verify the benchmark name matches the folder name in `queue_config.yml`
- Confirm the queue exists in the `queues` section

### GPU Jobs Failing
- Verify GPU agents are available and tagged with `queue: "gpu"`
- Check that CUDA dependencies are properly configured on GPU nodes
- Review environment variables for GPU-specific settings