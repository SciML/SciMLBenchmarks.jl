# Catalyst Chemical Reaction Network Benchmarks

This directory contains benchmarks for the Catalyst chemical reaction network modeling ecosystem, based on the comprehensive study published in *PLOS Computational Biology* (2023).

## Overview

The benchmark suite compares the performance of different chemical reaction network (CRN) modeling approaches:
- **Catalyst (Julia)** - The primary focus, with multiple ODE and SSA solvers
- **GillesPy2 (Python)** - Popular Python package for stochastic simulation
- **COPASI (Python)** - Comprehensive biochemical modeling tool
- **GillespieSSA2 (R)** - R package for stochastic simulation

## Files

### Main Benchmark Files
- `catalyst_suite.jmd` - Main benchmark file (Weave.jl format)
- `model_utils.jl` - Utility functions for model loading and processing
- `test_benchmarks.jl` - Test suite to validate benchmark setup

### External Language Support
- `python_benchmarks.py` - Python benchmarking script for GillesPy2 and COPASI
- `r_benchmarks.R` - R benchmarking script for GillespieSSA2

### Configuration
- `Project.toml` - Julia package dependencies
- `README.md` - This file

## Models Tested

The benchmark suite includes five biological network models of increasing complexity:

1. **Multistate** (9 species, 18 reactions)
   - Simple multi-state protein system
   - Good for testing basic functionality

2. **Multisite2** (66 species, 288 reactions)
   - Multi-site protein phosphorylation
   - Medium complexity system

3. **EGFR Network** (356 species, 3,749 reactions)
   - Epidermal growth factor receptor signaling
   - Real biological complexity

4. **BCR Network** (1,122 species, 24,388 reactions)
   - B-cell receptor signaling pathway
   - Large-scale biological system

5. **FcεRI Network** (3,744 species, 58,276 reactions)
   - High-affinity IgE receptor signaling
   - Extremely large biological system

## Usage

### Prerequisites

**Julia packages** (automatically installed via Project.toml):
- Catalyst
- DifferentialEquations
- JumpProcesses
- BenchmarkTools
- PyCall
- RCall

**Python packages** (install manually):
```bash
pip install gillespy2 python-copasi python-libsbml numpy
```

**R packages** (install manually):
```r
install.packages(c("GillespieSSA2", "jsonlite"))
```

### Running the Benchmarks

1. **Test the setup** (recommended first step):
```julia
julia> include("test_benchmarks.jl")
```

2. **Run the full benchmark suite**:
```julia
julia> using Weave
julia> weave("catalyst_suite.jmd")
```

3. **Run specific components**:
```julia
julia> include("model_utils.jl")
julia> problems = create_catalyst_problems("multistate")
julia> # Run your own benchmarks...
```

### Python Benchmarks

The Python benchmarks can be run independently:
```bash
python python_benchmarks.py model_file.xml model_name --num-runs 10
```

### R Benchmarks

The R benchmarks can be run independently:
```bash
Rscript r_benchmarks.R model_name simple 10.0 10
```

## Expected Results

Based on the original study, you should expect:
- **Catalyst typically outperforms other tools by 1-2 orders of magnitude**
- **Performance scales well with model complexity**
- **Different solvers perform optimally for different model types**
- **ODE methods generally faster than SSA for deterministic simulations**

## Troubleshooting

### Common Issues

1. **Model files not found**: The benchmark automatically downloads model files from the GitHub repository. If this fails, check your internet connection.

2. **Python/R packages not found**: Install the required packages manually (see Prerequisites).

3. **PyCall/RCall issues**: These packages require proper Python/R installations. See their documentation for setup instructions.

4. **Memory issues with large models**: The largest models (BCR, FcεRI) require significant memory. Consider running on a machine with at least 8GB RAM.

### Testing

Always run the test suite first:
```julia
julia> include("test_benchmarks.jl")
```

This will verify:
- Model loading works correctly
- Required packages are available
- Basic solving functionality works
- Benchmarking infrastructure is operational

## Customization

### Adding New Models

To add a new model:
1. Place the model file (`.net`, `.xml`, or `.bngl`) in the `data/` directory
2. Update the model list in `catalyst_suite.jmd`
3. Add appropriate initial conditions in `model_utils.jl`

### Adding New Solvers

To benchmark new solvers:
1. Add the solver to the appropriate list in `catalyst_suite.jmd`
2. Ensure the solver is available in the current Julia environment
3. Test with a simple model first

### Modifying Benchmarking Parameters

Common parameters to adjust:
- `tspan`: Simulation time span
- `abstol`/`reltol`: Solver tolerances
- `num_runs`: Number of benchmark repetitions
- `seconds`: Benchmarking duration

## Citation

If you use these benchmarks in your research, please cite:

```bibtex
@article{loman2023catalyst,
  title={Catalyst: fast and flexible modeling of reaction networks},
  author={Loman, Torkel E and Ma, Yingbo and Ilin, Vasily and Gowda, Shashi and Korsbo, Niklas and Yewale, Nikhil and Rackauckas, Chris and Isaacson, Samuel A},
  journal={PLoS computational biology},
  volume={19},
  number={10},
  pages={e1011530},
  year={2023},
  publisher={Public Library of Science San Francisco, CA USA}
}
```

## Support

For issues with:
- **Catalyst**: [Catalyst.jl GitHub Issues](https://github.com/SciML/Catalyst.jl/issues)
- **SciMLBenchmarks**: [SciMLBenchmarks.jl GitHub Issues](https://github.com/SciML/SciMLBenchmarks.jl/issues)
- **Original benchmark code**: [Catalyst_PLOS_COMPBIO_2023 GitHub Issues](https://github.com/SciML/Catalyst_PLOS_COMPBIO_2023/issues)