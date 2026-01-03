#!/usr/bin/env python3
"""
Python benchmark script for chemical reaction networks using GillesPy2 and COPASI.
This script is designed to be called from Julia using PyCall.
"""

import os
import sys
import time
import json
import numpy as np
from pathlib import Path
import argparse

# Try to import required packages
try:
    import gillespy2
    from gillespy2.core import gillespyError
    GILLESPY2_AVAILABLE = True
except ImportError:
    GILLESPY2_AVAILABLE = False
    print("Warning: GillesPy2 not available. Install with: pip install gillespy2")

try:
    import basico
    COPASI_AVAILABLE = True
except ImportError:
    COPASI_AVAILABLE = False
    print("Warning: COPASI not available. Install with: pip install python-copasi")

try:
    import libsbml
    LIBSBML_AVAILABLE = True
except ImportError:
    LIBSBML_AVAILABLE = False
    print("Warning: libsbml not available. Install with: pip install python-libsbml")

def benchmark_gillespy2_ode(model_file, model_name, time_span=10.0, num_runs=10):
    """Benchmark GillesPy2 ODE solver"""
    if not GILLESPY2_AVAILABLE:
        return {
            "solver": "GillesPy2_ODE",
            "median_time": float('inf'),
            "success": False,
            "error": "GillesPy2 not available"
        }
    
    try:
        # Load model from SBML file
        model = gillespy2.import_SBML(model_file)
        
        # Set simulation parameters
        model.timespan = np.linspace(0, time_span, 101)
        
        # Warmup run
        try:
            results = model.run(solver=gillespy2.ODESolver, show_labels=False)
        except Exception as e:
            return {
                "solver": "GillesPy2_ODE",
                "median_time": float('inf'),
                "success": False,
                "error": f"Warmup failed: {str(e)}"
            }
        
        # Benchmark runs
        times = []
        for _ in range(num_runs):
            start_time = time.time()
            results = model.run(solver=gillespy2.ODESolver, show_labels=False)
            elapsed = (time.time() - start_time) * 1000  # Convert to ms
            times.append(elapsed)
        
        return {
            "solver": "GillesPy2_ODE",
            "median_time": np.median(times),
            "min_time": np.min(times),
            "max_time": np.max(times),
            "std_time": np.std(times),
            "success": True,
            "model": model_name,
            "num_runs": num_runs
        }
    except Exception as e:
        return {
            "solver": "GillesPy2_ODE",
            "median_time": float('inf'),
            "success": False,
            "error": str(e),
            "model": model_name
        }

def benchmark_gillespy2_ssa(model_file, model_name, time_span=10.0, num_runs=10):
    """Benchmark GillesPy2 SSA solver"""
    if not GILLESPY2_AVAILABLE:
        return {
            "solver": "GillesPy2_SSA",
            "median_time": float('inf'),
            "success": False,
            "error": "GillesPy2 not available"
        }
    
    try:
        # Load model from SBML file
        model = gillespy2.import_SBML(model_file)
        
        # Set simulation parameters
        model.timespan = np.linspace(0, time_span, 101)
        
        # Warmup run
        try:
            results = model.run(solver=gillespy2.SSACSolver, show_labels=False)
        except Exception as e:
            return {
                "solver": "GillesPy2_SSA",
                "median_time": float('inf'),
                "success": False,
                "error": f"Warmup failed: {str(e)}"
            }
        
        # Benchmark runs
        times = []
        for _ in range(num_runs):
            start_time = time.time()
            results = model.run(solver=gillespy2.SSACSolver, show_labels=False)
            elapsed = (time.time() - start_time) * 1000  # Convert to ms
            times.append(elapsed)
        
        return {
            "solver": "GillesPy2_SSA",
            "median_time": np.median(times),
            "min_time": np.min(times),
            "max_time": np.max(times),
            "std_time": np.std(times),
            "success": True,
            "model": model_name,
            "num_runs": num_runs
        }
    except Exception as e:
        return {
            "solver": "GillesPy2_SSA",
            "median_time": float('inf'),
            "success": False,
            "error": str(e),
            "model": model_name
        }

def benchmark_copasi_ode(model_file, model_name, time_span=10.0, num_runs=10):
    """Benchmark COPASI ODE solver"""
    if not COPASI_AVAILABLE:
        return {
            "solver": "COPASI_ODE",
            "median_time": float('inf'),
            "success": False,
            "error": "COPASI not available"
        }
    
    try:
        # Load model
        model = basico.load_model(model_file)
        
        # Set simulation parameters
        basico.set_task_settings(basico.T.TIME_COURSE, {
            'duration': time_span,
            'intervals': 100,
            'output_event': False
        })
        
        # Warmup run
        try:
            result = basico.run_time_course()
        except Exception as e:
            return {
                "solver": "COPASI_ODE",
                "median_time": float('inf'),
                "success": False,
                "error": f"Warmup failed: {str(e)}"
            }
        
        # Benchmark runs
        times = []
        for _ in range(num_runs):
            start_time = time.time()
            result = basico.run_time_course()
            elapsed = (time.time() - start_time) * 1000  # Convert to ms
            times.append(elapsed)
        
        return {
            "solver": "COPASI_ODE",
            "median_time": np.median(times),
            "min_time": np.min(times),
            "max_time": np.max(times),
            "std_time": np.std(times),
            "success": True,
            "model": model_name,
            "num_runs": num_runs
        }
    except Exception as e:
        return {
            "solver": "COPASI_ODE",
            "median_time": float('inf'),
            "success": False,
            "error": str(e),
            "model": model_name
        }

def benchmark_copasi_ssa(model_file, model_name, time_span=10.0, num_runs=10):
    """Benchmark COPASI SSA solver"""
    if not COPASI_AVAILABLE:
        return {
            "solver": "COPASI_SSA",
            "median_time": float('inf'),
            "success": False,
            "error": "COPASI not available"
        }
    
    try:
        # Load model
        model = basico.load_model(model_file)
        
        # Set simulation parameters for stochastic simulation
        basico.set_task_settings(basico.T.TIME_COURSE, {
            'duration': time_span,
            'intervals': 100,
            'output_event': False,
            'method': {
                'name': 'Stochastic (Gibson & Bruck)',
                'use_random_seed': True,
                'random_seed': 1
            }
        })
        
        # Warmup run
        try:
            result = basico.run_time_course()
        except Exception as e:
            return {
                "solver": "COPASI_SSA",
                "median_time": float('inf'),
                "success": False,
                "error": f"Warmup failed: {str(e)}"
            }
        
        # Benchmark runs
        times = []
        for _ in range(num_runs):
            start_time = time.time()
            result = basico.run_time_course()
            elapsed = (time.time() - start_time) * 1000  # Convert to ms
            times.append(elapsed)
        
        return {
            "solver": "COPASI_SSA",
            "median_time": np.median(times),
            "min_time": np.min(times),
            "max_time": np.max(times),
            "std_time": np.std(times),
            "success": True,
            "model": model_name,
            "num_runs": num_runs
        }
    except Exception as e:
        return {
            "solver": "COPASI_SSA",
            "median_time": float('inf'),
            "success": False,
            "error": str(e),
            "model": model_name
        }

def run_all_benchmarks(model_file, model_name, time_span=10.0, num_runs=10):
    """Run all available Python benchmarks for a model"""
    results = []
    
    # GillesPy2 benchmarks
    results.append(benchmark_gillespy2_ode(model_file, model_name, time_span, num_runs))
    results.append(benchmark_gillespy2_ssa(model_file, model_name, time_span, num_runs))
    
    # COPASI benchmarks
    results.append(benchmark_copasi_ode(model_file, model_name, time_span, num_runs))
    results.append(benchmark_copasi_ssa(model_file, model_name, time_span, num_runs))
    
    return results

def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(description="Benchmark chemical reaction networks with Python tools")
    parser.add_argument("model_file", help="Path to SBML model file")
    parser.add_argument("model_name", help="Name of the model")
    parser.add_argument("--time-span", type=float, default=10.0, help="Simulation time span")
    parser.add_argument("--num-runs", type=int, default=10, help="Number of benchmark runs")
    parser.add_argument("--output", help="Output JSON file path")
    
    args = parser.parse_args()
    
    # Run benchmarks
    results = run_all_benchmarks(args.model_file, args.model_name, args.time_span, args.num_runs)
    
    # Print results
    print(f"\\nBenchmark results for {args.model_name}:")
    for result in results:
        if result["success"]:
            print(f"  {result['solver']}: {result['median_time']:.2f} ms")
        else:
            print(f"  {result['solver']}: FAILED - {result['error']}")
    
    # Save results if requested
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\\nResults saved to {args.output}")
    
    return results

if __name__ == "__main__":
    main()