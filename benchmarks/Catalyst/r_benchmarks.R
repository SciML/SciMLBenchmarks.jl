#!/usr/bin/env Rscript
#
# R benchmark script for chemical reaction networks using GillespieSSA2
# This script is designed to be called from Julia using RCall
#

# Load required libraries
library(GillespieSSA2)
library(jsonlite)

# Check if required packages are available
check_packages <- function() {
  packages_available <- list(
    GillespieSSA2 = requireNamespace("GillespieSSA2", quietly = TRUE),
    jsonlite = requireNamespace("jsonlite", quietly = TRUE)
  )
  
  if (!all(unlist(packages_available))) {
    missing <- names(packages_available)[!unlist(packages_available)]
    stop(paste("Missing required packages:", paste(missing, collapse = ", ")))
  }
  
  return(packages_available)
}

# Helper function to create simple reaction network for demonstration
create_simple_network <- function() {
  # Simple birth-death process
  reactions <- list(
    reaction("birth", c(), c("X"), function(x, p) p$alpha),
    reaction("death", c("X"), c(), function(x, p) p$beta * x["X"])
  )
  
  initial_state <- c(X = 10)
  params <- list(alpha = 1.0, beta = 0.1)
  
  return(list(reactions = reactions, initial_state = initial_state, params = params))
}

# More complex reaction network (multistate-like)
create_multistate_network <- function() {
  # Multi-state protein system
  reactions <- list(
    reaction("A_on", c(), c("A"), function(x, p) p$kon_A),
    reaction("A_off", c("A"), c(), function(x, p) p$koff_A * x["A"]),
    reaction("B_on", c("A"), c("A", "B"), function(x, p) p$kon_B * x["A"]),
    reaction("B_off", c("B"), c(), function(x, p) p$koff_B * x["B"]),
    reaction("C_form", c("A", "B"), c("C"), function(x, p) p$kform * x["A"] * x["B"]),
    reaction("C_degrade", c("C"), c(), function(x, p) p$kdeg * x["C"])
  )
  
  initial_state <- c(A = 100, B = 0, C = 0)
  params <- list(
    kon_A = 1.0,
    koff_A = 0.1,
    kon_B = 0.5,
    koff_B = 0.2,
    kform = 0.01,
    kdeg = 0.05
  )
  
  return(list(reactions = reactions, initial_state = initial_state, params = params))
}

# Benchmark GillespieSSA2
benchmark_gillespie_ssa2 <- function(model_name, model_complexity = "simple", 
                                   time_span = 10.0, num_runs = 10) {
  tryCatch({
    # Check packages
    check_packages()
    
    # Create reaction network based on complexity
    if (model_complexity == "simple") {
      network <- create_simple_network()
    } else if (model_complexity == "multistate") {
      network <- create_multistate_network()
    } else {
      stop("Unknown model complexity")
    }
    
    # Warmup run
    tryCatch({
      out <- ssa(
        initial_state = network$initial_state,
        reactions = network$reactions,
        params = network$params,
        final_time = time_span,
        method = ssa_exact(),
        census_interval = 0.1,
        verbose = FALSE
      )
    }, error = function(e) {
      stop(paste("Warmup failed:", e$message))
    })
    
    # Benchmark runs
    times <- numeric(num_runs)
    for (i in 1:num_runs) {
      start_time <- Sys.time()
      out <- ssa(
        initial_state = network$initial_state,
        reactions = network$reactions,
        params = network$params,
        final_time = time_span,
        method = ssa_exact(),
        census_interval = 0.1,
        verbose = FALSE
      )
      times[i] <- as.numeric(difftime(Sys.time(), start_time, units = "secs")) * 1000
    }
    
    result <- list(
      solver = "GillespieSSA2",
      median_time = median(times),
      min_time = min(times),
      max_time = max(times),
      std_time = sd(times),
      success = TRUE,
      model = model_name,
      complexity = model_complexity,
      num_runs = num_runs
    )
    
    return(result)
    
  }, error = function(e) {
    result <- list(
      solver = "GillespieSSA2",
      median_time = Inf,
      success = FALSE,
      error = as.character(e),
      model = model_name
    )
    return(result)
  })
}

# Benchmark different SSA methods
benchmark_ssa_methods <- function(model_name, model_complexity = "simple", 
                                time_span = 10.0, num_runs = 10) {
  methods <- list(
    "SSA_Exact" = ssa_exact(),
    "SSA_ETL" = ssa_etl(),
    "SSA_BTL" = ssa_btl()
  )
  
  results <- list()
  
  for (method_name in names(methods)) {
    tryCatch({
      # Check packages
      check_packages()
      
      # Create reaction network
      if (model_complexity == "simple") {
        network <- create_simple_network()
      } else if (model_complexity == "multistate") {
        network <- create_multistate_network()
      } else {
        stop("Unknown model complexity")
      }
      
      # Warmup run
      tryCatch({
        out <- ssa(
          initial_state = network$initial_state,
          reactions = network$reactions,
          params = network$params,
          final_time = time_span,
          method = methods[[method_name]],
          census_interval = 0.1,
          verbose = FALSE
        )
      }, error = function(e) {
        stop(paste("Warmup failed:", e$message))
      })
      
      # Benchmark runs
      times <- numeric(num_runs)
      for (i in 1:num_runs) {
        start_time <- Sys.time()
        out <- ssa(
          initial_state = network$initial_state,
          reactions = network$reactions,
          params = network$params,
          final_time = time_span,
          method = methods[[method_name]],
          census_interval = 0.1,
          verbose = FALSE
        )
        times[i] <- as.numeric(difftime(Sys.time(), start_time, units = "secs")) * 1000
      }
      
      result <- list(
        solver = paste0("GillespieSSA2_", method_name),
        median_time = median(times),
        min_time = min(times),
        max_time = max(times),
        std_time = sd(times),
        success = TRUE,
        model = model_name,
        complexity = model_complexity,
        num_runs = num_runs
      )
      
      results[[method_name]] <- result
      
    }, error = function(e) {
      result <- list(
        solver = paste0("GillespieSSA2_", method_name),
        median_time = Inf,
        success = FALSE,
        error = as.character(e),
        model = model_name
      )
      results[[method_name]] <- result
    })
  }
  
  return(results)
}

# Main function for command-line usage
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript r_benchmarks.R <model_name> <complexity> [time_span] [num_runs] [output_file]\\n")
    cat("  model_name: Name of the model to benchmark\\n")
    cat("  complexity: 'simple' or 'multistate'\\n")
    cat("  time_span: Simulation time span (default: 10.0)\\n")
    cat("  num_runs: Number of benchmark runs (default: 10)\\n")
    cat("  output_file: JSON output file path (optional)\\n")
    quit(status = 1)
  }
  
  model_name <- args[1]
  complexity <- args[2]
  time_span <- if (length(args) >= 3) as.numeric(args[3]) else 10.0
  num_runs <- if (length(args) >= 4) as.integer(args[4]) else 10
  output_file <- if (length(args) >= 5) args[5] else NULL
  
  # Run benchmarks
  cat(sprintf("\\nRunning R benchmarks for %s (%s complexity)\\n", model_name, complexity))
  
  # Simple benchmark
  simple_result <- benchmark_gillespie_ssa2(model_name, complexity, time_span, num_runs)
  
  # Method comparison
  method_results <- benchmark_ssa_methods(model_name, complexity, time_span, num_runs)
  
  # Combine results
  all_results <- c(list(simple_result), method_results)
  
  # Print results
  cat("\\nBenchmark results:\\n")
  for (result in all_results) {
    if (result$success) {
      cat(sprintf("  %s: %.2f ms\\n", result$solver, result$median_time))
    } else {
      cat(sprintf("  %s: FAILED - %s\\n", result$solver, result$error))
    }
  }
  
  # Save results if requested
  if (!is.null(output_file)) {
    write_json(all_results, output_file, pretty = TRUE)
    cat(sprintf("\\nResults saved to %s\\n", output_file))
  }
  
  return(all_results)
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}