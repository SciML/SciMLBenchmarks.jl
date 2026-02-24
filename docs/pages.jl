# This file assumes `dir` is the directory for the package! dir = @__DIR__() * "/.."

dir = @__DIR__() * "/.."

cp(joinpath(dir, "markdown"), joinpath(dir, "docs", "src"), force=true)
cp(joinpath(dir, "docs", "extrasrc", "assets"), joinpath(dir, "docs", "src", "assets"), force=true)
cp(joinpath(dir, "README.md"), joinpath(dir, "docs", "src", "index.md"), force=true)
benchmarksdir = joinpath(dir, "docs", "src")

@show readdir(benchmarksdir)

pages = Any["SciMLBenchmarks.jl: Benchmarks for Scientific Machine Learning (SciML), Equation Solvers, and AI for Science"=>"index.md"]

for folder in readdir(benchmarksdir)
    newpages = Any[]
    if folder[end-2:end] != ".md" && folder != "Testing" && folder != "figures" && folder != "assets"
        for file in filter(x -> x[end-2:end] == ".md", readdir(
            joinpath(benchmarksdir, folder)))
            try
                filecontents = readlines(joinpath(benchmarksdir, folder, file))
                title = filecontents[3][9:end-1]

                # Cut out the first 5 lines from the file to remove the Weave header stuff
                open(joinpath(benchmarksdir, folder, file), "w") do output
                    println(output, "# $title")
                    for line in Iterators.drop(filecontents, 4)
                        println(output, line)
                    end
                end
                push!(newpages, title => joinpath(folder, file))
            catch e
                @show folder, file, e
            end
        end
        push!(pages, folder => newpages)
    end
end


# The result is in alphabetical order, change to the wanted order

section_titles = [
    "MultiLanguage" => "Multi-Language Wrapper Benchmarks",
    "LinearSolve" => "Linear Solvers",
    "IntervalNonlinearProblem" => "Interval Rootfinding",
    "NonlinearProblem" => "Nonlinear Solvers",
    "AutomaticDifferentiation" => "Automatic Differentiation",
    "AutomaticDifferentiationSparse" => "Sparse Automatic Differentiation",
    "NonStiffODE" => "Non-Stiff Ordinary Differential Equations (ODEs)",
    "StiffODE" => "Stiff Ordinary Differential Equations (ODEs)",
    "Bio" => "Biological Differential Equations",
    "AstroChem" => "Astrochemistry Differential Equations",
    "DAE" => "Differential-Algebraic Equations (DAEs)",
    "NonStiffBVP" => "Non-Stiff Boundary Value Problems (BVPs)",
    "StiffBVP" => "Stiff Boundary Value Problems (BVPs)",
    "ModelingToolkit" => "ModelingToolkit Acausal Modeling / Symbolic-Numeric Benchmarks",
    "SimpleHandwrittenPDE" => "Simple Handwritten Partial Differential Equations (PDEs) as ODEs",
    "ComplicatedPDE" => "Complicated Partial Differential Equations (PDEs)",
    "DynamicalODE" => "Dynamical ODEs (Hamiltonian and Second Order)",
    "NBodySimulator" => "N-Body Problem Benchmarks",
    "NonStiffSDE" => "Non-Stiff Stochastic Differential Equations (SDEs)",
    "StiffSDE" => "Stiff Stochastic Differential Equations (SDEs)",
    "NonStiffDDE" => "Non-Stiff Delay Differential Equations (DDEs)",
    "StiffDDE" => "Stiff Delay Differential equations (DDEs)",
    "Jumps" => "Jump Process Equations (Gillespie Benchmarks)",
    "HybridJumps" => "Hybrid (Time-Dependent) Jump Processes",
    "Optimization" => "Nonlinear Optimization Solver Benchmarks",
    "OptimizationCUTEst" => "CUTEst Optimization Solver Benchmarks",
    "GlobalOptimization" => "Global Optimization Benchmarks",
    "OptimizationFrameworks" => "Optimization Framework Benchmarks",
    "ParameterEstimation" => "Parameter Estimation and Inverse Problem Benchmarks",
    "BayesianInference" => "Bayesian Inference and Probabilistic Inverse Problem Benchmarks",
    "MethodOfLinesPDE" => "MethodOfLines.jl Partial Differential Equation (PDE) Formulations",
    "PINNErrorsVsTime" => "Physics-Informed Neural Network (Neural Network PDE Solver) Cost Function Benchmarks",
    "PINNOptimizers" => "Physics-Informed Neural Network (Neural Network PDE Solver) Optimizer Benchmarks",
    "AdaptiveSDE" => "SDE Adaptivity Benchmarks",
    "Surrogates" => "Surrogate Benchmarks",
    "Symbolics" => "Symbolic Manipulation Benchmarks"
]

renamed_index = "SciMLBenchmarks.jl: Benchmarks for Scientific Machine Learning (SciML) and Equation Solvers" =>
                pages[1][2]
remaining_pages = Dict{String,Any}(pages[2:end])
ordered_pages = Any[renamed_index]

for (folder, title) in section_titles
    if haskey(remaining_pages, folder)
        push!(ordered_pages, title => remaining_pages[folder])
        delete!(remaining_pages, folder)
    end
end

# Keep docs generation robust when new benchmark folders are added.
for folder in sort!(collect(keys(remaining_pages)))
    push!(ordered_pages, folder => remaining_pages[folder])
end

pages = ordered_pages
