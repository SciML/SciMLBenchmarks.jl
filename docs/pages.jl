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

permute!(pages,
    [1, 16, 13, 11, 22, 3, 4, 20, 30, 6, 8, 18, 15, 28, 7, 9, 17, 21, 31, 19, 29, 12, 10, 24, 23, 27, 5, 14, 25, 26, 2, 32]
)

names = [
    "SciMLBenchmarks.jl: Benchmarks for Scientific Machine Learning (SciML) and Equation Solvers",
    "Multi-Language Wrapper Benchmarks",
    "Linear Solvers",
    "Interval Rootfinding",
    "Nonlinear Solvers",
    "Automatic Differentiation",
    "Sparse Automatic Differentiation",
    "Non-Stiff Ordinary Differential Equations (ODEs)",
    "Stiff Ordinary Differential Equations (ODEs)",
    "Biological Differential Equations",
    "Differential-Algebraic Equations (DAEs)",
    "Boundary Value Problems (BVPs)",
    "ModelingToolkit Acausal Modeling / Symbolic-Numeric Benchmarks",
    "Simple Handwritten Partial Differential Equations (PDEs) as ODEs",
    "Complicated Partial Differential Equations (PDEs)",
    "Dynamical ODEs (Hamiltonian and Second Order)",
    "N-Body Problem Benchmarks",
    "Non-Stiff Stochastic Differential Equations (SDEs)",
    "Stiff Stochastic Differential Equations (SDEs)",
    "Non-Stiff Delay Differential Equations (DDEs)",
    "Stiff Delay Differential equations (DDEs)",
    "Jump Process Equations (Gillespie Benchmarks)",
    "Hybrid (Time-Dependent) Jump Processes",
    "Nonlinear Optimization Solver Benchmarks",
    "Optimization Framework Benchmarks",
    "Parameter Estimation and Inverse Problem Benchmarks",
    "Bayesian Inference and Probabilistic Inverse Problem Benchmarks",
    "MethodOfLines.jl Partial Differential Equation (PDE) Formulations",
    "Physics-Informed Neural Network (Neural Network PDE Solver) Cost Function Benchmarks",
    "Physics-Informed Neural Network (Neural Network PDE Solver) Optimizer Benchmarks",
    "SDE Adaptivity Benchmarks",
    "Surrogate Benchmarks"]

for i in 1:length(pages)
    pages[i] = names[i] => pages[i][2]
end
