pages = Any["SciMLBenchmarks"=>"index.md"]

for folder in readdir(benchmarksdir)
    newpages = Any[]
    if folder[end-2:end] != ".md" && folder != "Testing" && folder != "figures"
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
    [1, 7, 10, 16, 3, 4, 6, 8, 11, 17, 9, 15, 5, 14, 12, 13, 2]
)

names = [
    "SciMLBenchmarks",
    "Multi-Language Wrapper Benchmarks",
    "Non-Stiff Ordinary Differential Equations",
    "Stiff Ordinary Differential Equations",
    "Biological Differential Equations",
    "Differential-Algebraic Equations (DAEs)",
    "Method of Lines Partial Differential Equations (PDEs)",
    "Dynamical ODEs (Hamiltonian and Second Order)",
    "N-Body Problem Benchmarks",
    "Non-Stiff Stochastic Differential Equations",
    "Stiff Stochastic Differential Equations",
    "Non-Stiff Delay Differential Equations",
    "Stiff Delay Differential equations",
    "Jump Process Equations (Gillespie Benchmarks)",
    "Parameter Estimation and Inverse Problem Benchmarks",
    "Physics-Informed Neural Network (Neural Network PDE Solver) Cost Function Benchmarks",
    "Physics-Informed Neural Network (Neural Network PDE Solver) Optimizer Benchmarks",
    "SDE Adaptivity Benchmarks"]

for i in 1:length(pages)
    pages[i] = names[i] => pages[i][2]
end
