using Documenter

makedocs(
    sitename="The SciML Benchmarks: Differential Equations, Inverse Problems, Physics-Informed ML, Science-Guided AI",
    authors="Chris Rackauckas",
    modules=[],
    clean=true,doctest=false,
    format = Documenter.HTML(#analytics = "UA-90474609-3",
                             assets = ["assets/favicon.ico"],
                             canonical="https://benchmarks.sciml.ai/stable/"),
    pages=[
        "Home" => "../README.md",
        "Non-Stiff Ordinary Differential Equation (ODE) Solver Benchmarks" => Any[
            "../NonStiffODE/linear_wpd.md"
        ],
    ]
)

deploydocs(;
    repo="github.com/SciML/SciMLBenchmarksOutput",
    devbranch="main",
)
