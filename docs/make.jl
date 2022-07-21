using Documenter, SciMLBenchmarksOutput

pkgdir(SciMLBenchmarksOutput)
cp(joinpath(mod, "README.md"), joinpath(dir, "docs", "src", "index.md"), force=true)
cp(joinpath(mod, "markdown"), joinpath(dir, "docs", "src"), force=true)

makedocs(
    sitename="The SciML Benchmarks: Differential Equations, Inverse Problems, Physics-Informed ML, Science-Guided AI",
    authors="Chris Rackauckas",
    modules=[SciMLBenchmarksOutput],
    clean=true,doctest=false,
    format = Documenter.HTML(#analytics = "UA-90474609-3",
                             assets = ["assets/favicon.ico"],
                             canonical="https://benchmarks.sciml.ai/stable/"),
    pages=[
        "Home" => "README.md",
        "Non-Stiff Ordinary Differential Equation (ODE) Solver Benchmarks" => Any[
            "markdown/MultiLanguage/wrapper_packages.md"
        ],
    ]
)

deploydocs(;
    repo="github.com/SciML/SciMLBenchmarksOutput",
    devbranch="main",
    branch = "main"
)
