using Documenter, SciMLBenchmarksOutput

dir = pkgdir(SciMLBenchmarksOutput)
cp(joinpath(dir, "markdown"), joinpath(dir, "docs", "src"), force=true)
cp(joinpath(dir, "README.md"), joinpath(dir, "docs", "src", "index.md"), force=true)
benchmarksdir = joinpath(dir, "docs", "src")

include("pages.jl")

makedocs(
    sitename="The SciML Benchmarks",
    authors="Chris Rackauckas",
    modules=[SciMLBenchmarksOutput],
    clean=true, doctest=false,
    format=Documenter.HTML(#analytics = "UA-90474609-3",
        assets=["assets/favicon.ico"],
        canonical="https://benchmarks.sciml.ai/stable/"),
    pages=pages
)

deploydocs(;
    repo="github.com/SciML/SciMLBenchmarksOutput",
    devbranch="main",
    branch="main"
)
