using Documenter, SciMLBenchmarksOutput

dir = @__DIR__() * "/.."

@show dir
@show readdir(dir)

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
