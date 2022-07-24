using Documenter, SciMLBenchmarksOutput

dir = pkgdir(SciMLBenchmarksOutput)
cp(joinpath(dir, "README.md"), joinpath(dir, "docs", "src", "index.md"), force=true)
cp(joinpath(dir, "markdown"), joinpath(dir, "docs", "src"), force=true)

pages = Any["SciMLBenchmarks"=>"index.md"]

for folder in readdir(joinpath(dir, "docs", "src"))
    @show folder[end-2:end]
    if folder[end-2:end] != ".md" && folder != "Testing" && folder != "figures"
        newpages = [joinpath(folder, file) for file in
                    filter(x -> x[end-2:end] == ".md", readdir(
            joinpath(dir, "docs", "src", folder)))]
        push!(pages, [folder => newpages])
    end
end

makedocs(
    sitename="The SciML Benchmarks",
    authors="Chris Rackauckas",
    modules=[SciMLBenchmarksOutput],
    clean=true, doctest=false,
    format=Documenter.HTML(#analytics = "UA-90474609-3",
        assets=["assets/favicon.ico"],
        canonical="https://benchmarks.sciml.ai/stable/"),
    pages=[
        "Home" => "index.md",
        "Non-Stiff Ordinary Differential Equation (ODE) Solver Benchmarks" => Any[
            "markdown/MultiLanguage/wrapper_packages.md"
        ],
    ]
)

deploydocs(;
    repo="github.com/SciML/SciMLBenchmarksOutput",
    devbranch="main",
    branch="main"
)
