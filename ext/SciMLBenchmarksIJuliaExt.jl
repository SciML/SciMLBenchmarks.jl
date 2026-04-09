module SciMLBenchmarksIJuliaExt

using SciMLBenchmarks: SciMLBenchmarks, repo_directory, weave_all
using IJulia: IJulia

function SciMLBenchmarks.open_notebooks()
    weave_all((:notebook,))
    path = joinpath(repo_directory, "notebook")
    IJulia.notebook(; dir = path)
    newpath = joinpath(pwd(), "generated_notebooks")
    mv(path, newpath)
    return IJulia.notebook(; dir = newpath)
end

end # module SciMLBenchmarksIJuliaExt
