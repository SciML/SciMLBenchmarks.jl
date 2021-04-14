using Pkg
import Git

const git = Git.git()

benchpath = joinpath(@__DIR__, "..", "..", "benchmarks")

for dir in readdir(benchpath)
    model_dir = joinpath(benchpath, dir)
    isdir(model_dir) || continue
    println("--- Inspecting $dir ---")
    cd(model_dir)
    Pkg.activate(".") do
        Pkg.update()
    end
    manpath = (joinpath(benchpath, "Manifest.toml"))
    if length(readlines(`git status . --porcelain`)) > 1
        run(`git checkout -B $dir`)
        run(`$git add -A . :!$(manpath)`)
        run(`$git commit -m "Updated $(dir)"`)
        run(`$git push --set-upstream origin $(dir)`)
    end
end