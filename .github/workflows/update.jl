using Pkg
Pkg.add(["Git", "GitHub", "Dates"])
using Git, GitHub, Dates

gh_token = ARGS[1]
myauth = GitHub.authenticate(gh_token)

(@isdefined myauth) ? @info("Authentication token is found...") : @info("Coudn't find the authentication token")

const git = Git.git()
date = Dates.format(now(), "yyyy-mm-dd")
benchpath = joinpath(@__DIR__, "..", "..", "benchmarks")

# Get all the open PRs and their number
gh_prs = GitHub.pull_requests("SciML/SciMLBenchmarks.jl"; auth=myauth)
prs = Dict{String, Int64}()
for i in 1:length(gh_prs[1])
    prs[gh_prs[1][i].head.ref] = gh_prs[1][i].number
end

# Get all the branches from the repo
gh_branches = GitHub.branches("SciML/SciMLBenchmarks.jl"; auth=myauth)
branches = [gh_branches[1][i].name for i in 1:length(gh_branches[1])]

@info("PRs and branches", prs, branches)

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
        if dir ∉ branches
            run(`git checkout -b $(dir) master`)
            run(`$git add -A . :!$(manpath)`)
            run(`$git commit -m "Updated $(dir) on $(date)"`)
            run(`$git push --set-upstream origin $(dir)`)
        else
            run(`$git fetch origin`)
            run(`$git checkout $(dir)`)
            run(`$git pull -Xours`)
            run(`$git add -A . :!$(manpath)`)
            run(`$git commit -m "Updated $(dir) on $(date)"`)
            run(`$git push`)
        end
        if dir ∉ keys(prs)
            params = Dict(
                "title" => "Updated $(dir) for benchmarks",
                "head"  => "$(dir)",
                "base"  => "master"
            )
            @info("Creating a pull request from head: ", dir) 
            GitHub.create_pull_request("SciML/SciMLBenchmarks.jl"; params=params, auth=myauth)
        else
            @info("Updating the pull request numbered: ", prs[dir])
            GitHub.update_pull_request("SciML/SciMLBenchmarks.jl", prs[dir]; auth=myauth)
        end
    end
end
