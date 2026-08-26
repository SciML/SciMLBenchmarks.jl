module WeeklyUpdateHelpers

function git_cmd(repo::AbstractString, args::AbstractString...)
    return Cmd(Cmd(["git", args...]); dir = repo)
end

function has_changes(repo::AbstractString, path::AbstractString)
    return !isempty(readlines(git_cmd(repo, "status", "--porcelain", "--", path)))
end

function update_branch_name(dir::AbstractString, date::AbstractString, run_id::AbstractString)
    return "update/$dir/$date-$run_id"
end

function find_update_pull_request(prs, dir::AbstractString)
    legacy = findfirst(pr -> pr.head.ref == dir, prs)
    legacy !== nothing && return prs[legacy]

    prefix = "update/$dir/"
    current = findfirst(pr -> startswith(pr.head.ref, prefix), prs)
    return current === nothing ? nothing : prs[current]
end

function add_update_worktree(
        repo::AbstractString, worktree::AbstractString, branch::AbstractString;
        base_branch::AbstractString = "master", remote_branch::Bool = false
    )
    if remote_branch
        run(git_cmd(repo, "fetch", "origin", branch))
        run(git_cmd(repo, "worktree", "add", "-B", branch, worktree, "origin/$branch"))
        run(git_cmd(worktree, "merge", "--no-edit", base_branch))
    else
        run(git_cmd(repo, "worktree", "add", "-b", branch, worktree, base_branch))
    end
    return nothing
end

end
