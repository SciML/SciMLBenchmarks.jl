module WeeklyUpdateHelpers

using Dates

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

"""
    select_update_directory(dirs, requested, date) -> String

Pick the single benchmark directory to refresh on this run.

Refreshing every directory in one run opened a pull request per directory (~40 at once), which
is more review load than the weekly cadence absorbs. A non-empty `requested` (the
`workflow_dispatch` input) selects that directory explicitly; otherwise the choice rotates by
week, so each directory comes up on a predictable cycle rather than at random.
"""
function select_update_directory(dirs, requested::AbstractString, date::Date)
    isempty(dirs) && error("No benchmark directories found")
    req = strip(requested)
    if !isempty(req)
        req in dirs ||
            error("Requested directory $req not found. Available: $(join(dirs, ", "))")
        return String(req)
    end
    # Count cron weeks, not calendar weeks. 1970-01-03 was a Saturday, and the schedule is
    # `0 6 * * 6`, so anchoring there makes Saturday-through-Friday one bucket: a manual
    # re-run later in the same week picks the same directory as that week's scheduled run.
    epoch_saturday = Date(1970, 1, 3)
    weeks = Dates.value(date - epoch_saturday) ÷ 7
    return String(dirs[mod1(weeks + 1, length(dirs))])
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
