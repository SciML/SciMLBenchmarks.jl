using Pkg

Pkg.activate(; temp = true)
Pkg.add("GitHub")

using Dates, GitHub

include("update_helpers.jl")
using .WeeklyUpdateHelpers

const REPOSITORY = "SciML/SciMLBenchmarks.jl"
const BASE_BRANCH = "master"
const SESSION_URL = "https://chatgpt.com/codex/tasks/01a03a17-ad6f-7131-82fc-d0fd57ea6512"

gh_token = only(ARGS)
auth = GitHub.authenticate(gh_token)
repo = normpath(joinpath(@__DIR__, "..", ".."))
benchpath = joinpath(repo, "benchmarks")
date = Dates.format(now(UTC), "yyyy-mm-dd")
run_id = get(ENV, "GITHUB_RUN_ID", string(getpid()))
scratch_parent = get(ENV, "RUNNER_TEMP", get(ENV, "TMPDIR", repo))

open_prs, _ = GitHub.pull_requests(
    REPOSITORY;
    auth,
    params = Dict("state" => "open", "per_page" => 100),
)
filter!(pr -> startswith(something(pr.head.label, ""), "SciML:"), open_prs)

failures = String[]
for dir in readdir(benchpath)
    model_dir = joinpath(benchpath, dir)
    isdir(model_dir) || continue
    println("--- Inspecting $dir ---")

    existing_pr = WeeklyUpdateHelpers.find_update_pull_request(open_prs, dir)
    branch = if existing_pr === nothing
        WeeklyUpdateHelpers.update_branch_name(dir, date, run_id)
    else
        existing_pr.head.ref
    end

    mktempdir(scratch_parent) do scratch
        worktree = joinpath(scratch, "repo")
        try
            WeeklyUpdateHelpers.add_update_worktree(
                repo, worktree, branch;
                base_branch = BASE_BRANCH,
                remote_branch = existing_pr !== nothing,
            )

            model_path = joinpath(worktree, "benchmarks", dir)
            Pkg.activate(model_path) do
                Pkg.update()
                Pkg.precompile(; strict = true)
            end

            relative_model_path = joinpath("benchmarks", dir)
            changed = WeeklyUpdateHelpers.has_changes(worktree, relative_model_path)
            if changed
                run(WeeklyUpdateHelpers.git_cmd(worktree, "add", "-A", "--", relative_model_path))
                commit_message = """Updated $dir on $date

                Co-Authored-By: Chris Rackauckas <accounts@chrisrackauckas.com>
                Co-Authored-By: Claude <noreply@anthropic.com>
                Claude-Session: $SESSION_URL"""
                run(WeeklyUpdateHelpers.git_cmd(worktree, "commit", "-m", commit_message))
            end

            remote_head = existing_pr === nothing ? nothing : readchomp(
                    WeeklyUpdateHelpers.git_cmd(worktree, "rev-parse", "origin/$branch")
                )
            local_head = readchomp(WeeklyUpdateHelpers.git_cmd(worktree, "rev-parse", "HEAD"))
            should_push = changed || (remote_head !== nothing && local_head != remote_head)
            if should_push
                if existing_pr === nothing
                    run(WeeklyUpdateHelpers.git_cmd(worktree, "push", "--set-upstream", "origin", branch))
                    body = """This automated dependency refresh should be ignored until reviewed by @ChrisRackauckas.

                    Verification:

                    - `Pkg.update()` completed for `benchmarks/$dir`.
                    - `Pkg.precompile(; strict = true)` completed for the updated environment.

                    🤖 Generated with [Claude Code](https://claude.com/claude-code)
                    $SESSION_URL"""
                    params = Dict(
                        "title" => "Updated $dir for benchmarks",
                        "head" => branch,
                        "base" => BASE_BRANCH,
                        "body" => body,
                        "draft" => true,
                    )
                    created_pr = GitHub.create_pull_request(REPOSITORY; params, auth)
                    @info "Created draft pull request" dir branch url = created_pr.html_url
                else
                    run(WeeklyUpdateHelpers.git_cmd(worktree, "push", "origin", branch))
                    @info "Updated pull request" dir branch number = existing_pr.number
                end
            end
        catch error
            push!(failures, dir)
            @error "Failed to update benchmark environment" dir branch exception = (
                error, catch_backtrace(),
            )
        finally
            isdir(worktree) && run(
                WeeklyUpdateHelpers.git_cmd(repo, "worktree", "remove", "--force", worktree)
            )
        end
    end
end

isempty(failures) || error("Failed to update: $(join(failures, ", "))")
