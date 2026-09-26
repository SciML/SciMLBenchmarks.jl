using Dates
using Test

include("update_helpers.jl")
using .WeeklyUpdateHelpers

function git(repo, args...)
    return run(WeeklyUpdateHelpers.git_cmd(repo, args...))
end

@testset "weekly update helpers" begin
    @test WeeklyUpdateHelpers.update_branch_name("Example", "2026-08-25", "123") ==
        "update/Example/2026-08-25-123"

    prs = [
        (head = (ref = "unrelated",), number = 1),
        (head = (ref = "update/Example/2026-08-18-99",), number = 2),
    ]
    @test WeeklyUpdateHelpers.find_update_pull_request(prs, "Example").number == 2
    @test isnothing(WeeklyUpdateHelpers.find_update_pull_request(prs, "Other"))
    legacy_prs = [(head = (ref = "Example",), number = 3)]
    @test WeeklyUpdateHelpers.find_update_pull_request(legacy_prs, "Example").number == 3

    # One directory per run: explicit request wins, otherwise rotate by week.
    dirs = ["Alpha", "Beta", "Gamma"]
    @test WeeklyUpdateHelpers.select_update_directory(dirs, "Beta", Date(2026, 8, 29)) == "Beta"
    @test WeeklyUpdateHelpers.select_update_directory(dirs, "  Beta  ", Date(2026, 8, 29)) == "Beta"
    @test_throws ErrorException WeeklyUpdateHelpers.select_update_directory(
        dirs, "Missing", Date(2026, 8, 29)
    )
    @test_throws ErrorException WeeklyUpdateHelpers.select_update_directory(
        String[], "", Date(2026, 8, 29)
    )

    # The rotation advances exactly once per week and stays put within a week.
    base = Date(2026, 8, 29)
    picks = [WeeklyUpdateHelpers.select_update_directory(dirs, "", base + Day(7i)) for i in 0:5]
    @test all(p -> p in dirs, picks)
    @test picks[1:3] != picks[2:4]                       # consecutive weeks differ
    @test picks[1] == picks[4] && picks[2] == picks[5]   # wraps with the directory count
    # A manual re-run any day before the next Saturday picks that week's directory.
    @test all(
        WeeklyUpdateHelpers.select_update_directory(dirs, "", base + Day(d)) == picks[1]
            for d in 0:6
    )
    @test WeeklyUpdateHelpers.select_update_directory(dirs, "", base + Day(7)) != picks[1]

    # Every directory is reached over a full cycle, so nothing is starved.
    cycle = Set(
        WeeklyUpdateHelpers.select_update_directory(dirs, "", base + Day(7i))
            for i in 0:(length(dirs) - 1)
    )
    @test cycle == Set(dirs)

    scratch_parent = get(ENV, "TMPDIR", tempdir())
    mktempdir(scratch_parent) do scratch
        repo = joinpath(scratch, "repo")
        remote = joinpath(scratch, "remote.git")
        worktree = joinpath(scratch, "worktree")
        mkpath(joinpath(repo, "benchmarks", "Example"))
        git(repo, "init", "-b", "master")
        git(repo, "config", "user.name", "Weekly Update Test")
        git(repo, "config", "user.email", "weekly-update@example.com")
        write(joinpath(repo, "benchmarks", "Example", "Manifest.toml"), "version = 1\n")
        git(repo, "add", ".")
        git(repo, "commit", "-m", "Initial state")
        git(scratch, "init", "--bare", remote)
        git(repo, "remote", "add", "origin", remote)
        git(repo, "push", "-u", "origin", "master")

        try
            WeeklyUpdateHelpers.add_update_worktree(
                repo, worktree, "update/Example/2026-08-25-123";
                base_branch = "master"
            )
            manifest = joinpath(worktree, "benchmarks", "Example", "Manifest.toml")
            @test !WeeklyUpdateHelpers.has_changes(worktree, "benchmarks/Example")
            write(manifest, "version = 2\n")
            @test WeeklyUpdateHelpers.has_changes(worktree, "benchmarks/Example")
            @test !WeeklyUpdateHelpers.has_changes(repo, "benchmarks/Example")
        finally
            isdir(worktree) && git(repo, "worktree", "remove", "--force", worktree)
        end

        git(repo, "switch", "-c", "Example")
        write(joinpath(repo, "branch.txt"), "existing update\n")
        git(repo, "add", "branch.txt")
        git(repo, "commit", "-m", "Existing update")
        git(repo, "push", "-u", "origin", "Example")
        git(repo, "switch", "master")
        write(joinpath(repo, "base.txt"), "new base\n")
        git(repo, "add", "base.txt")
        git(repo, "commit", "-m", "Advance base")

        existing_worktree = joinpath(scratch, "existing-worktree")
        try
            WeeklyUpdateHelpers.add_update_worktree(
                repo, existing_worktree, "Example";
                base_branch = "master", remote_branch = true
            )
            @test isfile(joinpath(existing_worktree, "branch.txt"))
            @test isfile(joinpath(existing_worktree, "base.txt"))
            @test !WeeklyUpdateHelpers.has_changes(repo, ".")
        finally
            isdir(existing_worktree) &&
                git(repo, "worktree", "remove", "--force", existing_worktree)
        end
    end
end
