@static if VERSION < v"1.4"
    @testset "legacy footer process context" begin
        package_project = Base.active_project()
        mktempdir() do directory
            project = joinpath(directory, "Project.toml")
            write(
                project, """
                [deps]
                Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
                """
            )
            write(
                joinpath(directory, "Manifest.toml"), """
                [[Base64]]
                uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

                [[Markdown]]
                deps = ["Base64"]
                uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
                """
            )
            julia = joinpath(Sys.BINDIR, Sys.iswindows() ? "julia.exe" : "julia")
            script = joinpath(@__DIR__, "footer_process_child.jl")
            separator = Sys.iswindows() ? ';' : ':'
            for (case, load_path, active) in (
                    ("no-active", join((directory, dirname(package_project), "@stdlib"), separator), package_project),
                    ("restricted", "@", project),
                )
                environment = copy(ENV)
                environment["JULIA_LOAD_PATH"] = load_path
                environment["FOOTER_PACKAGE_PROJECT"] = dirname(package_project)
                delete!(environment, "JULIA_PROJECT")
                cmd = case == "no-active" ? `$julia --startup-file=no $script $case` :
                    `$julia --startup-file=no --project=$active $script $case`
                cmd = setenv(cmd, environment)
                output = joinpath(directory, "$case.log")
                ok = open(output, "w") do io
                    success(pipeline(cmd; stdout = io, stderr = io))
                end
                ok || @error "Footer process regression failed" case output = read(output, String)
                @test ok
            end
        end
    end
end
