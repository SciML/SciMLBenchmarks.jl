using SciMLBenchmarks
using Pkg
using Test

struct FooterCapture <: AbstractDisplay
    markdown::Vector{String}
end
Base.display(::FooterCapture, x) = nothing
Base.display(capture::FooterCapture, ::MIME"text/markdown", x) =
    push!(capture.markdown, string(x))

@testset "rendered package footer" begin
    previous_project = Base.active_project()
    mktempdir() do directory
        write(
            joinpath(directory, "Project.toml"), """
            [deps]
            Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
            """
        )
        write(
            joinpath(directory, "JuliaProject.toml"), """
            [deps]
            Base64 = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
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
        try
            for (filename, direct, transitive) in (
                    ("Project.toml", "Markdown", "Base64"),
                    ("JuliaProject.toml", "Base64", "Markdown"),
                )
                project_path = abspath(directory, filename)
                Pkg.activate(project_path)
                @test Base.active_project() == project_path
                capture = FooterCapture(String[])
                output = joinpath(directory, "stdout")
                open(output, "w") do io
                    redirect_stdout(io) do
                        pushdisplay(capture)
                        try
                            SciMLBenchmarks.bench_footer()
                        finally
                            popdisplay(capture)
                        end
                    end
                end
                @test isempty(read(output, String))
                footer = last(capture.markdown)
                sections = split(footer, "And the full manifest:")
                @test length(sections) == 2
                if length(sections) == 2
                    project, manifest = sections
                    @test occursin(direct, project)
                    @test !occursin(transitive, project)
                    @test occursin(direct, manifest)
                    if filename == "Project.toml"
                        @test occursin(transitive, manifest)
                    end
                end
            end
        finally
            if previous_project === nothing
                Pkg.activate()
            else
                Pkg.activate(previous_project)
            end
        end
        @test Base.active_project() == previous_project
    end
end
