case = ARGS[1]
if case == "restricted"
    push!(LOAD_PATH, ENV["FOOTER_PACKAGE_PROJECT"])
    push!(LOAD_PATH, "@stdlib")
end
using SciMLBenchmarks
using Pkg
using Test

if case == "no-active"
    empty!(LOAD_PATH)
end

struct ProcessFooterCapture <: AbstractDisplay
    markdown::Vector{String}
end
Base.display(::ProcessFooterCapture, x) = nothing
Base.display(capture::ProcessFooterCapture, ::MIME"text/markdown", x) =
    push!(capture.markdown, string(x))

@testset "footer $case" begin
    capture = ProcessFooterCapture(String[])
    pushdisplay(capture)
    try
        if case == "no-active"
            @test Base.active_project() === nothing
            native_error = try
                Pkg.status()
                nothing
            catch e
                e
            end
            footer_error = try
                SciMLBenchmarks.bench_footer()
                nothing
            catch e
                e
            end
            @test native_error isa Exception
            @test footer_error isa Exception
            if native_error isa Exception && footer_error isa Exception
                @test occursin("no active project", sprint(showerror, native_error))
                @test typeof(footer_error) === typeof(native_error)
                @test sprint(showerror, footer_error) == sprint(showerror, native_error)
            end
            @test all(text -> !occursin("Markdown", text), capture.markdown)
        else
            SciMLBenchmarks.bench_footer()
            sections = split(last(capture.markdown), "And the full manifest:")
            @test length(sections) == 2
            if length(sections) == 2
                @test occursin("Markdown", sections[1])
                @test !occursin("Base64", sections[1])
                @test occursin("Markdown", sections[2])
                @test occursin("Base64", sections[2])
            end
        end
    finally
        popdisplay(capture)
    end
end
