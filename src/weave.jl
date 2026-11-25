repo_directory = joinpath(@__DIR__, "..")

function _notebook_path(folder, file)
    return joinpath(repo_directory, "notebook", basename(folder), file[1:(end - 4)] * ".ipynb")
end

function _render_file(folder, file, build_list)
    Weave.set_chunk_defaults!(:error => false)
    target = joinpath(folder, file)
    args = Dict{Symbol, String}(:folder => folder, :file => file)
    if :script ∈ build_list
        println("Building Script")
        dir = joinpath(repo_directory, "script", basename(folder))
        mkpath(dir)
        tangle(target; out_path = dir)
    end
    if :html ∈ build_list
        println("Building HTML")
        dir = joinpath(repo_directory, "html", basename(folder))
        mkpath(dir)
        weave(target, doctype = "md2html", out_path = dir, args = args, fig_ext = ".svg")
    end
    if :pdf ∈ build_list
        println("Building PDF")
        dir = joinpath(repo_directory, "pdf", basename(folder))
        mkpath(dir)
        try
            weave(target, doctype = "md2pdf", out_path = dir, args = args)
        catch ex
            @warn "PDF generation failed" exception = (ex, catch_backtrace())
        end
    end
    if :github ∈ build_list
        println("Building Github Markdown")
        dir = joinpath(repo_directory, "markdown", basename(folder))
        mkpath(dir)
        weave(target, doctype = "github", out_path = dir, args = args)
    end
    return if :notebook ∈ build_list
        println("Building Notebook")
        output = _notebook_path(folder, file)
        mkpath(dirname(output))
        Weave.convert_doc(target, output)
    end
end
