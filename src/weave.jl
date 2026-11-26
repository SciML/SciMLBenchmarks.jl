repo_directory = joinpath(@__DIR__, "..")

function _notebook_path(folder, file)
    return joinpath(repo_directory, "notebook", basename(folder), file[1:(end - 4)] * ".ipynb")
end


function _literal_include_path(argument)
    argument isa String && return argument
    Meta.isexpr(argument, :call) && length(argument.args) > 2 || return nothing
    callee, directory, parts... = argument.args
    callee === :joinpath && Meta.isexpr(directory, :macrocall) &&
        directory.args[1] === Symbol("@__DIR__") && all(part -> part isa String, parts) ||
        return nothing
    return joinpath(parts...)
end

function _collect_include_paths!(paths, ex)
    ex isa Expr && ex.head !== :quote || return paths
    if Meta.isexpr(ex, :call, 2) && ex.args[1] === :include
        path = _literal_include_path(ex.args[2])
        path === nothing || push!(paths, path)
    end
    foreach(argument -> _collect_include_paths!(paths, argument), ex.args)
    return paths
end

function _script_include_paths(script)
    source = read(script, String)
    paths = String[]
    position = 1
    while position <= ncodeunits(source)
        expression, position = Meta.parse(source, position; raise = false)
        _collect_include_paths!(paths, expression)
    end
    return unique!(paths)
end

function _directory_resolves_inside(directory, root)
    while !islink(directory) && !ispath(directory)
        directory = dirname(directory)
    end
    ispath(directory) || return false
    return first(splitpath(relpath(realpath(directory), realpath(root)))) != ".."
end

function _stage_script_includes(folder, script, script_root)
    for path in _script_include_paths(script)
        source = normpath(folder, path)
        isfile(source) || continue
        destination = normpath(dirname(script), path)
        if first(splitpath(relpath(destination, script_root))) in (".", "..") ||
                !_directory_resolves_inside(dirname(destination), script_root)
            @warn "Not staging an include that resolves outside the script tree" script path destination
            continue
        end
        mkpath(dirname(destination))
        cp(source, destination; force = true)
    end
    return nothing
end

function _render_file(folder, file, build_list)
    Weave.set_chunk_defaults!(:error => false)
    target = joinpath(folder, file)
    args = Dict{Symbol, String}(:folder => folder, :file => file)
    if :script ∈ build_list
        println("Building Script")
        script_root = joinpath(repo_directory, "script")
        dir = joinpath(script_root, basename(folder))
        mkpath(dir)
        tangle(target; out_path = dir)
        script = joinpath(dir, string(first(splitext(file)), ".jl"))
        _stage_script_includes(folder, script, script_root)
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
