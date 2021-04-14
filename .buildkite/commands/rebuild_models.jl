#!/usr/bin/env julia

using Pkg
repo_dir = dirname(dirname(@__DIR__))

model = ARGS[1]
folder, file = split(model, "/")

project_abspath = joinpath(repo_dir, "library", model)
if !isfile(project_abspath) || !occursin(".jmd", file)
    error("Invalid file $file")
end

@info("Rebuilding $folder/$file")
weave_file(folder, file)
