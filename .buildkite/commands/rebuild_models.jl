#!/usr/bin/env julia

using Pkg
Pkg.instantiate()
using SciMLBenchmarks

repo_dir = dirname(dirname(@__DIR__))

model = ARGS[1]
split_path = split(model, "/")
file = split_path[end]
folder = split_path[end-1]

project_abspath = joinpath(repo_dir, model)

@info("debugging!", model, file, folder, project_abspath, isfile(project_abspath), occursin(".jmd", file))

if !isfile(project_abspath) || !occursin(".jmd", file)
    error("Invalid file $file")
end

@info("Rebuilding $folder/$file")
SciMLBenchmarks.weave_file(folder, file)
