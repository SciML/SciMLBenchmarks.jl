#!/usr/bin/env julia

using Pkg

# Set to clone everything via SSH (if it can't be found in the package server)
# so we authenticate via SSH key to GitHub for our private packages/registries
Pkg.setprotocol!(protocol="ssh")

# Add our registries
Pkg.Registry.add([
    RegistrySpec(name="General"),
    RegistrySpec(url="https://github.com/JuliaComputing/JuliaSimRegistry")
])

Pkg.instantiate()
using SciMLBenchmarks

repo_dir = dirname(dirname(@__DIR__))

model = ARGS[1]
split_path = split(model, "/")
file = split_path[end]
folder = split_path[end-1]

project_abspath = joinpath(repo_dir, model)

@info("debugging!", model, file, folder, project_abspath, isfile(project_abspath), occursin(".jmd", file))

if !isfile(project_abspath)
    error("Invalid file $file")
end

if occursin(".jmd", file)
    @info("Rebuilding $folder/$file")
    SciMLBenchmarks.weave_file(folder, file)
elseif occursin(".toml", file)
    @info("Rebuilding $folder")
    SciMLBenchmarks.weave_folder(folder)
else
    error("Invalid file $file")
end
