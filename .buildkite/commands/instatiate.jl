#!/usr/bin/env julia
using Pkg

# Set to clone everything via SSH (if it can't be found in the package server)
# so we authenticate via SSH key to GitHub for our private packages/registries
Pkg.setprotocol!(protocol="ssh")

# Add our registries
Pkg.Registry.add([
    RegistrySpec(name="General"),
    RegistrySpec(url="https://github.com/JuliaRegistries/General")
])

# Instantiate the current project
Pkg.instantiate()

# Precompile in parallel
Pkg.precompile()