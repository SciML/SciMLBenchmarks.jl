using Franklin
using Weave


"""
    weaveall()

Weave all lecture notes in the `_weave` directory. Run from site root.
"""
function weaveall()
    for (root, _, files) in walkdir("_weave")
        for file in files
            if endswith(file, "jmd")
                @info "Weaving Document: $(joinpath(root, file))"
                weave(joinpath(root, file); out_path=:doc)
            end
        end
    end
end


"""
    cleanall()

Cleanup all Weave generated subdirectories. Run from site root.
"""
function cleanall()
    for (root, dirs, _) in walkdir("_weave")
        for dir in dirs
            if startswith(dir, "jl_")
                @info "Removing Directory: $(joinpath(root, dir))"
                rm(joinpath(root, dir); recursive=true, force=true)
            end
        end
    end
end
