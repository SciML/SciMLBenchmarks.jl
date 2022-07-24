pages = Any["SciMLBenchmarks"=>"index.md"]

for folder in readdir(benchmarksdir)
    newpages = Any[]
    if folder[end-2:end] != ".md" && folder != "Testing" && folder != "figures"
        for file in filter(x -> x[end-2:end] == ".md", readdir(
            joinpath(benchmarksdir, folder)))
            try
                filecontents = readlines(joinpath(benchmarksdir, folder, file))
                title = filecontents[3][9:end-1]

                # Cut out the first 5 lines from the file to remove the Weave header stuff
                open(joinpath(benchmarksdir, folder, file), "w") do output
                    for line in Iterators.drop(filecontents, 4)
                        println(output, line)
                    end
                end
                push!(newpages, title => joinpath(folder, file))
            catch e
                @show folder, file, e
            end
        end
        push!(pages, Any[folder=>newpages])
    end
end
