module SciMLBenchmarks

using Weave, Pkg, IJulia, InteractiveUtils, Markdown

repo_directory = joinpath(@__DIR__,"..")

macro subprocess(ex, wait=true)
  quote
    local project = Pkg.project().path
    local ex_str = $(esc(sprint(Base.show_unquoted, ex)))
    run(`$(Base.julia_cmd()) --project=$(project) -e "$(ex_str)"`; wait=$(wait))
  end
end

function weave_file(folder,file,build_list=(:script,:github))
  target = joinpath(folder, file)
  @info("Weaving $(target)")

  if isfile(joinpath(folder, "Project.toml")) && build_list != (:notebook,)
    @info("Instantiating", folder)
    Pkg.activate(folder)
    Pkg.instantiate()
    Pkg.build()
  end

  args = Dict{Symbol,String}(:folder=>folder,:file=>file)
  if :script ∈ build_list
    println("Building Script")
    dir = joinpath(repo_directory,"script",basename(folder))
    mkpath(dir)
    tangle(target; out_path=dir)
  end
  if :html ∈ build_list
    println("Building HTML")
    dir = joinpath(repo_directory,"html",basename(folder))
    mkpath(dir)
    weave(target,doctype = "md2html",out_path=dir,args=args,fig_ext=".svg")
  end
  if :pdf ∈ build_list
    println("Building PDF")
    dir = joinpath(repo_directory,"pdf",basename(folder))
    mkpath(dir)
    try
      weave(target,doctype="md2pdf",out_path=dir,args=args)
    catch ex
      @warn "PDF generation failed" exception=(ex, catch_backtrace())
    end
  end
  if :github ∈ build_list
    println("Building Github Markdown")
    dir = joinpath(repo_directory,"markdown",basename(folder))
    mkpath(dir)
    weave(target,doctype = "github",out_path=dir,args=args)
  end
  if :notebook ∈ build_list
    println("Building Notebook")
    dir = joinpath(repo_directory,"notebook",basename(folder))
    mkpath(dir)
    Weave.convert_doc(target,joinpath(dir,file[1:end-4]*".ipynb"))
  end
end

function weave_all(build_list=(:script,:github))
  for folder in readdir(joinpath(repo_directory,"benchmarks"))
    folder == "test.jmd" && continue
    weave_folder(joinpath(repo_directory,"benchmarks",folder),build_list)
  end
end

function weave_folder(folder, build_list=(:script,:github))
  for file in readdir(folder)
    # Skip non-`.jmd` files
    if !endswith(file, ".jmd")
      continue
    end

    try
      @eval @subprocess begin
        using SciMLBenchmarks
        SciMLBenchmarks.weave_file($folder, $file, $build_list)
      end
    catch e
      @show folder, file
      @error(e)
    end
  end
end

function bench_footer(folder=nothing, file=nothing)
    display(md"""
    ## Appendix

    These benchmarks are a part of the SciMLBenchmarks.jl repository, found at: <https://github.com/SciML/SciMLBenchmarks.jl>.
    For more information on high-performance scientific machine learning, check out the SciML Open Source Software Organization <https://sciml.ai>.

    """)
    if folder !== nothing && file !== nothing
        display(Markdown.parse("""
        To locally run this benchmark, do the following commands:
        ```
        using SciMLBenchmarks
        SciMLBenchmarks.weave_file("$folder","$file")
        ```
        """))
    end
    display(md"Computer Information:")
    vinfo = sprint(InteractiveUtils.versioninfo)
    display(Markdown.parse("""
    ```
    $(vinfo)
    ```
    """))

    display(md"""
    Package Information:
    """)

    proj = sprint(io -> Pkg.status(io=io))
    mani = sprint(io -> Pkg.status(io=io, mode = Pkg.PKGMODE_MANIFEST))

    md = """
    ```
    $(chomp(proj))
    ```

    And the full manifest:

    ```
    $(chomp(mani))
    ```
    """
    display(Markdown.parse(md))
end

function open_notebooks()
    Base.eval(Main, Meta.parse("import IJulia"))
    weave_all((:notebook,))
    path = joinpath(repo_directory,"notebook")
    IJulia.notebook(;dir=path)
    newpath = joinpath(pwd(),"generated_notebooks")
    mv(path, newpath)
    IJulia.notebook(;dir=newpath)
  end

end # module SciMLBenchmarks
