module DiffEqBenchmarks

using Weave, Pkg

repo_directory = joinpath(@__DIR__,"..")

function dir_exists(dir)
  exists = true
  try
    readdir(dir)
  catch err
    exists = false
  end
  exists
end

function weave_file(folder,file,build_list=(:script,:html,:pdf))
  println("File: $file")
  tmp = joinpath(repo_directory,"benchmarks",folder,file)
  args = Dict{Symbol,String}(:folder=>folder,:file=>file)
  if :script ∈ build_list
    println("Building Script")
    dir = joinpath(repo_directory,"script",folder)
    dir_exists(dir) || mkdir(dir)
    tangle(tmp;out_path=dir)
  end
  if :html ∈ build_list
    println("Building HTML")
    dir = joinpath(repo_directory,"html",folder)
    dir_exists(dir) || mkdir(dir)
    weave(tmp,doctype = "md2html",out_path=dir,args=args)
  end
  if :pdf ∈ build_list
    println("Building PDF")
    dir = joinpath(repo_directory,"pdf",folder)
    dir_exists(dir) || mkdir(dir)
    weave(tmp,doctype="md2pdf",out_path=dir,args=args)
  end
  if :github ∈ build_list
    println("Building Github Markdown")
    dir = joinpath(repo_directory,"markdown",folder)
    dir_exists(dir) || mkdir(dir)
    weave(tmp,doctype = "github",out_path=dir,args=args)
  end
end


#=
# Needs two arg form
function weave_all()
  foreach(weave_file,
          file for file in readdir("tutorials") if endswith(file, ".jmd"))
end
=#

function bench_footer(folder,file)
  println("""
  These benchmarks are part of the DiffEqBenchmarks.jl repository, found at:

  https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl

  To locally run this benchmark, do the following commands:

  ```julia
  using DiffEqBenchmarks
  DiffEqBenchmarks.weave_file($folder,$file)
  ```

  """)
  println("Computer Information:")
  versioninfo()
  println()
  println("Package Information:")
  Pkg.status()
end

end
