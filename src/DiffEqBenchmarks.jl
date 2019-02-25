module DiffEqBenchmarks

using Weave, Pkg, IJulia

repo_directory = joinpath(@__DIR__,"..")

function weave_file(folder,file,build_list=(:script,:html,:pdf,:notebook))
  println("File: $file")
  tmp = joinpath(repo_directory,"benchmarks",folder,file)
  args = Dict{Symbol,String}(:folder=>folder,:file=>file)
  if :script ∈ build_list
    println("Building Script")
    dir = joinpath(repo_directory,"script",folder)
    isdir(dir) || mkdir(dir)
    tangle(tmp;out_path=dir)
  end
  if :html ∈ build_list
    println("Building HTML")
    dir = joinpath(repo_directory,"html",folder)
    isdir(dir) || mkdir(dir)
    weave(tmp,doctype = "md2html",out_path=dir,args=args)
  end
  if :pdf ∈ build_list
    println("Building PDF")
    dir = joinpath(repo_directory,"pdf",folder)
    isdir(dir) || mkdir(dir)
    weave(tmp,doctype="md2pdf",out_path=dir,args=args)
  end
  if :github ∈ build_list
    println("Building Github Markdown")
    dir = joinpath(repo_directory,"markdown",folder)
    isdir(dir) || mkdir(dir)
    weave(tmp,doctype = "github",out_path=dir,args=args)
  end
  if :notebook ∈ build_list
    println("Building Notebook")
    dir = joinpath(repo_directory,"notebook",folder)
    isdir(dir) || mkdir(dir)
    Weave.notebook(tmp,dir)
  end
end

function weave_all()
  for folder in readdir(joinpath(repo_directory,"benchmarks"))
    folder == "test.jmd" && continue
    for file in readdir(joinpath(repo_directory,"benchmarks",folder))
      println("Building $(joinpath(folder,file)))")
      try
        weave_file(folder,file)
      catch
      end
    end
  end
end

function bench_footer(folder,file)
  println("""
  These benchmarks are part of the DiffEqBenchmarks.jl repository, found at:

  https://github.com/JuliaDiffEq/DiffEqBenchmarks.jl

  To locally run this benchmark, do the following commands:

  using DiffEqBenchmarks
  DiffEqBenchmarks.weave_file("$folder","$file")

  """)
  println("Computer Information:\n")
  Main.versioninfo()
  println()
  println("Package Information:\n")
  DiffEqBenchmarks.Pkg.status()
end

end
