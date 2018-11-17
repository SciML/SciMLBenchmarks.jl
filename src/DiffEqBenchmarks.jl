module DiffEqBenchmarks

using Weave

repo_directory = joinpath(@__DIR__,"..")

function weave_file(folder,file)
  println("File: $file")
  tmp = joinpath(repo_directory,folder,file)
  #weave(tmp,doctype="pandoc",out_path=:pwd)
  println("Building Script")
  tangle(tmp;out_path=joinpath(repo_directory,folder))
  println("Building HTML")
  weave(tmp,doctype = "md2html",out_path=joinpath(repo_directory,folder))
  println("Building PDF")
  weave(tmp,doctype="md2pdf",out_path=joinpath(repo_directory,folder))
end

function weave_all()
  foreach(weave_file,
          file for file in readdir("tutorials") if endswith(file, ".jmd"))
end

end
