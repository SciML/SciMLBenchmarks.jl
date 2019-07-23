using Pkg; Pkg.add.(["Plots", "DSP"])
using Weave
#HTML
weave(joinpath(dirname(pathof(Weave)), "$(homedir())/Desktop/adff.jmd"),
  out_path=:pwd,
  doctype = "md2html")
