folder, file, caller_project, benchmark_project, output_count = ARGS[1:5]
use_benchmark_project = benchmark_project == "true"
first_load_path_arg = 6 + parse(Int, output_count)
build_list = Symbol.(ARGS[6:(first_load_path_arg - 1)])
append!(empty!(LOAD_PATH), ARGS[first_load_path_arg:end])

include("weave_bootstrap.jl")
_bootstrap_weave(folder, file, build_list, caller_project, use_benchmark_project)
