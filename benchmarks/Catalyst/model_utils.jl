# Model utility functions for Catalyst benchmarks
# Provides helper functions for loading and processing reaction network models

using Catalyst, ModelingToolkit, JSON
using Downloads

"""
    download_catalyst_data(; force=false)

Download the Catalyst benchmark data files from the repository.
"""
function download_catalyst_data(; force=false)
    data_dir = joinpath(@__DIR__, "data")
    mkpath(data_dir)
    
    base_url = "https://raw.githubusercontent.com/SciML/Catalyst_PLOS_COMPBIO_2023/duplication_fix/Benchmarks/Data/"
    
    # Files to download
    files = [
        "multistate.net", "multistate.xml", "multistate.bngl",
        "multisite2.net", "multisite2.xml", "multisite2.bngl",
        "egfr_net.net", "egfr_net.xml", "egfr_net.bngl",
        "BCR.net", "BCR.xml", "BCR.bngl", 
        "fceri_gamma2.net", "fceri_gamma2.xml", "fceri_gamma2.bngl"
    ]
    
    downloaded_files = []
    
    for file in files
        filepath = joinpath(data_dir, file)
        if !isfile(filepath) || force
            try
                url = base_url * file
                Downloads.download(url, filepath)
                push!(downloaded_files, filepath)
                println("Downloaded: $(file)")
            catch e
                println("Failed to download $(file): $(e)")
            end
        else
            println("File already exists: $(file)")
        end
    end
    
    return downloaded_files
end

"""
    load_catalyst_model(model_name::String)

Load a Catalyst reaction network model from the benchmark data.
"""
function load_catalyst_model(model_name::String)
    data_dir = joinpath(@__DIR__, "data")
    net_file = joinpath(data_dir, "$(model_name).net")
    
    if !isfile(net_file)
        println("Model file not found: $(net_file)")
        println("Attempting to download...")
        download_catalyst_data()
    end
    
    if !isfile(net_file)
        error("Could not find or download model file: $(net_file)")
    end
    
    # Load the reaction network
    try
        rn = loadrxnetwork(BNGNetwork(), net_file)
        return rn
    catch e
        println("Error loading model $(model_name): $(e)")
        rethrow(e)
    end
end

"""
    create_catalyst_problems(model_name::String; tspan=(0.0, 10.0))

Create ODE and Jump problems for a Catalyst model.
"""
function create_catalyst_problems(model_name::String; tspan=(0.0, 10.0))
    # Load the reaction network
    rn = load_catalyst_model(model_name)
    
    # Get species and parameters
    species_list = species(rn)
    params_list = parameters(rn)
    
    # Create initial conditions
    u0 = []
    for sp in species_list
        # Default to zero, but set some initial concentrations for specific models
        if model_name == "multistate" && string(sp) in ["A", "B"]
            push!(u0, sp => 100.0)
        elseif model_name == "multisite2" && string(sp) == "RecMon"
            push!(u0, sp => 10.0)
        elseif model_name == "egfr_net" && string(sp) in ["EGF", "EGFR"]
            push!(u0, sp => 680.0)
        elseif model_name == "BCR" && string(sp) == "Lyn"
            push!(u0, sp => 28000.0)
        elseif model_name == "fceri_gamma2" && string(sp) == "Lyn"
            push!(u0, sp => 28000.0)
        else
            push!(u0, sp => 0.0)
        end
    end
    
    # Create parameter values (use default values or common biochemical values)
    p = []
    for param in params_list
        param_name = string(param)
        if occursin("kon", param_name) || occursin("kf", param_name)
            push!(p, param => 1e-3)  # Forward rate constants
        elseif occursin("koff", param_name) || occursin("kr", param_name)
            push!(p, param => 1e-2)  # Reverse rate constants
        elseif occursin("kcat", param_name) || occursin("k", param_name)
            push!(p, param => 1e-1)  # Catalytic rate constants
        else
            push!(p, param => 1.0)   # Default value
        end
    end
    
    # Create ODE problem
    ode_prob = ODEProblem(rn, u0, tspan, p)
    
    # Create Jump problem  
    jump_prob = JumpProblem(rn, u0, tspan, p)
    
    return (ode_prob=ode_prob, jump_prob=jump_prob, rn=rn)
end

"""
    get_model_info(model_name::String)

Get information about a benchmark model.
"""
function get_model_info(model_name::String)
    model_specs = Dict(
        "multistate" => (species=9, reactions=18, description="Simple multi-state system"),
        "multisite2" => (species=66, reactions=288, description="Multi-site protein interactions"),
        "egfr_net" => (species=356, reactions=3749, description="Epidermal growth factor receptor"),
        "BCR" => (species=1122, reactions=24388, description="B-cell receptor signaling"),
        "fceri_gamma2" => (species=3744, reactions=58276, description="High-affinity IgE receptor")
    )
    
    if haskey(model_specs, model_name)
        return model_specs[model_name]
    else
        return (species=0, reactions=0, description="Unknown model")
    end
end

"""
    save_benchmark_results(results::Dict, filename::String)

Save benchmark results to a JSON file.
"""
function save_benchmark_results(results::Dict, filename::String)
    # Convert any non-serializable types to strings
    serializable_results = deepcopy(results)
    
    # Recursively convert problematic types
    function make_serializable(obj)
        if obj isa Dict
            return Dict(string(k) => make_serializable(v) for (k, v) in obj)
        elseif obj isa Array
            return [make_serializable(item) for item in obj]
        elseif obj isa Number && (isnan(obj) || isinf(obj))
            return string(obj)
        else
            return obj
        end
    end
    
    serializable_results = make_serializable(serializable_results)
    
    # Save to JSON
    open(filename, "w") do f
        JSON.print(f, serializable_results, 4)
    end
    
    println("Results saved to: $(filename)")
end

"""
    load_benchmark_results(filename::String)

Load benchmark results from a JSON file.
"""
function load_benchmark_results(filename::String)
    if !isfile(filename)
        error("File not found: $(filename)")
    end
    
    return JSON.parsefile(filename)
end

"""
    create_model_summary()

Create a summary of all benchmark models.
"""
function create_model_summary()
    models = ["multistate", "multisite2", "egfr_net", "BCR", "fceri_gamma2"]
    
    summary = []
    for model in models
        info = get_model_info(model)
        push!(summary, (
            name=model,
            species=info.species,
            reactions=info.reactions,
            description=info.description
        ))
    end
    
    return summary
end

"""
    print_model_summary()

Print a formatted summary of all benchmark models.
"""
function print_model_summary()
    summary = create_model_summary()
    
    println("\\n=== Catalyst Benchmark Models ===")
    println("Model Name       | Species | Reactions | Description")
    println("-"^70)
    
    for model in summary
        println("$(rpad(model.name, 15)) | $(lpad(model.species, 7)) | $(lpad(model.reactions, 9)) | $(model.description)")
    end
    
    println("-"^70)
end