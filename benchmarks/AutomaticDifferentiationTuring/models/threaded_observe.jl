#=
Note: this model is run with 4 threads
=#

@model function threaded_observe(x)
    a ~ Normal()
    Threads.@threads for i in eachindex(x)
        x[i] ~ Normal(a)
    end
end

x = randn(100)
model = setthreadsafe(threaded_observe(x), true)
