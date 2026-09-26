#=
Note: this example is run with 4 threads
=#

@model function threaded_assume(x)
    a = Vector{Float64}(undef, length(x))
    Threads.@threads for i in eachindex(x)
        a[i] ~ Normal()
        x[i] ~ Normal(a[i])
    end
end

x = randn(50)
model = setthreadsafe(threaded_assume(x), true)
