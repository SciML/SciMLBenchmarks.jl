@model function demo_dot_assume_observe_index(
        x = [1.5, 2.0],
        ::Type{TV} = Vector{Float64},
    ) where {TV}
    # `dot_assume` and `observe` with indexing
    s = TV(undef, length(x))
    s .~ InverseGamma(2, 3)
    m = TV(undef, length(x))
    m ~ product_distribution(Normal.(0, sqrt.(s)))
    for i in eachindex(x)
        x[i] ~ Normal(m[i], sqrt(s[i]))
    end
end

model = demo_dot_assume_observe_index()
