using LinearAlgebra: Diagonal

@model function demo_assume_index_observe(
        x = [1.5, 2.0],
        ::Type{TV} = Vector{Float64},
    ) where {TV}
    # `assume` with indexing and `observe`
    s = TV(undef, length(x))
    for i in eachindex(s)
        s[i] ~ InverseGamma(2, 3)
    end
    m = TV(undef, length(x))
    for i in eachindex(m)
        m[i] ~ Normal(0, sqrt(s[i]))
    end
    x ~ MvNormal(m, Diagonal(s))
end

model = demo_assume_index_observe()
