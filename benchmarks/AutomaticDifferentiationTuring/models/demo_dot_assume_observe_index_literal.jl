@model function demo_dot_assume_observe_index_literal(
        ::Type{TV} = Vector{Float64},
    ) where {TV}
    # `dot_assume` and literal `observe` with indexing
    s = TV(undef, 2)
    m = TV(undef, 2)
    s .~ InverseGamma(2, 3)
    m ~ product_distribution(Normal.(0, sqrt.(s)))

    1.5 ~ Normal(m[1], sqrt(s[1]))
    2.0 ~ Normal(m[2], sqrt(s[2]))
end

model = demo_dot_assume_observe_index_literal()
