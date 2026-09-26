using LinearAlgebra: Diagonal

@model function demo_dot_assume_observe_matrix_index(
        x = transpose([1.5 2.0;]),
        ::Type{TV} = Vector{Float64},
    ) where {TV}
    s = TV(undef, length(x))
    s .~ InverseGamma(2, 3)
    m = TV(undef, length(x))
    m ~ product_distribution(Normal.(0, sqrt.(s)))
    x[:, 1] ~ MvNormal(m, Diagonal(s))
end

model = demo_dot_assume_observe_matrix_index()
