using LinearAlgebra: Diagonal

@model function demo_assume_multivariate_observe_literal()
    # multivariate `assume` and literal `observe`
    s ~ product_distribution([InverseGamma(2, 3), InverseGamma(2, 3)])
    m ~ MvNormal(zeros(2), Diagonal(s))
    [1.5, 2.0] ~ MvNormal(m, Diagonal(s))
end

model = demo_assume_multivariate_observe_literal()
