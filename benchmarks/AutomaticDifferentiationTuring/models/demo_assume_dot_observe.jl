@model function demo_assume_dot_observe(x = [1.5, 2.0])
    # `assume` and `dot_observe`
    s ~ InverseGamma(2, 3)
    m ~ Normal(0, sqrt(s))
    x .~ Normal(m, sqrt(s))
end

model = demo_assume_dot_observe()
