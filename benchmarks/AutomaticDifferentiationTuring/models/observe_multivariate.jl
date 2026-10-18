@model function observe_multivariate(
        x = [1.5, 2.0, 2.5],
        ::Type{TV} = Vector{Float64},
    ) where {TV}
    a = TV(undef, length(x))
    a .~ Normal()
    x ~ MvNormal(a, I)
end

model = observe_multivariate()
