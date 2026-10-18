# PosteriorDB: garch-garch11
using PosteriorDB

pdb = PosteriorDB.database()
post = PosteriorDB.posterior(pdb, "garch-garch11")
data = PosteriorDB.load(PosteriorDB.dataset(post))

T = data["T"]
y = Float64.(data["y"])
sigma1 = Float64(data["sigma1"])

@model function pdb_garch11(T, y, sigma1)
    mu ~ Flat()
    alpha0 ~ FlatPos(0.0)
    alpha1 ~ Uniform(0, 1)
    beta1 ~ Uniform(0, 1 - alpha1 + eps())

    sigma = Vector{typeof(mu)}(undef, T)
    sigma[1] = sigma1
    for t in 2:T
        sigma[t] = sqrt(alpha0 + alpha1 * (y[t - 1] - mu)^2 + beta1 * sigma[t - 1]^2)
    end
    for t in 1:T
        y[t] ~ Normal(mu, sigma[t])
    end
end

model = pdb_garch11(T, y, sigma1)
