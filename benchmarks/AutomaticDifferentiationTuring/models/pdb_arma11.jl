# PosteriorDB: arma-arma11
using PosteriorDB

pdb = PosteriorDB.database()
post = PosteriorDB.posterior(pdb, "arma-arma11")
data = PosteriorDB.load(PosteriorDB.dataset(post))

T = data["T"]
y = Float64.(data["y"])

@model function pdb_arma11(T, y)
    mu ~ Normal(0, 10)
    phi ~ Normal(0, 2)
    theta ~ Normal(0, 2)
    sigma ~ truncated(Cauchy(0, 2.5); lower = 0)

    nu = Vector{typeof(mu)}(undef, T)
    err = Vector{typeof(mu)}(undef, T)
    nu[1] = mu + (phi * mu)
    err[1] = y[1] - nu[1]
    for t in 2:T
        nu[t] = mu + (phi * y[t - 1]) + (theta * err[t - 1])
        err[t] = y[t] - nu[t]
    end
    y ~ MvNormal(nu, sigma^2 * I)
end

model = pdb_arma11(T, y)
