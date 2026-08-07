# PosteriorDB: sblri-blr
using PosteriorDB

pdb = PosteriorDB.database()
post = PosteriorDB.posterior(pdb, "sblri-blr")
data = PosteriorDB.load(PosteriorDB.dataset(post))

X = Float64.(data["X"])
y = Float64.(data["y"])
D = data["D"]

@model function pdb_sblri(X, y, D)
    beta ~ filldist(Normal(0, 10), D)
    sigma ~ truncated(Normal(0, 10); lower = 0)
    y ~ MvNormal(X * beta, sigma^2 * I)
end

model = pdb_sblri(X, y, D)
