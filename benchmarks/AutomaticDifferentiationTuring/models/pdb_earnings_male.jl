# PosteriorDB: earnings-logearn_height_male
using PosteriorDB

pdb = PosteriorDB.database()
post = PosteriorDB.posterior(pdb, "earnings-logearn_height_male")
data = PosteriorDB.load(PosteriorDB.dataset(post))

log_earn = log.(Float64.(data["earn"]))
height = Float64.(data["height"])
male = Bool.(data["male"])
N = data["N"]

@model function pdb_earnings_male(log_earn, height, male, N)
    beta ~ filldist(Flat(), 3)
    sigma ~ FlatPos(0.0)
    log_earn ~ MvNormal(beta[1] .+ beta[2] .* height .+ beta[3] .* male, sigma^2 * I)
end

model = pdb_earnings_male(log_earn, height, male, N)
