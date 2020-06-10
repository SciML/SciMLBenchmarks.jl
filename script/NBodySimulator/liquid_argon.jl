
using ProgressLogging
using NBodySimulator, OrdinaryDiffEq, StaticArrays
using Plots, DataFrames, StatsPlots

function setup(t)
    T = 120.0 # K
    kb = 1.38e-23 # J/K
    ϵ = T * kb # J
    σ = 3.4e-10 # m
    ρ = 1374 # kg/m^3
    m = 39.95 * 1.6747 * 1e-27 # kg
    N = 216
    L = (m*N/ρ)^(1/3)
    R = 2.25σ
    v_dev = sqrt(kb * T / m) # m/s

    _L = L / σ
    _σ = 1.0
    _ϵ = 1.0
    _m = 1.0
    _v = v_dev / sqrt(ϵ / m)
    _R = R / σ

    bodies = generate_bodies_in_cell_nodes(N, _m, _v, _L)
    lj_parameters = LennardJonesParameters(_ϵ, _σ, _R)
    pbc = CubicPeriodicBoundaryConditions(_L)
    lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => lj_parameters));
    simulation = NBodySimulation(lj_system, (0.0, t), pbc, _ϵ/T)

    return simulation
end


config(τ, at, rt) = [
    # symplectic
    (alg=VelocityVerlet, dt=τ),
    (alg=McAte2, dt=τ),
    # (alg=CalvoSanz4, dt=τ),
    # (alg=McAte5, dt=τ),
    # (alg=Yoshida6, dt=τ),
    # (alg=KahanLi8, dt=τ),
    # DPRKN
    (alg=DPRKN6, abstol=at, rtol=rt),
    # (alg=DPRKN8, abstol=at, rtol=rt),
    (alg=DPRKN12, abstol=at, rtol=rt),
    # others
    (alg=Tsit5, abstol=at, rtol=rt),
    (alg=Vern7, abstol=at, rtol=rt),
    (alg=Vern9, abstol=at, rtol=rt)
]


function benchmark(energyerr, rts, ts, t, configs)
    simulation = setup(t)
    prob = SecondOrderODEProblem(simulation)
    for config in configs
        alg = config.alg
        sol, rt = @timed solve(prob, alg(); progress=true, progress_name="$alg", config...)
        result = NBodySimulator.SimulationResult(sol, simulation)
        ΔE(t) = total_energy(result, t) - total_energy(result, 0)
        energyerr[alg] = [ΔE(t) for t in sol.t[2:end]]
        rts[alg] = rt
        ts[alg] = sol.t[2:end]
    end
end

ΔE = Dict()
rt = Dict()
ts = Dict()
configs = config(2.3e-4, 1e-20, 1e-20)
benchmark(ΔE, rt, ts, 40., configs)

plt = plot(xlabel="Rescaled Time", ylabel="Energy error", legend=:topleft);
for c in configs
    plot!(plt, ts[c.alg], abs.(ΔE[c.alg]), label="$(c.alg), $(rt[c.alg])s", xscale=:log10, yscale=:log10)
end
plt


t = 35.0
results = DataFrame(:integrator=>String[], :runtime=>Float64[], :τ=>Float64[], :abstol=>Float64[],
    :EnergyError=>Float64[], :timesteps=>Int[], :f_evals=>Int[]);

function benchmark(energyerr, rts, bytes, allocs, nt, nf, t, configs)
    simulation = setup(t)
    prob = SecondOrderODEProblem(simulation)
    for config in configs
        alg = config.alg
        sol, rt, b, gc, memalloc = @timed solve(prob, alg();
            save_everystep=false, progress=true, progress_name="$alg", config...)
        result = NBodySimulator.SimulationResult(sol, simulation)
        ΔE = total_energy(result, t) - total_energy(result, 0)
        energyerr[alg] = ΔE
        rts[alg] = rt
        bytes[alg] = b
        allocs[alg] = memalloc
        nt[alg] = sol.destats.naccept
        nf[alg] = sol.destats.nf + sol.destats.nf2
    end
end


τs = 10 .^range(-4, -3, length=5)
ats = 10 .^range(-20, -14, length=5)
rts = 10 .^range(-20, -14, length=5)


@progress "Variable dt" for (τ, at, rt) in zip(τs, ats, rts)
    runtime = Dict()
    ΔE = Dict()
    nt = Dict()
    nf = Dict()
    b = Dict()
    allocs = Dict()

    GC.gc()
    benchmark(ΔE, runtime, b, allocs, nt, nf, t, config(τ, at, rt))

    for (k,v) in ΔE
        push!(results, [string(k), runtime[k], at, rt, abs(v), nt[k], nf[k]])
    end
end


results


@df results plot(:EnergyError, :runtime, group=:integrator,
    xscale=:log10, yscale=:log10, xlabel="Energy error", ylabel="Runtime (s)")


@df results plot(:EnergyError, :timesteps, group=:integrator,
    xscale=:log10, yscale=:log10, xlabel="Energy error", ylabel="Number of timesteps")


@df results plot(:f_evals, :runtime, group=:integrator,
    xscale=:log10, yscale=:log10, xlabel="Number of f calls", ylabel="Runtime (s)")

