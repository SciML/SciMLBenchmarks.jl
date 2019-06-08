
using OrdinaryDiffEq, Plots, DiffEqCallbacks, LinearAlgebra
using TaylorIntegration

T(p) = A / 2 * norm(p)^2
V(q) = A / 2 * (q[1]^2 + q[2]^2) + B / √2 * q[1] * (3 * q[2]^2 - q[1]^2) + D / 4 * (q[1]^2 + q[2]^2)^2
H(p, q, params) = T(p) + V(q)

const A, B, D = 1., 0.55, 0.4

module InPlace
    using ParameterizedFunctions

    function q̇(dq, p, q, params, t)
        dq[1] = A * p[1]
        dq[2] = A * p[2]
    end

    function ṗ(dp, p, q, params, t)
        dp[1] = -A * q[1] - 3 * B / √2 * (q[2]^2 - q[1]^2) - D * q[1] * (q[1]^2 + q[2]^2)
        dp[2] = -q[2] * (A + 3 * √2 * B * q[1] + D * (q[1]^2 + q[2]^2))
    end

    const A, B, D = 1., 0.55, 0.4

    const q0 = [4.919080920016389, 2.836942666663649]
    const p0 = [0., 0.]
    const u0 = vcat(p0,q0)

    h_eqs = @ode_def HamiltEqs begin
      dp₀ = -A * q₀ - 3 * B / √2 * (q₂^2 - q₀^2) - D * q₀ * (q₀^2 + q₂^2)
      dp₂ = -q₂ * (A + 3 * √2 * B * q₀ + D * (q₀^2 + q₂^2))
      dq₀ = A * p₀
      dq₂ = A * p₂
    end A B D

    p = [1.0,0.55,0.4]

end

module OutOfPlace
    using StaticArrays

    function q̇(p, q, params, t)
        p
    end

    function ṗ(p, q, params, t)
        dp1 = -A * q[1] - 3 * B / √2 * (q[2]^2 - q[1]^2) - D * q[1] * (q[1]^2 + q[2]^2)
        dp2 = -q[2] * (A + 3 * √2 * B * q[1] + D * (q[1]^2 + q[2]^2))
        @SVector [dp1, dp2]
    end

    const A, B, D = 1., 0.55, 0.4

    const q0 = @SVector [4.919080920016389, 2.836942666663649]
    const p0 = @SVector [0., 0.]
    const u0 = vcat(p0,q0)

    h_eqs(z, params, t) = SVector(
        -A * z[3] - 3 * B / √2 * (z[4]^2 - z[3]^2) - D * z[3] * (z[3]^2 + z[4]^2),
        -z[4] * (A + 3 * √2 * B * z[3] + D * (z[3]^2 + z[4]^2)),
        z[1],
        z[2]
    )

    p = nothing

end

function g(resid,u)
    resid[1] = H([u[1],u[2]],[u[3],u[4]],nothing) - E
    resid[2:4] .= 0
end

const E = H(InPlace.p0, InPlace.q0, InPlace.p)
const cb = ManifoldProjection(g, nlopts=Dict(:ftol=>1e-13))


energy_err(sol) = map(i->H([sol[1,i], sol[2,i]], [sol[3,i], sol[4,i]],nothing)-E, 1:length(sol.u))
abs_energy_err(sol) = [abs.(H([sol[1,j], sol[2,j]], [sol[3,j], sol[4,j]],nothing) - E) for j=1:length(sol.u)]

function compare(mode=InPlace, all=true, plt=nothing; tmax=1e2)
    prob1 = DynamicalODEProblem(mode.ṗ,  mode.q̇, mode.p0, mode.q0, (0., tmax))
    prob2 = ODEProblem(mode.h_eqs, mode.u0, (0., tmax), mode.p)

    GC.gc()
    (mode == InPlace  && all) && @time sol1 = solve(prob2, Vern9(), callback=cb, abstol=1e-14, reltol=1e-14)
    GC.gc()
    @time sol2 = solve(prob1, KahanLi8(), dt=1e-2, maxiters=1e10)
    GC.gc()
    @time sol3 = solve(prob1, SofSpa10(), dt=1e-2, maxiters=1e8)
    GC.gc()
    @time sol4 = solve(prob2, Vern9(), abstol=1e-14, reltol=1e-14)
    GC.gc()
    @time sol5 = solve(prob1, DPRKN12(), abstol=1e-14, reltol=1e-14)
    GC.gc()
    (mode == InPlace && all) && @time sol6 = solve(prob2, TaylorMethod(50), abstol=1e-20)

    (mode == InPlace && all) && println("Vern9 + ManifoldProjection max energy error:\t"*
        "$(maximum(abs_energy_err(sol1)))\tin\t$(length(sol1.u))\tsteps.")
    println("KahanLi8 max energy error:\t\t\t$(maximum(abs_energy_err(sol2)))\tin\t$(length(sol2.u))\tsteps.")
    println("SofSpa10 max energy error:\t\t\t$(maximum(abs_energy_err(sol3)))\tin\t$(length(sol3.u))\tsteps.")
    println("Vern9 max energy error:\t\t\t\t$(maximum(abs_energy_err(sol4)))\tin\t$(length(sol4.u))\tsteps.")
    println("DPRKN12 max energy error:\t\t\t$(maximum(abs_energy_err(sol5)))\tin\t$(length(sol5.u))\tsteps.")
    (mode == InPlace && all) && println("TaylorMethod max energy error:\t\t\t$(maximum(abs_energy_err(sol6)))"*
        "\tin\t$(length(sol6.u))\tsteps.")

    if plt == nothing
        plt = plot(xlabel="t", ylabel="Energy error")
    end

    (mode == InPlace && all) && plot!(sol1.t, energy_err(sol1), label="Vern9 + ManifoldProjection")
    plot!(sol2.t, energy_err(sol2), label="KahanLi8", ls=mode==InPlace ? :solid : :dash)
    plot!(sol3.t, energy_err(sol3), label="SofSpa10", ls=mode==InPlace ? :solid : :dash)
    plot!(sol4.t, energy_err(sol4), label="Vern9", ls=mode==InPlace ? :solid : :dash)
    plot!(sol5.t, energy_err(sol5), label="DPRKN12", ls=mode==InPlace ? :solid : :dash)
    (mode == InPlace && all) && plot!(sol6.t, energy_err(sol6), label="TaylorMethod")

    return plt
end


compare(tmax=1e2)


compare(tmax=1e3)


compare(tmax=1e4)


compare(tmax=2e4)


function in_vs_out(;all=false, tmax=1e2)
    println("In place versions:")
    plt = compare(InPlace, all, tmax=tmax)
    println("\nOut of place versions:")
    plt = compare(OutOfPlace, false, plt; tmax=tmax)
end


in_vs_out(all=true, tmax=1e3)


in_vs_out(tmax=1e2)


in_vs_out(tmax=1e3)


in_vs_out(tmax=1e4)


in_vs_out(tmax=2e4)


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

