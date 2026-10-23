using DiffEqDevTools, PETScDiffEq, Plots, Printf, SciMLBase

const COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]

is_petsc(alg) = parentmodule(typeof(alg)) === PETScDiffEq
uses_lu(alg) = is_petsc(alg) && "lu" in alg.petsc_options

function line_style(alg)
    uses_lu(alg) && return (:dot, :utriangle)
    return is_petsc(alg) ? (:solid, :circle) : (:dash, :rect)
end

function reference(prob, alg; tol)
    return solve(prob, alg; abstol = tol, reltol = tol, save_everystep = false, maxiters = 10^8)
end

function reference_margin(wp, gap)
    smallest = minimum(e for w in wp.wps for e in w.errors.final if isfinite(e))
    @printf("smallest error %.2e, reference gap %.2e, ratio %.0f\n", smallest, gap, smallest / gap)
    return nothing
end

function wpd_plot(wp, colors; title)
    styles = [line_style(s[:alg]) for s in wp.setups]
    return plot(
        wp; title, label = permutedims(wp.names), color = permutedims(COLORS[colors]),
        linestyle = permutedims(first.(styles)), markershape = permutedims(last.(styles)),
        xticks = exp10.(-16:0), linewidth = 2, markersize = 4, markerstrokewidth = 0,
        legend = :outerright, size = (900, 480), titlefontsize = 11, legendfontsize = 8,
        left_margin = 4Plots.mm, bottom_margin = 4Plots.mm,
    )
end

function print_points(wp)
    pad = maximum(textwidth, wp.names) + 2
    @printf(
        "%s%7s %7s %10s %10s %6s %7s %6s %7s\n", rpad("solver", pad), "abstol", "reltol",
        "error", "time (ms)", "steps", "rejects", "jacs", "newton",
    )
    for (name, w) in zip(wp.names, wp.wps), i in 1:(w.N)
        st = w.stats[i]
        @printf(
            "%s%7.0e %7.0e %10.3e %10.4g %6d %7d %6d %7d\n", rpad(name, pad), w.abstols[i],
            w.reltols[i], w.errors.final[i], 1.0e3 * w.times[i], st.naccept, st.nreject,
            st.njacs, st.nnonliniter,
        )
    end
    return nothing
end

function time_at(points, target)
    front = Tuple{Float64, Float64}[]
    for p in sort(points; by = last)
        (isempty(front) || p[1] < front[end][1]) && push!(front, p)
    end
    isempty(front) && return NaN, false
    front[1][1] <= target && return front[1][2], true
    for ((e1, t1), (e2, t2)) in zip(front, front[2:end])
        e1 >= target >= e2 || continue
        w = log(e1 / target) / log(e1 / e2)
        return exp((1 - w) * log(t1) + w * log(t2)), false
    end
    return NaN, false
end

function target_times(wp, target)
    return map(wp.wps) do w
        points = [(e, t) for (e, t) in zip(w.errors.final, w.times) if isfinite(e * t)]
        time_at(points, target)
    end
end

function time_to_error(wp, target)
    times = target_times(wp, target)
    best = minimum(t for (t, _) in times if isfinite(t))
    pad = maximum(textwidth, wp.names) + 2
    @printf("Time to an error of %.0e\n\n", target)
    @printf("%s%14s %10s\n", rpad("solver", pad), "time (ms)", "relative")
    for (name, (t, bound)) in zip(wp.names, times)
        if isfinite(t)
            cell = (bound ? "at most " : "") * @sprintf("%.3g", 1.0e3 * t)
            @printf("%s%14s %10.1f\n", rpad(name, pad), cell, t / best)
        else
            @printf("%s%14s\n", rpad(name, pad), "not reached")
        end
    end
    return nothing
end

function ratio_text(a, b)
    qualifier = a[2] == b[2] ? "" : a[2] ? "at most " : "at least "
    return qualifier * @sprintf("%.1f", a[1] / b[1])
end

function summarize(wp, target, pairs)
    times = target_times(wp, target)
    reached = findall(t -> isfinite(t[1]), times)
    best = reached[argmin([times[i][1] for i in reached])]
    petsc = filter(i -> is_petsc(wp.setups[i][:alg]), reached)
    fastest = petsc[argmin([times[i][1] for i in petsc])]
    @printf("At an error of %.0e the fastest method is %s. ", target, wp.names[best])
    @printf(
        "The fastest PETScDiffEq method, %s, takes %s times as long.\n", wp.names[fastest],
        ratio_text(times[fastest], times[best]),
    )
    for (a, b) in pairs
        i, j = findfirst(==(a), wp.names), findfirst(==(b), wp.names)
        @printf("%s takes %s times as long as %s.\n", a, ratio_text(times[i], times[j]), b)
    end
    return nothing
end

function compare_pairs(wp, pairs)
    w0 = wp.wps[1]
    @printf("At abstol %.0e, reltol %.0e\n\n", w0.abstols[end], w0.reltols[end])
    pad = maximum(textwidth, wp.names) + 2
    @printf(
        "%s%9s %10s %10s %15s\n", rpad("solver", pad), "attempts", "error", "time (ms)",
        "µs per attempt",
    )
    for pair in pairs, name in pair
        w = wp.wps[findfirst(==(name), wp.names)]
        st = w.stats[w.N]
        n = st.naccept + st.nreject
        @printf(
            "%s%9d %10.2e %10.3g %15.3g\n", rpad(name, pad), n, w.errors.final[w.N],
            1.0e3 * w.times[w.N], 1.0e6 * w.times[w.N] / n,
        )
    end
    return nothing
end

function fixed_step_gap(prob, a, b; dt)
    sa = solve(prob, a; dt, adaptive = false, save_everystep = false)
    sb = solve(prob, b; dt, adaptive = false, save_everystep = false)
    return maximum(abs, sa.u[end] - sb.u[end])
end

function one_step_time(prob, alg; numruns = 20)
    h = 1.0e-6 * (prob.tspan[2] - prob.tspan[1])
    short = remake(prob; tspan = (prob.tspan[1], prob.tspan[1] + h))
    run() = solve(short, alg; dt = h, adaptive = false, save_everystep = false)
    run()
    return minimum(@elapsed(run()) for _ in 1:numruns)
end

function one_step_table(probs, entries)
    pad = maximum(textwidth, first.(entries)) + 2
    @printf("%s%20s\n", rpad("solver", pad), "one-step solve (ms)")
    for (name, alg, choice) in entries
        @printf("%s%20.3g\n", rpad(name, pad), 1.0e3 * one_step_time(probs[choice], alg))
    end
    return nothing
end
