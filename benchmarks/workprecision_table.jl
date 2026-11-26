using Printf: @sprintf

function usable_points(wp)
    return findall(eachindex(wp.times)) do i
        error = wp.errors.final[i]
        return isfinite(error) && error > 0 && isfinite(wp.times[i]) && wp.times[i] > 0
    end
end

function refinement_rate(steps, errors, error_floor)
    for i in (length(steps) - 2):-1:1
        h = steps[i:(i + 2)]
        e = errors[i:(i + 2)]
        all(x -> isfinite(x) && x > error_floor, e) || continue
        h[1] > 0 && 0 < h[2] <= 0.8 * h[1] && 0 < h[3] <= 0.8 * h[2] || continue
        log_steps, log_errors = log.(h), log.(e)
        slopes = diff(log_errors) ./ diff(log_steps)
        all(>(0), slopes) && abs(slopes[1] - slopes[2]) <= max(0.25, 0.15 * minimum(slopes)) || continue
        return @sprintf "%.2f" (last(log_errors) - first(log_errors)) / (last(log_steps) - first(log_steps))
    end
    return "—"
end

function measured_order(wp, points, error_floor)
    isnothing(wp.dts) && return "—"
    interval = last(wp.prob.tspan) - first(wp.prob.tspan)
    errors = fill(NaN, length(wp.dts))
    for i in points
        expected = round(Int, interval / wp.dts[i])
        stats = wp.stats[i]
        if !isnothing(stats) && stats.naccept == expected && stats.nreject == 0
            errors[i] = wp.errors.final[i]
        end
    end
    return refinement_rate(wp.dts, errors, error_floor)
end

function adaptive_step_rate(wp, points, error_floor)
    :adaptive in wp.tags || return "—"
    interval = last(wp.prob.tspan) - first(wp.prob.tspan)
    steps = fill(NaN, wp.N)
    for i in points
        stats = wp.stats[i]
        isnothing(stats) && continue
        stats.naccept > 0 || continue
        steps[i] = interval / stats.naccept
    end
    return refinement_rate(steps, wp.errors.final, error_floor)
end

function controlled_order(prob, alg, counts, endpoint, error_floor)
    interval = last(prob.tspan) - first(prob.tspan)
    steps = interval ./ counts
    errors = map(counts, steps) do count, dt
        sol = solve(
            prob, alg; adaptive = false, dt, abstol = 1.0e-14, reltol = 1.0e-12,
            save_everystep = false, dense = false, verbose = None(),
        )
        if successful_retcode(sol) && sol.t[end] == last(prob.tspan) &&
                sol.stats.naccept == count && sol.stats.nreject == 0
            return sum(abs, sol.u[end] - endpoint) / length(endpoint)
        end
        return NaN
    end
    return refinement_rate(steps, errors, error_floor)
end

function work_table_markdown(wps, error_floor; orders = Dict{String, String}())
    io = IOBuffer()
    println(io, "| Algorithm | Tags | Usable points | Observed order | Adaptive-step rate | Best error | Time (s) | RHS calls | Accepted steps |")
    println(io, "|:--|:--|--:|--:|--:|--:|--:|--:|--:|")
    for wp in wps[:]
        points = usable_points(wp)
        tags = join(wp.tags, ", ")
        order = get(() -> measured_order(wp, points, error_floor), orders, wp.name)
        rate = adaptive_step_rate(wp, points, error_floor)
        if isempty(points)
            println(io, "| $(wp.name) | $tags | 0/$(wp.N) | $order | $rate | — | — | — | — |")
            continue
        end
        i = points[argmin(wp.errors.final[points])]
        stats = wp.stats[i]
        rhs_calls = isnothing(stats) || wp.name == "PDIRK44" ? "—" :
            iszero(stats.nf + stats.nf2) ? "—" : string(stats.nf + stats.nf2)
        accepted = isnothing(stats) ? "—" : string(stats.naccept)
        println(
            io,
            "| $(wp.name) | $tags | $(length(points))/$(wp.N) | $order | $rate | ",
            @sprintf("%.3g | %.3g | ", wp.errors.final[i], wp.times[i]),
            rhs_calls, " | ", accepted, " |",
        )
    end
    return String(take!(io))
end

function print_work_table(wps, error_floor; kwargs...)
    return display("text/markdown", work_table_markdown(wps, error_floor; kwargs...))
end

function assert_work_table_markdown_fixture()
    fixtures = (
        ("PDIRK44", 1384, 32, "—"), ("KuttaPRK2p5", 0, 128, "—"),
        ("Tsit5", 441, 73, "441"), ("lsoda", nothing, nothing, "—"),
    )
    wps = map(fixtures) do (name, calls, steps, _)
        stats = isnothing(calls) ? nothing : (nf = calls, nf2 = 0, naccept = steps)
        return (
            name = name, times = [1.0], errors = (final = [0.01],), tags = (),
            dts = nothing, N = 1, stats = [stats],
        )
    end
    rows = split(chomp(work_table_markdown(collect(wps), 0.0)), '\n')[3:end]
    @assert length(rows) == length(fixtures)
    for (row, (name, _, steps, expected_rhs)) in zip(rows, fixtures)
        columns = strip.(split(row, '|'; keepempty = false))
        @assert columns[1] == name
        rhs, accepted = columns[8:9]
        expected_accepted = isnothing(steps) ? "—" : string(steps)
        @assert rhs == expected_rhs "$name RHS calls should be $expected_rhs, got $rhs"
        @assert accepted == expected_accepted "$name accepted steps should be $expected_accepted, got $accepted"
    end
    return nothing
end
