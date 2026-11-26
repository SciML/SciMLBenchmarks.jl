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

function _setup_value(setup, key)
    setup === nothing && return nothing
    if setup isa AbstractDict
        return haskey(setup, key) ? setup[key] : nothing
    end
    return hasproperty(setup, key) ? getproperty(setup, key) : nothing
end

function _effective_problem(wp, setup = nothing)
    prob = wp.prob
    prob isa AbstractArray || return prob
    choice = something(_setup_value(setup, :prob_choice), 1)
    return prob[choice]
end

function measured_order(wp, points, error_floor, prob = _effective_problem(wp))
    isnothing(wp.dts) && return "—"
    interval = last(prob.tspan) - first(prob.tspan)
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

function adaptive_step_rate(wp, points, error_floor, prob = _effective_problem(wp))
    :adaptive in wp.tags || return "—"
    :accepted_steps_unavailable in wp.tags && return "—"
    interval = last(prob.tspan) - first(prob.tspan)
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

function _work_setup(wps, setups, index)
    source = setups
    if source === nothing && hasproperty(wps, :setups)
        source = wps.setups
    end
    source isa AbstractVector && length(source) >= index || return nothing
    return source[index]
end

"""
    work_table_markdown(wps, error_floor; orders = Dict(), setups = nothing)

Render endpoint work-precision measurements as a Markdown table. `orders` overrides
measured fixed-step orders by series name. Array-valued problems use the corresponding
setup's `:prob_choice`, defaulting to the first problem as `WorkPrecisionSet` does.
Tag a setup with `:accepted_steps_unavailable` when its solver wrapper does not report
accepted-step counts; the table then displays an em dash while preserving genuine
measured zero counts.
"""
function work_table_markdown(wps, error_floor; orders = Dict{String, String}(), setups = nothing)
    io = IOBuffer()
    println(io, "| Algorithm | Tags | Usable points | Observed order | Adaptive-step rate | Best error | Time (s) | RHS calls | Accepted steps |")
    println(io, "|:--|:--|--:|--:|--:|--:|--:|--:|--:|")
    for (index, wp) in enumerate(wps[:])
        setup = _work_setup(wps, setups, index)
        prob = _effective_problem(wp, setup)
        points = usable_points(wp)
        tags = join(wp.tags, ", ")
        order = get(() -> measured_order(wp, points, error_floor, prob), orders, wp.name)
        rate = adaptive_step_rate(wp, points, error_floor, prob)
        if isempty(points)
            println(io, "| $(wp.name) | $tags | 0/$(wp.N) | $order | $rate | — | — | — | — |")
            continue
        end
        i = points[argmin(wp.errors.final[points])]
        stats = wp.stats[i]
        rhs_calls = isnothing(stats) || wp.name == "PDIRK44" ? "—" :
            iszero(stats.nf + stats.nf2) ? "—" : string(stats.nf + stats.nf2)
        accepted = isnothing(stats) || :accepted_steps_unavailable in wp.tags ? "—" : string(stats.naccept)
        println(
            io,
            "| $(wp.name) | $tags | $(length(points))/$(wp.N) | $order | $rate | ",
            @sprintf("%.3g | %.3g | ", wp.errors.final[i], wp.times[i]),
            rhs_calls, " | ", accepted, " |",
        )
    end
    return String(take!(io))
end

struct MarkdownSource
    source::String
end

Base.show(io::IO, ::MIME"text/markdown", markdown::MarkdownSource) = print(io, markdown.source)

# `Markdown.parse` would turn underscores in solver tags into emphasis.
display_markdown(source) = display("text/markdown", MarkdownSource(source))

function print_work_table(wps, error_floor; kwargs...)
    return display_markdown(work_table_markdown(wps, error_floor; kwargs...))
end

function assert_work_table_markdown_fixture()
    prob = (tspan = (0.0, 1.0),)
    fixtures = (
        (name = "PDIRK44", calls = 1384, steps = 32, rhs = "—", accepted = "32", tags = ()),
        (name = "KuttaPRK2p5", calls = 0, steps = 128, rhs = "—", accepted = "128", tags = ()),
        (name = "Tsit5", calls = 441, steps = 73, rhs = "441", accepted = "73", tags = ()),
        (name = "lsoda", calls = nothing, steps = nothing, rhs = "—", accepted = "—", tags = [:adaptive]),
        (
            name = "UnavailableSteps", calls = 210, steps = 0, rhs = "210", accepted = "—",
            tags = [:adaptive, :accepted_steps_unavailable],
        ),
        (
            name = "MeasuredZero", calls = 210, steps = 0, rhs = "210", accepted = "0",
            tags = [:adaptive],
        ),
    )
    wps = map(fixtures) do fixture
        stats = isnothing(fixture.calls) ? nothing :
            (nf = fixture.calls, nf2 = 0, naccept = fixture.steps)
        return (
            name = fixture.name, prob, times = [1.0], errors = (final = [0.01],),
            tags = fixture.tags, dts = nothing, N = 1, stats = [stats],
        )
    end
    rows = split(chomp(work_table_markdown(collect(wps), 0.0)), '\n')[3:end]
    @assert length(rows) == length(fixtures)
    for (row, fixture) in zip(rows, fixtures)
        columns = strip.(split(row, '|'; keepempty = false))
        @assert columns[1] == fixture.name
        rhs, accepted = columns[8:9]
        @assert rhs == fixture.rhs "$(fixture.name) RHS calls should be $(fixture.rhs), got $rhs"
        @assert accepted == fixture.accepted "$(fixture.name) accepted steps should be $(fixture.accepted), got $accepted"
    end
    lsoda = only(row for row in rows if startswith(row, "| lsoda |"))
    lsoda_columns = strip.(split(lsoda, '|'; keepempty = false))
    @assert lsoda_columns[5] == "—" "lsoda adaptive-step rate should be —, got $(lsoda_columns[5])"

    array_probs = [(tspan = (0.0, 2.0),), (tspan = (0.0, 8.0),)]
    dts = [0.5, 0.25, 0.125]
    errors = [0.08, 0.04, 0.02]
    function fixed_array_wp(name, choice)
        interval = last(array_probs[choice].tspan) - first(array_probs[choice].tspan)
        stats = map(dts) do dt
            return (nf = 1, nf2 = 0, naccept = round(Int, interval / dt), nreject = 0)
        end
        return (
            name, prob = array_probs, times = ones(3), errors = (final = copy(errors),),
            tags = [:fixed_step], dts, N = 3, stats,
        )
    end
    function adaptive_array_wp(name)
        stats = map((20, 40, 80)) do naccept
            return (nf = naccept, nf2 = 0, naccept, nreject = 0)
        end
        return (
            name, prob = array_probs, times = ones(3), errors = (final = copy(errors),),
            tags = [:adaptive], dts = nothing, N = 3, stats,
        )
    end
    array_wps = [
        fixed_array_wp("ArrayFixedDefault", 1), fixed_array_wp("ArrayFixedExplicit", 2),
        adaptive_array_wp("ArrayAdaptiveDefault"), adaptive_array_wp("ArrayAdaptiveExplicit"),
    ]
    array_setups = [Dict{Symbol, Any}(), Dict(:prob_choice => 2), Dict{Symbol, Any}(), Dict(:prob_choice => 2)]
    array_rows = split(
        chomp(work_table_markdown(array_wps, 0.0; setups = array_setups)), '\n'
    )[3:end]
    @assert length(array_rows) == length(array_wps)
    @assert _effective_problem(array_wps[1], array_setups[1]).tspan == (0.0, 2.0)
    @assert _effective_problem(array_wps[2], array_setups[2]).tspan == (0.0, 8.0)
    for (row, order, rate) in zip(array_rows, ("1.00", "1.00", "—", "—"), ("—", "—", "1.00", "1.00"))
        columns = strip.(split(row, '|'; keepempty = false))
        @assert columns[4] == order "array-problem order should be $order, got $(columns[4])"
        @assert columns[5] == rate "array-problem adaptive rate should be $rate, got $(columns[5])"
    end
    return nothing
end
