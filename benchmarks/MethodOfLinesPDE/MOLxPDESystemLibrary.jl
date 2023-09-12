
using MethodOfLines, DomainSets, OrdinaryDiffEq, ModelingToolkit, DiffEqDevTools, LinearAlgebra,
      LinearSolve, Plots
gr()
using PDESystemLibrary

function center_uniform_grid(ex, ivs, N)
    map(ivs) do x
        xdomain = ex.domain[findfirst(d -> isequal(x, d.variables), ex.domain)]
        x => trunc(Int, N^(1 / length(ivs)))
    end
end

function edge_uniform_grid(ex, ivs, N)
    map(ivs) do x
        xdomain = ex.domain[findfirst(d -> isequal(x, d.variables), ex.domain)]
        x => trunc(Int, N^(1 / length(ivs)))
    end
end

function center_chebygrid(ex, ivs, N)
    map(ivs) do x
        xdomain = ex.domain[findfirst(d -> isequal(x, d.variables), ex.domain)]
        chebyspace(trunc(Int, N^(1 / length(ivs))), xdomain)
    end
end

function edge_chebygrid(ex, ivs, N)
    map(ivs) do x
        xdomain = ex.domain[findfirst(d -> isequal(x, d.variables), ex.domain)]
        chebyspace(trunc(Int, N^(1 / length(ivs))) - 1, xdomain)
    end
end

function uniformupwind1(ex, ivs, t, N)
    dxs = center_uniform_grid(ex, ivs, N)

    MOLFiniteDifference(dxs, t, advection_scheme=UpwindScheme())
end

function uniformupwind2(ex, ivs, t, N)
    dxs = edge_uniform_grid(ex, ivs, N)

    MOLFiniteDifference(dxs, t, advection_scheme=UpwindScheme(), grid_align=edge_align)
end

function chebyupwind1(ex, ivs, t, N)
    dxs = center_chebygrid(ex, ivs, N)

    MOLFiniteDifference(dxs, t, advection_scheme=UpwindScheme())
end

function chebyupwind2(ex, ivs, t, N)
    dxs = edge_chebygrid(ex, ivs, N)

    MOLFiniteDifference(dxs, t, advection_scheme=UpwindScheme(), grid_align=edge_align)
end


function discweno1(ex, ivs, t, N)
    dxs = center_uniform_grid(ex, ivs, N)

    MOLFiniteDifference(dxs, t, advection_scheme=WENOScheme())
end

function discweno2(ex, ivs, t, N)
    dxs = edge_uniform_grid(ex, ivs, N)

    MOLFiniteDifference(dxs, t, advection_scheme=WENOScheme(), grid_align=edge_align)
end

N = 100

for ex in PDESystemLibrary.all_systems
    try
        if ex.analytic_func === nothing
            continue
        end
        ivs = filter(x -> !isequal(Symbol(x), :t), ex.ivs)
        if length(ivs) == 0
            continue
        elseif length(ivs) == length(ex.ivs)
            # Skip nonlinear systems until I know the syntax for that
            continue

            # advection = false
            # discuu1 = uniformupwind1(ex, ivs, nothing, N)
            # discuu2 = uniformupwind2(ex, ivs, nothing, N)
            # discnu1 = chebyupwind1(ex, ivs, nothing, N)
            # discnu2 = chebyupwind2(ex, ivs, nothing, N)
            # discs = [discuu1, discuu2, discnu1, discnu2]
            # if "Advection" in ex.metadata
            #     advection = true
            #     discw1 = discweno1(ex, ivs, nothing, N)
            #     discw2 = discweno2(ex, ivs, nothing, N)
            #     push!(discs, discw1, discw2)
            # end

            # probs = map(discs) do disc
            #     discretize(ex, disc, analytic = ex.analytic_func)
            # end

            # title = "Work Precision Diagram for $(ex.name), Tags: $(ex.metadata)"
            # println("Running $title")
            # if advection
            #     dummy_appxsol = [nothing for i in 1:length(probs1)]
            #     abstols = 1.0 ./ 10.0 .^ (5:8)
            #     reltols = 1.0 ./ 10.0 .^ (1:4);
            #     setups = [Dict(:alg => solver, :prob_choice => 1),
            #         Dict(:alg => solver, :prob_choice => 2),
            #         Dict(:alg => solver, :prob_choice => 3),
            #         Dict(:alg => solver, :prob_choice => 4),
            #         Dict(:alg => solver, :prob_choice => 5),
            #         Dict(:alg => solver, :prob_choice => 6),]
            #     names = ["Uniform Upwind, center_align", "Uniform Upwind, edge_align",
            #              "Chebyshev Upwind, center_align", "Chebyshev Upwind, edge_align",
            #              "Uniform WENO, center_align", "Uniform WENO, edge_align"];

            #     wp = WorkPrecisionSet(probs, abstols, reltols, setups; names=names,
            #         save_everystep=false, appxsol = dummy_appxsol, maxiters=Int(1e5),
            #         numruns=10, wrap=Val(false))
            #     plot(wp, title=title)
            # else
            #     dummy_appxsol = [nothing for i in 1:length(probs)]
            #     abstols = 1.0 ./ 10.0 .^ (5:8)
            #     reltols = 1.0 ./ 10.0 .^ (1:4);
            #     setups = [Dict(:alg => solver, :prob_choice => 1),
            #         Dict(:alg => solver, :prob_choice => 2),
            #         Dict(:alg => solver, :prob_choice => 3),
            #         Dict(:alg => solver, :prob_choice => 4),]
            #     names = ["Uniform Upwind, center_align", "Uniform Upwind, edge_align",
            #              "Chebyshev Upwind, center_align", "Chebyshev Upwind, edge_align"];

            #     wp = WorkPrecisionSet(probs1, abstols, reltols, setups; names=names,
            #         save_everystep=false, appxsol = dummy_appxsol, maxiters=Int(1e5),
            #         numruns=10, wrap=Val(false))
            #     plot(wp, title=title)
            # end

        else
            @parameters t
            # Create discretizations
            advection = false
            discuu1 = uniformupwind1(ex, ivs, t, N)
            discuu2 = uniformupwind2(ex, ivs, t, N)
            discs = [discuu1, discuu2]
            if !("Periodic" in ex.metadata)
                discnu1 = chebyupwind1(ex, ivs, t, N)
                discnu2 = chebyupwind2(ex, ivs, t, N)
                push!(discs, discnu1, discnu2)
            end
            if "Advection" in ex.metadata
                discw1 = discweno1(ex, ivs, t, N)
                discw2 = discweno2(ex, ivs, t, N)
                push!(discs, discw1, discw2)
            end

            # Create problems
            probs = map(discs) do disc
                discretize(ex, disc, analytic = ex.analytic_func)
            end

            title = "Work Precision Diagram for $(ex.name), Tags: $(ex.metadata)"
            println("Running $title")
            dummy_appxsol = [nothing for i in 1:length(probs)]
            abstols = 1.0 ./ 10.0 .^ (5:8)
            reltols = 1.0 ./ 10.0 .^ (1:4);
            if "Advection" in ex.metadata
                solver = FBDF()
                if "Periodic" in ex.metadata
                    setups = [Dict(:alg => solver, :prob_choice => 1),
                    Dict(:alg => solver, :prob_choice => 2),
                    Dict(:alg => solver, :prob_choice => 3),
                    Dict(:alg => solver, :prob_choice => 4),]
                    names = ["Uniform Upwind, center_align", "Uniform Upwind, edge_align",
                            "Uniform WENO, center_align", "Uniform WENO, edge_align"];
                else
                    setups = [Dict(:alg => solver, :prob_choice => 1),
                        Dict(:alg => solver, :prob_choice => 2),
                        Dict(:alg => solver, :prob_choice => 3),
                        Dict(:alg => solver, :prob_choice => 4),
                        Dict(:alg => solver, :prob_choice => 5),
                        Dict(:alg => solver, :prob_choice => 6),]
                    names = ["Uniform Upwind, center_align", "Uniform Upwind, edge_align",
                            "Chebyshev Upwind, center_align", "Chebyshev Upwind, edge_align",
                            "Uniform WENO, center_align", "Uniform WENO, edge_align"];
                end
            else
                solver = TRBDF2()
                if "Periodic" in ex.metadata
                    setups = [Dict(:alg => solver, :prob_choice => 1),
                    Dict(:alg => solver, :prob_choice => 2)]
                names = ["Uniform, center_align", "Uniform, edge_align"];
                else
                    setups = [Dict(:alg => solver, :prob_choice => 1),
                        Dict(:alg => solver, :prob_choice => 2),
                        Dict(:alg => solver, :prob_choice => 3),
                        Dict(:alg => solver, :prob_choice => 4),]
                    names = ["Uniform, center_align", "Uniform, edge_align",
                             "Chebyshev, center_align", "Chebyshev, edge_align"];
                end
            end
            wp = WorkPrecisionSet(probs, abstols, reltols, setups; names=names,
                save_everystep=false, appxsol = dummy_appxsol, maxiters=Int(1e5),
                numruns=10, wrap=Val(false))
            p = plot(wp, title=title)
            display(ex.name)
            display(ex)
            display(p)

        end
    catch e
        println("Failed on $(ex.name):")
        println(e)
    end
end 
