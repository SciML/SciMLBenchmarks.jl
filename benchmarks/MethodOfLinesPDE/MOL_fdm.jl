"""
---
title: Burgers FDM Work-Precision Diagrams with Various MethodOfLines Methods
author: Alex Jones
---
"""
# This benchmark is for the MethodOfLines package, which is an automatic PDE discretization package.
# It is concerned with comparing the performance of various discretization methods for the Burgers equation.

using MethodOfLines, DomainSets, OrdinaryDiffEq, ModelingToolkit, DiffEqDevTools, LinearAlgebra,
      LinearSolve, Plots
# Note: ModelingToolkit is held at v8.22.1, as somewhere between this and v8.29.1 a bug was introduced that breaks this benchmark.
# See MethodOfLines/#194 for more information.

#Here is the burgers equation with a Dirichlet and Neumann boundary conditions,
# pdesys1 has Dirichlet BCs, pdesys2 has Neumann BCs


const N = 30

@parameters x t
@variables u(..)
Dx = Differential(x)
Dt = Differential(t)
x_min = 0.0
x_max = 1.0
t_min = 0.0
t_max = 6.0

solver = FBDF()

analytic_u(t, x) = x / (t + 1)


eq = Dt(u(t, x)) ~ -u(t, x) * Dx(u(t, x))

bcs1 = [u(0, x) ~ x,
    u(t, x_min) ~ analytic_u(t, x_min),
    u(t, x_max) ~ analytic_u(t, x_max)]

bcs2 = [u(0, x) ~ x,
    Dx(u(t, x_min)) ~ 1 / (t + 1),
    Dx(u(t, x_max)) ~ 1 / (t + 1)]

domains = [t ∈ Interval(t_min, t_max),
    x ∈ Interval(x_min, x_max)]

@named pdesys1 = PDESystem(eq, bcs1, domains, [t, x], [u(t, x)])
@named pdesys2 = PDESystem(eq, bcs2, domains, [t, x], [u(t, x)])

# Here is a uniform discretization with the Upwind scheme:

discupwind1 = MOLFiniteDifference([x => N], t, advection_scheme=UpwindScheme())
discupwind2 = MOLFiniteDifference([x => N-1], t, advection_scheme=UpwindScheme(), grid_align=edge_align)


# Here is a uniform discretization with the WENO scheme:

discweno1 = MOLFiniteDifference([x => N], t, advection_scheme=WENOScheme())
discweno2 = MOLFiniteDifference([x => N-1], t, advection_scheme=WENOScheme(), grid_align=edge_align)

# Grid is the same for the above two types of discretization
griduniform1 = get_discrete(pdesys1, discupwind1)
griduniform2 = get_discrete(pdesys1, discupwind2)

griduniform1 = griduniform1[x]
griduniform2 = griduniform2[x]

# Here is a non-uniform discretization with the Upwind scheme, using tanh (nonuniform WENO is not implemented yet):
gridf(x) = tanh.(x) ./ 2 .+ 0.5
gridnu1 = gridf(vcat(-Inf, range(-3.0, 3.0, length=N-2), Inf))
gridnu2 = gridf(vcat(-Inf, range(-3.0, 3.0, length=N - 3), Inf))

discnu1 = MOLFiniteDifference([x => gridnu1], t, advection_scheme=UpwindScheme())
discnu2 = MOLFiniteDifference([x => gridnu2], t, advection_scheme=UpwindScheme(), grid_align=edge_align)

gridnu1 = get_discrete(pdesys1, discnu1)
gridnu2 = get_discrete(pdesys1, discnu2)

gridnu1 = gridnu1[x]
gridnu2 = gridnu2[x]

# Here are the problems for pdesys1:
probupwind1 = discretize(pdesys1, discupwind1)
probupwind2 = discretize(pdesys1, discupwind2)

probweno1 = discretize(pdesys1, discweno1)
probweno2 = discretize(pdesys1, discweno2)

probnu1 = discretize(pdesys1, discnu1)
probnu2 = discretize(pdesys1, discnu2)

probs1 = [probupwind1, probupwind2, probnu1, probnu2, probweno1, probweno2]

# Here are the reference solutions for pdesys1:
solupwind1 = solve(probupwind1, solver, abstol=1 / 10^14, reltol=1 / 10^14)
solupwind2 = solve(probupwind2, solver, abstol=1 / 10^14, reltol=1 / 10^14)

solweno1 = solve(probweno1, solver, abstol=1 / 10^14, reltol=1 / 10^14)
solweno2 = solve(probweno2, solver, abstol=1 / 10^14, reltol=1 / 10^14)

solnu1 = solve(probnu1, solver, abstol=1 / 10^14, reltol=1 / 10^14)
solnu2 = solve(probnu2, solver, abstol=1 / 10^14, reltol=1 / 10^14)

test_sols1 = getfield.([solupwind1, solupwind2, solnu1, solnu2, solweno1, solweno2], (:original_sol,))

# Here are the analytic solutions for pdesys1:
solupwind1_analytic = [analytic_u(t_, x_) for t_ in solupwind1.t, x_ in griduniform1]
solupwind2_analytic = [analytic_u(t_, x_) for t_ in solupwind2.t, x_ in griduniform2]

solweno1_analytic = [analytic_u(t_, x_) for t_ in solweno1.t, x_ in griduniform1]
solweno2_analytic = [analytic_u(t_, x_) for t_ in solweno2.t, x_ in griduniform2]

solnu1_analytic = [analytic_u(t_, x_) for t_ in solnu1.t, x_ in gridnu1]
solnu2_analytic = [analytic_u(t_, x_) for t_ in solnu2.t, x_ in gridnu2]

# Here are the absolute errors for pdesys1:
errupwind1 = norm(solupwind1_analytic .- solupwind1[u(t, x)])
errupwind2 = norm(solupwind2_analytic .- solupwind2[u(t, x)])

errweno1 = norm(solweno1_analytic .- solweno1[u(t, x)])
errweno2 = norm(solweno2_analytic .- solweno2[u(t, x)])

errnu1 = norm(solnu1_analytic .- solnu1[u(t, x)])
errnu2 = norm(solnu2_analytic .- solnu2[u(t, x)])

atols1 = [errupwind1, errupwind2, errnu1, errnu2, errweno1, errweno2]

# Here are the relative errors for pdesys1:
errupwind1 = norm(solupwind1_analytic .- solupwind1[u(t, x)]) / norm(solupwind1_analytic)
errupwind2 = norm(solupwind2_analytic .- solupwind2[u(t, x)]) / norm(solupwind2_analytic)

errweno1 = norm(solweno1_analytic .- solweno1[u(t, x)]) / norm(solweno1_analytic)
errweno2 = norm(solweno2_analytic .- solweno2[u(t, x)]) / norm(solweno2_analytic)

errnu1 = norm(solnu1_analytic .- solnu1[u(t, x)]) / norm(solnu1_analytic)
errnu2 = norm(solnu2_analytic .- solnu2[u(t, x)]) / norm(solnu2_analytic)

rtols1 = [errupwind1, errupwind2, errnu1, errnu2, errweno1, errweno2]

# Work-Precision Plot for Burgers Equation, Dirichlet BCs

abstols = 1.0 ./ 10.0 .^ (5:8)
reltols = 1.0 ./ 10.0 .^ (1:4);
setups = [Dict(:alg => solver),
    Dict(:alg => solver, :prob_choice => 2),
    Dict(:alg => solver, :prob_choice => 3),
    Dict(:alg => solver, :prob_choice => 4),
    Dict(:alg => solver, :prob_choice => 5),
    Dict(:alg => solver, :prob_choice => 6),]
names = ["Uniform Upwind, center_align", "Uniform Upwind, edge_align", "Nonuniform Upwind, center_align",
         "Nonuniform Upwind, edge_align", "WENO, center_align", "WENO, edge_align"];

wp = WorkPrecisionSet(probs1, abstols, reltols, setups; names=names,
    save_everystep=false, appxsol=test_sols1, maxiters=Int(1e5),
    numruns=10, wrap=Val(false))
plot(wp)

# Here are the problems for pdesys2:
probupwind1 = discretize(pdesys2, discupwind1)
probupwind2 = discretize(pdesys2, discupwind2)

probweno1 = discretize(pdesys2, discweno1)
probweno2 = discretize(pdesys2, discweno2)

probnu1 = discretize(pdesys2, discnu1)
probnu2 = discretize(pdesys2, discnu2)

probs2 = [probupwind1, probupwind2, probnu1, probnu2, probweno1, probweno2]
# Here are the reference solutions for pdesys2:
solupwind1 = solve(probupwind1, solver)
solupwind2 = solve(probupwind2, solver)

solweno1 = solve(probweno1, solver)
solweno2 = solve(probweno2, solver)

solnu1 = solve(probnu1, solver)
solnu2 = solve(probnu2, solver)

test_sols2 = getfield.([solupwind1, solupwind2, solnu1, solnu2, solweno1, solweno2], (:original_sol,))

# Here are the analytic solutions for pdesys2:
solupwind1_analytic = [analytic_u(t_, x_) for t_ in solupwind1.t, x_ in griduniform1]
solupwind2_analytic = [analytic_u(t_, x_) for t_ in solupwind2.t, x_ in griduniform2]

solweno1_analytic = [analytic_u(t_, x_) for t_ in solweno1.t, x_ in griduniform1]
solweno2_analytic = [analytic_u(t_, x_) for t_ in solweno2.t, x_ in griduniform2]

solnu1_analytic = [analytic_u(t_, x_) for t_ in solnu1.t, x_ in gridnu1]
solnu2_analytic = [analytic_u(t_, x_) for t_ in solnu2.t, x_ in gridnu2]

# Here are the absolute errors for pdesys2:
errupwind1 = norm(solupwind1_analytic .- solupwind1[u(t, x)])
errupwind2 = norm(solupwind2_analytic .- solupwind2[u(t, x)])

errweno1 = norm(solweno1_analytic .- solweno1[u(t, x)])
errweno2 = norm(solweno2_analytic .- solweno2[u(t, x)])

errnu1 = norm(solnu1_analytic .- solnu1[u(t, x)])
errnu2 = norm(solnu2_analytic .- solnu2[u(t, x)])

atols2 = [errupwind1, errupwind2, errnu1, errnu2, errweno1, errweno2]

# Here are the relative errors for pdesys2:
errupwind1 = norm(solupwind1_analytic .- solupwind1[u(t, x)]) / norm(solupwind1_analytic)
errupwind2 = norm(solupwind2_analytic .- solupwind2[u(t, x)]) / norm(solupwind2_analytic)

errweno1 = norm(solweno1_analytic .- solweno1[u(t, x)]) / norm(solweno1_analytic)
errweno2 = norm(solweno2_analytic .- solweno2[u(t, x)]) / norm(solweno2_analytic)

errnu1 = norm(solnu1_analytic .- solnu1[u(t, x)]) / norm(solnu1_analytic)
errnu2 = norm(solnu2_analytic .- solnu2[u(t, x)]) / norm(solnu2_analytic)

rtols2 = [errupwind1, errupwind2, errnu1, errnu2, errweno1, errweno2]

# Work-Precision Plot for Burgers Equation, Neumann BCs

abstols = 1.0 ./ 10.0 .^ (5:8)
reltols = 1.0 ./ 10.0 .^ (1:4);
setups = [Dict(:alg => solver),
          Dict(:alg => solver, :prob_choice => 2),
          Dict(:alg => solver, :prob_choice => 3),
          Dict(:alg => solver, :prob_choice => 4),
          Dict(:alg => solver, :prob_choice => 5),
          Dict(:alg => solver, :prob_choice => 6),]
names = ["Uniform Upwind, center_align", "Uniform Upwind, edge_align", "Nonuniform Upwind, center_align",
         "Nonuniform Upwind, edge_align", "WENO, center_align", "WENO, edge_align"];

wp = WorkPrecisionSet(probs2, abstols, reltols, setups; names=names,
                      save_everystep=false, appxsol=test_sols2, maxiters=Int(1e5),
                      numruns=10, wrap=Val(false))
plot(wp)
