
using MethodOfLines, DomainSets, OrdinaryDiffEq, ModelingToolkit, DiffEqDevTools, LinearAlgebra,
      LinearSolve, Plots

# pdesys1 has Dirichlet BCs, pdesys2 has Neumann BCs
const N = 30

@parameters x t
@variables u(..)
Dx = Differential(x)
Dt = Differential(t)
x_min = 0.0
x_max = 1.0
t_min = 0.0
t_max = 20.0

solver = FBDF()

analytic_u(p, t, x) = x / (t + 1)

analytic = [u(t, x) ~ analytic_u([], t, x)]

eq = Dt(u(t, x)) ~ -u(t, x) * Dx(u(t, x))

bcs1 = [u(0, x) ~ x,
    u(t, x_min) ~ analytic_u([], t, x_min),
    u(t, x_max) ~ analytic_u([], t, x_max)]

bcs2 = [u(0, x) ~ x,
    Dx(u(t, x_min)) ~ 1 / (t + 1),
    Dx(u(t, x_max)) ~ 1 / (t + 1)]

domains = [t ∈ Interval(t_min, t_max),
    x ∈ Interval(x_min, x_max)]

@named pdesys1 = PDESystem(eq, bcs1, domains, [t, x], [u(t, x)], analytic=analytic)
@named pdesys2 = PDESystem(eq, bcs2, domains, [t, x], [u(t, x)], analytic=analytic)

discupwind1 = MOLFiniteDifference([x => N], t, advection_scheme=UpwindScheme())
discupwind2 = MOLFiniteDifference([x => N-1], t, advection_scheme=UpwindScheme(), grid_align=edge_align)

discweno1 = MOLFiniteDifference([x => N], t, advection_scheme=WENOScheme())
discweno2 = MOLFiniteDifference([x => N-1], t, advection_scheme=WENOScheme(), grid_align=edge_align)

gridnu1 = chebyspace(N, domains[2])
gridnu2 = chebyspace(N - 1, domains[2])

discnu1 = MOLFiniteDifference([gridnu1], t, advection_scheme=UpwindScheme())
discnu2 = MOLFiniteDifference([gridnu2], t, advection_scheme=UpwindScheme(), grid_align=edge_align)

probupwind1 = discretize(pdesys1, discupwind1; analytic=pdesys1.analytic_func)
probupwind2 = discretize(pdesys1, discupwind2; analytic=pdesys1.analytic_func)

probweno1 = discretize(pdesys1, discweno1; analytic=pdesys1.analytic_func)
probweno2 = discretize(pdesys1, discweno2; analytic=pdesys1.analytic_func)

probnu1 = discretize(pdesys1, discnu1; analytic=pdesys1.analytic_func)
probnu2 = discretize(pdesys1, discnu2; analytic=pdesys1.analytic_func)

probs1 = [probupwind1, probupwind2, probnu1, probnu2, probweno1, probweno2]

dummy_appxsol = [nothing for i in 1:length(probs1)]
abstols = 1.0 ./ 10.0 .^ (5:8)
reltols = 1.0 ./ 10.0 .^ (1:4);
setups = [Dict(:alg => solver, :prob_choice => 1),
    Dict(:alg => solver, :prob_choice => 2),
    Dict(:alg => solver, :prob_choice => 3),
    Dict(:alg => solver, :prob_choice => 4),
    Dict(:alg => solver, :prob_choice => 5),
    Dict(:alg => solver, :prob_choice => 6),]
names = ["Uniform Upwind, center_align", "Uniform Upwind, edge_align", "Nonuniform Upwind, center_align",
         "Nonuniform Upwind, edge_align", "WENO, center_align", "WENO, edge_align"];

wp = WorkPrecisionSet(probs1, abstols, reltols, setups; names=names,
    save_everystep=false, appxsol = dummy_appxsol, maxiters=Int(1e5),
    numruns=10, wrap=Val(false))
plot(wp)

probupwind1 = discretize(pdesys2, discupwind1; analytic=pdesys2.analytic_func)
probupwind2 = discretize(pdesys2, discupwind2; analytic=pdesys2.analytic_func)

probweno1 = discretize(pdesys2, discweno1; analytic=pdesys2.analytic_func)
probweno2 = discretize(pdesys2, discweno2; analytic=pdesys2.analytic_func)

probnu1 = discretize(pdesys2, discnu1; analytic=pdesys2.analytic_func)
probnu2 = discretize(pdesys2, discnu2; analytic=pdesys2.analytic_func)

probs2 = [probupwind1, probupwind2, probnu1, probnu2, probweno1, probweno2]

abstols = 1.0 ./ 10.0 .^ (5:8)
reltols = 1.0 ./ 10.0 .^ (1:4);
setups = [Dict(:alg => solver, :prob_choice => 1),
          Dict(:alg => solver, :prob_choice => 2),
          Dict(:alg => solver, :prob_choice => 3),
          Dict(:alg => solver, :prob_choice => 4),
          Dict(:alg => solver, :prob_choice => 5),
          Dict(:alg => solver, :prob_choice => 6),]
names = ["Uniform Upwind, center_align", "Uniform Upwind, edge_align", "Nonuniform Upwind, center_align",
         "Nonuniform Upwind, edge_align", "WENO, center_align", "WENO, edge_align"];

wp = WorkPrecisionSet(probs2, abstols, reltols, setups; names=names,
                      save_everystep=false, appxsol = dummy_appxsol, maxiters=Int(1e5),
                      numruns=10, wrap=Val(false))
plot(wp)
