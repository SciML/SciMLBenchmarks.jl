
using ReactionNetworkImporters, OrdinaryDiffEq, DiffEqBiological,
      Sundials, Plots, DiffEqDevTools, ODEInterface, ODEInterfaceDiffEq,
      LSODA, TimerOutputs

gr()
prnbng = loadrxnetwork(BNGNetwork(), "BNGRepressilator",
                       joinpath(dirname(pathof(ReactionNetworkImporters)),"..","data","bcr","bcr.net"))

rn = deepcopy(prnbng.rn)
addodes!(rn; build_jac=false, build_symfuncs=false, build_paramjac=false)
tf = 100000.0
oprob = ODEProblem(rn, prnbng.u₀, (0.,tf), prnbng.p);

densejac_rn = deepcopy(prnbng.rn)
# zeroout_jac=true is needed to keep the generated expressions from being too big for the compiler
addodes!(densejac_rn; build_jac=true, zeroout_jac = true, sparse_jac = false, build_symfuncs=false, build_paramjac=false)
densejacprob = ODEProblem(densejac_rn, prnbng.u₀, (0.,tf), prnbng.p);

sparsejac_rn = deepcopy(prnbng.rn)
addodes!(sparsejac_rn; build_jac=true, sparse_jac = true, build_symfuncs=false, build_paramjac=false)
sparsejacprob = ODEProblem(sparsejac_rn, prnbng.u₀, (0.,tf), prnbng.p);


@show numspecies(rn) # Number of ODEs
@show numreactions(rn) # Apprx. number of terms in the ODE
@show numparams(rn) # Number of Parameters


const to = TimerOutput()
u₀ = prnbng.u₀
u = copy(u₀);
du = similar(u);
p = prnbng.p
@timeit to "ODERHS Eval1" rn.f(du,u,p,0.)
@timeit to "ODERHS Eval2" rn.f(du,u,p,0.)
sparsejac_rn.f(du,u,p,0.)

J = zeros(length(u),length(u))
@timeit to "DenseJac Eval1" densejac_rn.jac(J,u,p,0.)
@timeit to "DenseJac Eval2" densejac_rn.jac(J,u,p,0.)

Js = similar(sparsejac_rn.odefun.jac_prototype)
@timeit to "SparseJac Eval1" sparsejac_rn.jac(Js,u,p,0.)
@timeit to "SparseJac Eval2" sparsejac_rn.jac(Js,u,p,0.)
show(to)


sol = solve(oprob, CVODE_BDF(), saveat=tf/1000., reltol=1e-5, abstol=1e-5)
plot(sol,legend=false, fmt=:png)


@time sol = solve(oprob,CVODE_BDF(),abstol=1/10^12,reltol=1/10^12)
test_sol = TestSolution(sol)


abstols = 1.0 ./ 10.0 .^ (5:8)
reltols = 1.0 ./ 10.0 .^ (5:8);
setups = [
          #Dict(:alg=>Rosenbrock23(autodiff=false)),
          Dict(:alg=>TRBDF2(autodiff=false)),
          Dict(:alg=>CVODE_BDF()),
          #Dict(:alg=>rodas()),
          #Dict(:alg=>radau()),
          #Dict(:alg=>Rodas4(autodiff=false)),
          #Dict(:alg=>Rodas5(autodiff=false)),
          Dict(:alg=>KenCarp4(autodiff=false)),
          #Dict(:alg=>RadauIIA5(autodiff=false)),
          #Dict(:alg=>lsoda()),
          ]


wp = WorkPrecisionSet(oprob,abstols,reltols,setups;error_estimate=:l2,
                      saveat=tf/10000.,appxsol=test_sol,maxiters=Int(1e5),numruns=1)
plot(wp)


wp = WorkPrecisionSet(densejacprob,abstols,reltols,setups;error_estimate=:l2,
                      saveat=tf/10000.,appxsol=test_sol,maxiters=Int(1e5),numruns=1)
plot(wp)


setups = [
          #Dict(:alg=>Rosenbrock23(autodiff=false)),
          Dict(:alg=>TRBDF2(autodiff=false)),
          #Dict(:alg=>CVODE_BDF()),
          #Dict(:alg=>rodas()),
          #Dict(:alg=>radau()),
          #Dict(:alg=>Rodas4(autodiff=false)),
          #Dict(:alg=>Rodas5(autodiff=false)),
          Dict(:alg=>KenCarp4(autodiff=false)),
          #Dict(:alg=>RadauIIA5(autodiff=false)),
          #Dict(:alg=>lsoda()),
          ]
wp = WorkPrecisionSet(sparsejacprob,abstols,reltols,setups;error_estimate=:l2,
                      saveat=tf/10000.,appxsol=test_sol,maxiters=Int(1e5),numruns=1)
plot(wp)

