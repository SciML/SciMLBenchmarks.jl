
using StochasticDiffEq, DiffEqNoiseProcess, LinearAlgebra, Statistics

function generate_stiff_stoch_heat(D=1,k=1;N = 100, t_end = 3.0, adaptivealg = :RSwM3)
    A = Array(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
    dx = 1/N
    A = D/(dx^2) * A
    function f(du,u,p,t)
        mul!(du,A,u)
    end
    #=
    function f(::Type{Val{:analytic}},u0,p,t,W)
        exp((A-k/2)*t+W*I)*u0 # no -k/2 for Strat
    end
    =#
    function g(du,u,p,t)
        @. du = k*u
    end
    SDEProblem(f,g,ones(N),(0.0,t_end),noise=WienerProcess(0.0,0.0,0.0,rswm=RSWM(adaptivealg=adaptivealg)))
end

N = 100
D = 1; k = 1
    A = Array(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
    dx = 1/N
    A = D/(dx^2) * A;


prob = generate_stiff_stoch_heat(1.0,1.0)
@time sol = solve(prob,SRIW1(),progress=true,abstol=1e-6,reltol=1e-6);


@time sol = solve(generate_stiff_stoch_heat(1.0,1.0),SRIW1());


@time sol = solve(generate_stiff_stoch_heat(1.0,1.0),SRIW1(),progress=true,adaptive=false,dt=0.00005);


@time sol = solve(generate_stiff_stoch_heat(1.0,1.0),EM(),progress=true,adaptive=false,dt=0.00005);


@time sol = solve(generate_stiff_stoch_heat(1.0,1.0),ImplicitRKMil(),progress=true,dt=0.1);


@time sol = solve(generate_stiff_stoch_heat(1.0,1.0),ImplicitRKMil(),progress=true,dt=0.01);


@time sol = solve(generate_stiff_stoch_heat(1.0,1.0),ImplicitRKMil(),progress=true,dt=0.001);


@time sol = solve(generate_stiff_stoch_heat(1.0,1.0),ImplicitEM(),progress=true,dt=0.001);


function simple_error(alg;kwargs...)
    sol = solve(generate_stiff_stoch_heat(1.0,1.0,t_end=0.25),alg;kwargs...);
    sum(abs2,sol[end] - exp(A*sol.t[end]+sol.W[end]*I)*prob.u0)
end


mean(simple_error(EulerHeun(),dt=0.00005) for i in 1:400)


mean(simple_error(ImplicitRKMil(interpretation=:Stratanovich),dt=0.1) for i in 1:400)


mean(simple_error(ImplicitRKMil(interpretation=:Stratanovich),dt=0.01) for i in 1:400)


mean(simple_error(ImplicitRKMil(interpretation=:Stratanovich),dt=0.001) for i in 1:400)


mean(simple_error(ImplicitEulerHeun(),dt=0.001) for i in 1:400)


mean(simple_error(ImplicitEulerHeun(),dt=0.01) for i in 1:400)


mean(simple_error(ImplicitEulerHeun(),dt=0.1) for i in 1:400)


sol = solve(generate_stiff_stoch_heat(1.0,1.0,adaptivealg=:RSwM1),SRIW1());


using DiffEqBenchmarks
DiffEqBenchmarks.bench_footer(WEAVE_ARGS[:folder],WEAVE_ARGS[:file])

