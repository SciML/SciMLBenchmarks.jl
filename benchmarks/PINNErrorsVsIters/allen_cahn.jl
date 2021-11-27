using NeuralPDE
using Quadrature, Cubature, Cuba
using Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using QuasiMonteCarlo

print("Precompiling Done")


# 4 spatial dimensions
function allen_cahn_time(strategy, minimizer, maxIters)

    ##  DECLARATIONS
    @parameters  t x1 x2 x3 x4
    @variables   u(..)

    Dt = Differential(t)

    Dxx1 = Differential(x1)^2
    Dxx2 = Differential(x2)^2
    Dxx3 = Differential(x3)^2
    Dxx4 = Differential(x4)^2


    # Discretization
    tmax         = 1.0
    x1width      = 1.0
    x2width      = 1.0
    x3width      = 1.0
    x4width      = 1.0

    tMeshNum     = 10
    x1MeshNum    = 10
    x2MeshNum    = 10
    x3MeshNum    = 10
    x4MeshNum    = 10

    dt   = tmax/tMeshNum
    dx1  = x1width/x1MeshNum
    dx2  = x2width/x2MeshNum
    dx3  = x3width/x3MeshNum
    dx4  = x4width/x4MeshNum

    domains = [t ∈ IntervalDomain(0.0,tmax),
               x1 ∈ IntervalDomain(0.0,x1width),
               x2 ∈ IntervalDomain(0.0,x2width),
               x3 ∈ IntervalDomain(0.0,x3width),
               x4 ∈ IntervalDomain(0.0,x4width)]

    ts  = 0.0 : dt : tmax
    x1s = 0.0 : dx1 : x1width
    x2s = 0.0 : dx2 : x2width
    x3s = 0.0 : dx3 : x3width
    x4s = 0.0 : dx4 : x4width

    # Operators
    Δu = Dxx1(u(t,x1,x2,x3,x4)) + Dxx2(u(t,x1,x2,x3,x4)) + Dxx3(u(t,x1,x2,x3,x4)) + Dxx4(u(t,x1,x2,x3,x4)) # Laplacian


    # Equation
    eq = Dt(u(t,x1,x2,x3,x4)) - Δu - u(t,x1,x2,x3,x4) + u(t,x1,x2,x3,x4)*u(t,x1,x2,x3,x4)*u(t,x1,x2,x3,x4) ~ 0  #ALLEN CAHN EQUATION

    initialCondition =  1/(2 + 0.4 * (x1*x1 + x2*x2 + x3*x3 + x4*x4)) # see PNAS paper

    bcs = [u(0,x1,x2,x3,x4) ~ initialCondition]  #from literature

    ## NEURAL NETWORK
    n = 20   #neuron number

    chain = FastChain(FastDense(5,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1))   #Neural network from Flux library

    discretization = NeuralPDE.PhysicsInformedNN(chain, strategy)

    indvars = [t,x1,x2,x3,x4]   #phisically independent variables
    depvars = [u]       #dependent (target) variable

    dim = length(domains)

    losses = []

    cb = function (p,l)
        append!(losses, l)    #loss function handling
        println(length(losses), " Current loss is: $l")
        return false
    end

    pde_system = PDESystem(eq, bcs, domains, indvars, depvars)
    prob = discretize(pde_system, discretization)

    res = GalacticOptim.solve(prob, minimizer; cb = cb, maxiters = 1) #allow_f_increase = false,

    init_θ = res.minimizer

    discretization2 = NeuralPDE.PhysicsInformedNN(chain, strategy; init_params = init_θ)   #Second learning phase, lower learning parameter
    init_θ == discretization2.init_params
    prob2 = NeuralPDE.discretize(pde_system,discretization2)

    t_0 = time_ns()

    res2 = GalacticOptim.solve(prob2, minimizer, cb = cb, maxiters=maxIters)

    t_f = time_ns()
    training_time = (t_f - t_0)/10^9

    return [training_time, init_θ]
end


function allen_cahn(strategy, minimizer, maxIters, params)

    ##  DECLARATIONS
    @parameters  t x1 x2 x3 x4
    @variables   u(..)

    Dt = Differential(t)

    Dxx1 = Differential(x1)^2
    Dxx2 = Differential(x2)^2
    Dxx3 = Differential(x3)^2
    Dxx4 = Differential(x4)^2


    # Discretization
    tmax         = 1.0
    x1width      = 1.0
    x2width      = 1.0
    x3width      = 1.0
    x4width      = 1.0

    tMeshNum     = 10
    x1MeshNum    = 10
    x2MeshNum    = 10
    x3MeshNum    = 10
    x4MeshNum    = 10

    dt   = tmax/tMeshNum
    dx1  = x1width/x1MeshNum
    dx2  = x2width/x2MeshNum
    dx3  = x3width/x3MeshNum
    dx4  = x4width/x4MeshNum

    domains = [t ∈ IntervalDomain(0.0,tmax),
               x1 ∈ IntervalDomain(0.0,x1width),
               x2 ∈ IntervalDomain(0.0,x2width),
               x3 ∈ IntervalDomain(0.0,x3width),
               x4 ∈ IntervalDomain(0.0,x4width)]

    ts  = 0.0 : dt : tmax
    x1s = 0.0 : dx1 : x1width
    x2s = 0.0 : dx2 : x2width
    x3s = 0.0 : dx3 : x3width
    x4s = 0.0 : dx4 : x4width

    # Operators
    Δu = Dxx1(u(t,x1,x2,x3,x4)) + Dxx2(u(t,x1,x2,x3,x4)) + Dxx3(u(t,x1,x2,x3,x4)) + Dxx4(u(t,x1,x2,x3,x4)) # Laplacian


    # Equation
    eq = Dt(u(t,x1,x2,x3,x4)) - Δu - u(t,x1,x2,x3,x4) + u(t,x1,x2,x3,x4)*u(t,x1,x2,x3,x4)*u(t,x1,x2,x3,x4) ~ 0  #ALLEN CAHN EQUATION

    initialCondition =  1/(2 + 0.4 * (x1*x1 + x2*x2 + x3*x3 + x4*x4)) # see PNAS paper

    bcs = [u(0,x1,x2,x3,x4) ~ initialCondition]  #from literature

    ## NEURAL NETWORK
    n = 20   #neuron number

    chain = FastChain(FastDense(5,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1))   #Neural network from Flux library

    indvars = [t,x1,x2,x3,x4]   #phisically independent variables
    depvars = [u]       #dependent (target) variable

    dim = length(domains)

    losses = []
    error = []

    phi = NeuralPDE.get_phi(chain)
    derivative = NeuralPDE.get_numeric_derivative()
    initθ = DiffEqFlux.initial_params(chain)

    dim = length(domains)

    error_strategy = NeuralPDE.QuadratureTraining(quadrature_alg=CubatureJLh(),reltol= 1e-4,abstol= 1e-3,maxiters=10, batch=10)
    _pde_loss_function = NeuralPDE.build_loss_function(eq,indvars,depvars,phi,derivative,chain,initθ,error_strategy)
    _pde_loss_function(rand(dim,10), initθ)

    bc_indvars = NeuralPDE.get_argument(bcs,indvars,depvars)
    _bc_loss_functions = [NeuralPDE.build_loss_function(bc,indvars,depvars, phi, derivative,chain,initθ,error_strategy,
                                                  bc_indvars = bc_indvar) for (bc,bc_indvar) in zip(bcs,bc_indvars)]
    map(loss_f -> loss_f(rand(dim-1,10), initθ),_bc_loss_functions)

    dx = 0.1
    train_sets = NeuralPDE.generate_training_sets(domains,dx,[eq],bcs,indvars,depvars)
    pde_train_set,bcs_train_set = train_sets
    pde_bounds, bcs_bounds = NeuralPDE.get_bounds(domains,bcs,indvars,depvars,error_strategy)



    pde_loss_function = NeuralPDE.get_loss_function([_pde_loss_function],
                                                    pde_bounds,
                                                    error_strategy;
                                                    τ = 1/100)

    pde_loss_function(initθ)
    error_strategy = NeuralPDE.QuadratureTraining(quadrature_alg=CubatureJLh(),reltol= 1e-2,abstol= 1e-1,maxiters=5, batch=100)
    bc_loss_function = NeuralPDE.get_loss_function(_bc_loss_functions,
                                                   bcs_bounds,
                                                   error_strategy;
                                                   τ = 1/40)
    bc_loss_function(initθ)

    function loss_function_(θ,p)
        return pde_loss_function(θ) + bc_loss_function(θ)
    end

    cb_ = function (p,l)
        append!(losses, l)
        append!(error, pde_loss_function(p) + bc_loss_function(p))
        println(length(losses), " Current loss is: ", l, " uniform error is, ",  pde_loss_function(p) + bc_loss_function(p))
        return false
    end

    pde_system = PDESystem(eq, bcs, domains, indvars, depvars)

    discretization = NeuralPDE.PhysicsInformedNN(chain, strategy; init_params = params)
    prob = NeuralPDE.discretize(pde_system,discretization)
    res = GalacticOptim.solve(prob, minimizer, cb = cb_, maxiters=maxIters)

    phi = discretization.phi

    final_params = res.minimizer

    ##SAVE PARAMETERS
    outfile = "final_params_allen_cahn.txt"
    open(outfile, "w") do f
      for i in final_params
        println(f, i)
      end
    end


    # Model prediction
    domain = [ts,x1s,x2s,x3s,x4s]

    u_predict = [reshape([first(phi([t,x1,x2,x3,x4],res.minimizer)) for x1 in x1s for x2 in x2s for x3 in x3s for x4 in x4s], (length(x1s),length(x2s), length(x3s),length(x4s))) for t in ts]  #matrix of model's prediction

    return [error, u_predict, u_predict, domain]
end

#losses, u_predict, u_predict, domain, training_time = allen_cahn(NeuralPDE.QuadratureTraining(algorithm = CubaCuhre(), reltol = 1e-8, abstol = 1e-8, maxiters = 100), GalacticOptim.ADAM(0.01), 500)

## Numerical Part
#=
"""
# Parameters, variables, and derivatives
@parameters  t x1 x2 x3 x4
@variables   u(..)
@derivatives Dt'~t
@derivatives Dxx1''~x1
@derivatives Dxx2''~x2
@derivatives Dxx3''~x3
@derivatives Dxx4''~x4
# Operators
Δu = Dxx1(u(t,x1,x2,x3,x4)) + Dxx2(u(t,x1,x2,x3,x4)) + Dxx3(u(t,x1,x2,x3,x4)) + Dxx4(u(t,x1,x2,x3,x4)) # Laplacian
# Equation
eq = Dt(u(t,x1,x2,x3,x4)) - Δu - u(t,x1,x2,x3,x4) + u(t,x1,x2,x3,x4)*u(t,x1,x2,x3,x4)*u(t,x1,x2,x3,x4) ~ 0  #LEVEL SET EQUATION
initialCondition =  1/(2 + 0.4 * (x1*x1 + x2*x2 + x3*x3 + x4*x4)) # see PNAS paper
bcs = [u(0,x1,x2,x3,x4) ~ initialCondition,
       u(t,0,0,0,0) ~ 0.5,
       u(t,1,1,1,1) ~ 1/2.4]  #from literature
# Discretization
tmax         = 1.0
x1width      = 1.0
x2width      = 1.0
x3width      = 1.0
x4width      = 1.0
tMeshNum     = 10
x1MeshNum    = 10
x2MeshNum    = 10
x3MeshNum    = 10
x4MeshNum    = 10
dt   = tmax/tMeshNum
dx1  = x1width/x1MeshNum
dx2  = x2width/x2MeshNum
dx3  = x3width/x3MeshNum
dx4  = x4width/x4MeshNum
# Space and time domains
domains = [t ∈ IntervalDomain(0.0,tmax),
           x1 ∈ IntervalDomain(0.0,x1width),
           x2 ∈ IntervalDomain(0.0,x2width),
           x3 ∈ IntervalDomain(0.0,x3width),
           x4 ∈ IntervalDomain(0.0,x4width)]
# PDE system
pdesys = PDESystem(eq,bcs,domains,[t,x1,x2,x3,x4],[u])
# Method of lines discretization
order = 2
discretization = MOLFiniteDifference([dx1,dx2,dx3,dx4],order)
# Convert the PDE problem into an ODE problem
prob = DiffEqOperators.discretize(pdesys,discretization)
# Solve ODE problem
using OrdinaryDiffEq
sol = solve(prob,Tsit5(),saveat=0.2)
# Plot results and compare with exact solution
x = prob.space[2]
t = sol.t
"""
## Manually Constructing 4D Laplacian
## ref: https://github.com/SciML/DiffEqOperators.jl/blob/master/test/3D_laplacian.jl
"""
s = x1, x2, x3, x4 = (-5:0.2:5, -5:0.2:5, -5:0.2:5, -5:0.2:5)
dx1 = dx2 = dx3 = dx4 = x[2] - x[1]
ricker(x1::T, x2::T, x3::T, x4 ::T) where T = (4*(x1^2+x2^2+x3^2 +x4^2) - 6)*exp(-(x1^2+x2^2+x3^2 +x4^2))
u0 = [ricker(X1, X2, X3, X4) for X4 in x4, X3 in x3, X2 in x2, X1 in x1]
Dxx1 = CenteredDifference{1}(2, 4, dx1, length(x1))
Dxx2 = CenteredDifference{1}(2, 4, dx2, length(x2))
Dxx3 = CenteredDifference{1}(2, 4, dx3, length(x3))
Dxx4 = CenteredDifference{1}(2, 4, dx4, length(x4))
A = Dxx1 + Dxx2 + Dxx3 + Dxx4
Q = compose(Dirichlet0BC(Float64, length.(s))...)
dt = dx/(sqrt(3)*3e8)
t = 0.0:dt:10/3e8
f(u,p,t) = (3e8)^2 .*(A*Q*u) #.+u .-u^3
using OrdinaryDiffEq
prob = ODEProblem(f, u0, (0., 0.5))
solve(prob,Tsit5(),abstol=1e-6,reltol=1e-6);
"""
=#
