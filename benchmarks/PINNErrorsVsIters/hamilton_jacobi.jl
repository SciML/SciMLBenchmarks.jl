using NeuralPDE
using Quadrature, Cubature, Cuba
using Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Plots
using PyPlot
using QuasiMonteCarlo

print("Precompiling Done")

function hamilton_jacobi_time(strategy, minimizer, maxIters)

    ##  DECLARATIONS
    @parameters  t x1 x2 x3 x4
    @variables   u(..)

    Dt = Differential(t)

    Dx1 = Differential(x1)
    Dx2 = Differential(x2)
    Dx3 = Differential(x3)
    Dx4 = Differential(x4)

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

    λ = 1.0f0

    # Operators
    Δu = Dxx1(u(t,x1,x2,x3,x4)) + Dxx2(u(t,x1,x2,x3,x4)) + Dxx3(u(t,x1,x2,x3,x4)) + Dxx4(u(t,x1,x2,x3,x4)) # Laplacian
    ∇u = [Dx1(u(t,x1,x2,x3,x4)), Dx2(u(t,x1,x2,x3,x4)),Dx3(u(t,x1,x2,x3,x4)),Dx4(u(t,x1,x2,x3,x4))]

    # Equation
    eq = Dt(u(t,x1,x2,x3,x4)) + Δu - λ*sum(∇u.^2) ~ 0  #HAMILTON-JACOBI-BELLMAN EQUATION

    terminalCondition =  log((1 + x1*x1 + x2*x2 + x3*x3 + x4*x4)/2) # see PNAS paper

    bcs = [u(tmax,x1,x2,x3,x4) ~ terminalCondition]  #PNAS paper again

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


function hamilton_jacobi(strategy, minimizer, maxIters, params)

    ##  DECLARATIONS
    @parameters  t x1 x2 x3 x4
    @variables   u(..)

    Dt = Differential(t)

    Dx1 = Differential(x1)
    Dx2 = Differential(x2)
    Dx3 = Differential(x3)
    Dx4 = Differential(x4)

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

    λ = 1.0f0

    # Operators
    Δu = Dxx1(u(t,x1,x2,x3,x4)) + Dxx2(u(t,x1,x2,x3,x4)) + Dxx3(u(t,x1,x2,x3,x4)) + Dxx4(u(t,x1,x2,x3,x4)) # Laplacian
    ∇u = [Dx1(u(t,x1,x2,x3,x4)), Dx2(u(t,x1,x2,x3,x4)),Dx3(u(t,x1,x2,x3,x4)),Dx4(u(t,x1,x2,x3,x4))]

    # Equation
    eq = Dt(u(t,x1,x2,x3,x4)) + Δu - λ*sum(∇u.^2) ~ 0  #HAMILTON-JACOBI-BELLMAN EQUATION

    terminalCondition =  log((1 + x1*x1 + x2*x2 + x3*x3 + x4*x4)/2) # see PNAS paper

    bcs = [u(tmax,x1,x2,x3,x4) ~ terminalCondition]  #PNAS paper again

    ## NEURAL NETWORK
    n = 20   #neuron number

    chain = FastChain(FastDense(5,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1))   #Neural network from Flux library

    indvars = [t,x1,x2,x3,x4]   #phisically independent variables
    depvars = [u]       #dependent (target) variable

    dim = length(domains)

    losses = []
    error  = []

    phi = NeuralPDE.get_phi(chain)
    derivative = NeuralPDE.get_numeric_derivative()
    initθ = DiffEqFlux.initial_params(chain)

    dim = length(domains)

    dx = 0.1
    
    error_strategy = GridTraining(dx) #NeuralPDE.QuadratureTraining(quadrature_alg=CubatureJLh(),reltol= 1e-4,abstol= 1e-3,maxiters=10, batch=10)
    _pde_loss_function = NeuralPDE.build_loss_function(eq,indvars,depvars,
                                         phi,derivative,chain,initθ,error_strategy)

    bc_indvars = NeuralPDE.get_variables(bcs,indvars,depvars)
    _bc_loss_functions = [NeuralPDE.build_loss_function(bc,indvars,depvars,
                                              phi,derivative,chain,initθ,error_strategy,
                                              bc_indvars = bc_indvar) for (bc,bc_indvar) in zip(bcs,bc_indvars)]

    train_sets = NeuralPDE.generate_training_sets(domains,dx,[eq],bcs,indvars,depvars)
    train_domain_set, train_bound_set = train_sets


    pde_loss_function = NeuralPDE.get_loss_function([_pde_loss_function],
                                          train_domain_set,
                                          error_strategy)

    bc_loss_function = NeuralPDE.get_loss_function(_bc_loss_functions,
                                         train_bound_set,
                                         error_strategy)

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
    outfile = "final_params_hamilton_jacobi.txt"
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
