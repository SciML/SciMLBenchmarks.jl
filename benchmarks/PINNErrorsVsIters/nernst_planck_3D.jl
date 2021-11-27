using NeuralPDE
using Quadrature, Cubature, Cuba
using Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Plots
using PyPlot
using DelimitedFiles
using QuasiMonteCarlo

print("Precompiling Done")

function nernst_planck_time(strategy, minimizer, maxIters)

    ##  DECLARATIONS
    @parameters t x y z
    @variables c(..)

    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)
    Dz = Differential(z)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2
    Dzz = Differential(z)^2



    ## DOMAINS AND OPERATORS

    # Discretization
    xwidth      = 1.0
    ywidth      = 1.0
    zwidth      = 1.0
    tmax        = 1.0
    xMeshNum    = 10
    yMeshNum    = 10
    zMeshNum    = 10
    tMeshNum    = 10

    dx  = xwidth/xMeshNum
    dy  = ywidth/yMeshNum
    dz  = zwidth/zMeshNum
    dt  = tmax/tMeshNum

    domains = [t ∈ IntervalDomain(0.0,tmax),
               x ∈ IntervalDomain(0.0,xwidth),
               y ∈ IntervalDomain(0.0,ywidth),
               z ∈ IntervalDomain(0.0,zwidth)]

    xs = 0.0 : dx : xwidth
    ys = 0.0 : dy : ywidth
    zs = 0.0 : dz : zwidth
    ts = 0.0 : dt : tmax

    # Constants
    D = 1  #dummy
    ux = 10 #dummy
    uy = 10 #dummy
    uz = 10 #dummy

    # Operators
    div = - D*(Dxx(c(t,x,y,z)) + Dyy(c(t,x,y,z)) + Dzz(c(t,x,y,z)))
          + (ux*Dx(c(t,x,y,z)) + uy*Dy(c(t,x,y,z)) + uz*Dz(c(t,x,y,z)))

    # Equation
    eq = Dt(c(t,x,y,z)) + div ~ 0      #NERNST-PLANCK EQUATION

    # Boundary conditions
    bcs = [c(0,x,y,z) ~ 0]

    ## NEURAL NETWORK
    n = 16   #neuron number

    chain = FastChain(FastDense(4,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1))   #Neural network from Flux library

    discretization = NeuralPDE.PhysicsInformedNN(chain,strategy)

    indvars = [t,x,y,z]   #independent variables
    depvars = [c]       #dependent (target) variable

    dim = length(domains)

    losses = []

    cb = function (p,l)     #loss function handling
        append!(losses, l)
        println(length(losses), " Current loss is: $l")
        return false
    end

    pde_system = PDESystem(eq, bcs, domains, indvars, depvars)
    prob = discretize(pde_system, discretization)

    res = GalacticOptim.solve(prob, minimizer; cb = cb, maxiters=1) #allow_f_increase = false,

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



function nernst_planck(strategy, minimizer, maxIters, params)

    ##  DECLARATIONS
    @parameters t x y z
    @variables c(..)

    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)
    Dz = Differential(z)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2
    Dzz = Differential(z)^2



    ## DOMAINS AND OPERATORS

    # Discretization
    xwidth      = 1.0
    ywidth      = 1.0
    zwidth      = 1.0
    tmax        = 1.0
    xMeshNum    = 10
    yMeshNum    = 10
    zMeshNum    = 10
    tMeshNum    = 10

    dx  = xwidth/xMeshNum
    dy  = ywidth/yMeshNum
    dz  = zwidth/zMeshNum
    dt  = tmax/tMeshNum

    domains = [t ∈ IntervalDomain(0.0,tmax),
               x ∈ IntervalDomain(0.0,xwidth),
               y ∈ IntervalDomain(0.0,ywidth),
               z ∈ IntervalDomain(0.0,zwidth)]

    xs = 0.0 : dx : xwidth
    ys = 0.0 : dy : ywidth
    zs = 0.0 : dz : zwidth
    ts = 0.0 : dt : tmax

    # Constants
    D = 1  #dummy
    ux = 10 #dummy
    uy = 10 #dummy
    uz = 10 #dummy

    # Operators
    div = - D*(Dxx(c(t,x,y,z)) + Dyy(c(t,x,y,z)) + Dzz(c(t,x,y,z)))
          + (ux*Dx(c(t,x,y,z)) + uy*Dy(c(t,x,y,z)) + uz*Dz(c(t,x,y,z)))

    # Equation
    eq = Dt(c(t,x,y,z)) + div ~ 0      #NERNST-PLANCK EQUATION

    # Boundary conditions
    bcs = [c(0,x,y,z) ~ 0]

    ## NEURAL NETWORK
    n = 16   #neuron number

    chain = FastChain(FastDense(4,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1))   #Neural network from Flux library

    discretization = NeuralPDE.PhysicsInformedNN(chain,strategy)

    indvars = [t,x,y,z]   #independent variables
    depvars = [c]       #dependent (target) variable

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
    outfile = "final_params_nernst_planck.txt"
    open(outfile, "w") do f
      for i in final_params
        println(f, i)
      end
    end


    # Model prediction
    domain = [ts, xs, ys, zs]

    u_predict  = [reshape([phi([t,x,y,z],res.minimizer) for x in xs for y in ys for z in zs],
                 (length(xs),length(ys),length(zs))) for t in ts]


    return [error, u_predict, u_predict, domain]
end
