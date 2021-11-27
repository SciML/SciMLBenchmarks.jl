using NeuralPDE
using Quadrature, Cubature, Cuba
using Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Plots
using PyPlot
using DelimitedFiles
using QuasiMonteCarlo

print("Precompiling Done")

function level_set_time(strategy, minimizer, maxIters)

    ##  DECLARATIONS
    @parameters  t x y
    @variables   u(..)

    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)

    # Discretization
    xwidth      = 1.0      #ft
    ywidth      = 1.0
    tmax        = 1.0      #min
    xScale      = 1.0
    yScale      = 1.0
    xMeshNum    = 10
    yMeshNum    = 10
    tMeshNum    = 10
    dx  = xwidth/xMeshNum
    dy  = ywidth/yMeshNum
    dt  = tmax/tMeshNum


    domains = [t ∈ IntervalDomain(0.0,tmax),
               x ∈ IntervalDomain(0.0,xwidth),
               y ∈ IntervalDomain(0.0,ywidth)]

    xs = 0.0 : dx : xwidth
    ys = 0.0 : dy : ywidth
    ts = 0.0 : dt : tmax

    # Definitions
    x0    = 0.5
    y0    = 0.5
    Uwind = [0.0, 2.0]  #wind vector

    # Operators
    gn   = (Dx(u(t,x,y))^2 + Dy(u(t,x,y))^2)^0.5  #gradient's norm
    ∇u   = [Dx(u(t,x,y)), Dy(u(t,x,y))]
    n    = ∇u/gn              #normal versor
    #U    = ((Uwind[1]*n[1] + Uwind[2]*n[2])^2)^0.5 #inner product between wind and normal vector

    R0 = 0.112471
    ϕw = 0#0.156927*max((0.44*U)^0.04086,1.447799)
    ϕs = 0
    S  = R0*(1 + ϕw + ϕs)

    # Equation
    eq = Dt(u(t,x,y)) + S*gn ~ 0  #LEVEL SET EQUATION

    initialCondition = (xScale*(x - x0)^2 + (yScale*(y - y0)^2))^0.5 - 0.2   #Distance from ignition

    bcs = [u(0,x,y) ~ initialCondition]  #from literature


    ## NEURAL NETWORK
    n = 16   #neuron number

    chain = FastChain(FastDense(3,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1))   #Neural network from Flux library

    discretization = NeuralPDE.PhysicsInformedNN(chain, strategy)

    indvars = [t,x,y]   #phisically independent variables
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


function level_set(strategy, minimizer, maxIters, params)

    ##  DECLARATIONS
    @parameters  t x y
    @variables   u(..)

    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)

    # Discretization
    xwidth      = 1.0      #ft
    ywidth      = 1.0
    tmax        = 1.0      #min
    xScale      = 1.0
    yScale      = 1.0
    xMeshNum    = 10
    yMeshNum    = 10
    tMeshNum    = 10
    dx  = xwidth/xMeshNum
    dy  = ywidth/yMeshNum
    dt  = tmax/tMeshNum


    domains = [t ∈ IntervalDomain(0.0,tmax),
               x ∈ IntervalDomain(0.0,xwidth),
               y ∈ IntervalDomain(0.0,ywidth)]

    xs = 0.0 : dx : xwidth
    ys = 0.0 : dy : ywidth
    ts = 0.0 : dt : tmax

    # Definitions
    x0    = 0.5
    y0    = 0.5
    Uwind = [0.0, 2.0]  #wind vector

    # Operators
    gn   = (Dx(u(t,x,y))^2 + Dy(u(t,x,y))^2)^0.5  #gradient's norm
    ∇u   = [Dx(u(t,x,y)), Dy(u(t,x,y))]
    n    = ∇u/gn              #normal versor
    #U    = ((Uwind[1]*n[1] + Uwind[2]*n[2])^2)^0.5 #inner product between wind and normal vector

    R0 = 0.112471
    ϕw = 0#0.156927*max((0.44*U)^0.04086,1.447799)
    ϕs = 0
    S  = R0*(1 + ϕw + ϕs)

    # Equation
    eq = Dt(u(t,x,y)) + S*gn ~ 0  #LEVEL SET EQUATION

    initialCondition = (xScale*(x - x0)^2 + (yScale*(y - y0)^2))^0.5 - 0.2   #Distance from ignition

    bcs = [u(0,x,y) ~ initialCondition]  #from literature


    ## NEURAL NETWORK
    n = 16   #neuron number

    chain = FastChain(FastDense(3,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1))   #Neural network from Flux library

    indvars = [t,x,y]   #phisically independent variables
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
    outfile = "final_params_level_set.txt"
    open(outfile, "w") do f
      for i in final_params
        println(f, i)
      end
    end


    # Model prediction
    domain = [ts, xs, ys]

    u_predict = [reshape([first(phi([t,x,y],res.minimizer)) for x in xs for y in ys], (length(xs),length(ys))) for t in ts]  #matrix of model's prediction

    return [error, u_predict, u_predict, domain] #add numeric solution
end

#level_set(NeuralPDE.QuadratureTraining(algorithm = CubaCuhre(), reltol = 1e-8, abstol = 1e-8, maxiters = 100), GalacticOptim.ADAM(0.01), 500)
