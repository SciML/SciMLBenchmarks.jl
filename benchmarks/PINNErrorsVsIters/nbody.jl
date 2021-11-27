using NeuralPDE
using Quadrature, Cubature, Cuba
using Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Plots
using PyPlot
using DelimitedFiles
using QuasiMonteCarlo

print("Precompiling Done")

function schrodinger(strategy, minimizer, maxIters)

    ##  DECLARATIONS
    @parameters t x y z
    @variables u(..)
    @derivatives Dt'~t
    @derivatives Dx'~x
    @derivatives Dy'~y
    @derivatives Dz'~z
    @derivatives Dxx''~x
    @derivatives Dyy''~y
    @derivatives Dzz''~z

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
    dy  = ywidth/xMeshNum
    dz  = zwidth/xMeshNum
    dt  = tmax/tMeshNum

    domains = [t ∈ IntervalDomain(0.0,tmax),
               x ∈ IntervalDomain(0.0,xwidth),
               y ∈ IntervalDomain(0.0,ywidth),
               z ∈ IntervalDomain(0.0,zwidth)]

    xs = 0.0 : dx : xwidth
    ys = 0.0 : dy : ywidth
    zs = 0.0 : dz : zwidth

    # Definitions
    lapl = Dxx(u(t,x,y,z)) + Dyy(u(t,x,y,z)) + Dzz(u(t,x,y,z))  #Laplace's operator
    V = 0       #Potential
    ħ = 1       #reduced Planck's Constant
    m = 1       #particle's mass

    initialCondition = exp(- x^2 - y^2 - z^2)

    bcs = [u(0,x,y,z) ~ initialCondition]

    eq = Dt(u(t,x,y,z)) + lapl - V*u(t,x,y,z) ~ 0
    #=
    DISCLAIMER: This is NOT the Schrodingere equations, it only resembles it.
    In orde to make it work properly, we first need to find the enrgy eigenvalue
    =#


    ## NEURAL NETWORK
    n = 16   #neuron number

    chain = FastChain(FastDense(4,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1))   #Neural network from Flux library

    discretization = NeuralPDE.PhysicsInformedNN(chain, strategy = strategy)

    indvars = [t,x,y,z]   #phisically independent variables
    depvars = [u]       #dependent (target) variable

    dim = length(domains)

    losses = []
    cb = function (p,l)     #loss function handling
        println("Current loss is: $l")
        append!(losses, l)
        return false
    end

    pde_system = PDESystem(eq, bcs, domains, indvars, depvars)
    prob = discretize(pde_system, discretization)

    t_0 = time_ns()

    res = GalacticOptim.solve(prob, minimizer; cb = cb, maxiters=maxIters) #allow_f_increase = false,

    t_f = time_ns()
    print(string("Training time = ",(t_f - t_0)/10^9))

    phi = discretization.phi

    printBCSComp = true     #prints initial condition comparison and training loss plot

    xs = 0.0 : dx : xwidth
    ts = 1 : dt : tmax

    u_predict = [reshape([first(phi([t,x,y,z],res.minimizer)) for x in xs for y in ys for z in zs],
                    (length(xs),length(ys),length(zs))) for t in ts]  #matrix of model's prediction

    maxlim = maximum(maximum(u_predict[t]) for t = 1:length(ts))
    minlim = minimum(minimum(u_predict[t]) for t = 1:length(ts))

    trainingPlot = Plots.plot(1:(maxIters + 1), losses, yaxis=:log, title = string("Training time = ",(t_f - t_0)/10^9," s",
    "\\n Iterations: ", maxIters), ylabel = "log(loss)", legend = false) #loss plot

    #GLMakie.contour(u_predict[9],colormap=:balance,algorithm = :iso)
    #GLMakie.volume(u_predict[2])
end

schrodinger(NeuralPDE.QuadratureTraining(algorithm = CubaCuhre(), reltol = 1e-8, abstol = 1e-8, maxiters = 100),GalacticOptim.ADAM(0.01), 100)
