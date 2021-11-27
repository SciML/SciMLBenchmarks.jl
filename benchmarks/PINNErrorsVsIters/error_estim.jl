using NeuralPDE
using Quadrature, Cubature, Cuba
using Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Plots
using BenchmarkTools, CUDA

CUDA.allowscalar(false)

strategy = NeuralPDE.QuadratureTraining(quadrature_alg = CubaCuhre(), reltol = 1e-8, abstol = 1e-8, maxiters = 100, batch=100) #use algorithm for older versions
minimizer = GalacticOptim.ADAM(0.01)
maxIters = 1000

#function hamilton_jacobi(strategy, minimizer, maxIters)

    ##  DECLARATIONS
@parameters  t x1 x2 x3 x4
@variables   u(..)

#=@derivatives Dt'~t

@derivatives Dx1'~x1
@derivatives Dx2'~x2
@derivatives Dx3'~x3
@derivatives Dx4'~x4

@derivatives Dxx1''~x1
@derivatives Dxx2''~x2
@derivatives Dxx3''~x3
@derivatives Dxx4''~x4
=#
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

chain = FastChain(FastDense(5,n,Flux.σ),FastDense(n,n,Flux.σ),FastDense(n,1)) |> gpu  #Neural network from Flux library
initθ = DiffEqFlux.initial_params(chain) |> gpu

discretization = NeuralPDE.PhysicsInformedNN(chain, strategy) #Discretization used for training
#dx = 0.1
#ref_dscr = NeuralPDE.PhysicsInformedNN(chain, NeuralPDE.GridTraining(dx)) #Discretization used for error estimation

indvars = [t,x1,x2,x3,x4]   #phisically independent variables

depvars = [u]       #dependent (target) variable

dim = length(domains)

losses = []

#grid_strategy = NeuralPDE.GridTraining(dx)

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
    println(length(losses), " Current loss is: ", l, " uniform error is, ",  pde_loss_function(p) + bc_loss_function(p))
    return false
end


pde_system = PDESystem(eq, bcs, domains, [t,x1,x2,x3,x4], [u]) #substituted with indvars and depvars
prob = discretize(pde_system, discretization)

t_0 = time_ns()

res = GalacticOptim.solve(prob, minimizer; cb = cb_, maxiters=maxIters) #allow_f_increase = false,

t_f = time_ns()

training_time = (t_f - t_0)/10^9
#print(string("Training time = ",(t_f - t_0)/10^9))

phi = discretization.phi

# Model prediction
domain = [ts,x1s,x2s,x3s,x4s]

u_predict = [reshape([first(phi([t,x1,x2,x3,x4],res.minimizer)) for x1 in x1s for x2 in x2s for x3 in x3s for x4 in x4s], (length(x1s),length(x2s), length(x3s),length(x4s))) for t in ts]  #matrix of model's prediction

#return [losses, u_predict, u_predict, domain, training_time] #add numeric solution


losses, u_predict, u_predict, domain, training_time = hamilton_jacobi(NeuralPDE.QuadratureTraining(algorithm = CubaCuhre(), reltol = 1e-8, abstol = 1e-8, maxiters = 100),
                                                      GalacticOptim.ADAM(0.01), 500)
