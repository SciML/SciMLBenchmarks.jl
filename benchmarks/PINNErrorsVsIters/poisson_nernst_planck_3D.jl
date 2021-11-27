## THIS CODE SOLVES THE 3D PNP SYSTEM, DOESN'T CONVERGE (RETURNS NaNs)

using NeuralPDE, Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Quadrature, Cubature, Cuba, QuasiMonteCarlo
using Parameters
using Plots, LaTeXStrings

include("./nernst_planck_1d_num.jl")


# Parameters and variables

@parameters t x y z
@variables Phi(..) Na(..) Cl(..)
@derivatives Dt'~t
@derivatives Dx'~x
@derivatives Dy'~y
@derivatives Dz'~z
@derivatives Dxx''~x

t_ref = 1.0       # s
x_ref = 0.38      # dm
y_ref = 0.014     # dm
z_ref = 0.014     # dm
C_ref = 0.16      # mol/dm^3
Phi_ref = 1.0     # V

epsilon = 78.5    # K
F = 96485.3415    # A s mol^-1
R = 831.0         # kg dm^2 s^-2 K^-1 mol^-1
T = 298.0         # K

z_Na = 1.0        # non-dim
z_Cl = -1.0       # non-dim

D_Na = 0.89e-7    # dm^2 s^−1
D_Cl = 1.36e-7    # dm^2 s^−1

u_Na = D_Na * abs(z_Na) * F / (R * T)
u_Cl = D_Cl * abs(z_Cl) * F / (R * T)

t_max = 1.0  / t_ref   # non-dim
x_max = 0.38  / x_ref   # non-dim
y_max = 0.014 / y_ref   # non-dim
z_max = 0.014 / z_ref   # non-dim

Na_0 = 0.16 / C_ref     # non-dim
Cl_0 = 0.16 / C_ref     # non-dim
Phi_0 = 4.0 / Phi_ref   # non-dim

Na_anode = 0.0            # non-dim
Na_cathode = 2.0 * Na_0   # non-dim
Cl_anode = 1.37 * Cl_0    # non-dim
Cl_cathode = 0.0          # non-dim

Pe_Na = x_ref^2 / ( t_ref * D_Na )  # non-dim
Pe_Cl = x_ref^2 / ( t_ref * D_Cl )  # non-dim

M_Na = x_ref^2 / ( t_ref * Phi_ref * u_Na )  # non-dim
M_Cl = x_ref^2 / ( t_ref * Phi_ref * u_Cl )  # non-dim

Po_1 = (epsilon * Phi_ref) / (F * x_ref * C_ref) # non-dim


# Definitions

grad_Phi = [Dx(Phi(t,x,y,z)),Dy(Phi(t,x,y,z)),Dz(Phi(t,x,y,z))]
grad_Na  = [Dx(Na(t,x,y,z)),Dy(Na(t,x,y,z)),Dz(Na(t,x,y,z))]
grad_Cl  = [Dx(Cl(t,x,y,z)),Dy(Cl(t,x,y,z)),Dz(Cl(t,x,y,z))]
norm_y_pos  = [0.,1.,0.]
norm_y_neg  = [0.,-1.,0.]
norm_z_pos  = [0.,0.,1.]
norm_z_neg  = [0.,0.,-1.]


# Space and time domains ###################################################

domains = [t ∈ IntervalDomain(0.0, t_max),
           x ∈ IntervalDomain(0.0, x_max),
           y ∈ IntervalDomain(0.0, y_max),
           z ∈ IntervalDomain(0.0, z_max)]

dt = 0.1
dx = 0.1 # non-dim
dy = 0.1
dz = 0.1

xs = 0.0 : dx : x_max
ys = 0.0 : dy : y_max
zs = 0.0 : dz : z_max
ts = 0.0 : dt : t_max


# Equations, initial and boundary conditions ###############################

eq1 = Dxx(Phi(t,x,y,z)) ~ (1.0/Po_1)*(z_Na*Na(t,x,y,z) + z_Cl*Cl(t,x,y,z))
eq2 = Dt(Na(t,x,y,z)) ~ (1.0/Pe_Na)*Dxx(Na(t,x,y,z)) + z_Na/(abs(z_Na)*M_Na)*(Dx(Na(t,x,y,z))*Dx(Phi(t,x,y,z))
      + Na(t,x,y,z)*Dxx(Phi(t,x,y,z)))
eq3 = Dt(Cl(t,x,y,z)) ~ (1.0/Pe_Cl)*Dxx(Cl(t,x,y,z)) + z_Cl/(abs(z_Cl)*M_Cl)*(Dx(Cl(t,x,y,z))*Dx(Phi(t,x,y,z))
      + Cl(t,x,y,z)*Dxx(Phi(t,x,y,z)))
eqs = [eq1, eq2, eq3]


bcs = [
        Phi(t,0.0,y,z) ~ Phi_0,
        Phi(t,x_max,y,z) ~ 0.0,
        Dy(Phi(t,x,y,z)) ~ 0,  #scalar product between grad_Phi and norm_y_pos
        #- Dy(Phi(t,x,y,z)) ~ 0, #equivalent to the above condition
        Dz(Phi(t,x,y,z)) ~ 0,  #scalar product between grad_Phi and norm_z_pos
        #- Dz(Phi(t,x,y,z)) ~ 0, #equivalent to the above condition

        Na(0.0,x,y,z) ~ Na_0,
        Na(t,0.0,y,z) ~ Na_anode,
        Na(t,x_max,y,z) ~ Na_cathode,
        Dy(Na(t,x,y,z)) ~ 0,
        #- Dy(Na(t,x,y,z)) ~ 0,
        Dz(Na(t,x,y,z)) ~ 0,
        #- Dz(Na(t,x,y,z)) ~ 0,

        Cl(0.0,x,y,z) ~ Cl_0,
        Cl(t,0.0,y,z) ~ Cl_anode,
        Cl(t,x_max,y,z) ~ Cl_cathode,
        Dy(Cl(t,x,y,z)) ~ 0,
        #- Dy(Cl(t,x,y,z)) ~ 0,
        Dz(Cl(t,x,y,z)) ~ 0
        #- Dz(Cl(t,x,y,z)) ~ 0
      ]


# Neural network, Discretization ###########################################

dim = length(domains)
output = length(eqs)
neurons = 16
#chain = FastChain(FastDense(dim, neurons, Flux.σ),FastDense(neurons, neurons, Flux.σ),FastDense(neurons, output))

chain1 = FastChain( FastDense(dim, neurons, Flux.σ),
                    FastDense(neurons, neurons, Flux.σ),
                    FastDense(neurons, neurons, Flux.σ),
                    FastDense(neurons, 1))
chain2 = FastChain( FastDense(dim, neurons, Flux.σ),
                    FastDense(neurons, neurons, Flux.σ),
                    FastDense(neurons, neurons, Flux.σ),
                    FastDense(neurons, 1))
chain3 = FastChain( FastDense(dim, neurons, Flux.σ),
                    FastDense(neurons, neurons, Flux.σ),
                    FastDense(neurons, neurons, Flux.σ),
                    FastDense(neurons, 1))


q_strategy = NeuralPDE.QuadratureTraining(algorithm = CubaCuhre(),reltol = 1e-8, abstol = 1e-8, maxiters = 100)
#discretization = NeuralPDE.PhysicsInformedNN(chain, q_strategy)
discretization = PhysicsInformedNN([chain1, chain2, chain3], strategy = NeuralPDE.GridTraining())

pde_system = PDESystem(eqs, bcs, domains, [t,x,y,z], [Phi, Na, Cl])
prob = discretize(pde_system, discretization)

loss = []
cb = function (p, l)
    println("Current loss is: $l")
    push!(loss, l)
    return false
end

t_0 = time_ns()
res = GalacticOptim.solve(prob, GalacticOptim.ADAM(0.01); cb=cb, maxiters = 100)
t_f = time_ns()

training_time = (t_f - t_0)/10^9

phi = discretization.phi
initθ = discretization.initθ
acum = [0;accumulate(+, length.(initθ))]
sep = [acum[i]+1 : acum[i+1] for i in 1:length(acum)-1]
minimizers = [res.minimizer[s] for s in sep]

dt = t_max / 4.0
ts = 0.0:dt:t_max
xs = 0.0:dx:x_max
ys = 0.0:dy:y_max
zs = 0.0:dz:z_max

Phi_predict = [ [ phi[1]([t,x,y,z],minimizers[1])[1] for x in xs for y in ys for z in zs] for t in ts ]
Na_predict  = [ [ phi[2]([t,x,y,z],minimizers[2])[1] for x in xs for y in ys for z in zs] for t in ts ]
Cl_predict  = [ [ phi[3]([t,x,y,z],minimizers[3])[1] for x in xs for y in ys for z in zs] for t in ts ]


u_predict = [Phi_predict, Na_predict, Cl_predict]
u_num = nernst_planck_1d_num()

return [losses, u_predict, u_num, domain, training_time]

#end

function plot_PNP(res, loss, discretization, pars)

    discretization = discretization

    labels = permutedims([ "$t s" for t in ts ])
    p1 = plot(xs * x_ref, Phi_predict * Phi_ref,
              xlabel = "dm", ylabel = "V",title = L"$\Phi$",labels=labels);
    savefig("Phi.svg")
    p2 = plot(xs * x_ref, Na_predict * C_ref,
              xlabel = "dm", ylabel = "M",title = L"$Na^+$",labels=labels);
    savefig("Na.svg")
    p3 = plot(xs * x_ref, Cl_predict * C_ref,
              xlabel = "dm", ylabel = "M",title = L"$Cl^-$",labels=labels);
    savefig("Cl.svg")

end

function solve_PNP()
    strategy = NeuralPDE.GridTraining(dx=1e-2)
    return solve_PNP(strategy)
end
