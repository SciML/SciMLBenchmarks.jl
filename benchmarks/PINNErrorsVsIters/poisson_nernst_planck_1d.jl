################################################################################
#   This code solves a simplified 1D Poisson-Nernst-Planck system using NeuralPDE
#
#   Equations:
#        d2Phi/dx2 =    ( 1.0 / Po_1 ) * ( z_Na^+ * Na^+ + z_Cl^- * Cl^- )
#        dNa^+/dt =     ( 1.0 / Pe_Na^+ ) * d2Na^+/dx2
#                     + z_Na^+ / ( abs(z_Na^+) * M_Na^+ ) *
#                     ( dNa^+/dx * dPhi/dx + Na^+ * d2Phi/dx2 )
#        dCl^-/dt =     ( 1.0 / Pe_Cl^- ) * d2Cl^-/dx2
#                     + z_Cl^- / ( abs(z_Cl^-) * M_Cl^- )
#                     * ( dCl^-/dx * dPhi/dx + Cl^- * d2Phi/dx2 )
#
#   Initial conditions:
#        Phi(0,x) = 0.0
#        Na^+(0,x) = Na^+_0
#        Cl^-(0,x) = Cl^-_0
#
#   Boundary conditions (simplified):
#        Phi(t,0) = Phi_0
#        Phi(t,n) = 0.0
#        Na^+(t,0) = 0.0
#        Na^+(t,n) = 2.0 * Na^+_0
#        Cl^-(t,0) = 1.37 * Cl^-_0
#        Cl^-(t,n) = 0.0
#
#   How to run:
#        $ julia solve-pnp.jl
#
#        or
#
#        $ julia
#        julia> include("1d-poisson-nerst-planck-Cl-Na-adim-neuralpde.jl")
#        julia> res, loss, discretization, pars = solve_PNP()
#        julia> plot_PNP(res, loss, discretization, pars)
#
#################################################################################

using NeuralPDE, Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Quadrature, Cubature, Cuba, QuasiMonteCarlo
using Parameters
using Plots, LaTeXStrings

include("./nernst_planck_1d_num.jl")

@with_kw struct Params
    t_ref = 1.0       # s
    x_ref = 0.38      # dm
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

    t_max = 0.01 / t_ref    # non-dim
    x_max = 0.38 / x_ref    # non-dim
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

    Po_1 = (epsilon * Phi_ref) / (F * x_ref * C_ref)  # non-dim

    dx = 0.01 # non-dim
end

function nernst_planck_1d(strategy, minimizer, maxIters)

    # Parameters and variables

    @parameters t,x
    @variables Phi(..),Na(..),Cl(..)
    @derivatives Dt'~t
    @derivatives Dx'~x
    @derivatives Dxx''~x

    pars = Params()

    # Equations, initial and boundary conditions ###############################

    eqs = [
            ( Dxx(Phi(t,x)) ~ ( 1.0 / pars.Po_1 ) *
                              ( pars.z_Na * Na(t,x) + pars.z_Cl * Cl(t,x) ) )
            ,
            ( Dt(Na(t,x)) ~ ( 1.0 / pars.Pe_Na ) * Dxx(Na(t,x))
                          +   pars.z_Na / ( abs(pars.z_Na) * pars.M_Na )
                          * ( Dx(Na(t,x)) * Dx(Phi(t,x)) + Na(t,x) * Dxx(Phi(t,x)) ) )
            ,
            ( Dt(Cl(t,x)) ~ ( 1.0 / pars.Pe_Cl ) * Dxx(Cl(t,x))
                          +   pars.z_Cl / ( abs(pars.z_Cl) * pars.M_Cl )
                          * ( Dx(Cl(t,x)) * Dx(Phi(t,x)) + Cl(t,x) * Dxx(Phi(t,x)) ) )
          ]

    bcs = [
            Phi(t,0.0) ~ pars.Phi_0,
            Phi(t,pars.x_max) ~ 0.0
            ,
            Na(0.0,x) ~ pars.Na_0,
            Na(t,0.0) ~ pars.Na_anode,
            Na(t,pars.x_max) ~ pars.Na_cathode
            ,
            Cl(0.0,x) ~ pars.Cl_0,
            Cl(t,0.0) ~ pars.Cl_anode,
            Cl(t,pars.x_max) ~ pars.Cl_cathode
          ]

    # Space and time domains ###################################################

    domains = [
                t ∈ IntervalDomain(0.0, pars.t_max),
                x ∈ IntervalDomain(0.0, pars.x_max)
              ]

    # Neural network, Discretization ###########################################

    dim = length(domains)
    output = length(eqs)
    neurons = 16
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

    discretization = PhysicsInformedNN([chain1, chain2, chain3], strategy=strategy)

    pde_system = PDESystem(eqs, bcs, domains, [t,x], [Phi, Na, Cl])
    prob = discretize(pde_system, discretization)

    loss = []
    cb = function (p, l)
        println("Current loss is: $l")
        push!(loss, l)
        return false
    end

    t_0 = time_ns()
    res = GalacticOptim.solve(prob, minimizer ; cb=cb, maxiters = maxIters)
    t_f = time_ns()

    training_time = (t_f - t_0)/10^9

    phi = discretization.phi
    initθ = discretization.initθ

    Phi_predict = [ [ phi[1]([t,x],minimizers[1])[1] for x in xs ] for t in ts ]
    Na_predict  = [ [ phi[2]([t,x],minimizers[2])[1] for x in xs ] for t in ts ]
    Cl_predict  = [ [ phi[3]([t,x],minimizers[3])[1] for x in xs ] for t in ts ]

    u_predict = [Phi_predict, Na_predict, Cl_predict]
    u_num = nernst_planck_1d_num()
    #return res, loss, discretization, pars
    return [losses, u_predict, u_num, domain, training_time]

end

function plot_PNP(res, loss, discretization, pars)

    discretization = discretization
    t_max = pars.t_max
    x_max = pars.x_max
    dx = pars.dx
    x_ref = pars.x_ref
    Phi_ref = pars.Phi_ref
    C_ref = pars.C_ref

    phi = discretization.phi
    initθ = discretization.initθ
    acum = [0;accumulate(+, length.(initθ))]
    sep = [acum[i]+1 : acum[i+1] for i in 1:length(acum)-1]
    minimizers = [res.minimizer[s] for s in sep]

    dt = t_max / 4.0
    ts = 0.0:dt:t_max
    xs = 0.0:dx:x_max

    Phi_predict = [ [ phi[1]([t,x],minimizers[1])[1] for x in xs ] for t in ts ]
    Na_predict  = [ [ phi[2]([t,x],minimizers[2])[1] for x in xs ] for t in ts ]
    Cl_predict  = [ [ phi[3]([t,x],minimizers[3])[1] for x in xs ] for t in ts ]

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
    strategy  = NeuralPDE.QuadratureTraining(algorithm = HCubatureJL(), reltol = 1e-8, abstol = 1e-8, maxiters = 100)
    minimizer = Optim.BFGS()
    maxIters  = 200
    return nernst_planck_1d(strategy, minimizer, maxIters)
end

solve_PNP()
