using OrdinaryDiffEq, DiffEqDevTools, Sundials, Plots
using LinearAlgebra, StaticArrays, FiniteDiff

# System parameters
const N = 40  # Number of heated units

# Physical parameters
const Cu = N == 1 ? [2e7] : (ones(N) .+ range(0,1.348,length=N))*1e7  # Heat capacity of heated units
const Cd = 2e6*N  # Heat capacity of distribution circuit
const Gh = 200    # Thermal conductance of heating elements
const Gu = 150    # Thermal conductance of heated units to the atmosphere
const Qmax = N*3000  # Maximum power output of heat generation unit
const Teps = 0.5     # Threshold of heated unit temperature controllers
const Td0 = 343.15   # Set point of distribution circuit temperature
const Tu0 = 293.15   # Heated unit temperature set point
const Kp = Qmax/4    # Proportional gain of heat generation unit temperature controller
const a = 50         # Gain of the hysteresis function
const b = 15         # Slope of the saturation function at the origin

# Helper functions for parameter and variable organization
function set_up_params(p)
    Cu   = @view p[1:N]
    Cd   = p[N+1]
    Gh   = p[N+2]
    Gu   = p[N+3]
    Qmax = p[N+4]
    Teps = p[N+5]
    Td0  = p[N+6]
    Tu0  = p[N+7]
    Kp   = p[N+8]
    a    = p[N+9]
    b    = p[N+10]
    
    Cu, Cd, Gh, Gu, Qmax, Teps, Td0, Tu0, Kp, a, b
end

function set_up_vars(u)
    Td = u[1]
    Tu = @view u[2:N+1]
    x  = @view u[N+2:2N+1]
    
    Td, Tu, x
end

function combine_vars!(du, dTd, dTu, dx)
    du[1] = dTd
    du[2:N+1] .= dTu[1:N]
    du[N+2:2N+1] .= dx[1:N]
    
    nothing
end

# Hysteresis function for temperature controllers
hist(x, p, e = 1) = -(x + 0.5)*(x - 0.5) * x * 1/(0.0474)*e + p

# Saturation function for control outputs
sat(x, xmin, xmax) = tanh(2*(x-xmin)/(xmax-xmin)-1) * (xmax-xmin)/2 + (xmax+xmin)/2

# Main ODE function for the heating system
function heating(du, u, (p, Qh, Que), t)
    (Td, Tu, x) = set_up_vars(u)
    (dTd, dTu, dx) = set_up_vars(du)
    (Cu, Cd, Gh, Gu, Qmax, Teps, Td0, Tu0, Kp, a, b) = set_up_params(p)
    
    # External temperature with daily variation
    Text = 278.15 + 8*sin(2π*t/86400)
    
    # Heat transfer from distribution to units (with control)
    @inbounds for i in 1:length(Qh)
        Qh[i] = Gh * (Td - Tu[i]) * (sat(b * x[i], -0.5, 0.5) + 0.5)
    end
    
    # Heat loss from units to environment
    @inbounds for i in 1:length(Que)
        Que[i] = Gu * (Tu[i] - Text)
    end
    
    # Heat input to distribution circuit (with saturation)
    Qd = sat(Kp*(Td0 - Td), 0, Qmax)
    
    # Temperature dynamics
    dTd = (Qd - sum(Qh))/Cd
    @inbounds for i in 1:length(dTu)
        dTu[i] = (Qh[i] - Que[i]) / Cu[i]
    end
    
    # Controller states (hysteretic behavior)
    @inbounds for i in 1:length(dx)
        dx[i] = a * hist(x[i], Tu0 - Tu[i], Teps)
    end
    
    combine_vars!(du, dTd, dTu, dx)
end

# Initial conditions and parameters
u0 = [Td0, fill(Tu0, N)..., fill(-0.5, N)...]
Qh0 = zeros(N)
Que0 = zeros(N)
p = [Cu..., Cd, Gh, Gu, Qmax, Teps, Td0, Tu0, Kp, a, b]

println("Testing heating system benchmark...")
println("Number of variables: ", length(u0))
println("Initial distribution temperature: ", u0[1])
println("Initial unit temperatures: ", u0[2])
println("Number of parameters: ", length(p))

# Test with shorter time span first
tspan = (0., 1000.)
prob = ODEProblem(heating, u0, tspan, (p, Qh0, Que0))

println("\nSolving with short time span...")
@time sol = solve(prob, Rodas4P(autodiff=false), reltol=1e-6, abstol=1e-8)

println("Solution length: ", length(sol))
println("Final distribution temperature: ", sol[end][1])
println("Final unit temperature (first): ", sol[end][2])

# Test a simple work-precision comparison  
println("\nTesting work-precision setup...")
test_sol = TestSolution(sol)

abstols = [1e-5, 1e-6]
reltols = [1e-2, 1e-3]

setups = [Dict(:alg=>Rodas4P(autodiff=false)),
          Dict(:alg=>CVODE_BDF())]

wp = WorkPrecisionSet(prob, abstols, reltols, setups;
                      appxsol=test_sol, maxiters=Int(1e4), error_estimate=:final)

println("Work-precision set created successfully!")
println("Number of algorithms tested: ", length(setups))

# Test plotting
try
    p1 = plot(sol, vars=1, title="Distribution Temperature", legend=false)
    println("Plotting successful!")
catch e
    println("Plotting failed: ", e)
end

println("\nHeating system benchmark test completed successfully!")