using LinearAlgebra, OrdinaryDiffEq, Sundials

# MINIMAL RIGID-BODY SLIDER-CRANK (no elasticity)

const L1 = 0.15
const L2 = 0.30
const M1 = 0.36
const M2 = 0.151104
const M3 = 0.075552
const J1 = 0.002727
const J2 = 0.0045339259
const GRAV = 0.0
const OMEGA = 150.0

# Rigid only: 3 positions (φ₁, φ₂, x₃) + 3 velocities + 3 accelerations + 3 multipliers = 12 vars

function slider_crank_rigid_resmbs!(residual, y, yprime, t)
    """
    Minimal rigid-body index-2 DAE
    y = [φ₁, φ₂, x₃, ω₁, ω₂, v_x₃, ϕ̈₁, ϕ̈₂, a_x₃, λ₁, λ₂, λ₃]
    """
    
    # Extract positions and velocities
    p1, p2, x3 = y[1], y[2], y[3]
    v1, v2, vx3 = y[4], y[5], y[6]
    w1, w2, wx3 = y[7], y[8], y[9]
    lam1, lam2, lam3 = y[10], y[11], y[12]
    
    cosp1 = cos(p1)
    cosp2 = cos(p2)
    sinp1 = sin(p1)
    sinp2 = sin(p2)
    cosp12 = cos(p1 - p2)
    sinp12 = sin(p1 - p2)
    
    # Mass matrix (rigid only)
    M = zeros(3, 3)
    M[1,1] = J1 + M2*L1^2
    M[1,2] = 0.5*L1*L2*M2*cosp12
    M[2,1] = 0.5*L1*L2*M2*cosp12
    M[2,2] = J2
    M[3,3] = M3
    
    # Constraint Jacobian
    G = zeros(3, 3)
    G[1,1] = L1*cosp1
    G[1,2] = L2*cosp2
    G[2,1] = L1*sinp1
    G[2,2] = L2*sinp2
    G[2,3] = 1.0
    G[3,1] = 1.0
    
    # Forces
    F = zeros(3)
    F[1] = -0.5*L1*GRAV*(M1 + 2*M2)*cosp1 - 0.5*L1*L2*M2*v2^2*sinp12
    F[2] = -0.5*L2*GRAV*M2*cosp2 + 0.5*L1*L2*M2*v1^2*sinp12
    F[3] = 0.0
    
    # Residual rows 1-3: dp/dt = v
    residual[1:3] = yprime[1:3] - y[4:6]
    
    # Residual rows 4-6: dv/dt = w
    residual[4:6] = yprime[4:6] - y[7:9]
    
    # Residual rows 7-9: M*w = F + G'*λ (rearranged: M*w - F - G'*λ = 0)
    Mw = M * y[7:9]
    GTlam = G' * y[10:12]
    residual[7:9] = Mw - F - GTlam
    
    # Residual rows 10-12: G*v + r'(t) = 0 (velocity constraints)
    # r'₁, r'₂ = 0 (holonomic)
    # r'₃ = -Ω (prescribed crank angular velocity)
    Gv = G * y[4:6]
    residual[10] = Gv[1]
    residual[11] = Gv[2]
    residual[12] = Gv[3] - OMEGA
end

# Initial conditions
y0 = zeros(12)
y0[1] = 0.0      # φ₁(0)
y0[2] = 0.0      # φ₂(0)
y0[3] = L1 + L2  # x₃(0) = L₁ + L₂ (positions at equilibrium)
y0[4] = 150.0    # ω₁(0) = Ω
y0[5] = -74.9957670  # ω₂(0)
y0[6] = -0.268938672e-05  # v_x₃(0)
y0[7] = 0.0      # ẇ₁(0) will be computed
y0[8] = 0.0      # ẇ₂(0) will be computed
y0[9] = 0.0      # ȧ_x₃(0) will be computed
y0[10] = 0.0     # λ₁(0) will be computed
y0[11] = 0.0     # λ₂(0) will be computed
y0[12] = 0.0     # λ₃(0) will be computed

# Problem setup
tspan = (0.0, 0.01)

# Initial derivative: yp0
# yp0[1:3] = v0 (= y0[4:6])
# yp0[4:6] = w0 (= y0[7:9])
# yp0[7:12] = 0 (algebraic vars, no time derivative)
yp0 = zeros(12)
yp0[1:3] = y0[4:6]  # dp/dt = v
yp0[4:6] = y0[7:9]  # dv/dt = w (accelerations)
# yp0[7:12] = 0 (algebraic)

differential_vars = [true, true, true, true, true, true, false, false, false, false, false, false]

function residual_wrapper!(residual, yprime, y, p, t)
    slider_crank_rigid_resmbs!(residual, y, yprime, t)
end

prob_dae = DAEProblem(residual_wrapper!, yp0, y0, tspan;
                      differential_vars=differential_vars)

println("="^70)
println("MINIMAL RIGID-BODY SLIDER-CRANK (IDA test)")
println("="^70)
println("12 variables: 3 pos + 3 vel + 3 accel + 3 multipliers")
println("Index-2 semi-explicit DAE")
println()

try
    println("Solving with IDA...")
    sol = solve(prob_dae, IDA(); abstol=1e-6, reltol=1e-4)
    println("✓ Solution computed!")
    println("  Steps: $(length(sol))")
    println("  φ₁(0.01) = $(sol.u[end][1])")
    println("  x₃(0.01) = $(sol.u[end][3])")
catch e
    println("✗ Failed: $(e)")
end
