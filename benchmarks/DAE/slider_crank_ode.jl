"""
    SliderCrankODEForm

ODE mass-matrix form of Simeon (1998) slider-crank mechanism.
Derived from the index-2 DAE by solving for accelerations and multipliers at each step.

The 14-variable ODE form is:
    dp/dt = v
    dv/dt = M(p,q)^{-1} * [f(p,v,q) + G(p)'*λ(p,v,q)]

where λ is solved from the constraint projection:
    λ = (G*M^{-1}*G')^{-1} * (G*M^{-1}*f + G*v + r'(t))
"""

using OrdinaryDiffEq, LinearAlgebra

# Physical constants
const M1 = 0.36
const M2 = 0.151104
const M3 = 0.075552
const L1 = 0.15
const L2 = 0.30
const J1 = 0.002727
const J2 = 0.0045339259
const PI = 3.1415927
const EE = 0.20e12
const BB = 0.0080
const HH = 0.0080
const RHO = 7870.0
const GRAV = 0.0
const OMEGA = 150.0

const NQ = 4
const NP = 7
const NV = 2 * NP  # 14 variables: positions + velocities

# Pre-compute FE matrices
function initialize_fe_matrices()
    FACM = RHO * BB * HH * L2
    FACK = EE * BB * HH / L2
    FACB = BB * HH * L2
    
    MQ = zeros(NQ, NQ)
    MQ[1,1] = FACM * 0.5
    MQ[2,2] = FACM * 0.5
    MQ[3,3] = FACM * 8.0
    MQ[3,4] = FACM * 1.0
    MQ[4,3] = FACM * 1.0
    MQ[4,4] = FACM * 2.0
    
    KQ = zeros(NQ, NQ)
    KQ[1,1] = FACK * PI^4 / 24.0 * (HH/L2)^2
    KQ[2,2] = FACK * PI^4 * 2.0 / 3.0 * (HH/L2)^2
    KQ[3,3] = FACK * 16.0 / 3.0
    KQ[3,4] = -FACK * 8.0 / 3.0
    KQ[4,3] = -FACK * 8.0 / 3.0
    KQ[4,4] = FACK * 7.0 / 3.0
    
    BQ = zeros(NQ, NQ)
    DQ = zeros(NQ, NQ)
    
    c1 = zeros(NQ)
    c2 = zeros(NQ)
    c12 = zeros(NQ)
    c21 = zeros(NQ)
    
    c1[3] = FACB * 2.0 / 3.0
    c1[4] = FACB * 1.0 / 6.0
    c2[1] = FACB * 2.0 / PI
    c12[3] = L2 * FACB * 1.0 / 3.0
    c12[4] = L2 * FACB * 1.0 / 6.0
    c21[1] = L2 * FACB * 1.0 / PI
    c21[2] = -L2 * FACB * 0.5 / PI
    
    return MQ, KQ, BQ, DQ, c1, c2, c12, c21
end

const MQ, KQ, BQ, DQ, c1, c2, c12, c21 = initialize_fe_matrices()

"""
Compute mass matrix AM(7×7) and forces F(7) for the ODE solver.
"""
function compute_mass_and_forces(y::Vector{Float64}, v::Vector{Float64})
    # Extract positions and elastic modes
    p1, p2, x3 = y[1], y[2], y[3]
    q = y[4:7]
    
    # Precompute trig
    cosp1 = cos(p1)
    cosp2 = cos(p2)
    sinp1 = sin(p1)
    sinp2 = sin(p2)
    cosp12 = cos(p1 - p2)
    sinp12 = sin(p1 - p2)
    
    # Extract velocities
    v1, v2, vx3 = v[1], v[2], v[3]
    vq = v[4:7]
    
    # Scalar products
    c1Tq = dot(c1, q)
    c1Tqd = dot(c1, vq)
    c2Tq = dot(c2, q)
    c2Tqd = dot(c2, vq)
    c12Tq = dot(c12, q)
    c12Tqd = dot(c12, vq)
    
    MQq = MQ * q
    KQq = KQ * q
    BQqd = BQ * vq
    DQqd = DQ * vq
    
    qtmqq = dot(q, MQq)
    qdtmqq = dot(vq, MQq)
    qdtbqqd = dot(vq, BQqd)
    
    QtBQ = zeros(NQ)
    for i in 1:NQ
        QtBQ[i] = dot(q, BQ[:, i])
    end
    
    # Mass matrix
    AM = zeros(NP, NP)
    
    AM[1,1] = J1 + M2 * L1^2
    AM[1,2] = 0.5 * L1 * L2 * M2 * cosp12
    AM[2,2] = J2
    AM[3,3] = M3
    
    AM[1,2] += RHO * L1 * (sinp12 * c2Tq + cosp12 * c1Tq)
    AM[2,2] += qtmqq + 2.0 * RHO * c12Tq
    
    for i in 1:NQ
        AM[1, 3+i] = RHO * L1 * (-sinp12 * c1[i] + cosp12 * c2[i])
        AM[2, 3+i] = RHO * c21[i] + RHO * QtBQ[i]
        AM[3, 3+i] = 0.0
    end
    
    for i in 1:NQ
        for j in 1:i
            AM[3+j, 3+i] = MQ[j, i]
        end
    end
    
    for i in 1:NP
        for j in i+1:NP
            AM[j, i] = AM[i, j]
        end
    end
    
    # Forces
    F = zeros(NP)
    
    F[1] = -0.5 * L1 * GRAV * (M1 + 2.0*M2) * cosp1 - 0.5 * L1 * L2 * M2 * v2^2 * sinp12
    F[2] = -0.5 * L2 * GRAV * M2 * cosp2 + 0.5 * L1 * L2 * M2 * v1^2 * sinp12
    F[3] = 0.0
    
    F[1] += RHO * L1 * v2^2 * (-sinp12 * c1Tq + cosp12 * c2Tq) - 2.0 * RHO * L1 * v2 * (cosp12 * c1Tqd + sinp12 * c2Tqd)
    
    F[2] += RHO * L1 * v1^2 * (sinp12 * c1Tq - cosp12 * c2Tq) - 2.0 * RHO * v2 * c12Tqd - 2.0 * v2 * qdtmqq - RHO * qdtbqqd - RHO * GRAV * (cosp2 * c1Tq - sinp2 * c2Tq)
    
    for i in 1:NQ
        F[3+i] = v2^2 * MQq[i]
        F[3+i] += RHO * (v2^2 * c12[i] + L1 * v1^2 * (cosp12 * c1[i] + sinp12 * c2[i]) + 2.0 * v2 * BQqd[i])
        F[3+i] -= RHO * GRAV * (sinp2 * c1[i] + cosp2 * c2[i])
    end
    
    for i in 1:NQ
        F[3+i] -= KQq[i] + DQqd[i]
    end
    
    return AM, F
end

"""
Compute constraint Jacobian GP(3×7)
"""
function compute_constraint_jacobian(y::Vector{Float64})
    p1, p2, x3 = y[1], y[2], y[3]
    q = y[4:7]
    
    cosp1 = cos(p1)
    cosp2 = cos(p2)
    sinp1 = sin(p1)
    sinp2 = sin(p2)
    
    qku = q[4]
    
    GP = zeros(3, NP)
    
    GP[1, 1] = L1 * cosp1
    GP[1, 2] = L2 * cosp2 + qku * cosp2
    GP[1, 3] = 0.0
    GP[1, 7] = sinp2
    
    GP[2, 1] = L1 * sinp1
    GP[2, 2] = L2 * sinp2 + qku * sinp2
    GP[2, 3] = 1.0
    GP[2, 7] = -cosp2
    
    GP[3, 1] = 1.0
    
    return GP
end

"""
ODE RHS: compute accelerations by solving the algebraic system
"""
function slider_crank_ode!(du, u, p, t)
    # u = [p1, p2, x3, q1, q2, q3, q4, v1, v2, vx3, vq1, vq2, vq3, vq4]
    y = u[1:7]      # positions
    v = u[8:14]     # velocities
    
    # Compute mass matrix and forces
    AM, F = compute_mass_and_forces(y, v)
    
    # Compute constraint Jacobian
    GP = compute_constraint_jacobian(y)
    
    # Solve for accelerations and multipliers
    # System: [AM   -GP'] [w]   [F]
    #         [GP    0  ] [λ] = [RHS_constraint]
    
    # RHS for constraints: G*w + r''(t) = 0
    # r''(t) = 0 (no explicit time derivative for holonomic constraints)
    rhs_constraint = zeros(3)
    rhs_constraint[3] = -OMEGA * 0  # No acceleration term, velocity constraint only
    
    # Form augmented system
    sys = zeros(10, 10)
    sys[1:7, 1:7] = AM
    sys[1:7, 8:10] = -GP'
    sys[8:10, 1:7] = GP
    
    rhs = zeros(10)
    rhs[1:7] = F
    rhs[8:10] = rhs_constraint
    
    sol = sys \ rhs
    
    w = sol[1:7]  # accelerations
    
    # Fill derivatives
    du[1:7] = v               # dp/dt = v
    du[8:14] = w              # dv/dt = w
end

"""
Get initial conditions for ODE form
"""
function get_ode_initial_conditions()
    u0 = zeros(14)
    u0[1] = 0.0           # p1
    u0[2] = 0.0           # p2
    u0[3] = L1 + L2       # x3
    u0[4:7] .= 0.0        # q
    
    u0[8] = OMEGA         # v1
    u0[9] = -(L1/L2) * OMEGA  # v2
    u0[10:14] .= 0.0      # vq
    
    return u0
end

# For testing
if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    
    u0 = get_ode_initial_conditions()
    tspan = (0.0, 0.1)
    
    prob = ODEProblem(slider_crank_ode!, u0, tspan)
    
    println("="^70)
    println("SLIDER-CRANK ODE FORM (14 variables)")
    println("="^70)
    println("Integrating with Tsit5 (non-stiff ODE solver)...")
    
    try
        sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-6)
        
        println("✓ Solution computed!")
        println("  Steps: $(length(sol))")
        println("  φ₁(0.1) = $(sol.u[end][1]) rad (expected ≈ 15.0)")
        println("  x₃(0.1) = $(sol.u[end][3]) m")
        println()
        
        # Plot
        p = plot(sol, vars=(1,8), label="φ₁ vs v₁", xlabel="t", legend=false, size=(900, 300))
        savefig(p, "slider_crank_ode_solution.png")
        println("  Saved plot to slider_crank_ode_solution.png")
        
    catch e
        println("✗ Failed: $e")
    end
end
