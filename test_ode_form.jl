#!/usr/bin/env julia
"""
Build ODE mass-matrix form for slider-crank system.

Instead of solving the DAE:
  ẏ₁ = y₂
  ẏ₂ = y₃
  M*ẏ₃ = f - G'*λ
  0 = G*y₂ + r'(t)

We solve the ODE mass-matrix form:
  du/dt = (1/M) * (f - G' * (G*M^{-1}*G')^{-1} * (G*M^{-1}*f + G*u_v + r'))

This avoids DAE initialization issues and is often more robust for constrained systems.
"""

using LinearAlgebra, OrdinaryDiffEq, Sundials

# Physical Constants
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
const NL = 3
const NX = 3*NP + NL
const KU = 4
const KV = 0

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

MQ, KQ, BQ, DQ, c1, c2, c12, c21 = initialize_fe_matrices()

function compute_mass_matrix_and_forces(y::Vector, t::Float64)
    """Compute mass matrix M(y) and force vector f(y,t) for current state."""
    
    p1, p2, x3 = y[1], y[2], y[3]
    q = y[4:7]
    v1, v2, vx3 = y[8], y[9], y[10]
    vq = y[11:14]
    
    cosp1 = cos(p1)
    cosp2 = cos(p2)
    sinp1 = sin(p1)
    sinp2 = sin(p2)
    cosp12 = cos(p1 - p2)
    sinp12 = sin(p1 - p2)
    
    qku = (KU == 0) ? 0.0 : q[KU]
    qkv = (KV == 0) ? 0.0 : q[KV]
    
    # Scalar products
    c1Tq = dot(c1, q)
    c1Tqd = dot(c1, vq)
    c2Tq = dot(c2, q)
    c2Tqd = dot(c2, vq)
    c12Tq = dot(c12, q)
    c12Tqd = dot(c12, vq)
    
    MQq = MQ * q
    KQq = KQ * q
    DQqd = DQ * vq
    BQqd = BQ * vq
    
    qtmqq = dot(q, MQq)
    qdtmqq = dot(vq, MQq)
    qdtbqqd = dot(vq, BQqd)
    
    QtBQ = zeros(NQ)
    for i in 1:NQ
        QtBQ[i] = dot(q, BQ[:, i])
    end
    
    # Mass matrix M(y) - 7x7
    M = zeros(7, 7)
    
    M[1,1] = J1 + M2 * L1^2
    M[1,2] = 0.5 * L1 * L2 * M2 * cosp12
    M[2,2] = J2
    M[3,3] = M3
    
    M[1,2] += RHO * L1 * (sinp12 * c2Tq + cosp12 * c1Tq)
    M[2,2] += qtmqq + 2.0 * RHO * c12Tq
    
    for i in 1:NQ
        M[1, 3+i] = RHO * L1 * (-sinp12 * c1[i] + cosp12 * c2[i])
        M[2, 3+i] = RHO * c21[i] + RHO * QtBQ[i]
        M[3, 3+i] = 0.0
    end
    
    for i in 1:NQ
        for j in 1:i
            M[3+j, 3+i] = MQ[j, i]
        end
    end
    
    # Symmetrize
    for i in 1:7
        for j in i+1:7
            M[j, i] = M[i, j]
        end
    end
    
    # Constraint Jacobian G(y) - 3x7
    G = zeros(3, 7)
    
    G[1, 1] = L1 * cosp1
    G[1, 2] = L2 * cosp2 + qku * cosp2 - qkv * sinp2
    G[2, 1] = L1 * sinp1
    G[2, 2] = L2 * sinp2 + qku * sinp2 + qkv * cosp2
    G[2, 3] = 1.0
    G[3, 1] = 1.0
    
    if KU != 0
        G[1, 3+KU] = sinp2
        G[2, 3+KU] = -cosp2
    end
    if KV != 0
        G[1, 3+KV] = cosp2
        G[2, 3+KV] = sinp2
    end
    
    # Force vector f(y,t) - 7x1
    f = zeros(7)
    
    f[1] = -0.5 * L1 * GRAV * (M1 + 2.0*M2) * cosp1 - 0.5 * L1 * L2 * M2 * v2^2 * sinp12
    f[2] = -0.5 * L2 * GRAV * M2 * cosp2 + 0.5 * L1 * L2 * M2 * v1^2 * sinp12
    f[3] = 0.0
    
    f[1] += RHO * L1 * v2^2 * (-sinp12 * c1Tq + cosp12 * c2Tq) - 2.0 * RHO * L1 * v2 * (cosp12 * c1Tqd + sinp12 * c2Tqd)
    
    f[2] += RHO * L1 * v1^2 * (sinp12 * c1Tq - cosp12 * c2Tq) - 2.0 * RHO * v2 * c12Tqd - 2.0 * v2 * qdtmqq - RHO * qdtbqqd - RHO * GRAV * (cosp2 * c1Tq - sinp2 * c2Tq)
    
    for i in 1:NQ
        f[3+i] = v2^2 * MQq[i]
        f[3+i] += RHO * (v2^2 * c12[i] + L1 * v1^2 * (cosp12 * c1[i] + sinp12 * c2[i]) + 2.0 * v2 * BQqd[i])
        f[3+i] -= RHO * GRAV * (sinp2 * c1[i] + cosp2 * c2[i])
    end
    
    for i in 1:NQ
        f[3+i] -= KQq[i] + dot(DQ[i,:], vq)
    end
    
    # Time derivative of constraint: r'(t) = [0, 0, -OMEGA]'
    r_prime = [0.0, 0.0, -OMEGA]
    
    return M, G, f, r_prime
end

function slider_crank_ode_mass_matrix!(du, u, p, t)
    """
    ODE mass-matrix form RHS.
    
    State: u = [p, v] (14 variables)
    - p: 7 positions/coordinates
    - v: 7 velocities
    
    Solves: M(p)*dv/dt = f(p,v) - G(p)'*(G*M^{-1}*G')^{-1}*(G*M^{-1}*f + G*v + r'(t))
    """
    
    # Extract state
    p_vec = u[1:7]
    v_vec = u[8:14]
    
    # Extend to full 24-var state for mass/force computation
    y_full = zeros(24)
    y_full[1:7] = p_vec
    y_full[8:14] = v_vec
    y_full[15:21] .= 0.0  # Not needed for M,G,f
    y_full[22:24] .= 0.0  # Not needed for M,G,f
    
    # Compute mass matrix, constraint Jacobian, forces
    M, G, f, r_prime = compute_mass_matrix_and_forces(y_full, t)
    
    # Compute the constraint force term
    # λ = (G*M^{-1}*G')^{-1} * (G*M^{-1}*f + G*v + r'(t))
    
    # Initialize acceleration
    a = zeros(7)
    
    try
        M_inv = inv(M)
        GG_prod = G * M_inv * G'
        GG_inv = inv(GG_prod)
        
        rhs_constraint = G * M_inv * f + G * v_vec + r_prime
        lambda = GG_inv * rhs_constraint
        
        # Compute acceleration: M*a = f - G'*lambda
        a = M_inv * (f - G' * lambda)
        
    catch err
        # Fallback if matrix inversion fails (shouldn't happen in normal operation)
        a = zeros(7)
        @warn "Matrix inversion failed at t=$t: $err"
    end
    
    # Fill RHS: dp/dt = v, dv/dt = a
    du[1:7] = v_vec
    du[8:14] = a
end

# Initial conditions
y0 = zeros(14)
y0[1] = 0.0
y0[2] = 0.0
y0[3] = L1 + L2 + 0.169327969e-04
y0[4:7] .= [0.0, 0.0, 0.103339863e-04, 0.169327969e-04]
y0[8]  = 150.0
y0[9]  = -74.9957670
y0[10:14] .= [-0.268938672e-05, 0.444896105, 0.463434311e-02, -0.178591076e-05, -0.268938672e-05]

tspan = (0.0, 0.1)

# Create and solve ODE
prob_ode = ODEProblem(slider_crank_ode_mass_matrix!, y0, tspan)

println("="^70)
println("SOLVING ODE MASS-MATRIX FORM")
println("="^70)
println("Initial state:")
println("  p = [$(y0[1:3])]")
println("  v = [$(y0[8:10])]")

try
    sol = solve(prob_ode, CVODE_BDF(); abstol=1e-8, reltol=1e-6)
    
    println("\n✓ SOLUTION SUCCESSFUL!")
    println("  Time steps: $(length(sol.t))")
    println("  Final time: $(sol.t[end])")
    
    println("\nFinal state (t=$(sol.t[end])):")
    println("  Positions:")
    println("    φ₁ = $(sol.u[end][1])")
    println("    φ₂ = $(sol.u[end][2])")
    println("    x₃ = $(sol.u[end][3])")
    println("  Velocities:")
    println("    ω₁ = $(sol.u[end][8])")
    println("    ω₂ = $(sol.u[end][9])")
    println("    v_x₃ = $(sol.u[end][10])")
    
    println("\n" * "="^70)
    println("COMPARISON TO REFERENCE (crank.f, MEBDFI)")
    println("="^70)
    println("x₃(0.1) computed = $(sol.u[end][3])")
    println("x₃(0.1) reference ≈ 0.169737")
    println("Error = $(abs(sol.u[end][3] - 0.169737))")
    
catch e
    println("\n✗ Solve failed: $e")
    println("\nTrying with looser tolerances...")
    try
        sol = solve(prob_ode, CVODE_BDF(); abstol=1e-4, reltol=1e-2)
        println("✓ Solved with loose tolerances")
        println("Final x₃(0.1) = $(sol.u[end][3])")
    catch e2
        println("✗ Also failed: $e2")
    end
end
