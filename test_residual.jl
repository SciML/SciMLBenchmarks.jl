#!/usr/bin/env julia
# Diagnostic: Check if residual function gives zero residual at initial conditions

using LinearAlgebra

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

function get_initial_conditions()
    y = zeros(NX)
    y[1] = 0.0
    y[2] = 0.0
    y[3] = L1 + L2 + 0.169327969e-04   # Position constraint satisfied exactly
    y[4:7] .= [0.0, 0.0, 0.103339863e-04, 0.169327969e-04]
    y[8]  = 150.0
    y[9]  = -74.9957670
    y[10:14] .= [-0.268938672e-05, 0.444896105, 0.463434311e-02, -0.178591076e-05, -0.268938672e-05]
    y[15] = 150.0
    y[16] = -74.99576703969453
    y[17] = -1065.7243315012481
    y[18] = -14.324753102521305
    y[19] = -7.161568115404762
    y[20] = 133.27763907453476
    y[21] = -1065.7243315012481
    y[22] = -0.049729013574148544
    y[23] = 80.5176046935823
    y[24] = -0.44988495614760715
    return y
end

function slider_crank_resmbs!(residual, y, yprime, t)
    p1, p2, x3 = y[1], y[2], y[3]
    q = y[4:7]
    v1, v2, vx3 = y[8], y[9], y[10]
    vq = y[11:14]
    w1, w2, wx3 = y[15], y[16], y[17]
    wq = y[18:21]
    lambda1, lambda2, lambda3 = y[22], y[23], y[24]
    
    cosp1 = cos(p1)
    cosp2 = cos(p2)
    sinp1 = sin(p1)
    sinp2 = sin(p2)
    cosp12 = cos(p1 - p2)
    sinp12 = sin(p1 - p2)
    
    qku = (KU == 0) ? 0.0 : q[KU]
    qkv = (KV == 0) ? 0.0 : q[KV]
    
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
    
    # Kinematic equations dp/dt = v
    residual[1:7] = yprime[1:7] - y[8:14]
    
    # Velocity equations dv/dt = w
    residual[8:14] = yprime[8:14] - y[15:21]
    
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
    
    # Constraint matrix
    GP = zeros(3, NP)
    
    GP[1, 1] = L1 * cosp1
    GP[1, 2] = L2 * cosp2 + qku * cosp2 - qkv * sinp2
    GP[2, 1] = L1 * sinp1
    GP[2, 2] = L2 * sinp2 + qku * sinp2 + qkv * cosp2
    GP[2, 3] = 1.0
    GP[3, 1] = 1.0
    
    if KU != 0
        GP[1, 3+KU] = sinp2
        GP[2, 3+KU] = -cosp2
    end
    if KV != 0
        GP[1, 3+KV] = cosp2
        GP[2, 3+KV] = sinp2
    end
    
    # Force vector
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
        F[3+i] -= KQq[i] + dot(DQ[i,:], vq)
    end
    
    # Dynamics: M*w - f + G'*lambda = 0
    AMw = AM * y[15:21]
    for i in 1:NP
        residual[14+i] = AMw[i] - F[i] + GP[1,i]*lambda1 + GP[2,i]*lambda2 + GP[3,i]*lambda3
    end
    
    # Constraints: G(p)*v + r'(t) = 0
    qdku = (KU == 0) ? 0.0 : vq[KU]
    qdkv = (KV == 0) ? 0.0 : vq[KV]
    
    vlc1 = 0.0
    vlc2 = 0.0
    vlc3 = -OMEGA
    
    for i in 1:NP
        vlc1 += GP[1, i] * y[7+i]
        vlc2 += GP[2, i] * y[7+i]
        vlc3 += GP[3, i] * y[7+i]
    end
    
    residual[22] = vlc1
    residual[23] = vlc2
    residual[24] = vlc3
end

# Get initial condition with computed consistent IC
y0 = get_initial_conditions()

# Compute yprime from dp/dt = v, dv/dt = w structure
yprime0 = zeros(24)
yprime0[1:7] = y0[8:14]     # dp/dt = v
yprime0[8:14] = y0[15:21]   # dv/dt = w
# yprime0[15:21] don't have standard meaning for algebraic vars
# yprime0[22:24] don't have standard meaning for algebraic vars

println("="^70)
println("RESIDUAL EVALUATION AT COMPUTED CONSISTENT INITIAL CONDITION")
println("="^70)

residual = zeros(24)
slider_crank_resmbs!(residual, y0, yprime0, 0.0)

println("\nRESIDUALS:")
println("\nRows 1-7 (kinematic dp/dt = v):")
for i in 1:7
    println("  r[$i] = $(residual[i])")
end

println("\nRows 8-14 (velocity dv/dt = w):")
for i in 8:14
    println("  r[$i] = $(residual[i])")
end

println("\nRows 15-21 (dynamics M*w = f + G'*λ):")
for i in 15:21
    println("  r[$i] = $(residual[i])")
end

println("\nRows 22-24 (constraints G*v + r' = 0):")
for i in 22:24
    println("  r[$i] = $(residual[i])")
end

max_res = maximum(abs.(residual))
println("\n" * "="^70)
println("MAX |RESIDUAL| = $max_res")
if max_res < 1e-12
    println("✓ Residual essentially zero - system is consistent!")
else
    println("✗ Residual NOT zero - indicates inconsistency or bug in formulation")
end
