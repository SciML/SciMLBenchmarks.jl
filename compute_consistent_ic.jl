#!/usr/bin/env julia
# Compute consistent initial accelerations by solving M*w = f + G'*lambda

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

# Initial conditions (p,q,v,u as given by crank.f, CORRECTED position)
y = zeros(24)
y[1] = 0.0
y[2] = 0.0
y[3] = L1 + L2 + 0.169327969e-04  # CORRECTED: position constraint satisfied exactly
y[4:7] .= [0.0, 0.0, 0.103339863e-04, 0.169327969e-04]
y[8] = 0.150000000e+03
y[9] = -0.749957670e+02
y[10:14] .= [-0.268938672e-05, 0.444896105e+00, 0.463434311e-02, -0.178591076e-05, -0.268938672e-05]
# Accelerations initially guess zero
y[15:21] .= 0.0
# Multipliers guess
y[22] = 0.0
y[23] = 0.0
y[24] = 0.0

p1, p2, x3 = y[1], y[2], y[3]
q = y[4:7]
v1, v2, vx3 = y[8], y[9], y[10]
vq = y[11:14]
w = y[15:21]

cosp1 = cos(p1)
cosp2 = cos(p2)
sinp1 = sin(p1)
sinp2 = sin(p2)
cosp12 = cos(p1 - p2)
sinp12 = sin(p1 - p2)

qku = (KU == 0) ? 0.0 : q[KU]
qkv = (KV == 0) ? 0.0 : q[KV]

# Compute masses
c1Tq = dot(c1, q)
c2Tq = dot(c2, q)
c12Tq = dot(c12, q)
MQq = MQ * q
qtmqq = dot(q, MQq)

QtBQ = zeros(NQ)
for i in 1:NQ
    QtBQ[i] = dot(q, BQ[:, i])
end

AM = zeros(7, 7)

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

for i in 1:7
    for j in i+1:7
        AM[j, i] = AM[i, j]
    end
end

# Compute constraint Jacobian
GP = zeros(3, 7)

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

# Compute force vector
KQq = KQ * q

F = zeros(7)
F[1] = -0.5 * L1 * GRAV * (M1 + 2.0*M2) * cosp1 - 0.5 * L1 * L2 * M2 * v2^2 * sinp12
F[2] = -0.5 * L2 * GRAV * M2 * cosp2 + 0.5 * L1 * L2 * M2 * v1^2 * sinp12
F[3] = 0.0

c1Tqd = dot(c1, vq)
c2Tqd = dot(c2, vq)
c12Tqd = dot(c12, vq)
qdtmqq = dot(vq, MQq)
BQqd = BQ * vq
qdtbqqd = dot(vq, BQqd)

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

println("="^70)
println("COMPUTING CONSISTENT INITIAL CONDITIONS")
println("="^70)

# At t=0, we need: G * v = -r'(t)
# For phi_1 = OMEGA * t, we have: phi_1 - OMEGA*t = 0, so d/dt = -OMEGA (at t=0, r' = -OMEGA)

r_prime = zeros(3)
r_prime[1] = 0.0  # x_h constraint
r_prime[2] = 0.0  # y_h constraint
r_prime[3] = -OMEGA  # phi_1 = OMEGA*t

println("\nConstraint right-hand side (r'(t)):")
for i in 1:3
    println("  r_prime[$i] = $(r_prime[i])")
end

# Velocity constraint: G(p) * v + r' = 0
constraint_velocity = zeros(3)
for i in 1:7
    constraint_velocity[1] += GP[1,i] * y[7+i]
    constraint_velocity[2] += GP[2,i] * y[7+i]
    constraint_velocity[3] += GP[3,i] * y[7+i]
end

println("\nVelocity constraint residual G(p)*v - r'(t):")
for i in 1:3
    println("  G[$i]*v - r_prime[$i] = $(constraint_velocity[i] - r_prime[i])")
end

# Augmented system: [AM, G'; G, 0] * [w; lambda] = [f + 0; -r_prime]
println("\n" * "="^70)
println("SOLVING AUGMENTED SYSTEM FOR CONSISTENT w AND lambda")
println("="^70)

# Build augmented matrix
Aug = zeros(10, 10)
Aug[1:7, 1:7] = AM
Aug[1:7, 8:10] = GP'
Aug[8:10, 1:7] = GP
Aug[8:10, 8:10] .= 0.0

# Build RHS
rhs = zeros(10)
rhs[1:7] = F
rhs[8:10] = -r_prime

println("\nAugmented matrix rank: $(rank(Aug))")
println("Aug determinant: $(det(Aug))")

# Solve
try
    sol = Aug \ rhs
    w_consistent = sol[1:7]
    lambda_consistent = sol[8:10]
    
    println("\n✓ Consistent initial conditions found!")
    println("\nConsistent accelerations w:")
    for i in 1:7
        println("  w[$i] = $(w_consistent[i])")
    end
    
    println("\nConsistent Lagrange multipliers λ:")
    for i in 1:3
        println("  λ[$i] = $(lambda_consistent[i])")
    end
    
    # Verify
    residual_dynamics = AM * w_consistent + GP' * lambda_consistent - F
    println("\n Residual in dynamics M*w + G'*λ - f:")
    for i in 1:7
        println("  res[$i] = $(residual_dynamics[i]) (should be ~0)")
    end
    
    # Check constraint 
    residual_constraint = GP * w_consistent + [0,0,-OMEGA]
    println("\nResidual in constraint G*w + r''(t):")
    for i in 1:3
        println("  res_constraint[$i] = $(residual_constraint[i]) (should be ~0)")
    end
    
catch e
    println("✗ Failed to solve: $e")
end
