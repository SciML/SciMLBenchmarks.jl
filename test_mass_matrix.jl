#!/usr/bin/env julia
# Diagnostic script to check AM matrix and initial consistency

using LinearAlgebra

# Physical Constants (EXACT from crank.f)
const M1 = 0.36
const M2 = 0.151104
const M3 = 0.075552
const L1 = 0.15
const L2 = 0.30
const J1 = 0.002727
const J2 = 0.0045339259
const PI = 3.1415927
const EE = 0.20e12
const NUE = 0.30
const BB = 0.0080
const HH = 0.0080
const RHO = 7870.0
const GRAV = 0.0
const OMEGA = 150.0

# Grid parameters
const NQ = 4
const NP = 7
const NL = 3
const NX = 3*NP + NL
const KU = 4
const KV = 0

# Initialize FE matrices
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

# Get initial conditions
y = zeros(24)
y[1] = 0.0
y[2] = 0.0
y[3] = 0.450016933e+00
y[4:7] .= [0.0, 0.0, 0.103339863e-04, 0.169327969e-04]
y[8] = 0.150000000e+03
y[9] = -0.749957670e+02
y[10:14] .= [-0.268938672e-05, 0.444896105e+00, 0.463434311e-02, -0.178591076e-05, -0.268938672e-05]
y[15] = 0.0
y[16] = -0.1344541576008661e-03
y[17] = -0.5062194923138079e+03
y[18:21] .= [-0.6833142732779555e-05, 0.1449382650173157e-08, -0.4268463211410861e+00, 0.2098334687947376e-01]
y[22] = -0.6397251492537153e-08
y[23] = 0.3824589508329281e+02
y[24] = -0.4376060460948886e-09

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

# Compute mass matrix at t=0
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

println("="^70)
println("MASS MATRIX AM @ t=0")
println("="^70)
println("\nFull 7×7 Matrix:")
show(stdout, "text/plain", AM)
println("\n")

println("\nRigid block (1:3, 1:3):")
show(stdout, "text/plain", AM[1:3, 1:3])
println("\n")

println("\nElastic block (4:7, 4:7):")
show(stdout, "text/plain", AM[4:7, 4:7])
println("\n")

println("\nRigid-Elastic coupling (1:3, 4:7):")
show(stdout, "text/plain", AM[1:3, 4:7])
println("\n")

println("Determinant: $(det(AM))")
println("Condition number: $(cond(AM))")
println("Min eigenvalue: $(minimum(eigvals(AM)))")
println("Max eigenvalue: $(maximum(eigvals(AM)))")

println("\n" * "="^70)
println("EIGENVALUE ANALYSIS")
println("="^70)
eigs = sort(eigvals(AM))
for i in 1:length(eigs)
    println("λ[$i] = $(eigs[i])")
end

# Constraint Jacobian at t=0
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

println("\n" * "="^70)
println("CONSTRAINT JACOBIAN GP @ t=0")
println("="^70)
show(stdout, "text/plain", GP)
println("\n")

println("Rank of GP: $(rank(GP))")
println("Full rank? $(rank(GP) == 3)")

# Augmented matrix for consistency
Aug = vcat(AM, GP)
println("\n" * "="^70)
println("AUGMENTED MATRIX [AM; GP] (10×7)")
println("="^70)
show(stdout, "text/plain", Aug)
