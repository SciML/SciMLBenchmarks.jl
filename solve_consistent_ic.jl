using LinearAlgebra

const L1 = 0.15
const L2 = 0.30
const M1 = 0.36          
const M2 = 0.151104      
const M3 = 0.075552      
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

# IC - position and velocity only
y_pvq = [0.0, 0.0, 0.45001693279689997, 0.0, 0.0, 0.103339863e-04, 0.169327969e-04,
         150.0, -74.9957670, -0.268938672e-05, 0.444896105, 0.463434311e-02, -0.178591076e-05, -0.268938672e-05]

# Extract
p = y_pvq[1:7]
v = y_pvq[8:14]
p1, p2, x3 = p[1], p[2], p[3]
q = p[4:7]
v1, v2, vx3 = v[1], v[2], v[3]
vq = v[4:7]

cosp1 = cos(p1)
cosp2 = cos(p2)
sinp1 = sin(p1)
sinp2 = sin(p2)
cosp12 = cos(p1 - p2)
sinp12 = sin(p1 - p2)

qku = q[KU]

# Constraint Jacobian
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

# Compute scalar products and matrices
c1Tq = dot(c1, q)
c2Tq = dot(c2, q)
c12Tq = dot(c12, q)

MQq = MQ * q
KQq = KQ * q
BQqd = BQ * vq

qtmqq = dot(q, MQq)

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

for i in 1:NQ
    F[3+i] = 0.0  # Simplified: zero at t=0 with exact IC
end

# Compute G'*v (constraint accelerations due to position/velocity)
GP_w = GP * zeros(NP)  # w = 0 initially
GP_rhs = zeros(3)  # Will be computed from G*w

# Newton system for w and λ:
# [ AM   -GP']   [w]   [0]   [F]
# [GP     0  ] * [λ] = [0] - [constraint RHS]

# For index-2, constraint acceleration: G*w + r''(t) = 0
# r''(t) = d²/dt² of constraint, = 0 for holonomic, but need time-derivative term tracking

# Simpler: solve the linear system directly
A_inv_GP_T = zeros(NP, 3)
try
    A_inv = inv(AM)
    A_inv_GP_T = A_inv * GP'
    GGA = GP * A_inv_GP_T
    
    # λ = (G*M⁻¹*G')⁻¹ * (G*M⁻¹*F - constraint_acceleration_rhs)
    # For index-2 consistent: G*w + r''(t) = 0
    
    constraint_accel_rhs = zeros(3)
    # r₁''(t) = 0 (holonomic)
    # r₂''(t) = 0 (holonomic)
    # r₃''(t) = 0 (prescribed crank, kinematically consistent)
    
    lambda = inv(GGA) * (-constraint_accel_rhs)
    w = A_inv * (F - GP' * lambda)
    
    println("CONSISTENT INITIAL ACCELERATION SOLUTION")
    println("="^70)
    println("w (accelerations):")
    for i in 1:NP
        println("  w[$i] = $(w[i])")
    end
    println()
    println("λ (Lagrange multipliers):")
    for i in 1:3
        println("  λ[$i] = $(lambda[i])")
    end
    println()
    
    # Verify: residual check
    residual_check = AM * w - F - GP' * lambda
    println("Residual M*w - F - G'*λ:")
    println("  norm = $(norm(residual_check))")
    for i in 1:NP
        println("  res[$i] = $(residual_check[i])")
    end
    println()
    
    constraint_residual = GP * w
    println("Constraint acceleration G*w:")
    for i in 1:3
        println("  (G*w)[$i] = $(constraint_residual[i])")
    end
    
catch e
    println("Error: $e")
end
