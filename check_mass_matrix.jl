using LinearAlgebra

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
const NQ = 4
const NP = 7

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
    
    c21 = zeros(NQ)
    c21[1] = L2 * FACB * 1.0 / PI
    c21[2] = -L2 * FACB * 0.5 / PI
    
    return MQ, c21
end

MQ, c21 = initialize_fe_matrices()

# At t=0
p1, p2 = 0.0, 0.0
cosp12 = cos(0.0)  # 1
sinp12 = sin(0.0)  # 0
q = zeros(4)
c1 = zeros(4)
c2 = zeros(4)
c12 = zeros(4)
QtBQ = zeros(4)

AM = zeros(NP, NP)

AM[1,1] = J1 + M2 * L1^2
AM[1,2] = 0.5 * L1 * L2 * M2 * cosp12
AM[2,2] = J2
AM[3,3] = M3

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

println("Mass matrix AM at t=0:")
println(AM)
println()

println("Eigenvalues:")
eigvals_am = eigvals(AM)
println(eigvals_am)
println()

println("Condition number: $(cond(AM))")
println("Rank: $(rank(AM))")
println("Determinant: $(det(AM))")
println()

if rank(AM) < NP
    println("⚠️ WARNING: Mass matrix is singular!")
    println("   Rank deficiency = $(NP - rank(AM))")
else
    println("✓ Mass matrix is non-singular")
end
