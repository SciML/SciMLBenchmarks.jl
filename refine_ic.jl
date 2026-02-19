using LinearAlgebra

const L1 = 0.15
const L2 = 0.30
const KU = 4

# Start IC
p1_0 = 0.0
p2_0 = 0.0
q_0 = [0.0, 0.0, 0.103339863e-04, 0.169327969e-04]

# Elastic coupling
q4_0 = q_0[4]  # axial deformation

# Constraint 1 (holonomic): g₁ = 0
# x₃ = L₁·cos(φ₁) + (L₂ + q₄)·cos(φ₂)
# At t=0: x₃ = L₁·1 + (L₂ + q₄)·1 = L₁ + L₂ + q₄

x3_correct = L1 + L2 + q4_0

println("Constraint geometry check:")
println("p1 = $p1_0, p2 = $p2_0")
println("q4 = $(q4_0)")
println()
println("Current IC x3 = $(L1 + L2 + 0.169327969e-04)")
println("Required IC x3 = $x3_correct for g₁(p₀) = 0")
println()
println("Difference: $(abs((L1 + L2 + 0.169327969e-04) - x3_correct))")
println()

# Check constraint residual with corrected x3
x3_0_new = x3_correct

cosp1 = cos(p1_0)
cosp2 = cos(p2_0)
sinp1 = sin(p1_0)
sinp2 = sin(p2_0)

# Constraint Jacobian rows at t=0
g1_row = [L1 * cosp1,  L2 * cosp2 + q4_0 * cosp2,  0.0,  0.0,  0.0,  0.0,  sinp2]
g2_row = [L1 * sinp1,  L2 * sinp2 + q4_0 * sinp2,  1.0,  0.0,  0.0,  0.0, -cosp2]
g3_row = [1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]

println("Constraint Jacobian rows at t=0:")
println("G[1,:] = $g1_row")
println("G[2,:] = $g2_row")
println("G[3,:] = $g3_row")
println()

# Test velocities
v1_0 = 150.0
v2_0 = -74.9957670
vx3_0 = -0.268938672e-05
vq_0 = [0.444896105, 0.463434311e-02, -0.178591076e-05, -0.268938672e-05]

v_full = [v1_0; v2_0; vx3_0; vq_0]

gv1 = dot(g1_row, v_full)
gv2 = dot(g2_row, v_full)
gv3 = dot(g3_row, v_full)

println("Constraint velocities G(p)*v at t=0:")
println("G[1,:]*v = $(gv1)")
println("G[2,:]*v = $(gv2)")
println("G[3,:]*v = $(gv3)")
println()

OMEGA = 150.0
println("Expected constraint residual (should be near 0):")
println("residual[22] = G[1,:]*v        = $(gv1)      (should be ≈ 0)")
println("residual[23] = G[2,:]*v        = $(gv2)      (should be ≈ 0)")
println("residual[24] = G[3,:]*v - Ω    = $(gv3 - OMEGA)  (should be ≈ 0)")
println()

# Summary
println("="^70)
println("CORRECTED INITIAL CONDITION")
println("="^70)
println("y[3] (x₃) should be: $x3_correct")
println("Instead of:        $(L1 + L2 + 0.169327969e-04)")
