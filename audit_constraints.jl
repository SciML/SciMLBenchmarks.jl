using LinearAlgebra

# Check constraint geometry at t=0
const L1 = 0.15
const L2 = 0.30

# From your IC
p1_0 = 0.0
p2_0 = 0.0
x3_0 = L1 + L2 + 0.169327969e-04
q_0 = [0.0, 0.0, 0.103339863e-04, 0.169327969e-04]

v1_0 = 150.0
v2_0 = -74.9957670
vx3_0 = -0.268938672e-05
vq_0 = [0.444896105, 0.463434311e-02, -0.178591076e-05, -0.268938672e-05]

# Constraint Jacobian at t=0
cosp1 = cos(p1_0)  # = 1
cosp2 = cos(p2_0)  # = 1
sinp1 = sin(p1_0)  # = 0
sinp2 = sin(p2_0)  # = 0

KU = 4
KV = 0
qku = q_0[KU]
qkv = (KV == 0) ? 0.0 : q_0[KV]

GP = zeros(3, 7)

# Row 1: ∂g₁/∂p
GP[1, 1] = L1 * cosp1              # = 0.15
GP[1, 2] = L2 * cosp2 + qku * cosp2 - qkv * sinp2  # = 0.30 + 0.169327969e-04 - 0
GP[1, 3] = 0.0

if KU != 0
    GP[1, 3+KU] = sinp2  # = 0
end

# Row 2: ∂g₂/∂p
GP[2, 1] = L1 * sinp1              # = 0
GP[2, 2] = L2 * sinp2 + qku * sinp2 + qkv * cosp2  # = 0 + 0 + 0
GP[2, 3] = 1.0

if KU != 0
    GP[2, 3+KU] = -cosp2  # = -1
end

# Row 3: ∂g₃/∂p (torque constraint - dummy)
GP[3, 1] = 1.0
GP[3, 2] = 0.0
GP[3, 3] = 0.0

println("Constraint Jacobian GP at t=0:")
println(GP)
println()

# Compute G*v
Gv = GP * [v1_0; v2_0; vx3_0; vq_0]
println("G(p₀) * v₀:")
println("  g₁: $(Gv[1])")
println("  g₂: $(Gv[2])")
println("  g₃: $(Gv[3])")
println()

# Constraint residuals
# For holonomic constraints g(p) = 0
# Velocity constraint is d/dt g(p) = G(p) v + r'(t) = 0

# r'(t) values (time derivatives of prescribed constraints)
# In crank.f: r'(t) should encode the prescribed velocity
# For crank angle: d(φ₁)/dt = Ω = 150 rad/s is PRESCRIBED

println("Expected velocity constraints:")
println("  Row 1: g₁(p) = 0 → d/dt[g₁] = G[1,·]·v + r₁'(t) = 0")
println("         $(Gv[1]) + r₁'(t) = 0  →  r₁'(t) = $(−Gv[1])")
println()
println("  Row 2: g₂(p) = 0 → d/dt[g₂] = G[2,·]·v + r₂'(t) = 0")
println("         $(Gv[2]) + r₂'(t) = 0  →  r₂'(t) = $(−Gv[2])")
println()
println("  Row 3: φ₁ is prescribed → dφ₁/dt = Ω = 150 rad/s = -r₃'(t)")
println("         φ₁·velocity_prescribed: v₁ = 150 (matches! ✓)")
println()

# In residual form, the constraint should be:
# residual[22:24] = G(p)*v + r'(t) = [g₁'(p)*v + r₁'(t), g₂'(p)*v + r₂'(t), g₃'(p)*v + r₃'(t)]
# 
# Given that v₁ = 150 = Ω and Gv[1] ≈ 0.45 (centrifugal coupling)
# There's a geometric inconsistency

println("="^70)
println("CONSTRAINT CONSISTENCY ANALYSIS")
println("="^70)
println()
println("Constraint equations (holonomic + 1 prescribed):")
println("  g₁(p): x₃ = L₁·cos(φ₁) + (L₂ + q₄)·cos(φ₂)")
println("  g₂(p): 0  = L₁·sin(φ₁) + (L₂ + q₄)·sin(φ₂) - x₃")
println("  g₃: φ₁ - Ω·t = 0 (prescribed crank)")
println()
println("Velocity constraints (from differentiation):")
println("  G[1,:] · v + r'₁(t) = 0")
println("  G[2,:] · v + r'₂(t) = 0")
println("  G[3,:] · v + r'₃(t) = 0")
println()
println("Where r'ᵢ(t) come from d/dt of implicit time dependence in constraints")
println()

# The issue: r'(t) term is NOT being included in the residual!
println("⚠️  CRITICAL ISSUE:")
println("Your constraint residual only computes G·v")
println("But should compute: G·v + r'(t)")
println()
println("From crank.f, the r'(t) term encodes:")
println("  r'₁(t) = 0  (holonomic, no explicit time dependence)")
println("  r'₂(t) = 0  (holonomic, no explicit time dependence)")
println("  r'₃(t) = -Ω = -150  (prescribed angular velocity, NOT holonomic)")
println()
println("So residual[24] should be: G[3,:]·v + r'₃(t)")
println("                         = $(Gv[3]) + (-150)")
println("                         = $(Gv[3] - 150)")
println()
println("Your current residual[24] computes:")
println("vlc3 = -OMEGA = -150")
println("Then residual[24] = vlc3 = -150")
println("But should = G[3,:]·v - OMEGA = $(Gv[3]) - 150 = $(Gv[3] - 150)")
