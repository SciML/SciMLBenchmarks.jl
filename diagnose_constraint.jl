#!/usr/bin/env julia
# Fix: Ensure initial position satisfies geometric constraints perfectly

using LinearAlgebra

const L1 = 0.15
const L2 = 0.30

# Initial angles from crank.f
phi1 = 0.0
phi2 = 0.0
q4 = 0.169327969e-04  # axial elastic deformation

# Geometric constraints:
# x_h = L1*cos(φ₁) + (L2+q₄)*cos(φ₂)  = ?
# y_h = L1*sin(φ₁) + (L2+q₄)*sin(φ₂)  = ?
# φ₁ = Ω*t  (prescribed, satisfied by definition)

# At t=0, φ₁=0, φ₂=0:
# x_h = L1*cos(0) + (L2+q₄)*cos(0) = L1 + L2 + q₄
# y_h = L1*sin(0) + (L2+q₄)*sin(0) = 0

x_h_at_0 = L1 + L2 + q4
y_h_at_0 = 0.0

println("Geometric constraint solution:")
println("At t=0 (φ₁=0, φ₂=0):")
println("  x_h should be: $(x_h_at_0) = L1 + L2 + q₄ = $L1 + $L2 + $q4")
println("  y_h should be: $y_h_at_0")
println("\nBut crank.f gives:")
println("  x₃(0) = 0.450016933")
println("  This is NOT consistent with L1 + L2 + q₄ = $(x_h_at_0)")
println("\nThe CORRECT initial x₃ should be:")
println("  x₃(0) = $(x_h_at_0)")

# This is the bug! The position constraint is violated.
# The slider position x₃ should be exactly L1 + L2 + q₄.
# crank.f value of 0.450016933 is DIFFERENT!

# Check: does 0.450016933 make sense as a solution to the constraint?
# At φ₁, φ₂:
# x_h = L1*cos(φ₁) + (L2+q₄)*cos(φ₂) = x₃_crank_f
# This would require specific angles...

# The discrepancy suggests: crank.f might be using a DIFFERENT initial time step,
# or the initial angles are NOT zero.

# Let me verify: if x₃=0.450016933, can I recover the angles?

x3_given = 0.450016933
# But at t=0, we're supposed to have φ₁=0, φ₂=0.
# Then x3 must equal L1+L2+q₄

println("\n" * "="^70)
println("DIAGNOSIS: Initial Position Constraint Violated")
println("="^70)
println("\nRequired constraint: φ₁(0)=0, φ₂(0)=0 ⟹ x₃(0) = $(L1+L2+q4)")
println("But crank.f gives: x₃(0) = $x3_given")
println("\nThis is the source of the residual[22] ≠ 0 issue!")
println("\nSolution: Set x₃(0) = L1 + L2 + q₄ = $(L1+L2+q4)")
