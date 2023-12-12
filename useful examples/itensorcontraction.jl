N = 10

s = siteinds("S=1/2", N; conserve_qns = true)

state = n -> isodd(n) ? "↑" : "↓"
ψ₁ = randomMPS(s, state, 2)
ψ₂ = randomMPS(s, state, 2)
ψ₃ = randomMPS(s, state, 2)

ψ = +(ψ₁, ψ₂; cutoff = 1e-8)

# Can use:
#
# ψ = ψ₁ + ψ₂
#
# but generally you want to set a custom `cutoff` and `maxdim`.

println()
@show inner(ψ, ψ)
@show inner(ψ₁, ψ₂) + inner(ψ₁, ψ₂) + inner(ψ₂, ψ₁) + inner(ψ₂, ψ₂)

# Computes ψ₁ + 2ψ₂
ψ = ψ₁ + 2ψ₂

println()
@show inner(ψ, ψ)
@show inner(ψ₁, ψ₁) + 2 * inner(ψ₁, ψ₂) + 2 * inner(ψ₂, ψ₁) + 4 * inner(ψ₂, ψ₂)

# Computes ψ₁ + 2ψ₂ + ψ₃
ψ = ψ₁ + 2ψ₂ + ψ₃

println()
@show inner(ψ, ψ)
@show inner(ψ₁, ψ₁) + 2 * inner(ψ₁, ψ₂) + inner(ψ₁, ψ₃) +
      2 * inner(ψ₂, ψ₁) + 4 * inner(ψ₂, ψ₂) + 2 * inner(ψ₂, ψ₃) +
      inner(ψ₃, ψ₁) + 2 * inner(ψ₃, ψ₂) + inner(ψ₃, ψ₃)