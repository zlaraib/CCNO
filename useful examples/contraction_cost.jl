#This example helps us learn that in the limit of large MPS bond dimension m,
# the first contraction sequence is faster, while in the limit of large MPO 
# bond dimension k, the second sequence is faster. This has practical
# implications for writing an efficient DMRG algorithm in both limits, 
# which we plan to incorporate into ITensors.jl.

using ITensors
using Symbolics

using ITensors: contraction_cost

@variables m, k, d

l = Index(m, "l")
r = Index(m, "r")
h₁ = Index(k, "h₁")
h₂ = Index(k, "h₂")
h₃ = Index(k, "h₃")
s₁ = Index(d, "s₁")
s₂ = Index(d, "s₂")

H₁ = emptyITensor(dag(s₁), s₁', dag(h₁), h₂)
H₂ = emptyITensor(dag(s₂), s₂', dag(h₂), h₃)
L = emptyITensor(dag(l), l', h₁)
R = emptyITensor(dag(r), r', h₃)
ψ = emptyITensor(l, s₁, s₂, r)

TN = [ψ, L, H₁, H₂, R]
sequence1 = Any[2, Any[3, Any[4, Any[1, 5]]]]
sequence2 = Any[Any[4, 5], Any[1, Any[2, 3]]]
cost1 = contraction_cost(TN; sequence = sequence1)
cost2 = contraction_cost(TN; sequence = sequence2)

println("First sequence")
display(sequence1)
display(cost1)
@show sum(cost1)
@show substitute(sum(cost1), Dict(d => 4))

println("\nSecond sequence")
display(sequence2)
display(cost2)
@show sum(cost2)
@show substitute(sum(cost2), Dict(d => 4))