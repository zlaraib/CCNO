#using ITensors



# # d = 2
# # N = 5
# # A = randn(d,d,d,d,d)

# #define indices 
# l= Index(2, "l")
# m= Index(2, "m")
# r= Index(2, "r")
# a= Index(2, "a")
# b= Index(2, "b")

# #define ITensors
# A = randomITensor(l,m,r)
# L = randomITensor(l,l',a)
# M = randomITensor(m,m',a,b)
# R = randomITensor(r,r',b)

# #contract
# B = A*L*M*R


# U,S,V = svd(A, (l,m))
# println("U", U)
# println("S", S)
# println("S",S)
# #norm(A - U*S*V) 

# i= Index(2, "i")
# j= Index(2, "j")
# k= Index(2, "k")
# y= Index(2, "y")
# z= Index(2, "z")
# cutoff = 1E-8
# maxdim =10
# T = randomITensor(i,j,k,y,z)
# #T = randomITensor(l,m,r,a,b)
# M = MPS(T,(i,j,k,y,z);cutoff=cutoff,maxdim=maxdim)
# #M = MPS(T,(l,m,r,a,b);cutoff=cutoff,maxdim=maxdim)
using ITensors


N = 10
s = siteinds(2,N)
chi = 4
psi = randomMPS(s;linkdims=chi)

# Make an array of integers of the element we
# want to obtain
el = [1,2,1,1,2,1,2,2,2,1]
println(state)
V = ITensor(1.)
for j=1:N
    global V *= (psi[j]*state(s[j],el[j]))
end
v = scalar(V)

# v is the element we wanted to obtain:
@show v