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