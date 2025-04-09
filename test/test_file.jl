using ITensors
using ITensorMPS

N = 4
chi = 4

#The line sites = siteinds("S=1/2", N) creates an array of Index objects that specify the quantum numbers 
#and the range of the sites in the MPS. In this case, we're creating a chain of N spin-1/2 particles, 
#where each site can have a quantum number of either +1/2 or -1/2, represented by the string "S=1/2". 
#The siteinds function returns an array of Index objects, one for each site in the chain.
sites = siteinds("S=1/2", N)

# Create a initial state with first index down and second index up
psi = productMPS(sites, n -> isodd(n) ? "Up" : "Dn")


# # Calculate the magnitude of avg Spin in the z direction
magz = expect(psi, "Sz")
for (j, mz) in enumerate(magz)
    println("$j $mz")
end

@assert magz[1] ==  0.5
@assert magz[2] == -0.5
@assert magz[3] ==  0.5
@assert magz[4] == -0.5
