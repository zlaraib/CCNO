using ITensors

N = 10
chi = 4

# Create a random MPS
#The line sites = siteinds("S=1/2", N) creates an array of Index objects that specify the quantum numbers 
#and the range of the sites in the MPS. In this case, we're creating a chain of N spin-1/2 particles, 
#where each site can have a quantum number of either +1/2 or -1/2, represented by the string "S=1/2". 
#The siteinds function returns an array of Index objects, one for each site in the chain.
sites = siteinds("S=1/2", N)
psi = randomMPS(sites, chi)


# # Calculate the magnitude of avg Spin in the z direction
magz = expect(psi, "Sz")
for (j, mz) in enumerate(magz)
    println("$j $mz")
end
