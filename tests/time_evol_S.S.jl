using ITensors
 
# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.

  N = 100 # number of sites
  cutoff = 1E-8 # specifies a truncation threshold for the SVD in MPS representation
  tau = 0.1 # time step 
  ttotal = 5.0 # total time of evolution 

  # Compute the number of steps to do
  Nsteps = Int(ttotal/tau)

  
  # Make an array of 'site' indices
  s = siteinds("S=1/2",N;conserve_qns=true) # conserve_qns=true specifies the quantum numbers (total spin in sys = quantum number "S") that are conserved by the system (as sys evolves)

  # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on neighboring pairs of sites in the chain
  gates = ITensor[] # empty array that will hold ITensors that will be our Trotter gates
  #gates is an array of ITensors, each representing a time evolution operator for a pair of neighboring sites in this quantum system.
  for j=1:N-1
    s1 = s[j]
    s2 = s[j+1]
    #hj = 0.5 * (op("S+",s1)*op("S-",s2) + op("S-",s1)*op("S+",s2))
    #hj = 0.25 *((op("S-",s1)*op("S-",s2)) + (op("S-",s1)*op("S+",s2)) + (op("S+",s1)*op("S-",s2)) + (op("S+",s1)*op("S+",s2)))
    #hj =       op("Sx",s1) * op("Sx",s2)
     hj =       op("Sz",s1) * op("Sz",s2) +
         1/2 * op("S+",s1) * op("S-",s2) +
         1/2 * op("S-",s1) * op("S+",s2)
    Gj = exp(-1.0im * tau/2 * hj) #making Trotter gate Gj that would correspond to each gate in the gate array of ITensors
    push!(gates,Gj) # The push! function adds an element to the end of an array, ! = performs an operation without creating a new object, (in a wy overwite the previous array in consideration) i.e. it appends a new element Gj (which is an ITensor object representing a gate) to the end of the gates array.
  # Include gates in reverse order too
  # (N,N-1),(N-1,N-2),...
  append!(gates,reverse(gates)) # append! adds all the elements of a collection to the end of an array.
end

  # Function that measures <Sz> on site n
  function measure_Sz(psi,n)
    # Without calling `orthogonalize`, this would not be
    # the local density matrix! We would have to trace out
    # the sites that are not `n`.
    psi = orthogonalize(psi,n) # orthogonalize function ensures that MPS is in its canonical form (i.e. where each tensor is orthogonal to its neighboring tensors upon contraction in a specific order) and thus has minimum possible number of parameters to work with. b/c when we apply the gates to the MPS, the MPS could deviate from the orthogonal form.
    sn = siteind(psi,n)
    # scalar = inner product , dag = dagger, prime creates a new index object with the same quantum number structure as the index of psi[n], but with a new label "Site",  op function returns Sz operator on siteindex of nth site  
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n]) 
    return real(Sz)
  end

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

  c = div(N,2) # c = N/2 because only half of the spins on the chain are utilized 

  # Compute and print initial <Sz> value
  t = 0.0
  Sz = measure_Sz(psi,c) # computing initial expectation value of Sz at the center of the chain (site c)
  println("t_initial=$t, Sz_initial= $Sz")
  println("Begin Evolution:\n")
  # Do the time evolution by applying the gates
  # for Nsteps steps
  for step=1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff) # apply each gate in gates successively to the wavefunction psi =  evolving the wavefunction in time according to the time-dependent Hamiltonian represented by the gates.
    # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. It automatically handles truncating the MPS and can even handle non-nearest-neighbor gates, though that feature is not used in this example.
    t += tau
    Sz = measure_Sz(psi,c) # After each gate is applied, we compute the expectation value of the Sz operator at the central site c again.
    println("t_new=$t, Sz_new=$Sz") #Sz measured at each time step.
  end
  # Thus Sz expectation value is being measured, at a specific site c, at each time step to monitor the time evolution of the system.
  return
end