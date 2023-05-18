using ITensors

# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.
let
  N = 100 # number of sites
  cutoff = 1E-8 # specifies a truncation threshold for the SVD in MPS representation
  tau = 0.1 # time step 
  ttotal = 5.0 # total time of evolution 

  # Make an array of 'site' indices and label as s 
  # conserve_qns=true conserves the total spin quantum number "S"(in z direction) in the system as it evolves
  s = siteinds("S=1/2", N; conserve_qns=false)  

  # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on neighboring pairs of sites in the chain.
  # Create an empty ITensors array that will be our Trotter gates
  # gates acts as a single-particle operator that does not actually entangle pairs
  gates = ITensor[] 

  for j in 1:(N - 1)
    #s1, s2 are neighbouring spin site/indices from the s array
    s1 = s[j] 
    s2 = s[j + 1]

    # Hamiltonian operator oriented in the x direction of neigbouring spin sites. 
    # "Sx" and "Id" are inbuilt ITensor operators in ITensors library and op function applies these ITensors to a set of neigboring site indices aka pairs s1, s2.
    # hj operator acts on the pair system such that h at site index j and j+1 is computed separately and then summed on this 1D chain.
    # Identity operator on adjacent index is added to form a pair that attaches to the gate operator, hence gates is effectively a one-particle operator.
    hj = pi * (
      op("Sx", s1)* op("Id", s2)  + op("Sx", s2) * op("Id", s1)
    )

    # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors
    Gj = exp(-im * tau / 2 * hj)  

    # The push! function adds (appends) an element to the end of an array;
    # ! performs an operation without creating a new object, (in a way overwites the previous array in consideration); 
    # i.e. we append a new element Gj (which is an ITensor object representing a gate) to the end of the gates array.
    push!(gates, Gj) 
  end

  # append! adds all the elements of a gates in reverse order (i.e. (N,N-1),(N-1,N-2),...) to the end of gates array.
  append!(gates, reverse(gates)) 

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

  # center site 
  c = div(N, 2) # c = N/2 as only half of the spins on the chain are utilized 

  # Create empty array to store sz values 
  Sz_array = Float64[] 

  # Compute and print <Sz> at each time step then apply the gates to go to the next time
  for t in 0.0:tau:ttotal
    # compute initial expectation value of Sz(inbuilt operator in ITensors library) at the center of the chain (site c)
    sz = expect(psi, "Sz"; sites=c)
    # add an element sz to the end of Sz array  
    push!(Sz_array, sz) 
    println("$t $sz")

    # Writing an if statement in a shorthand way that checks whether the current value of t is equal to ttotal, 
    # and if so, it executes the break statement, which causes the loop to terminate early.
    t ≈ ttotal && break  

    # apply each gate in gates successively to the wavefunction psi (it is equivalent to time evolving psi according to the time-dependent Hamiltonian represented by gates).
    # The apply function is smart enough to determine which site indices each gate has, and then figure out where to apply it to our MPS. 
    # It automatically handles truncating the MPS (and can even handle non-nearest-neighbor gates but this feature is not used in this example).
    psi = apply(gates, psi; cutoff)

    # The normalize! function is used to ensure that the MPS is properly normalized after each application of the time evolution gates. 
    # This is necessary to ensure that the MPS represents a valid quantum state.
    normalize!(psi) 
  end

  # Checking that the value of Sz at the center spin oscillates between -0.5 and 0.5 every 5 timesteps
  @assert Sz_array[1] == -0.5
  @assert Sz_array[6] == 0.5
  @assert Sz_array[11] == -0.5
  @assert Sz_array[16] == 0.5
  @assert Sz_array[21] == -0.5
  @assert Sz_array[26] == 0.5
  @assert Sz_array[31] == -0.5
  @assert Sz_array[36] == 0.5
  @assert Sz_array[41] == -0.5
  @assert Sz_array[46] == 0.5
  @assert Sz_array[51] == -0.5
end