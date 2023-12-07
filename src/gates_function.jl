using ITensors
include("constants.jl")

# Expected units of the quantities defined in the files in tests directory that are being used in the gates function.                                                                   
# s = site index array (dimensionless and unitless)          
# n = no.of neutrinos (dimensionless and unitless)
# ω = vacuum oscillation angular frequency (rad/s)
# B = Normalized vector related to mixing angle in vacuum oscillations (dimensionless constant)
# N = Total no.of sites (dimensionless and unitless)
# Δx = length of the box of interacting neutrinos at a site (cm) 
# τ = time step (sec)

# This test runs the create_gates function that holds ITensors Trotter gates and returns the dimensionless unitary 
# operators govered by the Hamiltonian which includes effects of the vacuum and self-interaction potential for each site.

function create_gates(s, n, ω, B, N, Δx, τ,outputlevel, use_splitblocks)
    # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on any (non-neighboring) pairs of sites in the chain.
    # Create an empty ITensors array that will be our Trotter gates
    gates = ITensor[]                                                              
    
    for i in 1:(N-1)
        for j in i+1:N
            #s_i, s_j are non-neighbouring spin site/indices from the s array
            s_i = s[i]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
            s_j = s[j]       
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B) == 1

            # Our neutrino system Hamiltonian of self-interaction term represents 1D Heisenberg model.
            # total Hamiltonian of the system is a sum of local terms hj, where hj acts on sites i and j which are paired for gates to latch onto.
            # op function returns these operators as ITensors and we tensor product and add them together to compute the operator hj.
            # ni and nj are the neutrions at site i and j respectively.
            # mu pairs divided by 2 to avoid double counting
            
            hj = 
            ((2.0* √2 * G_F * (n[i]+ n[j])/(2*((Δx)^3)) * 1/N) * 
            (op("Sz", s_i) * op("Sz", s_j) +
             1/2 * op("S+", s_i) * op("S-", s_j) +
             1/2 * op("S-", s_i) * op("S+", s_j)))
             
             if ω[i] != 0 && ω[j] != 0
                hj += (1/(N-1))* 
                ((ω[i] * B[1] * op("Sx", s_i)* op("Id", s_j))  + (ω[j] * op("Sx", s_j) * op("Id", s_i))) + 
                ((ω[i] * B[2] * op("Sy", s_i)* op("Id", s_j))  + (ω[j] * op("Sy", s_j) * op("Id", s_i))) +
                ((ω[i] * B[3] * op("Sz", s_i)* op("Id", s_j))  + (ω[j] * op("Sz", s_j) * op("Id", s_i))) 
                     
             end
            
            # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors             
            Gj = exp(-im * τ/2 * hj)

            # The push! function adds (appends) an element to the end of an array;
            # ! performs an operation without creating a new object, (in a way overwites the previous array in consideration); 
            # i.e. we append a new element Gj (which is an ITensor object representing a gate) to the end of the gates array.
            push!(gates, Gj)
        end
    end

    # append! adds all the elements of a gates in reverse order (i.e. (N,N-1),(N-1,N-2),...) to the end of gates array.
    # appending reverse gates to create a second-order Trotter-Suzuki integration
    append!(gates, reverse(gates))
    if outputlevel > 0
        @show use_splitblocks
    end

    # Number of structural nonzero elements in a bulk
    # Hamiltonian MPO tensor
    if outputlevel > 0
        @show nnz(gates[end÷2])
        @show nnzblocks(gates[end÷2])
    end
    return gates
end




function Hamiltonian_mpo(s, n, ω, B, N, Δx, τ,outputlevel, use_splitblocks)
    #input the operator terms
    os = OpSum()                                                            
    
    for i in 1:(N-1)
        for j in i+1:N
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B) == 1

            # Our neutrino system Hamiltonian of self-interaction term represents 1D Heisenberg model.
            # ni and nj are the neutrions at site i and j respectively.
            # mu pairs divided by 2 to avoid double counting
            interaction_strength = (2.0* √2 * G_F * (n[i]+ n[j])/(2*((Δx)^3)) * 1/N) 
            numerical_factor=(1/(N-1))
            os+= interaction_strength,"Sz",i,"Sz", j 
            os+= 1/2,interaction_strength,"S+",i,"S-",j
            os+=  1/2,interaction_strength,"S-",i,"S+",j
             
            if ω[i] != 0 && ω[j] != 0

                os+= numerical_factor,ω[i],B[1],"Sx",i,"Id",j  
                os+= numerical_factor,ω[j],B[1],"Sx",j,"Id",i 
                os+= numerical_factor,ω[i],B[2],"Sy",i,"Id",j  
                os+= numerical_factor,ω[j],B[2],"Sy",j,"Id",i
                os+= numerical_factor,ω[i],B[3],"Sz",i,"Id",j 
                os+= numerical_factor,ω[j],B[3],"Sz",j,"Id",i
                  
            end
            
        end
    end


    #Convert these terms to an MPO
    H = MPO(os,s; splitblocks=true)

    if outputlevel > 0
        @show use_splitblocks
    end
    #    # This step makes the MPO more sparse but also
    #    # introduces more blocks.
    #   #  It generally improves DMRG performance
    #    # at large bond dimensions.
    if use_splitblocks
        H = splitblocks(linkinds, H)
    end
    # Number of structural nonzero elements in a bulk
    # Hamiltonian MPO tensor
    if outputlevel > 0
        @show nnz(H[end÷2])
        @show nnzblocks(H[end÷2])
    end
    return H
end

