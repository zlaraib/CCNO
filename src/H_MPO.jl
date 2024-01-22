
include("constants.jl")
<<<<<<< HEAD

function Hamiltonian_mpo(s,N, ω, B, N_sites, Δx,energy_sign)
    #input the operator terms
    os = OpSum()                                                            

    for i in 1:(N_sites-1)
        for j in i+1:N_sites
=======
include("geometric_func.jl")
include("shape_func.jl")
include("momentum.jl")


function Hamiltonian_mpo(s, N, B, N_sites, Δx, Δm², p, x, Δp, shape_name,L, τ, energy_sign, periodic)
    #input the operator terms
    os = OpSum()                                                            
    
    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p̂ = momentum(p,N_sites)  

    # define an array of vacuum oscillation frequencies (units of ergs)
    if Δm² == 0 #specific to self-int only
        ω = zeros(N_sites)
    elseif Δm² == 2 * π
        global ω = fill(π, N_sites) # added global so we can access and use this global variable without the need to pass them as arguments to another function
    else
        ω = [Δm²/ (2 * p_i_mod) for p_i_mod in p_mod]
        # @assert Im(ω) == ((Δm²))/(2 *hbar * Eνₑ)
        # println("Im(ω)=",Im(ω))
        # println("p_i_mod= ",p_i_mod)
        # println("ω ", ω)
    end
    println("ω = ", ω)

    for i in 1:(N-1)
        for j in i+1:N
>>>>>>> c9f2089 (added src file for MPO Hamiltonian)
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B) == 1

            # Our neutrino system Hamiltonian of self-interaction term represents 1D Heisenberg model.
            # ni and nj are the neutrions at site i and j respectively.
            # mu pairs divided by 2 to avoid double counting
<<<<<<< HEAD
             
            if energy_sign[i]*energy_sign[j]>0
                interaction_strength = (2.0/N_sites * √2 * G_F * (N[i]+ N[j])/(2*((Δx)^3))) 
                os+= interaction_strength,"Sz",i,"Sz", j 
                os+= 1/2,interaction_strength,"S+",i,"S-",j
                os+=  1/2,interaction_strength,"S-",i,"S+",j
            else
                interaction_strength = (2.0/N_sites * √2 * G_F * (N[i]+ N[j])/(2*((Δx)^3))) 
=======
            # interaction_strength = (2.0* √2 * G_F * (N[i]+ N[j])/(2*((Δx)^3)) * 1/N_sites) 
            # numerical_factor=(1/(N_sites-1))
            # os+= interaction_strength,"Sz",i,"Sz", j 
            # os+= 1/2,interaction_strength,"S+",i,"S-",j
            # os+=  1/2,interaction_strength,"S-",i,"S+",j
             
            if energy_sign[i]*energy_sign[j]>0
                # Get the shape function result for each pair of i and j 
                shape_result = shape_func(x, Δp, i, j,L, shape_name, periodic)
                # Calculate the geometric factor for each pair of i and j within the loop
                geometric_factor = geometric_func(p, p̂, i, j)
                interaction_strength = (2.0* √2 * G_F * (N[i]+ N[j])/(2*((Δx)^3))) * shape_result * geometric_factor
                os+= interaction_strength,"Sz",i,"Sz", j 
                os+= 1/2,interaction_strength,"S+",i,"S-",j
                os+=  1/2,interaction_strength,"S-",i,"S+",j
                #println("hj= ", hj)
            else
                # Get the shape function result for each pair of i and j 
                shape_result = shape_func(x, Δp, i, j,L, shape_name, periodic)
                # Calculate the geometric factor for each pair of i and j within the loop
                geometric_factor = geometric_func(p, p̂, i, j)
                interaction_strength = (2.0* √2 * G_F * (N[i]+ N[j])/(2*((Δx)^3))) * shape_result * geometric_factor
>>>>>>> c9f2089 (added src file for MPO Hamiltonian)
                os+= 2,interaction_strength ,"Sz",i,"Sz", j 
                os+= -interaction_strength ,"S+",i,"S-",j
                os+=  -interaction_strength ,"S-",i,"S+",j
            end

            if ω[i] != 0 && ω[j] != 0
<<<<<<< HEAD
                numerical_factor = (1/(N_sites-1))
                os+= numerical_factor,energy_sign[i],ω[i],B[1],"Sx",i,"Id",j  
                os+= numerical_factor,energy_sign[i],ω[i],B[2],"Sy",i,"Id",j 
                os+= numerical_factor,energy_sign[i],ω[i],B[3],"Sz",i,"Id",j 
                os+= numerical_factor,energy_sign[j],ω[j],B[1],"Id",i,"Sx",j
                os+= numerical_factor,energy_sign[j],ω[j],B[2],"Id",i,"Sy",j
                os+= numerical_factor,energy_sign[j],ω[j],B[3],"Id",i,"Sz",j
=======
                numerical_factor = 1/(N_sites-1)
                os+= numerical_factor,ω[i],B[1],"Sx",i,"Id",j  
                os+= numerical_factor,ω[i],B[2],"Sy",i,"Id",j 
                os+= numerical_factor,ω[i],B[3],"Sz",i,"Id",j 
                os+= numerical_factor,ω[j],B[1],"Id",i,"Sx",j
                os+= numerical_factor,ω[j],B[2],"Id",i,"Sy",j
                os+= numerical_factor,ω[j],B[3],"Id",i,"Sz",j
                  
>>>>>>> c9f2089 (added src file for MPO Hamiltonian)
            end
            
        end
    end
<<<<<<< HEAD
    # return os
    H = MPO(os,s)
    return H
end
=======


    #Convert these terms to an MPO
    # H = MPO(os,s; splitblocks=true)
    H = MPO(os,s)
    # if outputlevel > 0
    #     @show use_splitblocks
    # end
    # #    # This step makes the MPO more sparse but also
    # #    # introduces more blocks.
    # #   #  It generally improves DMRG performance
    # #    # at large bond dimensions.
    # if use_splitblocks
    #     H = splitblocks(linkinds, H)
    # end
    # # Number of structural nonzero elements in a bulk
    # # Hamiltonian MPO tensor
    # if outputlevel > 0 && ω == fill(0, N_sites) 
    #     @show nnz(H[end÷2])
    #     @show nnzblocks(H[end÷2])
    # end
    return H
end

# function eigenvals(s, τ, n, ω, B, N, Δx, cutoff, tolerance, ttotal,outputlevel,
#     maxdim,use_splitblocks, sweeps,nsweeps,use_threaded_blocksparse)
    
#     H = Hamiltonian_mpo(s, n, ω, B, N, Δx, τ,outputlevel, use_splitblocks)
#     state = Vector{String}(undef, N)
#     ψ = nothing
#     for i in 1:N
#         # Assign "Up" or "Dn" based on even/odd index
#         if i % 2 == 0
#             state[i] = "Up"
#         else
#             state[i] = "Dn"
#         end
#     end
    
#     #create a MPS state with half spins pointing up and half pointing down
#     ψ0 = randomMPS(s, state,4)
#     # ψ0 = productMPS(s, n -> isodd(n) ? "Dn" : "Up")
#     #ψ0=randomMPS(s,linkdims=4)
#     apply_time_parallel = 0

#     #obs = EntanglementObserver()

#     Sz_observer = DMRGObserver(["Sz"],s,energy_tol=1E-7)
#     observer = EntanglementObserver()
    
#     ITensors.enable_combine_contract()
#     #apply_time_parallel += @elapsed energy, ψ =@time dmrg(H,ψ0; nsweeps, maxdim, cutoff, observer=Sz_observer, outputlevel=2)
#     apply_time_parallel += @elapsed energy, ψ = @time dmrg(H, ψ0, sweeps; outputlevel)
#     normalize!(ψ)
#     ITensors.disable_combine_contract()
    
#     for (sw,Szs) in enumerate(measurements(Sz_observer)["Sz"])
#         println("Total Sz after sweep $sw = ", sum(Szs)/N)
#     end

#     return energy,ψ,apply_time_parallel
# end
>>>>>>> c9f2089 (added src file for MPO Hamiltonian)
