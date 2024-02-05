
include("constants.jl")

function Hamiltonian_mpo(s,n, ω, B, N_sites, Δx,energy_sign)
    #input the operator terms
    os = OpSum()                                                            

    for i in 1:(N_sites-1)
        for j in i+1:N_sites
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B) == 1

            # Our neutrino system Hamiltonian of self-interaction term represents 1D Heisenberg model.
            # ni and nj are the neutrions at site i and j respectively.
            # mu pairs divided by 2 to avoid double counting
             
            if energy_sign[i]*energy_sign[j]>0
                interaction_strength = (2.0/N_sites * √2 * G_F * (n[i]+ n[j])/(2*((Δx)^3))) 
                os+= interaction_strength,"Sz",i,"Sz", j 
                os+= 1/2,interaction_strength,"S+",i,"S-",j
                os+=  1/2,interaction_strength,"S-",i,"S+",j
            else
                interaction_strength = (2.0/N_sites * √2 * G_F * (n[i]+ n[j])/(2*((Δx)^3))) 
                os+= 2,interaction_strength ,"Sz",i,"Sz", j 
                os+= -interaction_strength ,"S+",i,"S-",j
                os+=  -interaction_strength ,"S-",i,"S+",j
            end

            if ω[i] != 0 && ω[j] != 0
                if energy_sign[i]*energy_sign[j]>0
                    numerical_factor = (1/(N_sites-1))
                    os+= numerical_factor,ω[i],B[1],"Sx",i,"Id",j  
                    os+= numerical_factor,ω[i],B[2],"Sy",i,"Id",j 
                    os+= numerical_factor,ω[i],B[3],"Sz",i,"Id",j 
                    os+= numerical_factor,ω[j],B[1],"Id",i,"Sx",j
                    os+= numerical_factor,ω[j],B[2],"Id",i,"Sy",j
                    os+= numerical_factor,ω[j],B[3],"Id",i,"Sz",j
                else
                    numerical_factor = (1/(N_sites-1))
                    os+= numerical_factor,ω[i],B[1],"Sx",i,"Id",j  
                    os+= numerical_factor,ω[i],B[2],"Sy",i,"Id",j 
                    os+= numerical_factor,ω[i],B[3],"Sz",i,"Id",j 
                    os+= numerical_factor,-ω[j],B[1],"Id",i,"Sx",j
                    os+= numerical_factor,-ω[j],B[2],"Id",i,"Sy",j
                    os+= numerical_factor,-ω[j],B[3],"Id",i,"Sz",j
                end
                  
            end
            
        end
    end
    # return os
    H = MPO(os,s)
    return H
end
