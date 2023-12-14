""" 
  This function takes in a momentum vector (for all sites) along with their directions 
  and returns its modulus (for all sites) and unit vectors (for all sites and directions)
"""
function momentum(p, N_sites)
p_hat = []  # Initialize an empty array to collect p_i_hat vectors. This array will contain unit vector(hat) of all sites in all directions
p_mod = [] #Initialize array that contains mod of all sites 
  for i in 1:N_sites

    if p ==  ones(N_sites,3)
      p[i, :] /= norm(p[i, :]) # Normalize each site of p to have magnitude 1. # special addition for this omega = pi case
      # check if this needs to be obeyed for FFI test
    end 

    # mod = norm in julia. # The norm returns the magnitude of the vector (by default)
    p_i_mod = norm(p[i, :]) # Norm of each site 
    #calculate unit vector p for each site 
    p_i_hat = p[i, :] / p_i_mod

    # Append p_i_hat to the p_hat array
    push!(p_hat, p_i_hat)
    # Append p_i_mod to the p_mod array
    push!(p_mod, p_i_mod)
  end
  
  return p_mod, p_hat

end