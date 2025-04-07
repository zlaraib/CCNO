@doc """ 
  This function takes in a momentum vector (for all sites) along with their directions 
  and returns its modulus (for all sites) and unit vectors (for all sites and directions)
"""
Base.@pure function momentum(p::Array{Float64,2}, N_sites::Int)
    p_mod = sqrt.(sum(p.^2,dims=2))
    p_hat = p ./ p_mod
    return p_mod, p_hat
end
