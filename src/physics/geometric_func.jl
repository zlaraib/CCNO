""" 
    This function takes in a momentum vector (for all sites) along with their directions 
    and returns the geometric factor dependent on the unit momentum vectors for all
    i< j sites for a (non) zero p vector.

    params: CCNO.Parameters object defined in physics/Parameters.jl
    phat: [Nsites,3] array of momentum unit vectors, one unit vector for each site
    i,j: indices of the two interacting sites

    Returns: scalar geometric factor for the two sites
"""
@inline function geometric_func(params::CCNO.Parameters, phat::Array{Float64,2}, i::Int64, j::Int64)
    if params.geometric_name == "none"
        return 1
    elseif params.geometric_name == "physical"
        @views result = 1 - dot(phat[i, :], phat[j, :])
        return result
    else
        println("params.geometric_func not implemented")
        exit()
    end

end
