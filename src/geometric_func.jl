
""" 
    This function takes in a momentum vector (for all sites) along with their directions 
    and returns the geometric factor dependent on the unit momentum vectors for all
    i< j sites for a (non) zero p vector. 
"""
function geometric_func(p, p̂, i, j)

    if p[i] == 0  && p[j] == 0.0 
        return 1
    else
        return 1 - dot(p̂[i, :], p̂[j, :])
    end

end
