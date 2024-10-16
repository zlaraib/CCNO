include("heaviside.jl")

""" 
    This file introduces the different shapes for neutrinos and map those shapes to
    a dictionary to allow easy switching between shapes in test files. 
"""

function flat_top(Δp, ξ)
    return 1 / Δp * heaviside(1/2 - abs(ξ))
end

function triangular(Δp, ξ)
    return 1 / Δp * (1 - abs(ξ)) * heaviside(1 - abs(ξ))
end

function none(Δp, ξ)
    return 1 
end

function shape_func(x, Δp, i, j,L, shape_name, periodic)

    # Define a dictionary mapping shape names to functions
    shapes = Dict(
        "flat_top" => flat_top,
        "triangular" => triangular,
        "none" => none
    )
    
    if haskey(shapes, shape_name)
        shape_function = shapes[shape_name]  # assign the corresponding function (selected in shape_name in the test file) to shape_function.
        # ξ determines the overlap of neutrinos (i.e. their shapes) on top of each other to signify the neutrino-neutrino interaction strength.
        ξ = (x[i] - x[j]) / Δp # dependent upon difference in distance of i and j sites and the width of the shape function. 
        if periodic 
            if x[i]-x[j] > L/2
                ξ = ξ - L/Δp 
            end
            if x[i] - x[j] < -L/2
                ξ = ξ + L/Δp 
            end

            if shape_name == "none"
                @assert (Δp == L)
            else @assert (2*Δp <L)
            end
            @assert (abs(ξ)<= L/(2* Δp))
        end

        return shape_function(Δp, ξ)
    else
        error("Unknown shape: $shape_name")
    end

end