@doc """ 
    This file introduces the different shapes for neutrinos and map those shapes to
    a dictionary to allow easy switching between shapes in test files. 
"""

@inline function flat_top(Δp::Float64, ξ::Float64)
    return heaviside(1/2 - abs(ξ))
end

@inline function triangular(Δp::Float64, ξ::Float64)
    return (1 - abs(ξ)) * heaviside(1 - abs(ξ))
end

@inline function none(Δp::Float64, ξ::Float64)
    return 1 
end

Base.@pure function shape_func(params::CCNO.Parameters, x::Vector{Float64}, i::Int, j::Int,L::Float64)

    # Define a dictionary mapping shape names to functions
    shapes = Dict(
        "flat_top" => flat_top,
        "triangular" => triangular,
        "none" => none
    )
    
    if haskey(shapes, params.shape_name)
        shape_function = shapes[params.shape_name]  # assign the corresponding function (selected in shape_name in the test file) to shape_function.
        # ξ determines the overlap of neutrinos (i.e. their shapes) on top of each other to signify the neutrino-neutrino interaction strength.
        ξ = (x[i] - x[j]) / params.Δp # dependent upon difference in distance of i and j sites and the width of the shape function. 
        if params.periodic 
            if x[i]-x[j] > L/2
                ξ = ξ - L/params.Δp 
            end
            if x[i] - x[j] < -L/2
                ξ = ξ + L/params.Δp 
            end

            if params.shape_name == "none"
                @assert (params.Δp == L)
            else @assert (2*params.Δp <L)
            end
            @assert (abs(ξ)<= L/(2* params.Δp))
        end

        return shape_function(params.Δp, ξ)
    else
        error("Unknown shape: $params.shape_name")
    end

end
