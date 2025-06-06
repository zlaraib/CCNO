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

# Define a dictionary mapping shape names to functions
const shapes = Dict(
    "flat_top" => flat_top,
    "triangular" => triangular,
    "none" => none
)
    
@inline function shape_func(params::CCNO.Parameters, state::CCNO.SimulationState, d::Int64, i::Int64, j::Int64)

    shape_function = shapes[params.shape_name]  # assign the corresponding function (selected in shape_name in the test file) to shape_function.
    # ξ determines the overlap of neutrinos (i.e. their shapes) on top of each other to signify the neutrino-neutrino interaction strength.
    dx_ij::Float64 = state.xyz[i,d] - state.xyz[j,d]
    ξ::Float64 = dx_ij / params.Δp # dependent upon difference in distance of i and j sites and the width of the shape function. 
    if params.periodic 
        if dx_ij > params.L/2
            ξ = ξ - params.L/params.Δp 
        end
        if dx_ij < -params.L/2
            ξ = ξ + params.L/params.Δp 
        end
        
        if params.shape_name == "none"
            @assert (params.Δp == params.L)
        else @assert (2*params.Δp <params.L)
        end
        @assert (abs(ξ)<= params.L/(2* params.Δp))
    end
    
    return shapes[params.shape_name](params.Δp, ξ)
end
