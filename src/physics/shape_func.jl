@inline function flat_top(Delta_p::Float64, xi::Float64)
    return heaviside(1/2 - abs(xi))
end

@inline function triangular(Delta_p::Float64, xi::Float64)
    return (1 - abs(xi)) * heaviside(1 - abs(xi))
end

@inline function none(Delta_p::Float64, xi::Float64)
    return 1
end

# Define a dictionary mapping shape names to functions
const shapes =
    Dict("flat_top" => flat_top, "triangular" => triangular, "none" => none)

@doc """ 
    Return normalized interaction strength used to quantify inhomogeneity.

    params: CCNO.Parameters object defined in physics/Parameters.jl
    state: CCNO.SimulationState object defined in physics/SimulationState.jl
    d: [1,2,3] direction along which the shape function is to be evaluated (x/y/z)
    i,j: integers of the two interacting sites
"""
@inline function shape_func(
    params::CCNO.Parameters,
    state::CCNO.SimulationState,
    i::Int64,
    j::Int64,
)

    # Range of considered dimensions (x, y, z)
    shape_result::Float64 = 1.0
    
    for d = 1:3
        # Ignore dim (same as shape_func=none) if homogeneous
        if params.homogeneous[d]
            continue
        end
        shape_function = shapes[params.shape_name]  # assign the corresponding function (selected in shape_name in the test file) to shape_function.
        # xi determines the overlap of neutrinos (i.e. their shapes) on top of each other to signify the neutrino-neutrino interaction strength.
        dx_ij::Float64 = state.xyz[i, d] - state.xyz[j, d]
        xi::Float64 = dx_ij / params.Delta_p # dependent upon difference in distance of i and j sites and the width of the shape function. 
        if params.periodic
            if dx_ij > params.L/2
                xi = xi - params.L/params.Delta_p
            end
            if dx_ij < -params.L/2
                xi = xi + params.L/params.Delta_p
            end

            if params.shape_name == "none"
                @assert (params.Delta_p == params.L)
            else
                @assert (2*params.Delta_p < params.L)
            end
            @assert (abs(xi)<=params.L/(2 * params.Delta_p))
        end
        shape_result *= shapes[params.shape_name](params.Delta_p, xi)
    end

    return shape_result
end
