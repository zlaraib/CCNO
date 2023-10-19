include("heaviside.jl")
# Δp  = rand()  # width of the shape function # move it to the tests script outside

function flat_top(Δp, ξ)
    return 1 / Δp * heaviside(1/2 - abs(ξ))
end

function triangular(Δp, ξ)
    return 1 / Δp * (1 - abs(ξ)) * heaviside(1 - abs(ξ))
end

function shape_func(x, Δp, N)
    if Δp == Inf
        return 1  # using Δp infinite as Rog case for now 
    else
        # Define a dictionary mapping shape names to functions
        shapes = Dict(
            "flat_top" => flat_top,
            "triangular" => triangular
        )
        #Select a shape function based on the shape_name variable form the list defined in dictionary 
        shape_name = "flat_top"  # Change this to the desired shape name
        if haskey(shapes, shape_name)
            shape_function = shapes[shape_name] # assign the corresponding function (selected in shape_name above) to shape_function.
            result = zeros(N)  # Adjust the dimensions as needed rn just initialized to N 
            for i in 1:(N - 1)
                for j in (i+1):N
                    ξ = (x[i] - x[j]) / Δp
                    result = shape_function(Δp, ξ)
                end
            end
            return result
        else
            error("Unknown shape: $shape_name")
        end
    end
end