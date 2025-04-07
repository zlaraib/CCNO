@doc """ 
    This function takes in a scalar x  
    and returns its heaviside step-function result
"""
@inline function heaviside(x::Float64)
    if x < 0
        return 0
    else
        return 1
    end
end
