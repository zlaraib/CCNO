""" 
    This function takes in a scalar x  
    and returns its heaviside step-function result
"""
function heaviside(x)
    if x < 0
        return 0
    else
        return 1
    end
end