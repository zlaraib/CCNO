function geometric_func(N,p)

    for i in 1:(N-1)
        for j in (i+1):N
            if dot(p[i, :], p[j, :]) ≈ 0.0
                return 1
            else
                return 1 - dot(p[i, :], p[j, :])
            end
        end
    end

    return 0  # Add a default return value if no perpendicular pairs are found
end
