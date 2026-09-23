function uniform_sphere(nphi_at_equator)
    dtheta = pi * sqrt(3) / nphi_at_equator

    xyz = zeros(0)
    theta = 0
    phi0 = 0
    while theta < pi/2 
        if theta == 0
            nphi = nphi_at_equator
        else
            nphi = int(round(nphi_at_equator * cos(theta)))
        end
            
        dphi = 2*pi/nphi
        if nphi==1 
          theta = pi/2
        end
        
        for iphi in 1:nphi
            phi = phi0 = iphi*dphi
            x = cos(theta) * cos(phi)
            y = cos(theta) * sin(phi)
            z = sin(theta)
            if theta > 0
                append!(xyz, hcat(-x,-y,-z))
            else
                append!(xyz, hcat(x,y,z))
            end
        theta += dtheta
        phi0 += 0.5 * dphi # offset by half step so adjacent latitudes are not always aligned in longitude
      end
    end
        
    return xyz
end
