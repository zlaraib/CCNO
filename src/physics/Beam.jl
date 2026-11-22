Base.@kwdef struct Beam
    num_sites::Int64
    rotation::Float64
    energy::Float64
    flavor::Int64 # NOTE: 1 <- electron, -1 <- other
end

function init_beam_arr(num_beams::Int64)
    beam_arr = Array{CCNO.Beam, 1}(undef, num_beams) # Init (not zero init, for some reason?)
    # - If num_sites = 0, then the beam is nothing and we skip init and move on.
    for i = 1:num_beams
        beam_arr[i] = Beam(
                    num_sites = 0,
                    rotation = 0.0,
                    energy = 0.0,
                    flavor = 0,
                   )
    end
    return beam_arr
end

