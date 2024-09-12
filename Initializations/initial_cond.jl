
# BTW When you apply gates to an MPS, it will in general increase the bond dimension 
# (one exception is that single-site gates don’t change the bond dimension), 
# and if you continue applying gates without truncation the bond dimension will in general grow exponentially.
# There isn’t a way to set the truncation when you initialize the state
# this is all from julia itensors discourse answered by matthew fisherman


# We are simulating the time evolution of a 1D spin chain with N sites, where each site is a spin-1/2 particle. 
# The simulation is done by applying a sequence of unitary gates to an initial state of the system, 
# which is a product state where each site alternates between up and down.
# throughout this code the words particle and sites are used interchangeably.
# Where each site is occupied by either some neutrinos or some antineutrinos. 

# throughout this code I am assuming each site is occupied by a particle i.e. each site contains some number of neutrinos all of same flavor 
# so all neutrinos are electron flavored (at a site) which interact with electron flavored anti neutrinos (at a different site) in the opposing beam.


"""
PIC SETUP: 

big box (for 'all' sites/particles)
total vol = V 
length of each side = L 
V = L³ 
N_sites = total no.of sites or particles
N =  total no.of neutrinos at all sites
n = total no.density of neutrinos at all sites
n = N / V 

small box (for 'each' site/particle)
total vol = Vᵢ
length of each side = Δx 
Δx³ = Vᵢ
Nᵢ =  total no.of neutrinos at site i 
nᵢ  = total no.density of neutrinos at site i 
Nᵢ = nᵢ  * Vᵢ

Combining (where index i represent a site and runs from 1:N_sites)
∑ᵢ nᵢ = n # total no.density of neutrinos
∑ᵢ  Nᵢ = N # total no.of neutrinos
∑ᵢ Vᵢ = V # total volume of the grid
∑ᵢ Δxᵢ = L # domain size

"""
function Neutrino_number(s, τ, B,L, N_sites, N_sites_eachflavor, tolerance,
    n_νₑ,n_νₑ̄,Eνₑ,Eνₑ̄,Δx,Δm², p, x, Δp, ψ₀, shape_name, energy_sign, cutoff, maxdim, ttotal,periodic)
    V = L^3 # volume of the big box containing all sites/particles

    # Create an array of dimension N and fill it half with values of sites containing all electron neutrinos 
    # and other half with sites containing electron anti-neutrino. 
    N_νₑ  = n_νₑ * V 
    N_1 = fill(N_νₑ / (N_sites ÷ 2), N_sites ÷ 2)
    N_νₑ̄  = n_νₑ̄ * V 
    N_2 = fill(N_νₑ̄/ (N_sites ÷ 2), N_sites ÷ 2)
    # N_1 = fill(N_νₑ, N_sites ÷ 2)
    # N_νₑ̄  = n_νₑ̄ * V 
    # N_2 = fill(N_νₑ̄, N_sites ÷ 2)
    N = vcat(N_1, N_2) # This is the total number of neutrinos. 
    return N 

end
# generate x_array such that the first particle is at position L/(2*N_sites) while subsequent particles are at a position incremental by L/N_sites. # grid style
function generate_x_array(N_sites, L)
    return [(i - 0.5) * L / N_sites for i in 1:N_sites]
end

#generate a momentum array that depicts the energy of neutrinos and anti-neutrinos in opposing beams # for bipolar rogerro test 
function generate_p_array(N_sites)                                                                                                                                                                                   
    half_N_sites = div(N_sites, 2)
    return [fill(1, half_N_sites); fill(1, half_N_sites)]
end

function generate_px_array(N_sites, Eνₑ, Eνₑ̄)                                                                                                                                                                                   
    half_N_sites = div(N_sites, 2)
    return [fill(Eνₑ, half_N_sites); fill(Eνₑ̄, half_N_sites)]
end


function generate_py_array(N_sites)                                                                                                                                                                                   
    half_N_sites = div(N_sites, 2)
    return [fill(0, half_N_sites); fill(0, half_N_sites)]
end

function generate_pz_array(N_sites)                                                                                                                                                                                   
    half_N_sites = div(N_sites, 2)
    return [fill(0, half_N_sites); fill(0, half_N_sites)]
end