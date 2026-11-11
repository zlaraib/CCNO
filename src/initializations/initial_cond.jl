
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
length of each side = Delta_x 
Delta_x³ = Vᵢ
Nᵢ =  total no.of neutrinos at site i 
nᵢ  = total no.density of neutrinos at site i 
Nᵢ = nᵢ  * Vᵢ

Combining (where index i represent a site and runs from 1:N_sites)
∑ᵢ nᵢ = n # total no.density of neutrinos
∑ᵢ  Nᵢ = N # total no.of neutrinos
∑ᵢ Vᵢ = V # total volume of the grid
∑ᵢ Delta_xᵢ = L # domain size

"""
# generate x_array such that the first particle is at position L/(2*N_sites_eachflavor) while subsequent particles are at a position incremental by L/N_sites_eachflavor. 
function generate_mu_x_array(N_sites_eachflavor, L)
    return [(i - 0.5) * L / N_sites_eachflavor for i in 1:N_sites_eachflavor]
end

# Similar function for anti-mu particles for now, but  keeping it separate for easier generalization of different mu/anti mu and their initial config for later.
function generate_antimu_x_array(N_sites_eachflavor, L)
    return [(i - 0.5) * L / N_sites_eachflavor for i in 1:N_sites_eachflavor]
end
# concatenate the arrays for both mu and anti-mu particles into a single array.
function generate_x_array(N_sites_eachflavor, L) 
    x_array = vcat(generate_mu_x_array(N_sites_eachflavor, L), generate_antimu_x_array(N_sites_eachflavor, L))
    return x_array 
end 

#generate a momentum array that depicts the energy of neutrinos and anti-neutrinos in opposing beams # for bipolar rogerro test 
function generate_p_array(N_sites)
    half_N_sites = div(N_sites, 2)
    return [fill(1, half_N_sites); fill(1, half_N_sites)]
end

function generate_px_array(N_sites, Enu_e, Enu_ē)
    half_N_sites = div(N_sites, 2)
    return [fill(Enu_e, half_N_sites); fill(Enu_ē, half_N_sites)]
end

function generate_py_array(N_sites)
    half_N_sites = div(N_sites, 2)
    return [fill(0, half_N_sites); fill(0, half_N_sites)]
end

function generate_pz_array(N_sites)
    half_N_sites = div(N_sites, 2)
    return [fill(0, half_N_sites); fill(0, half_N_sites)]
end

function generate_params_from_beamarr(beam_arr, L, n_nu, bonddim)

    num_beams = size(beam_arr)[1]
    
    total_sites = 0
    for i = 1:num_beams
        total_sites = total_sites + beam_arr[i].num_sites
    end
    
    ds = (2 * L) / total_sites
    
    x = Float64[]
    y = Float64[]
    z = Float64[]
    p = Array{Float64}(undef, 0, 3)
    energy_sign = Int[]
    spin_vec = String[]
    
    for i = 1:num_beams
        num_sites = beam_arr[i].num_sites
        if num_sites == 0
            continue
        end
        
        spacing = (L) / beam_arr[i].num_sites

        x_tmp = fill(0.0, num_sites) 
        y_tmp = fill(0.0, num_sites)
        z_tmp = fill(0.0, num_sites)
        
        for j = 1:num_sites
            x_tmp[j] = (j - 0.5) * spacing
        end
        
        px_tmp = fill(beam_arr[i].energy, num_sites) 
        py_tmp = fill(0.0, num_sites)
        pz_tmp = fill(0.0, num_sites)
        p_tmp = hcat(px_tmp, py_tmp, pz_tmp)

        theta = beam_arr[i].rotation
        
        rotation = [
            cos(theta) -sin(theta) 0;
            sin(theta) cos(theta) 0;
            0 0 1
        ]
        
        for j = 1:num_sites
            p_tmp[j, :] = rotation * p_tmp[j, :]
        end
        
        if beam_arr[i].flavor == 1
            spin_vec_tmp = ["Up"]
            energy_sign_tmp = [1]
            for j = 2:num_sites
                append!(spin_vec_tmp, ["Up"])
                append!(energy_sign_tmp, [1])
            end
        else
            spin_vec_tmp = ["Dn"]
            energy_sign_tmp = [-1]
            for j = 2:num_sites
                append!(spin_vec_tmp, ["Dn"])
                append!(energy_sign_tmp, [-1])
            end
        end
        
        append!(x, x_tmp)
        append!(y, y_tmp)
        append!(z, z_tmp)
        p = vcat(p, p_tmp)
        append!(energy_sign, energy_sign_tmp)
        append!(spin_vec, spin_vec_tmp)
    end
    
    xyz = hcat(x, y, z)
    
    #fixed #switched conserve_qns to false to avoid fluxes error in expect function
    s = siteinds("S=1/2", total_sites; conserve_qns = false)
    xyz = hcat(x, y, z)

    Psi = productMPS(s, spin_vec)

    N = n_nu * (ds)^3
    
    return CCNO.Parameters(
        N_sites = total_sites,
        τ = 5E-13,
        ttotal = 9.0E-11, 
        m1 = 0.0,
        m2 = 0.0,
        maxdim = bonddim,
        cutoff = 1e-100,
        theta_nu = 1.74532925E-8,
        shape_name = "triangular",
        homogeneous = [false, true, true], # homogeneous in y, z dims
        geometric_name = "physical",
        Delta_x = ds,
        L = L,
        Delta_p = ds,
        periodic = true,
        checkpoint_every = 20,
        do_recover = false,
        recover_file = "",
        plotdir = joinpath(@__DIR__, "plots"),
        datadir = joinpath(@__DIR__, "datafiles"),
        chkptdir = joinpath(@__DIR__, "checkpoints"),
        save_plots_flag = false,
        alpha = 1e-6,
    ),  CCNO.SimulationState(
        Psi = Psi,
        s = s,
        s0 = s,
        p = p,
        energy_sign = energy_sign,
        N = N,
        xyz = xyz,
    )


end
