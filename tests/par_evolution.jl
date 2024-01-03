    
    using ITensors
    using Plots
    using Measures
    using LinearAlgebra
    using DelimitedFiles
    include("../src/momentum.jl")
    include("../src/constants.jl")


    """
    Expected (CGS) units of the quantities defined in the files in tests directory that are being used in the evolve function.                                                                            
    N = array of no.of neutrinos contained on each site (dimensionless and unitless)
    N_sites = Total no.of sites (dimensionless and unitless)
    L = domain size (cm)
    p = array of momentum vectors (erg)
    x = array of positions of sites (cm)
    τ = time step (sec)
    periodic = boolean indicating whether boundary conditions should be periodic
    (anti)neutrino_energy = energy of (anti)neutrinos in ergs
    """

    # This file generates the evolve function which uses particles confined in a domain and tracks the particles displacement in time for a certain boundary condition

    N_sites = 4  # number of sites # variable
    L = 1e7 # cm # domain size # (aka big box length)
    τ = 1.66e-4 # time step to include 50 steps every 10 picoseconds # sec # variable
    ttotal = 1.66e-2 # total time of evolution # sec #variable
    periodic = true # true = imposes periodic boundary conditions while false doesn't
    neutrino_energy =  50.0e6 # energy of all neutrinos (P.S the its negative is energy of all antineutrinos)
    antineutrino_energy = -1 * neutrino_energy # specific to my case only. Since all neutrinos have same energy, except in my case anti neutrinos are moving in opposite direction to give it a negative sign

    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "..","misc","datafiles","evol", "par_"*string(N_sites), "tt_"*string(ttotal))
    # Specify the relative directory path
    plotdir = joinpath(@__DIR__, "..","misc","plots","evol", "par_"*string(N_sites), "tt_"*string(ttotal))

    # check if a directory exists, and if it doesn't, create it using mkpath
    isdir(plotdir) || mkpath(plotdir)
    # check if a directory exists, and if it doesn't, create it using mkpath
    isdir(datadir) || mkpath(datadir)
    # generate x_array such that the first particle is at position L/(2*N_sites) while subsequent particles are at a position incremental by L/N_sites. # grid style
    function generate_x_array(N_sites, L)
        return [(i - 0.5) * L / N_sites for i in 1:N_sites]
    end
    x = generate_x_array(N_sites, L)
    println("Initial positions of the particles=",x)


    function generate_px_array(N_sites)
        half_N_sites = div(N_sites, 2)
        return [fill(neutrino_energy, half_N_sites); fill(antineutrino_energy, half_N_sites)]
    end

    function generate_py_array(N_sites)                                                                                                                                                                                   
        half_N_sites = div(N_sites, 2)
        return [fill(0, half_N_sites); fill(0, half_N_sites)]
    end
    
    function generate_pz_array(N_sites)                                                                                                                                                                                   
        half_N_sites = div(N_sites, 2)
        return [fill(0, half_N_sites); fill(0, half_N_sites)]
    end
    
    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(generate_px_array(N_sites), generate_py_array(N_sites), generate_pz_array(N_sites))
    println("Initial p_vector of all particles=",p)

function evolve(τ, L, N_sites, p, x, ttotal,periodic)
    x_values = []
    px_values = []
    p_mod, p_hat = momentum(p, N_sites)
    p_x_hat = [sub_array[1] for sub_array in p_hat]
    t_array = [] # to store t values 

    for t in 0.0:τ:ttotal
        push!(x_values, copy(x))
        px = p[:, 1] #p_x of all sites/particles
        push!(px_values, copy(px))

        # method 1 of evolving particles within the domain "with" the nested for loop on particles 
        for i in 1:N_sites
            println("particle $i's position at time $t (before evolution) = $(x[i])")
            x[i] += p_x_hat[i] * c * τ
            println("particle $i's position at time $t (after evolution) = $(x[i])")
            if periodic
                # wrap around position from 0 to domain size L
                x[i] = mod(x[i],L)
                println("rearranged particle $i's position at time $t within domain = $(x[i])")
                # Checking if the updated x[i] satisfies the boundary conditions
                @assert (x[i] >= 0 && x[i] <= L)
            end
        end

        # method 2 of evolving particles within the domain "without" the nested for loop on particles 
        # x .+=  ((p_x_hat.*c) .* τ)  # displacing particle's position at each timestep 
        # if periodic
        #     # Wrap around position from 0 to domain size L
        #     x = mod.(x, L)
        #     println(x)
        #     # Check if the updated x satisfies the boundary conditions
        #     @assert all(x .>= 0) && all(x .<= L)
        # end

        t ≈ ttotal && break
    end
    t_array = 0.0:τ:ttotal
    fname3 = joinpath(datadir, "t_xsiteval.dat")
    writedlm(fname3, [t_array x_values])
    fname4 = joinpath(datadir, "t_pxsiteval.dat")
    writedlm(fname4, [t_array px_values])
    return x_values, px_values
end
x_values,px_values = evolve(τ, L, N_sites, p, x, ttotal,periodic)
#println(x_values)

plot(title="Particle Position Evolution", xlabel= "Position (x)",ylabel="Time")
for site in 1:N_sites
    site_positions = [(x_values[t][site]) for t in 1:length(x_values)]
    plot!(site_positions, 0.0:τ:ttotal, label="Site $site",left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm)
end

savefig(joinpath(plotdir,"Particles position(x) evolution.pdf"))

plot(title="Particle Momentum Evolution", xlabel= "Momentum in x direction(p_x)",ylabel="Time")
for site in 1:N_sites
    site_momentum = [(px_values[t][site]) for t in 1:length(px_values)]
    plot!(site_momentum, 0.0:τ:ttotal, label="Site $site",left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm)
end

savefig(joinpath(plotdir,"Particles momentum(p_x) evolution.pdf"))