
include("constants.jl")
include("geometric_func.jl")
include("shape_func.jl")
include("momentum.jl")


"""
    Expected (CGS) units of the quantities defined in the files in tests directory that are being used in the gates function.                                                                   
    s = site index array (dimensionless and unitless)          
    N = array of no.of neutrinos contained on each site (dimensionless and unitless)
    B = array of normalized vector related to mixing angle in vacuum oscillations (dimensionless constant)
    N_sites = Total no.of sites (dimensionless and unitless)
    Δx = length of the box of interacting neutrinos at a site (cm)
    Δm² = difference in mass squared (erg^2)
    p = array of momentum vectors (erg)
    x = array of positions of sites (cm)
    Δp = width of shape function (cm)
    shape_name = name of the shape function (string) ["none","triangular","flat_top"]
    τ = time step (sec)
    energy_sign = array of sign of the energy (1 or -1): 1 for neutrinos and -1 for anti-neutrinos
    maxdim = max bond dimension in MPS truncation (unitless and dimensionless)
    cutoff = truncation threshold for the SVD in MPS representation (unitless and dimensionless)
    periodic = boolean indicating whether boundary conditions should be periodic
"""

# This file generates the create_gates function that holds ITensors Trotter gates and returns the dimensionless unitary 
# operators govered by the Hamiltonian which includes effects of the vacuum and self-interaction potential for each site.

function create_gates(s, N, B, N_sites, Δx, Δm², p, x, Δp, theta_nu, shape_name,L, τ, energy_sign, periodic)
    
    # Make gates (1,2),(2,3),(3,4),... i.e. unitary gates which act on any (non-neighboring) pairs of sites in the chain.
    # Create an empty ITensors array that will be our Trotter gates
    gates = ITensor[] 

    # extract output of p_hat and p_mod for the p vector defined above for all sites. 
    p_mod, p̂ = momentum(p,N_sites)  
    
    # define an array of vacuum oscillation frequencies (units of ergs)
    if Δm² == 0 # specific to self-int only
        ω = zeros(N_sites)
    elseif Δm² == 2 * π # specific to vac_osc only
       global ω = fill(π, N_sites) # added global so we can access and use this global variable without the need to pass them as arguments to another function
    elseif (Δm² == 0.5 || Δm² == 0.2 || Δm² == 2) &&  L==1 # addition for full Hamiltonian from main_Rogerro, Rog_bipolar, Rog_N_loop and t_p_vs_N_unsym tests 
        # Create arrays ω_a and ω_b
        global ω_a = fill(Δm², div(N_sites, 2))
        global ω_b = fill(0, div(N_sites, 2))
        # Concatenate ω_a and ω_b to form ω
        ω = vcat(ω_a, ω_b)
    elseif (Δm² == -0.5 || Δm² == 0.0 || Δm² ==0.05 || Δm² ==0.125 || Δm² ==0.25 || Δm² ==0.5 || Δm² ==1.0) && L==10 # addition for t_p_vs_N_sym and t_p_vs_sym_delta_w tests   
        Δω_array= fill(Δm², div(N_sites, 2))
        # Calculate ω_a and ω_b based on Δω
        global ω_a = Δω_array 
        global ω_b = -Δω_array 
        ω = vcat(ω_a, ω_b)
    else 
        ω = [Δm² / (2 * p_mod[i]) * energy_sign[i] for i in 1:N_sites]
    end
    println("ω = ", ω)

    for i in 1:(N_sites-1)
        for j in i+1:N_sites
            #s_i, s_j are non-neighbouring spin site/indices from the s array
            s_i = s[i]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
            s_j = s[j]
            # assert B vector to have a magnitude of 1 while preserving its direction.
            @assert norm(B) == 1

            # Our neutrino system Hamiltonian of self-interaction term represents 1D Heisenberg model.
            # total Hamiltonian of the system is a sum of local terms hj, where hj acts on sites i and j which are paired for gates to latch onto.
            # op function returns these operators as ITensors and we tensor product and add them together to compute the operator hj.
            # ni and nj are the neutrions at site i and j respectively.
            # mu pairs divided by 2 to avoid double counting
            
            # if energy_sign[i]*energy_sign[j]>0
                # Get the shape function result for each pair of i and j 
                shape_result = shape_func(x, Δp, i, j,L, shape_name, periodic)
                # Calculate the geometric factor for each pair of i and j within the loop
                geometric_factor = geometric_func(p, p̂, i, j, theta_nu)
                interaction_strength = (2.0* √2 * G_F * (N[i]+ N[j])/(2*((Δx)^3))) * shape_result * geometric_factor
                hj = interaction_strength *
                (op("Sz", s_i) * op("Sz", s_j) +
                1/2 * op("S+", s_i) * op("S-", s_j) +
                1/2 * op("S-", s_i) * op("S+", s_j))
            # else
                # # Get the shape function result for each pair of i and j 
                # shape_result = shape_func(x, Δp, i, j,L, shape_name, periodic)
                # # Calculate the geometric factor for each pair of i and j within the loop
                # geometric_factor = geometric_func(p, p̂, i, j, theta_nu)
                # interaction_strength = (2.0* √2 * G_F * (N[i]+ N[j])/(2*((Δx)^3))) * shape_result * geometric_factor
                # hj = - interaction_strength * 
                # ((-2 *op("Sz",s_i) * op("Sz",s_j)) + 
                # op("S+", s_i) * op("S-", s_j) +
                # op("S-", s_i) * op("S+", s_j))
            # end

            # add Vacuum Oscillation Hamiltonian 
            if ω[i] != 0 || ω[j] != 0
                
                hj += (1/(N_sites-1))*( 
                    (ω[i] * B[1] * op("Sx", s_i)* op("Id", s_j))  + (ω[i] * B[2] * op("Sy", s_i)* op("Id", s_j))  + (ω[i] * B[3] * op("Sz", s_i)* op("Id", s_j)) )
                hj += (1/(N_sites-1))*(
                    (ω[j] * B[1] * op("Id", s_i) * op("Sx", s_j)) + (ω[j] * B[2]  * op("Id", s_i)* op("Sy", s_j)) + (ω[j] * B[3]  * op("Id", s_i)* op("Sz", s_j)) )
            end
            
            # make Trotter gate Gj that would correspond to each gate in the gate array of ITensors             
            if theta_nu == 0 ||  theta_nu == π/4 || theta_nu == π/2 
                Gj = exp(-im * τ/2 * hj)
            elseif theta_nu == 0.1 && Δm²== 0.2 # for Rog_bipolar
                Gj = exp(-im * τ/2 * hj)
            elseif theta_nu == 0.01 # for Richers bipolar
                t_bipolar = 8.96e-4
                Gj = exp(-im * τ/2 * hj* t_bipolar/hbar)
            else 
                Gj = exp(-im * τ/2 * hj* 1/hbar)
            end
            # println(imag(hj/hbar))
            # println((Δm²))/(2 *hbar * p_mod[1])
            # @assert imag(hj) == ((Δm²))/(2 *hbar * p_mod[1])
            # println("Gj= ",Gj)            # has_fermion_string(hj) = true
            # The push! function adds (appends) an element to the end of an array;
            # ! performs an operation without creating a new object, (in a way overwites the previous array in consideration); 
            # i.e. we append a new element Gj (which is an ITensor object representing a gate) to the end of the gates array.
            push!(gates, Gj)
        end 
    end

    # append! adds all the elements of a gates in reverse order (i.e. (N_sites,N_sites-1),(N_sites-1,N_sites-2),...) to the end of gates array.
    # appending reverse gates to create a second-order Trotter-Suzuki integration
    append!(gates, reverse(gates))
    return gates
end