push!(LOAD_PATH, "..")
using CCNO

using ITensors
using ITensorMPS
using Plots
using Measures
using DelimitedFiles
using HDF5

save_data = false  # true = saves datafiles for science runs while false doesn't. So change it to false for jenkins test runs
save_plots_flag = false # true = saves plots for science runs while false doesn't. So change it to false for jenkins test runs

# this file plots δω on x axis while minimum time tp on y axis and compares
# the values of our calculations with Rogerros calculations in Table I.
# Here symmetric δω is used i.e. ω_a = - ω_b for a given δω plotted on x axis.

N_sites = 4
ttotal = 5
function main(Δω, N_sites, ttotal)
    cutoff = 1E-14
    τ = 0.005
    tolerance  = 5E-1
    Δx = 1E-3
    maxdim = 1000 #bond dimension
    L = 10 # cm # not being used in this test but defined to keep the evolve function arguments consistent.
    Δp = L # width of shape function # not being used in this test but defined to keep the evolve function arguments consistent. 
    t1 = 0.0084003052 #choose initial time for growth rate calculation #variable, not being used in this test
    t2 = 0.011700318 #choose final time for growth rate calculation #variable, not being used in this test 
    periodic = false 
    checkpoint_every = 4
    do_recover = false
    recover_file = "" 
    
    if Δω==-0.5
        a_t = 0
        b_t = 0
        c_t = 1.82
    end
    if Δω==0.0
        a_t = 0.0
        b_t = 2.105
        c_t = 0
    end
    if Δω==0.05
        a_t = 2.586
        b_t = 0
        c_t = 0
    end
    if Δω==0.125
        a_t = 1.656
        b_t = 0
        c_t = 1.42
    end
    if Δω==0.25
        a_t = 1.224
        b_t = 0
        c_t = 1.62
    end
    if Δω==0.5
        a_t = 1.064
        b_t = 0
        c_t = 1.42
    end
    if Δω==1.0
        a_t = 0.965
        b_t = 0
        c_t = 0
    end
    Δm² = Δω
    s = siteinds("S=1/2", N_sites; conserve_qns=false)
    mu = ones(N_sites)
    N = mu .* fill(((Δx)^3 )/(√2 * CCNO.G_F * N_sites), N_sites)
    theta_nu = 0 # mixing_angle #rad 
    B = [sin(2*theta_nu), 0, -cos(2*theta_nu)] # is equivalent to B = [0, 0, -1] # fixed for Rogerro's case
    B = B / norm(B)

    # p matrix with numbers generated from the p_array for all components (x, y, z)
    p = hcat(CCNO.generate_p_array(N_sites),fill(0, N_sites), fill(0, N_sites))
    x = fill(rand(), N_sites) # variable.
    y = fill(rand(), N_sites) # variable.
    z = fill(rand(), N_sites) # variable.

    ψ = productMPS(s, N -> N <= N_sites/2 ? "Dn" : "Up")
    energy_sign = [i <= N_sites ÷ 2 ? 1 : 1 for i in 1:N_sites]
    shape_name = "none"  # Change this to the desired shape name # variable.
    # Specify the relative directory path
    datadir = joinpath(@__DIR__, "datafiles")
    chkptdir = joinpath(@__DIR__, "checkpoints")

    #extract output for the survival probability values at each timestep
    Sz_array, Sy_array, Sx_array,  prob_surv_array, x_values, pₓ_values, ρₑₑ_array, ρ_μμ_array, ρₑμ_array, t_array, t_recover = CCNO.evolve(
        s, τ, N, B, L, N_sites, Δx, Δm², p, x, Δp, theta_nu, ψ, shape_name, energy_sign, cutoff, maxdim, datadir, t1, t2, ttotal,chkptdir, checkpoint_every,  do_recover, recover_file, save_data , periodic)
    
    # extract the prob_surv on the first site 
    prob_surv_array_site1= [row[1] for row in prob_surv_array]
    function find_first_local_minima_index(arr)
        N = length(arr)
        for i in 2:(N-1)
            if arr[i] < arr[i-1] && arr[i] < arr[i+1]
                return i
            end
        end
        return -1  
    end
    
    tmin_ifirstlocalmin_file = joinpath(datadir, "tmin_ifirstlocalmin.dat")
    if !do_recover 
        # Index of first minimum of the prob_surv_array (containing survival probability values at each time step)
        global i_first_local_min = find_first_local_minima_index(prob_surv_array_site1)
        
        # Writing if_else statement to communicate if local minima (not) found
        if i_first_local_min != -1
            println("Index of the first local minimum: ", i_first_local_min)
            global t_min = τ * i_first_local_min - τ
            println(t_min)
        else
            global t_min = nothing
            println(t_min)
            println("No local minimum found in the array.")
        end
        println("Corresponding time of first minimum index= ", t_min)


        end
    if save_data
        # Save t_min and i_first_local_min to file
        writedlm(tmin_ifirstlocalmin_file, [t_min i_first_local_min])
        println("Saved t_min to file: ", tmin_ifirstlocalmin_file)
    else global t_min = t_min
    end
    if do_recover  

        # Read the saved t_min and i_first_local_min from the file
        tmin_data = readdlm(tmin_ifirstlocalmin_file)

        # Extracting t_min and i_first_local_min from the read data
        t_min_value = tmin_data[1]  # First value
        i_first_local_min_value = tmin_data[2]  # Second value

        # Convert `t_min` and `i_first_local_min` to appropriate types if needed
        t_min = t_min_value isa AbstractString ? tryparse(Float64, t_min_value) : t_min_value
        i_first_local_min = i_first_local_min_value isa AbstractString ? tryparse(Int, i_first_local_min_value) : i_first_local_min_value

        println("Recover t_min=",t_min)
        println("Recover i_first_local_min=",i_first_local_min)

        if t_min !== nothing && i_first_local_min != -1 
            # t_min = readdlm(tmin_ifirstlocalmin_file)[1]  # assign the previous t_min
            println("Recovered first local t_min from previous run: ", t_min)
            println("Recovered index of first local t_min from previous run: ", i_first_local_min)
        end
        if t_min === nothing && i_first_local_min == -1 
            i_first_local_min = find_first_local_minima_index(prob_surv_array_site1)
            t_min = (τ * i_first_local_min) - τ + t_recover
            println("Recalculated t_min=",t_min)
            println("Recalculated i_first_local_min=",i_first_local_min)
            if i_first_local_min != -1
                println("Index of the first local minimum: ", i_first_local_min)
                println("Final t_min=",t_min)
            else
                t_min = nothing
                i_first_local_min = -1
                println(t_min, i_first_local_min)
                println("No local minimum found in the array.")
            end
            
        end
        # Save t_min and i_first_local_min to file (again, updated if required)
        writedlm(tmin_ifirstlocalmin_file, [t_min i_first_local_min])
        println("Saved t_min to file: ", tmin_ifirstlocalmin_file)
    end 

    println("Corresponding time of first minimum index= ", t_min)


    # Rogerro(2021)'s fit for the first minimum of the survival probability reached for a time t_p 
    t_p_Rog = a_t*log(N_sites) + b_t * sqrt(N_sites) + c_t
    println("t_p_Rog= ",t_p_Rog)

   
    if save_data
        t_min_value = readdlm(tmin_ifirstlocalmin_file)[1]  # read the previous t_min
        t_min = t_min_value isa AbstractString ? tryparse(Float64, t_min_value) : t_min_value
        println(t_min)
    end

    # Check that our time of first minimum survival probability compared to Rogerro(2021) remains within the timestep and tolerance.
    if t_min !== nothing
        @assert abs(t_min - t_p_Rog) < τ + tolerance
    elseif t_min === nothing
        if do_recover
            t_min = (τ * i_first_local_min) - τ + t_recover
        else
            t_min = τ * i_first_local_min - τ
        end
        @assert abs(t_min - t_p_Rog) <  τ + tolerance
    end
    return t_p_Rog, t_min
end

# Arrays to store t_p_Rog and t_min for each Δω for 1st site
Δω_values = [-0.5, 0.0, 0.125, 0.25, 0.5, 1.0]
t_p_Rog_array = Float64[]
t_min_array = Float64[]

datadir = joinpath(@__DIR__, "datafiles")

for Δω in Δω_values
    t_p_Rog, t_min = main(Δω, N_sites, ttotal)
    push!(t_p_Rog_array, t_p_Rog)
    push!(t_min_array, t_min)
end
if save_data
    save_data= isdir(datadir) || mkpath(datadir)
    fname1 = joinpath(datadir, "δω_tpRog_tpmine.dat")
    writedlm(fname1, [Δω_values t_p_Rog_array t_min_array])
end 
if save_plots_flag
    plotdir = joinpath(@__DIR__, "plots")
    # check if a directory exists, and if it doesn't, create it using mkpath
    save_plots_flag = isdir(plotdir) || mkpath(plotdir)
    # Create the plot
    plot(Δω_values, t_p_Rog_array, label="Rogerro(2021)", xlabel="δω", ylabel="Minimum Time(tₚ)", title="Table I. Rogerro(2021) ", aspect_ratio=:auto,margin= 10mm)
    plot!(Δω_values, t_min_array, label="Our results")
    savefig(joinpath(plotdir,"t_p_vs_symmetric del_omega for N_sites$(N_sites).pdf"))
end 
if !save_plots_flag
    # Create the plot
    plot(Δω_values, t_p_Rog_array, label="Rogerro(2021)", xlabel="δω", ylabel="Minimum Time(tₚ)", title="Table I. Rogerro(2021) ", aspect_ratio=:auto,margin= 10mm)
    plot!(Δω_values, t_min_array, label="Our results")
    savefig("t_p_vs_symmetric del_omega for N_sites$(N_sites).pdf")
end
