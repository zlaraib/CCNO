# Function to append data as new rows to the existing file after every `checkpoint_every` values
function append_data(filename, new_data)
    open(filename, "a") do f
        writedlm(f, new_data)
    end
end

function store_data(datadir, t, ψ, x, p)
    #================================#
    # calculate necessary quantities #
    #================================#

    # compute the avg expectation value of Sz at all sites
    sz_arr = expect(ψ, "Sz")  # Compute Sz for each site and store the values in sz_arr
    
    # compute expectation value of sy and sx (inbuilt operator in ITensors library) at all sites on the chain
    sy_arr = expect(complex(ψ), "Sy")
    sx_arr = expect(ψ, "Sx")
    
    #survival probability for all sites (neutrino) to be found in its initial flavor state
    prob_surv_arr = 0.5 * (1 .- 2 .* sz_arr)

    # recall that in our code sigma_z = 2*Sz so make sure these expressions are consistent with "Sz in ITensors" 
    ρₑₑ_arr = ((2 .* sz_arr) .+ 1) ./ 2
    
    ρ_μμ_arr = ((-2 .* sz_arr) .+ 1) ./ 2
    
    ρₑμ_arr = sqrt.(sx_arr.^2 .+ sy_arr.^2)
    
    #===============#
    # save the data #
    #===============#
    
    mkpath(datadir)
    # Writing data to files with corresponding headers
    # 1. Sz arrays
    fname1 = joinpath(datadir, "t_<Sz>.dat")
    append_data(fname1, [t transpose(sz_arr)])

    # 2. Sy arrays
    fname2 = joinpath(datadir, "t_<Sy>.dat")
    append_data(fname2, [t transpose(sy_arr)])

    # 3. Sx arrays
    fname3 = joinpath(datadir, "t_<Sx>.dat")
    append_data(fname3, [t transpose(sx_arr)])

    # 4. prob_surv arrays
    fname4 = joinpath(datadir, "t_probsurv.dat")
    append_data(fname4, [t transpose(prob_surv_arr)])

    # 5. xsite values
    fname5 = joinpath(datadir, "t_xsiteval.dat")
    append_data(fname5, [t transpose(x)])

    # 6. pxsite values
    fname6 = joinpath(datadir, "t_pxsiteval.dat")
    append_data(fname6, [t transpose(p[:,1])])
    
    # 7. ρₑₑ arrays
    fname7 = joinpath(datadir, "t_ρₑₑ.dat")
    append_data(fname7, [t transpose(ρₑₑ_arr)])

    # 8. ρ_μμ arrays
    fname8 = joinpath(datadir, "t_ρ_μμ.dat")
    append_data(fname8, [t transpose(ρ_μμ_arr)])

    # 9. ρₑμ arrays
    fname9 = joinpath(datadir, "t_ρₑμ.dat")
    append_data(fname9, [t transpose(ρₑμ_arr)])
    
end 
