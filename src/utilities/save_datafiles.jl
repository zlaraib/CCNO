# Function to append data as new rows to the existing file after every `checkpoint_every` values
function append_data(filename, new_data)
    open(filename, "a") do f
        writedlm(f, new_data)
    end
end

function store_data(datadir, t, sz_tot, sy_tot, sx_tot, prob_surv_tot,x, px, ρₑₑ_tot,ρ_μμ_tot, ρₑμ_tot)
    
    mkpath(datadir)
    # Writing data to files with corresponding headers
    # 1. Sz arrays
    fname1 = joinpath(datadir, "t_<Sz>.dat")
    append_data(fname1, [t transpose(sz_tot)])

    # 2. Sy arrays
    fname2 = joinpath(datadir, "t_<Sy>.dat")
    append_data(fname2, [t transpose(sy_tot)])

    # 3. Sx arrays
    fname3 = joinpath(datadir, "t_<Sx>.dat")
    append_data(fname3, [t transpose(sx_tot)])

    # 4. prob_surv arrays
    fname4 = joinpath(datadir, "t_probsurv.dat")
    append_data(fname4, [t transpose(prob_surv_tot)])

    # 5. xsite values
    fname5 = joinpath(datadir, "t_xsiteval.dat")
    append_data(fname5, [t transpose(x)])

    # 6. pxsite values
    fname6 = joinpath(datadir, "t_pxsiteval.dat")
    append_data(fname6, [t transpose(px)])
    
    # 7. ρₑₑ arrays
    fname7 = joinpath(datadir, "t_ρₑₑ.dat")
    append_data(fname7, [t transpose(ρₑₑ_tot)])

    # 8. ρ_μμ arrays
    fname8 = joinpath(datadir, "t_ρ_μμ.dat")
    append_data(fname8, [t transpose(ρ_μμ_tot)])

    # 9. ρₑμ arrays
    fname9 = joinpath(datadir, "t_ρₑμ.dat")
    append_data(fname9, [t transpose(ρₑμ_tot)])
    
end 
