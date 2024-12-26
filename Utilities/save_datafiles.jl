

function store_data(do_recover, datadir, iteration, checkpoint_every, t_array, Sz_array, Sy_array, Sx_array, prob_surv_array,x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array)
    
    if !do_recover
        save_data = isdir(datadir) || mkpath(datadir)
        # Writing data to files with corresponding headers
        # 1. Sz arrays
        fname1 = joinpath(datadir, "t_<Sz>.dat")
        fname1_1 = joinpath(datadir, "t_<Sz>.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname1, [t_array Sz_array])
        writedlm(fname1_1, [t_array Sz_array])

        # 2. Sy arrays
        fname2 = joinpath(datadir, "t_<Sy>.dat")
        fname2_1 = joinpath(datadir, "t_<Sy>.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname2, [t_array Sy_array])
        writedlm(fname2_1, [t_array Sy_array])

        # 3. Sx arrays
        fname3 = joinpath(datadir, "t_<Sx>.dat")
        fname3_1 = joinpath(datadir, "t_<Sx>.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname3, [t_array Sx_array])
        writedlm(fname3_1, [t_array Sx_array])

        # 4. prob_surv arrays
        fname4 = joinpath(datadir, "t_probsurv.dat")
        fname4_1 = joinpath(datadir, "t_probsurv.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname4, [t_array prob_surv_array])
        writedlm(fname4_1, [t_array prob_surv_array])

        # 5. xsite values
        fname5 = joinpath(datadir, "t_xsiteval.dat")
        fname5_1 = joinpath(datadir, "t_xsiteval.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname5, [t_array x_values])
        writedlm(fname5_1, [t_array x_values])

        # 6. pxsite values
        fname6 = joinpath(datadir, "t_pxsiteval.dat")
        fname6_1 = joinpath(datadir, "t_pxsiteval.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname6, [t_array pₓ_values])
        writedlm(fname6_1, [t_array pₓ_values])

        # 7. ρₑₑ arrays
        fname7 = joinpath(datadir, "t_ρₑₑ.dat")
        fname7_1 = joinpath(datadir, "t_ρₑₑ.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname7, [t_array ρₑₑ_array])
        writedlm(fname7_1, [t_array ρₑₑ_array])

        # 8. ρ_μμ arrays
        fname8 = joinpath(datadir, "t_ρ_μμ.dat")
        fname8_1 = joinpath(datadir, "t_ρ_μμ.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname8, [t_array ρ_μμ_array])
        writedlm(fname8_1, [t_array ρ_μμ_array])

        # 9. ρₑμ arrays
        fname9 = joinpath(datadir, "t_ρₑμ.dat")
        fname9_1 = joinpath(datadir, "t_ρₑμ.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
        writedlm(fname9, [t_array ρₑμ_array])
        writedlm(fname9_1, [t_array ρₑμ_array])

        # fname8 = joinpath(datadir, "Im_Ω.dat")
        # writedlm(fname8, [Im_Ω])


    elseif do_recover

        # Determine the start index for appending data
        start_index = iteration >= checkpoint_every ? checkpoint_every + 1 : 0

        # Function to append data as new rows to the existing file after every `checkpoint_every` values
        function append_data(filename, new_data)
            open(filename, "a") do f
                writedlm(f, new_data)
            end
        end

        # write separate data file from the recovered loop and Append new values for each of the saved files, starting from the appropriate index
        if start_index > 0
            # 1. Sz arrays
            fname1 = joinpath(datadir, "t_<Sz>.dat")
            fname1_1 = joinpath(datadir, "t_<Sz>_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname1_1, [t_array Sz_array])
            append_data(fname1, [t_array[start_index:end] Sz_array[start_index:end]])

            # 2. Sy arrays
            fname2 = joinpath(datadir, "t_<Sy>.dat")
            fname2_1 = joinpath(datadir, "t_<Sy>_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname2_1, [t_array Sy_array])
            append_data(fname2, [t_array[start_index:end] Sy_array[start_index:end]])
    
            # 3. Sx arrays
            fname3 = joinpath(datadir, "t_<Sx>.dat")
            fname3_1 = joinpath(datadir, "t_<Sx>_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname3_1, [t_array Sx_array])
            append_data(fname3, [t_array[start_index:end] Sx_array[start_index:end]])

            # 4. prob_surv arrays
            fname4 = joinpath(datadir, "t_probsurv.dat")
            fname4_1 = joinpath(datadir, "t_probsurv_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname4_1, [t_array prob_surv_array])
            append_data(fname4, [t_array[start_index:end] prob_surv_array[start_index:end]])

            # 5. xsite values
            fname5 = joinpath(datadir, "t_xsiteval.dat")
            fname5_1 = joinpath(datadir, "t_xsiteval_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname5_1, [t_array x_values])
            append_data(fname5, [t_array[start_index:end] x_values[start_index:end]])

            # 6. pxsite values
            fname6 = joinpath(datadir, "t_pxsiteval.dat")
            fname6_1 = joinpath(datadir, "t_pxsiteval_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname6_1, [t_array pₓ_values])
            append_data(fname6, [t_array[start_index:end] pₓ_values[start_index:end]])

            # 7. ρₑₑ arrays
            fname7 = joinpath(datadir, "t_ρₑₑ.dat")
            fname7_1 = joinpath(datadir, "t_ρₑₑ_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname7_1, [t_array ρₑₑ_array])
            append_data(fname7, [t_array[start_index:end] ρₑₑ_array[start_index:end]])

            # 8. ρ_μμ arrays
            fname8 = joinpath(datadir, "t_ρ_μμ.dat")
            fname8_1 = joinpath(datadir, "t_ρ_μμ_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname8_1, [t_array ρ_μμ_array])
            append_data(fname8, [t_array[start_index:end] ρ_μμ_array[start_index:end]])

            # 9. ρₑμ arrays
            fname9 = joinpath(datadir, "t_ρₑμ.dat")
            fname9_1 = joinpath(datadir, "t_ρₑμ_recover.chkpt.it_end" * lpad(iteration, 6, "0") * ".dat")
            writedlm(fname9_1, [t_array ρₑμ_array])
            append_data(fname9, [t_array[start_index:end] ρₑμ_array[start_index:end]])

            # # Saving just a single value here, appended as a new row
            # fname8 = joinpath(datadir, "Im_Ω.dat")
            # open(fname8, "a") do f
            #     writedlm(f, [Im_Ω])
            # end
        end
    end
end 
