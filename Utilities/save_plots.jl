
# This function reads data from output files and saves the plots
function save_plots(τ, N_sites,L, tolerance, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array, datadir, plotdir, save_plot_flag::Bool)
    
    if save_plot_flag 
        save_plot_flag = isdir(plotdir) || mkpath(plotdir)

        # Define file paths for the data files
        fname1 = joinpath(datadir, "t_<Sz>_<Sy>_<Sx>.dat")
        fname2 = joinpath(datadir, "t_probsurv.dat")
        fname3 = joinpath(datadir, "t_xsiteval.dat")
        fname4 = joinpath(datadir, "t_pxsiteval.dat")
        fname5 = joinpath(datadir, "t_ρₑₑ.dat")
        fname6 = joinpath(datadir, "t_ρ_μμ.dat")
        fname7 = joinpath(datadir, "t_ρₑμ.dat")

        # Read the data files
        t_Sz_Sy_Sx = readdlm(fname1)
        t_probsurv = readdlm(fname2)
        t_xsiteval = readdlm(fname3)
        t_pxsiteval = readdlm(fname4)
        t_ρₑₑ = readdlm(fname5)
        t_ρ_μμ = readdlm(fname6)
        t_ρₑμ = readdlm(fname7)

        # Extract time arrays and corresponding values for plotting
        t_array = t_Sz_Sy_Sx[:, 1]  # Assuming time is the first column in all files

        # Plotting ρ_μμ vs t
        plot(t_array, t_ρ_μμ[:, 2], xlabel = "t", ylabel = "<ρ_μμ>", legend = false, 
            left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρ_μμ>_vs_t.pdf"))

        # Plotting ρ_ee vs t
        plot(t_array, t_ρₑₑ[:, 2], xlabel = "t", ylabel = "<ρₑₑ>", legend = false, 
            left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρₑₑ>_vs_t.pdf"))
        
        # Plotting ρₑμ vs t
        plot(t_array, t_ρₑμ[:, 2], xlabel = "t", ylabel = "<ρₑμ>", legend = false, 
            left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρₑμ>_vs_t.pdf"))

        # Plotting P_surv vs t
        plot(t_array, t_probsurv[:, 2], xlabel = "t", ylabel = "Survival Probability p(t)",
            legend = false, left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm)
        savefig(joinpath(plotdir, "Survival_probability_vs_t_for_N_sites$(N_sites).pdf"))

        # Plotting Sz vs t
        plot(t_array, t_Sz_Sy_Sx[:, 2], xlabel = "t", ylabel = "<Sz>", legend = false,
            left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 20mm, margin = 10mm,
            ylims = (minimum(t_Sz_Sy_Sx[:, 2]), maximum(t_Sz_Sy_Sx[:, 2]) + 0.1 * abs(maximum(t_Sz_Sy_Sx[:, 2]) - minimum(t_Sz_Sy_Sx[:, 2]))))
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<Sz>_vs_t_for_N_sites$(N_sites).pdf"))

        # Plotting Sy vs t
        plot(t_array, t_Sz_Sy_Sx[:, 3], xlabel = "t", ylabel = "<Sy>", legend = false,
            left_margin = 40mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<Sy>_vs_t.pdf"))

        # Plotting Sx vs t
        plot(t_array, t_Sz_Sy_Sx[:, 4], xlabel = "t", ylabel = "<Sx>", legend = false,
            left_margin = 40mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<Sx>_vs_t.pdf"))

        # Plotting particles positional evolution
        plot(title = "Position Evolution for $N_sites particles", ylabel = "Time (t)", xlabel = "Position (x)")
        t_array = t_xsiteval[:, 1]  # All rows, first column
        x_values = t_xsiteval[:, 2:end]  # All rows, all columns except the first
        xlims = (0, L)
        for site in 1:N_sites
            site_positions = x_values[:, site]  # Get all positions for a specific site
            # Convert each entry to a floating-point number
            site_positions = [
                parse(Float64, replace(strip(position, ['[', ']', ',']), "," => "")) for position in site_positions
            ]
            plot!(site_positions,t_array, label = "Site $site"
            ,left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm)
        end

        savefig(joinpath(plotdir,"Particles position(x) evolution for N_sites$(N_sites).pdf"))

        # Plotting particles momentum evolution
        plot(title = "Particle Momentum Evolution", ylabel = "Time (t)", xlabel = "Momentum in x direction (pₓ)")
        t_array = t_xsiteval[:, 1]  # All rows, first column
        px_values = t_pxsiteval[:, 2:end]  # All rows, all columns except the first
        xlims = (0, L)
        for site in 1:N_sites
            site_momentum = px_values[:, site]
            site_momentum = [
                parse(Float64, replace(strip(momentum, ['[', ']', ',']), "," => "")) for momentum in site_momentum
            ]
            plot!(site_momentum,t_array, label = "Site $site",
            left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm)
        end
        savefig(joinpath(plotdir, "Particles_momentum_pₓ_evolution_for_N_sites$(N_sites).pdf"))
    
    end
end
