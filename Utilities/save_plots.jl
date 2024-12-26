
# This function reads data from output files and saves the plots
function save_plots(τ, N_sites,L, t_array, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array, datadir, plotdir, save_plot_flag::Bool)
    
        save_plot_flag = isdir(plotdir) || mkpath(plotdir)

        # Plotting ρ_μμ vs t
        plot(t_array, ρ_μμ_array, xlabel = "t", ylabel = "<ρ_μμ>", legend = false, 
            left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρ_μμ>_vs_t.pdf"))

        # Plotting ρ_ee vs t
        plot(t_array, ρₑₑ_array, xlabel = "t", ylabel = "<ρₑₑ>", legend = false, 
            left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρₑₑ>_vs_t.pdf"))
        
        # Plotting ρₑμ vs t
        plot(t_array, ρₑμ_array, xlabel = "t", ylabel = "<ρₑμ>", legend = false, 
            left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρₑμ>_vs_t.pdf"))

        # Plotting P_surv vs t
        plot(t_array, prob_surv_array, xlabel = "t", ylabel = "Survival Probability p(t)",
            legend = false, left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm)
        savefig(joinpath(plotdir, "Survival_probability_vs_t_for_N_sites$(N_sites).pdf"))

        # Plotting Sz vs t
        plot(t_array, Sz_array, xlabel = "t", ylabel = "<Sz>", legend = false,
            left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 20mm, margin = 10mm,
            # ylims = (minimum(t_Sz[:, 2]), maximum(t_Sz[:, 2]) + 0.1 * abs(maximum(t_Sz[:, 2]) - minimum(t_Sz[:, 2])))
            )
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<Sz>_vs_t_for_N_sites$(N_sites).pdf"))

        # Plotting Sy vs t
        plot(t_array, Sy_array, xlabel = "t", ylabel = "<Sy>", legend = false,
            left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 20mm, margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<Sy>_vs_t.pdf"))

        # Plotting Sx vs t
        plot(t_array, Sx_array, xlabel = "t", ylabel = "<Sx>", legend = false,
            left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 20mm, margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<Sx>_vs_t.pdf"))

        # Plotting particles positional evolution
        plot(title = "Position Evolution for $N_sites particles", ylabel = "Time (t)", xlabel = "Position (x)")
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
        xlims = (0, L)
        for site in 1:N_sites
            site_momentum = pₓ_values[:, site]
            site_momentum = [
                parse(Float64, replace(strip(momentum, ['[', ']', ',']), "," => "")) for momentum in site_momentum
            ]
            plot!(site_momentum,t_array, label = "Site $site",
            left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm)
        end
        savefig(joinpath(plotdir, "Particles_momentum_pₓ_evolution_for_N_sites$(N_sites).pdf"))

end
