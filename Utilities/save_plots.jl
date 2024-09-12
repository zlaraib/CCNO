

# This file saves the plots for the different output datafiles 

function save_plots(τ, N_sites, ttotal,Sz_array, Sy_array, Sx_array, prob_surv_array, x_values, pₓ_values, ρₑₑ_array,ρ_μμ_array, ρₑμ_array, plotdir, save_plot_flag::Bool)

    if save_plot_flag 
        save_plot_flag = isdir(plotdir) || mkpath(plotdir)
        # Plotting ρ_μμ vs t
        plot(0.0:τ:τ*(length(ρ_μμ_array)-1), ρ_μμ_array, xlabel = "t", ylabel = "<ρ_μμ>", legend = false, 
        left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρ_μμ>_vs_t.pdf"))

        # Plotting ρ_ee vs t
        plot(0.0:τ:τ*(length(ρₑₑ_array)-1), ρₑₑ_array, xlabel = "t", ylabel = "<ρₑₑ>", legend = false, 
        left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρₑₑ>_vs_t.pdf"))
        
        # Plotting ρₑμ vs t
        plot(0.0:τ:τ*(length(ρₑμ_array)-1), ρₑμ_array, xlabel = "t", ylabel = "<ρₑμ>", legend = false, 
        left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm) 
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<ρₑμ>_vs_t.pdf"))

        # Plotting P_surv vs t
        plot(0.0:τ:τ*(length(prob_surv_array)-1), prob_surv_array, xlabel = "t", ylabel = "Survival Probabillity p(t)",
        legend = false, left_margin = 20mm, right_margin = 10mm, top_margin = 5mm, bottom_margin = 10mm)
        savefig(joinpath(plotdir,"Survival probability vs t for N_sites$(N_sites).pdf"))

        # Plotting Sz vs t
        plot(0.0:τ:τ*(length(Sz_array)-1), Sz_array, xlabel = "t", ylabel = "<Sz>", legend = false,
            left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 20mm,margin= 10mm,
            ylims = (minimum(Sz_array), maximum(Sz_array) + 0.1 * abs(maximum(Sz_array) - minimum(Sz_array))) )
        # Save the plot as a PDF file
        savefig(joinpath(plotdir, "<Sz> vs t for N_sites$(N_sites).pdf"))

        # Plotting Sy vs t
        plot(0.0:τ:τ*(length(Sy_array)-1), Sy_array, xlabel = "t", ylabel = "<Sy>", legend = false,
        left_margin = 40mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin= 10mm) 
        #Save the plot as a PDF file
        savefig(joinpath(plotdir,"<Sy> vs t.pdf"))

        # Plotting Sx vs t
        plot(0.0:τ:τ*(length(Sx_array)-1), Sx_array, xlabel = "t", ylabel = "<Sx>", legend = false,
        left_margin = 40mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin= 10mm) 
        #Save the plot as a PDF file
        savefig(joinpath(plotdir,"<Sx> vs t.pdf"))

        # Plotting particles positional evolution
        plot(title="Position Evolution for $N_sites particles", xlabel= "Position (x)",ylabel="Time(s)")
        for site in 1:N_sites
            site_positions = [(x_values[t][site]) for t in 1:length(x_values)]
            plot!(site_positions, 0.0:τ:ttotal, label="Site $site",
            left_margin = 25mm, right_margin = 5mm, top_margin = 5mm, bottom_margin = 10mm, margin= 10mm)
        end
        savefig(joinpath(plotdir,"Particles position(x) evolution for N_sites$(N_sites).pdf"))
        
        # Plotting particles momentum evolution
        plot(title="Particle Momentum Evolution", xlabel= "Momentum in x direction(pₓ)",ylabel="Time")
        for site in 1:N_sites
            site_momentum = [(pₓ_values[t][site]) for t in 1:length(pₓ_values)]
            plot!(site_momentum, 0.0:τ:ttotal, label="Site $site",left_margin = 25mm, right_margin = 5mm, 
            top_margin = 5mm, bottom_margin = 10mm)
        end
        savefig(joinpath(plotdir,"Particles momentum(pₓ) evolution for N_sites$(N_sites).pdf"))
        
    end
end 


