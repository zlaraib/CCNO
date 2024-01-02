using Plots
using LaTeXStrings  # Import the LaTeXStrings package

# Define constants based on the document's description
Δm2 = 7.5e-5  # Δm^2 in eV^2
E = 0.001     # Neutrino energy in GeV
θ12 = 0.5986479  # θ_12 in radians

# The transition probability function
function P_αβ(L, Δm2, E, θ)
    sin2θ = sin(2θ)^2
    Δm2L_E = 1.27 * Δm2 * L / E
    return sin2θ * sin(Δm2L_E)^2
end

# Generate data for the plot
L = 0:1:100  # Distance in km
probabilities_eμ = [P_αβ(l, Δm2, E, θ12) for l in L]
probabilities_ee = 1 .- probabilities_eμ  # Survival probability for P_ee

# Plot the data with the corrected labels and legend properties
p = plot(L, probabilities_eμ, label=L"P_{\nu_e \rightarrow \nu_{\mu}}", xlabel="L [km]", ylabel=L"P_{\alpha\beta}(L)",
         legend=:topright, legend_font=font(15), legend_background_color=RGBA(1, 1, 1, 0.5), color=:silver)
plot!(L, probabilities_ee, label=L"P_{\nu_e \rightarrow \nu_e}", color=:gold)


# Save the plot as a PDF file
savefig(p, "vacuum_oscillation_probabilities.pdf")  
