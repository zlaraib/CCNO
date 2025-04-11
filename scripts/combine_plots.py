import os
import numpy as np
import matplotlib.pyplot as plt
import sys

# Get data directory from command-line argument
if len(sys.argv) < 2:
    print("Usage: python combine_plots.py <data_directory>")
    sys.exit(1)

combined_data_dir = sys.argv[1]  # Directory containing combined data files
parent_dir = os.path.dirname(combined_data_dir)  # Get parent directory
combined_plot_dir = os.path.join(parent_dir, "combined_plots")  # Save plots in "combined_plots/"
os.makedirs(combined_plot_dir, exist_ok=True)  # Ensure directory exists

# Function to read data and extract time + multi-site values
def read_data(filename):
    file_path = os.path.join(combined_data_dir, filename)
    if not os.path.exists(file_path):
        print(f"Warning: {filename} not found, skipping.")
        return None, None
    data = np.genfromtxt(file_path)
    if data.shape[1] < 3:  # Ensure file has at least 3 columns
        print(f"Error: {filename} does not have expected shape (200, 3+). Found {data.shape}.")
        return None, None
    return data[:, 0], data[:, 1:]  # Extract first column as time, rest as values

# Load data
t_array, Sz_array = read_data("t_<Sz>.dat")
_, Sy_array = read_data("t_<Sy>.dat")
_, Sx_array = read_data("t_<Sx>.dat")
_, prob_surv_array = read_data("t_probsurv.dat")
_, x_values = read_data("t_xsiteval.dat")
_, pₓ_values = read_data("t_pxsiteval.dat")
_, ρₑₑ_array = read_data("t_ρₑₑ.dat")
_, ρ_μμ_array = read_data("t_ρ_μμ.dat")
_, ρₑμ_array = read_data("t_ρₑμ.dat")

# Ensure t_array was loaded properly
if t_array is None:
    print("Error: t_<Sz>.dat not found or unreadable. Ensure the file exists in combined_data.")
    sys.exit(1)

# Compute domain-averaged values (mean of absolute values for each row)
def domain_avg(array):
    if array is None:
        return None
    half_size = array.shape[1] // 2  # Take only half of the columns
    return np.mean(np.abs(array[:, :half_size]), axis=1)  # Compute mean over the first half of columns

Sz_array_domain_avgd = domain_avg(Sz_array)
Sy_array_domain_avgd = domain_avg(Sy_array)
Sx_array_domain_avgd = domain_avg(Sx_array)
prob_surv_array_domain_avgd = domain_avg(prob_surv_array)
ρₑₑ_array_domain_avgd = domain_avg(ρₑₑ_array)
ρ_μμ_array_domain_avgd = domain_avg(ρ_μμ_array)
ρₑμ_array_domain_avgd = domain_avg(ρₑμ_array)

# Function to plot and save figures
def plot_and_save(x, y, xlabel, ylabel, filename):
    if y is None:
        print(f"Skipping plot: {filename} (Data not available)")
        return
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker="o", linestyle="-", markersize=2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.savefig(os.path.join(combined_plot_dir, filename))
    plt.close()

# Generate and save plots for original site-wise data
plot_and_save(t_array, ρ_μμ_array[:, 0], "t", "<ρ_μμ>", "rho_mumu_first_site_vs_t.pdf")
plot_and_save(t_array, ρₑₑ_array[:, 0], "t", "<ρₑₑ>", "rho_ee_first_site_vs_t.pdf")
plot_and_save(t_array, ρₑμ_array[:, 0], "t", "<ρₑμ>", "rho_emu_first_site_vs_t.pdf")
plot_and_save(t_array, prob_surv_array[:, 0], "t", "Survival Probability p(t)", "survival_probability_first_site_vs_t.pdf")
plot_and_save(t_array, Sz_array[:, 0], "t", "<Sz>", "Sz_first_site_vs_t.pdf")
plot_and_save(t_array, Sy_array[:, 0], "t", "<Sy>", "Sy_first_site_vs_t.pdf")
plot_and_save(t_array, Sx_array[:, 0], "t", "<Sx>", "Sx_first_site_vs_t.pdf")

# Generate and save domain-averaged plots
plot_and_save(t_array, ρ_μμ_array_domain_avgd, "t", "Domain-Averaged <ρ_μμ>", "rho_mumu_domain_avg_vs_t.pdf")
plot_and_save(t_array, ρₑₑ_array_domain_avgd, "t", "Domain-Averaged <ρₑₑ>", "rho_ee_domain_avg_vs_t.pdf")
plot_and_save(t_array, ρₑμ_array_domain_avgd, "t", "Domain-Averaged <ρₑμ>", "rho_emu_domain_avg_vs_t.pdf")
plot_and_save(t_array, prob_surv_array_domain_avgd, "t", "Domain-Averaged Survival Probability", "survival_probability_domain_avg_vs_t.pdf")
plot_and_save(t_array, Sz_array_domain_avgd, "t", "Domain-Averaged <Sz>", "Sz_domain_avg_vs_t.pdf")
plot_and_save(t_array, Sy_array_domain_avgd, "t", "Domain-Averaged <Sy>", "Sy_domain_avg_vs_t.pdf")
plot_and_save(t_array, Sx_array_domain_avgd, "t", "Domain-Averaged <Sx>", "Sx_domain_avg_vs_t.pdf")

# Plot particles' position evolution
if x_values is not None:
    plt.figure(figsize=(8, 6))
    plt.title("Position Evolution")
    plt.xlabel("Time (t)")
    plt.ylabel("Position (x)")
    plt.grid(True)
    for i in range(x_values.shape[1]):  # Loop through sites
        plt.plot(t_array, x_values[:, i], label=f"Site {i+1}")  
    plt.legend()
    plt.savefig(os.path.join(combined_plot_dir, "particle_position_evolution.pdf"))
    plt.close()

# Plot particles' momentum evolution
if pₓ_values is not None:
    plt.figure(figsize=(8, 6))
    plt.title("Momentum Evolution")
    plt.xlabel("Time (t)")
    plt.ylabel("Momentum (pₓ)")
    plt.grid(True)
    for i in range(pₓ_values.shape[1]):  # Loop through sites
        plt.plot(t_array, pₓ_values[:, i], label=f"Site {i+1}")  
    plt.legend()
    plt.savefig(os.path.join(combined_plot_dir, "particle_momentum_evolution.pdf"))
    plt.close()

print(f"Plots saved in: {combined_plot_dir}")
