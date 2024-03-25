import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

par = 2
tt= 0.01666

file_path1 = f"misc/datafiles/FFI/par_{par}/tt_{tt}/t_ρ_μμ.dat"

# Read the data from the file
data1 = np.loadtxt(file_path1)

# Assuming data file has two columns, we separate them into x and y
x = data1[:, 0]  # First column
y = data1[:, 1]  # Second column


file_path2 = f"misc/datafiles/FFI/par_{par}/tt_{tt}/t_ρₑₑ.dat"

# Read the data from the file
data2 = np.loadtxt(file_path2)

x2 = data2[:, 0]
y2 = data2[:, 1]

file_path3 = f"misc/datafiles/FFI/par_{par}/tt_{tt}/t_ρₑμ.dat"

# Read the data from the file
data3 = np.loadtxt(file_path3)

x3 = data3[:, 0]
y3 = data3[:, 1]

mpl.rc('text', usetex=True)
mpl.rcParams['font.size'] = 20
mpl.rcParams['font.family'] = 'serif'

mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.major.pad'] = 8
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['axes.linewidth'] = 2
fig, ax = plt.subplots(figsize=(10, 8))

ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.minorticks_on()
ax.plot(x, y, color="blue",label="$\\rho_{\mu \mu}$", linestyle='solid') 
ax.plot(x2, y2, color="red",label="$\\rho_{e e}$", linestyle='solid') 
ax.plot(x3, y3, color="purple",label="$\\rho_{e \mu}$", linestyle='solid') 

# Setting y-axis scale to logarithmic
ax.set_yscale('log')

# Manually setting y-axis tick labels
ytick_labels = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
ax.set_yticks(ytick_labels)
ax.set_yticklabels(["$10^{-9}$", "$10^{-8}$", "$10^{-7}$", "$10^{-6}$", "$10^{-5}$", "$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "$1.0$"])
# Setting y-axis limits
ax.set_ylim(1e-9, 1.5)
xtick_labels_positions = [0, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014, 0.016]
xtick_labels = ["0.000", "0.002", "0.004", "0.006", "0.008", "0.010", "0.012", "0.014", "0.016"]
ax.set_xticks(xtick_labels_positions)
ax.set_xticklabels(xtick_labels)
ax.set_xlim(0, 0.01666)
# Adding title and labels
plt.xlabel("$t(s)$")
plt.ylabel("$\\rho_{ab}$")
ax.set_title("Richers (2021) Homogenous test \n Fig. 14 Appendix A", fontsize=24) 
plt.legend(frameon=False)
ax.grid(False)
plt.legend(loc='lower center', frameon=False, bbox_to_anchor=(0.8, 0.02))
plt.savefig("Richers_Homogenous_Fig_14_MF.pdf")
plt.show()
# Show grid (optional)
plt.grid(True)

# Display the plot
plt.show()
