
# This file plots using the datafiles generated for Fig 7 bipolar test(Rog) 
# from julia testfile of main_Bipolar_Rog. 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

par = 2
tt= 11.160714285714286

file_path1 = f"misc/datafiles/FFI/par_{par}/tt_{tt}/t_ρ_μμ.dat"

# Read the data from the file
data1 = np.loadtxt(file_path1)


t1_array = data1[:, 0]
ρ_μμ_τ1 = data1[:, 1]


file_path2 = f"misc/datafiles/FFI/par_{par}/tt_{tt}/t_ρₑₑ.dat"

# Read the data from the file
data2 = np.loadtxt(file_path2)


t2_array = data2[:, 0]
ρₑₑ_τ2 = data2[:, 1]

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

plt.plot(t1_array, ρ_μμ_τ1, color="green", label="$\rho_{\mu \mu}$", linestyle='dotted')
ax.plot(t2_array, ρₑₑ_τ2, color="purple", label="$\rho_{\e \e}$", linestyle='dashed') 

plt.xlabel("Time $t/ \\tau_{bipolar}$")
plt.ylabel("$\\rho_{ab}$")
ax.set_title("Richers(2021) Bipolar test \n Fig. 13 Appendix A", fontsize=24) 
plt.legend(frameon=False)
ax.grid(False)
plt.legend(loc='lower center', frameon=False, bbox_to_anchor=(0.45, 0.01))
plt.savefig("Richers_biplolar_Fig_13_MF.pdf")
plt.show()