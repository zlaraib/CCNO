
# This file plots using the datafiles generated for table I(Rog) 
# from julia testfile of t_p vs N_sites. 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Define the directory and file name, N and ttotal are inputs that need to be formatted into the file path
par = 24
tt= 50
τ= 0.002
file_path = f"misc/datafiles/Rog_bipolar/par_{par}/tt_{tt}/τ_{τ}/t_probsurv.dat"

# Read the data from the file
data = np.loadtxt(file_path)


t_array = data[:, 0]
probsurv_τ1 = data[:, 1]


τ= 0.05

file_path = f"misc/datafiles/Rog_bipolar/par_{par}/tt_{tt}/τ_{τ}/t_probsurv.dat"

# Read the data from the file
data = np.loadtxt(file_path)


t_array = data[:, 0]
probsurv_τ2 = data[:, 1]

τ=  0.25

file_path = f"misc/datafiles/Rog_bipolar/par_{par}/tt_{tt}/τ_{τ}/t_probsurv.dat"

# Read the data from the file
data = np.loadtxt(file_path)


t_array = data[:, 0]
probsurv_τ3 = data[:, 1]

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

plt.plot(t_array, probsurv_τ1, color="orange", label="$\delta =0.002$", linestyle='dotted')
ax.plot(t_array, probsurv_τ2, color="blue", label="$\delta =0.05$", linestyle='dashed') 
ax.plot(t_array, probsurv_τ3, color="red", label="$\delta =0.25$", linestyle='solid') 

plt.xlabel("Time $t$ [$\mu^{-1}$]")
plt.ylabel("Survival Probability $p(t)$")
ax.set_title("Rogerro(2021) Bipolar test Appendix.B Fig. 7 Direct Mean-Field Results", fontsize=24) 
plt.legend(frameon=False)
ax.grid(False)
plt.savefig("Rog_biplolar_Fig_7_MF.pdf")
plt.show()