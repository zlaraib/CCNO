
# This file plots using the datafiles generated for Fig 7 bipolar test(Rog) 
# from julia testfile of main_Bipolar_Rog. 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

par = 24
tt= 50
τ= 0.002
file_path1 = f"misc/datafiles/Rog_bipolar/par_{par}/tt_{tt}/τ_{τ}/t_probsurv.dat"

# Read the data from the file
data1 = np.loadtxt(file_path1)


t1_array = data1[:, 0]
probsurv_τ1 = data1[:, 1]


τ= 0.05

file_path2 = f"misc/datafiles/Rog_bipolar/par_{par}/tt_{tt}/τ_{τ}/t_probsurv.dat"

# Read the data from the file
data2 = np.loadtxt(file_path2)


t2_array = data2[:, 0]
probsurv_τ2 = data2[:, 1]

τ=  0.25

file_path3 = f"misc/datafiles/Rog_bipolar/par_{par}/tt_{tt}/τ_{τ}/t_probsurv.dat"

# Read the data from the file
data3 = np.loadtxt(file_path3)


t3_array = data3[:, 0]
probsurv_τ3 = data3[:, 1]

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

plt.plot(t1_array, probsurv_τ1, color="orange", label="$\delta =0.002 \mu^{-1}$", linestyle='dotted')
ax.plot(t2_array, probsurv_τ2, color="blue", label="$\delta =0.05 \mu^{-1}$", linestyle='dashed') 
ax.plot(t3_array, probsurv_τ3, color="red", label="$\delta =0.25 \mu^{-1}$", linestyle='solid') 

plt.xlabel("Time $t$ [$\mu^{-1}$]")
plt.ylabel("Survival Probability $p(t)$")
ax.set_title("Rogerro(2021) Bipolar test Fig. 7 \n Appendix B : Direct Mean-Field Results", fontsize=24) 
plt.legend(frameon=False)
ax.grid(False)
plt.legend(loc='lower center', frameon=False, bbox_to_anchor=(0.45, 0.01))
plt.savefig("Rog_biplolar_Fig_7_MF.pdf")
plt.show()