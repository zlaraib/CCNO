import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib.colorizer as mcolorizer
import matplotlib.colors as mcolors


datafile = sys.argv[1]

data = pd.read_csv(datafile, delimiter=' ', header=None, dtype=np.float64)

# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 6))
fig, axs = plt.subplots(1, 2, sharey=True, sharex=True)


# NOTE: normalizing data...
grid = data.to_numpy()
grid = grid / np.max(grid)


axs[0].set_aspect(1)
axs[0].set_title(r"Numerical (normalized)")
axs[0].set_xlabel(r"$c \ ( = \cos \theta )$")
axs[0].set_ylabel(r"$a$ (Asymmetry parameter)")

theta = np.linspace(0, np.pi, 100)
cn = np.cos(theta)
an = np.linspace(-1, 1, 100) # asymmetric param
c, a = np.meshgrid(cn, an)

omega = 1.0
mu = .5

z = .5 * ( np.imag( 
                    np.sqrt( 0j + 

        (a**2 * c**2 * mu**2) - (6 * a**2 * c * mu**2) + (9 * a**2 * mu**2) + (4 * c * mu * omega) - (12 * mu * omega) + (4 * omega**2)
                            
                           )
                 )
          )


axs[1].invert_xaxis()

axs[1].set_aspect(1)
axs[1].set_title(r"Analytical")
axs[1].set_xlabel(r"$c \ ( = \cos \theta )$")

norm = mcolors.Normalize(vmin=np.min(grid), vmax=np.max(grid))
colorizer = mcolorizer.Colorizer(norm=norm, cmap="viridis")

images = []
images.append(axs[0].imshow(grid, colorizer=colorizer, interpolation="none", extent=[1, -1, -1, 1]))
images.append(axs[1].imshow(z, colorizer=colorizer, interpolation="none", extent=[1, -1, -1, 1]))

colorbar = fig.colorbar(images[0], ax=axs, orientation="vertical", fraction=.021)

colorbar.ax.set_ylabel("Growth Rate")

# fig.suptitle(r"Homogeneous Growth Rate")
plt.show()

