import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

workingDir = os.getcwd()
ghia = []
U = []

# Read files
for folders in os.listdir(workingDir):
    foldersDirection = os.path.join(workingDir, folders)
    if not os.path.isdir(foldersDirection):
        continue
    for file in os.listdir(foldersDirection):

        if file.startswith("ghia") and file.endswith(".txt"):
            print(file)
            fullPath = os.path.join(foldersDirection, file)

            ghia.append(np.loadtxt(fullPath))
        elif file.startswith("02_U_") and file.endswith(".txt"):
            fullPath = os.path.join(foldersDirection, file)
            U.append(np.loadtxt(fullPath, skiprows=8))

# Limit to 5 plots
num_plots = min(5, len(ghia), len(U))
headline=['(a)','(b)','(c)','(d)','(e)']
# Setup figure and custom grid
import matplotlib as mpl  # needed for mpl.gridspec
fig = plt.figure(figsize=(23,28))
spec = gridspec.GridSpec(nrows=3, ncols=4, figure=fig, height_ratios=[1, 1, 1], hspace=0.2, wspace=1)

axes = [
    fig.add_subplot(spec[0, 0:2]),  # Top left
    fig.add_subplot(spec[0, 2:4]),  # Top right
    fig.add_subplot(spec[1, 0:2]),  # Middle left
    fig.add_subplot(spec[1, 2:4]),  # Middle right
    fig.add_subplot(spec[2, 1:3]),  # Bottom row, spanning both columns
]

for i in range(num_plots):
    ghia_data = ghia[i]
    U_data = U[i]

    ghiaY = ghia_data[:, 0]
    ghiaU = ghia_data[:, 1]
    u = U_data[:, 1]
    y = U_data[:, 0]

    ax = axes[i]
    ax.plot(ghiaU, ghiaY, "o", label="Ghia [8]", color="green",markersize='11')
    ax.plot(u, y, "-", linewidth=2.0, label="Current Study", color="blue")
    ax.set_title(headline[i], fontsize='25')
    

    ax.legend()
    ax.set_xlabel('u', fontsize='25')
    ax.set_ylabel('y', fontsize='25')
    ax.grid(True)
    
    # Make axis tick labels bigger
    ax.legend(fontsize=20)

    ax.tick_params(axis='both', which='major', labelsize=24)  # you can adjust size
    plt.savefig("step01ghiaValidation.svg", bbox_inches='tight', pad_inches=0.3)


    
    