import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

workingDir = os.getcwd()
uMid=[]

# Read files
for folders in os.listdir(workingDir):
    foldersDirection = os.path.join(workingDir, folders)
    if not os.path.isdir(foldersDirection):
        continue
    for file in os.listdir(foldersDirection):
        fullPath = os.path.join(foldersDirection, file)
        if file.startswith("01_average") and file.endswith(".txt"):
            uMid.append(np.loadtxt(fullPath,skiprows=5))


# Limit to 5 plots
num_plots = len(uMid)
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
    averageChange=uMid[i]
    
    time=averageChange[:,0]
    u=averageChange[:,4]

    ax = axes[i]
    ax.plot(time,u, "-", linewidth=2.0, color="blue")
    ax.set_title(headline[i], fontsize='25')
    

    # ax.legend()
    ax.set_xlabel('Time', fontsize='25')
    ax.set_ylabel('u', fontsize='25')
    ax.grid(True)
    
    # Make axis tick labels bigger
    # ax.legend(fontsize=20)

    ax.tick_params(axis='both', which='major', labelsize=24)  # you can adjust size
    plt.savefig("step01midVelocChange.svg", bbox_inches='tight', pad_inches=0.3)

    
    