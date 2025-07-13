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
        fullPath = os.path.join(foldersDirection, file)
        if file.startswith("ghia") and file.endswith(".txt"):
            ghia.append(np.loadtxt(fullPath))
        elif file.startswith("02_U_") and file.endswith(".txt"):
            U.append(np.loadtxt(fullPath, skiprows=8))
legend=["30x30","60x60","100x100","150x150","180x180"]
color=["red","blue","green","magenta","black"]
marker=[".","o","v","+","x"]
linestyle=["solid","dotted","dashed","dashdot",(0, (5, 1))]
plt.figure(figsize=(5,5))
ghia_data = ghia[0]
ghiaY = ghia_data[:, 0]
ghiaU = ghia_data[:, 1]# Limit to 5 plots
plt.plot(ghiaU,ghiaY,"o",label="Ghia [8]",color="green")
for i in range(len(U)):

    U_data = U[i]


    u = U_data[:, 1]
    y = U_data[:, 0]

    plt.plot(u,y,linestyle=linestyle[i],linewidth=1.5,markersize=4,label=legend[i],color=color[i])
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    # indices = np.linspace(0, len(u)-1, 10, dtype=int)

# Add markers at those points
    # plt.scatter(u[indices], y[indices], marker=marker[i])  # 'b^' = blue triangle
    plt.xlabel('U [m/s]')
    plt.ylabel('y [m]')
    plt.legend()

plt.savefig("5000Re_meshIndependenceStudy.svg")
