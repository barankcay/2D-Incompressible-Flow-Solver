import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

workingDir = os.getcwd()
ghia_3200 = []
ghia_5000=[]
U_3200 = []
U_5000=[]
i=-1
# Read files
for folders in os.listdir(workingDir):

    if folders.endswith("5000") or folders.endswith("3200"):
        i=i+1
        caseFolder=os.path.join(workingDir,folders)
        for folders2 in os.listdir(caseFolder):
            if folders2.endswith("_ghia"):
                caseFolder2=os.path.join(caseFolder,folders2)
                for folders3 in os.listdir(caseFolder2):
                    if "x" in folders3:
                        caseFolder3=os.path.join(caseFolder2,folders3)
                        for textFiles in os.listdir(caseFolder3):
                            if textFiles.startswith("ghia"):
                                ghiaPath=os.path.join(caseFolder3,textFiles)
   
                            elif textFiles.startswith("02_U_") and textFiles.endswith(".txt"):
                                fullPath = os.path.join(caseFolder3,textFiles)
                                if i==0:
                                    U_3200.append(np.loadtxt(fullPath, skiprows=8))
                                if i==1:
                                    U_5000.append(np.loadtxt(fullPath, skiprows=8))
        if i==0:
            ghia_3200.append(np.loadtxt(ghiaPath))
        if i==1:
            ghia_5000.append(np.loadtxt(ghiaPath))
# Limit to 5 plots

# headline=['(a)','(b)','(c)','(d)','(e)']
# # Setup figure and custom grid
# import matplotlib as mpl  # needed for mpl.gridspec

fig = plt.figure(figsize=(23,23))
spec = gridspec.GridSpec(nrows=2, ncols=4, figure=fig, hspace=0.2, wspace=1)

axes = [
    fig.add_subplot(spec[0:1, 0:2]),  # Top left
    fig.add_subplot(spec[0:1, 2:4]),  # Top right
   
]
# U=[]

# for i in range(5):
#     U.append(U_3200[i][:,1])
#     U
legend=["30x30","60x60","100x100","150x150","180x180"]
color=["red","blue","green","magenta","black"]
marker=[".","o","v","+","x"]
linestyle=["solid","dotted","dashed","dashdot",(0, (5, 1))]
headline=['(a)','(b)']
for i in range(2):
    if i==0:
        ax=axes[i]
        ghiaY=ghia_3200[0][:,0]
        ghiaU=ghia_3200[0][:,1]
        ax.plot(ghiaU, ghiaY, "o", label="Ghia [8]", color="green",markersize='11')
        ax.plot()
        for m in range(5):
            ax.plot(U_3200[m][:,1], U_3200[m][:,0], linestyle=linestyle[m],linewidth=3,markersize=4,label=legend[m],color=color[m])
        ax.legend()
        ax.set_xlabel('u', fontsize='25')
        ax.set_ylabel('y', fontsize='25')
        ax.grid(True)
        ax.set_title(headline[i], fontsize='25')
        ax.legend(fontsize=22)
        ax.tick_params(axis='both', which='major', labelsize=20)  # you can adjust size

    else:
        ax=axes[i]
        ghiaY=ghia_5000[0][:,0]
        ghiaU=ghia_5000[0][:,1]
        ax.plot(ghiaU, ghiaY, "o", label="Ghia [8]", color="green",markersize='11')
        ax.plot()
        for m in range(5):
            ax.plot(U_5000[m][:,1], U_5000[m][:,0],linestyle=linestyle[m],linewidth=3,markersize=4,label=legend[m],color=color[m])
        ax.legend()
        ax.set_xlabel('u', fontsize='25')
        ax.set_ylabel('y', fontsize='25')
        ax.grid(True)
        ax.set_title(headline[i], fontsize='25')
        ax.legend(fontsize=22)
        ax.tick_params(axis='both', which='major', labelsize=20)  # you can adjust size
plt.savefig("gridIndependence3200_5000.svg", bbox_inches='tight', pad_inches=0.3)

            
#     ghia_data = ghia[i]
#     U_data = U[i]

#     ghiaY = ghia_data[:, 0]
#     ghiaU = ghia_data[:, 1]
#     u = U_data[:, 1]
#     y = U_data[:, 0]

#     ax = axes[i]
#     ax.plot(ghiaU, ghiaY, "o", label="Ghia [8]", color="green",markersize='11')
#     ax.plot(u, y, "-", linewidth=2.0, label="Current Study", color="blue")
#     ax.set_title(headline[i], fontsize='25')
    

#     ax.legend()
#     ax.set_xlabel('u', fontsize='25')
#     ax.set_ylabel('y', fontsize='25')
#     ax.grid(True)
    
#     # Make axis tick labels bigger
#     ax.legend(fontsize=20)

#     ax.tick_params(axis='both', which='major', labelsize=24)  # you can adjust size
#     plt.savefig("step01ghiaValidation.svg", bbox_inches='tight', pad_inches=0.3)


    
    