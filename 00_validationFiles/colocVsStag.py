import numpy as np
from matplotlib import pyplot as plt
stag=np.loadtxt("stag1000.txt")
ghia=np.loadtxt("ghia1000.txt")
coloc=np.loadtxt("coloc1000.txt")

# kodCS=np.loadtxt("kodCS.txt")
stagVeloc=stag[:,1]
stagCoord=stag[:,0]

colocVeloc=coloc[:,1]
colocCoord=coloc[:,0]

ghiaVeloc=ghia[:,1]
ghiaCoord=ghia[:,0]


# colocCSVeloc=kodCS[:,1]
# colocCSCoord=kodCS[:,0]

plt.figure(figsize=(5,5))
plt.grid()
# plt.plot(colocVeloc,colocCoord,linestyle="-",marker='d',linewidth=0.5,markersize=2,label="11x11coloc",color="red")
plt.plot(stagVeloc,stagCoord,linestyle="-",marker='o',linewidth=0.5,markersize=2,label="11x11stag",color="blue")
plt.plot(ghiaVeloc,ghiaCoord,"o",label="ghia",color="green")

# plt.plot(colocCSVeloc,colocCSCoord,linestyle="--",linewidth=1,label="coloc CS",color="blue")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('Pressure [Pa]')
plt.ylabel('y [m]')
plt.title('LDC Re = 1000')
# plt.xlim(-0.387,-0.375)
# plt.xlim(-0.4,1)
# plt.xlim(-0.22,-0.18)
# plt.xticks(np.arange(-0.035, 0.01, step=0.005))

# plt.ylim(0,1)
plt.legend()
plt.savefig("uniformColocVsGhia_1000.png")

plt.show()