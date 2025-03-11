import numpy as np
from matplotlib import pyplot as plt
ghia=np.loadtxt("ghia5000.txt")

kod=np.loadtxt("kod5000.txt")

kodCS=np.loadtxt("kodCS.txt")
ghiaVeloc=ghia[:,1]
ghiaCoord=ghia[:,0]

codeVeloc=kod[:,1]
codeCoord=kod[:,0]

codeCSVeloc=kodCS[:,1]
codeCSCoord=kodCS[:,0]

plt.figure(figsize=(4,4))
plt.grid()
plt.plot(codeVeloc,codeCoord,linestyle="-",linewidth=1,label="Code",color="red")
plt.plot(ghiaVeloc,ghiaCoord,"o",color="black")
# plt.plot(codeCSVeloc,codeCSCoord,linestyle="--",linewidth=1,label="Code CS",color="blue")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('u[m/s]')
plt.ylabel('y [m]')
plt.title('LDC Re = 5000')
# plt.xlim(-0.387,-0.375)
# plt.xlim(-0.4,1)
# plt.xlim(-0.22,-0.18)

plt.ylim(0,1)

plt.savefig("LDC_Re_1000.png")
plt.legend()
plt.show()
