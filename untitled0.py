import numpy as np
from matplotlib import pyplot as plt
ghia=np.loadtxt("ghia3200.txt")

kod=np.loadtxt("kod3200Cr.txt")
kodstabil=np.loadtxt("kod3200Pe.txt")
# kodCS=np.loadtxt("kodCS.txt")
ghiaVeloc=ghia[:,1]
ghiaCoord=ghia[:,0]

codeVeloc=kod[:,1]
codeCoord=kod[:,0]

codeStabilVeloc=kodstabil[:,1]
codeStabilCoord=kodstabil[:,0]

# codeCSVeloc=kodCS[:,1]
# codeCSCoord=kodCS[:,0]

plt.figure(figsize=(5,5))
plt.grid()
plt.plot(codeVeloc,codeCoord,linestyle="-",linewidth=1,label="Code",color="red")
plt.plot(codeStabilVeloc,codeStabilCoord,linestyle="--",linewidth=1,label="Code Stabil",color="green")
plt.plot(ghiaVeloc,ghiaCoord,"o",label="ghia",color="black")
# plt.plot(codeCSVeloc,codeCSCoord,linestyle="--",linewidth=1,label="Code CS",color="blue")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('u[m/s]')
plt.ylabel('y [m]')
plt.title('LDC Re = 3200')
# plt.xlim(-0.387,-0.375)
# plt.xlim(-0.4,1)
# plt.xlim(-0.22,-0.18)

plt.ylim(0,1)

plt.savefig("LDC_Re_3200.png")
plt.legend()
plt.show()
