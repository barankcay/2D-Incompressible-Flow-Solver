import numpy as np
from matplotlib import pyplot as plt
averageChange=np.loadtxt("01_average_change.txt",skiprows=5)
ghia=np.loadtxt("ghia100.txt")

pressureOutput=np.loadtxt("03_P_output.txt",skiprows=8)
velocityOutput=np.loadtxt("02_U_output.txt",skiprows=8)

# Average Change Calculation
time=averageChange[:,0]
uChange=averageChange[:,1]
vChange=averageChange[:,2]
pChange=averageChange[:,3]
uMid=averageChange[:,4]

#Pressure andvelocity output
y=pressureOutput[:,0]
p=pressureOutput[:,1]
u=velocityOutput[:,1]

ghiaY=ghia[:,0]
ghiaU=ghia[:,1]

plt.figure(figsize=(5,5))
plt.grid()
# plt.plot(colocVeloc,colocCoord,linestyle="-",marker='d',linewidth=0.5,markersize=2,label="11x11coloc",color="red")
# plt.plot(stagVeloc,stagCoord,linestyle="-",marker='o',linewidth=0.5,markersize=2,label="11x11stag",color="blue")
plt.plot(ghiaU,ghiaY,"o",label="Ghia [10]",color="green")
plt.plot(u,y,linestyle="-",linewidth=0.5,markersize=2,label="fdmStag",color="blue")

# plt.plot(colocCSVeloc,colocCSCoord,linestyle="--",linewidth=1,label="coloc CS",color="blue")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('U [m/s]')
plt.ylabel('y [m]')
plt.title('LDC 100Re / 180x180')
plt.legend()
plt.savefig("01_FDM_StaggeredVsGhia_100.svg")
plt.show()



plt.figure(figsize=(5,5))
plt.semilogy(time,uChange,linestyle="-",linewidth=1,markersize=2,label="uChange",color="red")
plt.semilogy(time,vChange,linestyle="-",linewidth=1,markersize=2,label="vChange",color="green")
plt.semilogy(time,pChange,linestyle="-",linewidth=1,markersize=2,label="pChange",color="blue")

plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('time [s]')
plt.ylabel('Average change')
plt.title('LDC 100Re / 180x180')

# plt.xlim(0,30)
# plt.xlim(-0.4,1)
# plt.xlim(-0.22,-0.18)
# plt.xticks(np.arange(-0.035, 0.01, step=0.005))

# plt.ylim(0,1)
plt.legend()
plt.savefig("02_averageChange.svg")

plt.show()
plt.figure(figsize=(6,6))

plt.xlabel('time [s]')
plt.ylabel('u [m/s]')
plt.title('LDC 100Re / 180x180')

plt.grid()
plt.plot(time,uMid,linestyle="-",linewidth=1,markersize=2,label="uMid",color="black")
plt.legend()
plt.savefig("03_uMid.svg")
