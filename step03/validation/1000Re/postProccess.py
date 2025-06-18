import numpy as np
from matplotlib import pyplot as plt
averageChange=np.loadtxt("01_average_change.txt",skiprows=5)
pressureOutput=np.loadtxt("03_P_output.txt",skiprows=8)
velocityOutput=np.loadtxt("02_U_output.txt",skiprows=8)

# Average Change Calculation
time=averageChange[:,0]
uChange=averageChange[:,1]
vChange=averageChange[:,2]
pChange=averageChange[:,3]
uMid=averageChange[:,4]

#Pressure andvelocity output
x=pressureOutput[:,0]
p=pressureOutput[:,1]
u=velocityOutput[:,1]

plt.figure(figsize=(5,5))
plt.plot(x,u,linestyle="-",linewidth=1,markersize=2,label="u",color="red")

plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('x [m]')
plt.ylabel('u [m/s]')
plt.title('ChannelFlow 1000Re 200x20 Centerline u')

# plt.xlim(0,30)
plt.xlim(0, 12)
# plt.xlim(-0.22,-0.18)
# plt.xticks(np.arange(-0.035, 0.01, step=0.005))

# plt.ylim(0,1)
plt.legend()
plt.savefig("02_FVM_stag_centerlineVelocity.svg")

plt.show()
plt.figure(figsize=(5,5))
plt.plot(x,p,linestyle="-",linewidth=1,markersize=2,label="pressure",color="red")

plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('x [m]')
plt.ylabel('Pressure [Pa]')
plt.title('ChannelFlow 1000Re 200x20 Centerline p')
plt.xlim(0, 12)

# plt.xlim(0,30)
# plt.xlim(-0.4,1)
# plt.xlim(-0.22,-0.18)
# plt.xticks(np.arange(-0.035, 0.01, step=0.005))

# plt.ylim(0,1)
plt.legend()
plt.savefig("03_FVM_stag_centerlinePressure.svg")

plt.show()


plt.figure(figsize=(5,5))
plt.semilogy(time,uChange,linestyle="-",linewidth=1,markersize=2,label="uChange",color="red")
plt.semilogy(time,vChange,linestyle="-",linewidth=1,markersize=2,label="vChange",color="green")
plt.semilogy(time,pChange,linestyle="-",linewidth=1,markersize=2,label="pChange",color="blue")

plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel('time [s]')
plt.ylabel('Average change')
plt.title('Channel Flow 1000Re / 200x20')

# plt.xlim(0,30)
# plt.xlim(-0.4,1)
# plt.xlim(-0.22,-0.18)
# plt.xticks(np.arange(-0.035, 0.01, step=0.005))

# plt.ylim(0,1)
plt.legend()
plt.savefig("01_FVM_stag_averageChange.svg")

plt.show()

