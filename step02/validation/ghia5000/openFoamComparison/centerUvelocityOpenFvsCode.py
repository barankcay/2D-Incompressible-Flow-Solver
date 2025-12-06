import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

# CSV dosyasını başlıksız olarak oku
df = pd.read_csv('5000openFoamCenterU.csv', header=None)
dff= pd.read_csv('3200openFoamCenterU.csv', header=None)
print(df)
averageChange=np.loadtxt("01_average_change.txt",skiprows=5)

time=averageChange[:,0]
uMid=averageChange[:,4]
plt.figure(figsize=(8,8))

# 1. sütundaki (index 0) verileri başlıksız al, ilk satırı atla
Veloc_5000 = df.iloc[1:, 0].astype(float).tolist()
Time_5000  = df.iloc[1:, 1].astype(float).tolist()

Veloc_3200 = dff.iloc[1:, 1].astype(float).tolist()
Time_3200  = dff.iloc[1:, 0].astype(float).tolist()

fig = plt.figure(figsize=(23,23))
spec = gridspec.GridSpec(nrows=2, ncols=4, figure=fig, hspace=0.2, wspace=1)

axes = [
    fig.add_subplot(spec[0:1, 0:2]),  # Top left
    fig.add_subplot(spec[0:1, 2:4]),  # Top right
   
]

Time=[]

Time.append(Time_3200)
Time.append(Time_5000)

Veloc=[]

Veloc.append(Veloc_3200)
Veloc.append(Veloc_5000)
headline=['(a)','(b)']

# plt.grid()
# plt.plot(time,uMid,linestyle="-",linewidth=1,markersize=2,color="red",label="Current Study")

# plt.legend()
# Grafiği çiz
for i in range(2):


    ax = axes[i]
    ax.plot(Time[i],Veloc[i], "-", linewidth=2.0, color="blue")
    ax.set_title(headline[i], fontsize='25')
    

    # ax.legend()
    ax.set_xlabel('Time', fontsize='25')
    ax.set_ylabel('u', fontsize='25')
    ax.grid(True)
    
    # Make axis tick labels bigger
    # ax.legend(fontsize=20)

    ax.tick_params(axis='both', which='major', labelsize=24)  # you can adjust size
    plt.savefig("OpenFOAMvsCode.png", bbox_inches='tight', pad_inches=0.3)


