import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

# CSV dosyasını başlıksız olarak oku
df = pd.read_csv('openFoamCenterU.csv', header=None)
averageChange=np.loadtxt("01_average_change.txt",skiprows=5)

time=averageChange[:,0]
uMid=averageChange[:,4]
plt.figure(figsize=(8,8))

# 1. sütundaki (index 0) verileri başlıksız al, ilk satırı atla
Veloc = df.iloc[1:, 0].astype(float).tolist()
Time  = df.iloc[1:, 1].astype(float).tolist()



# plt.grid()
# plt.plot(time,uMid,linestyle="-",linewidth=1,markersize=2,color="red",label="Current Study")

# plt.legend()
# Grafiği çiz
plt.plot(Time, Veloc,linestyle="-",linewidth=1,markersize=2,color="blue",label="OpenFOAM")
plt.xlabel("Time")
plt.ylabel("u")
# plt.legend()
plt.grid(True)
plt.savefig("OpenFOAMvsCode.svg")
plt.show()

