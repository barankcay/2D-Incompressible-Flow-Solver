import matplotlib.pyplot as plt
import numpy as np
import time

plt.ion()  # interactive mode
fig, ax = plt.subplots()
line, = ax.plot([], [], 'b-', label='Center U velocity')
ax.set_xlabel('Time')
ax.set_ylabel('uMid')
ax.legend()
ax.grid(True)

while True:
    try:
        # Dosyayı oku, ilk 4 satırı atla
        data = np.loadtxt('01_average_change.txt', skiprows=4)
        
        # Eğer dosya boşsa atla
        if data.size == 0:
            continue
        
        # Tek satır bile olsa doğru şekle sok
        if data.ndim == 1:
            data = data.reshape(1, -1)

        t = data[:, 0]  # Time
        uMid = data[:, 4]  # uMid (Center U velocity)

        # Grafiği güncelle
        line.set_xdata(t)
        line.set_ydata(uMid)
        ax.relim()
        ax.autoscale_view()
        plt.draw()
        plt.pause(0.5)

    except Exception as e:
        # Veri henüz oluşmamışsa veya dosya yazılırken hata çıkarsa
        pass

    time.sleep(0.5)
