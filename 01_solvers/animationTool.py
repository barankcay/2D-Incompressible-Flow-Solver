import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time


# def read_csv(filename):
#     # CSV dosyasını pandas ile oku, ';' ile ayırarak
#     data = pd.read_csv(filename, header=None, delimiter=';')  # Noktalı virgülle ayır
#     data = data.astype(float)  # Verileri float türüne dönüştür
#     return data.values  # Pandas DataFrame'ini numpy array'e dönüştür

max_steps = 500


# Set the color range constant
MIN = 0
MAX = 300

for t in range(max_steps):
    filename = f'dens_t{t}.csv'
    print(f"Processing time step: {t}")  # Debug mesajı

    plt.ion()

    data = pd.read_csv(filename, header=None, delimiter=',')  # Noktalı virgülle ayır
    data = data.astype(float)  # Verileri float türüne dönüştür
    x = np.linspace(0, data.shape[0], data.shape[0])
    y = np.linspace(0, data.shape[1], data.shape[1])

        
    X, Y = np.meshgrid(x, y)
    fig, ax = plt.subplots()
    cax = ax.imshow(data, cmap='viridis', vmin=MIN, vmax=MAX)
    cbar = plt.colorbar(cax, ax=ax)  # Add colorbar
    cbar.set_label('Density', rotation=270, labelpad=20)  # Set label for colorbar


    ax.imshow(data, cmap='viridis', vmin=MIN, vmax=MAX)
    # plt.pause(0.5)




