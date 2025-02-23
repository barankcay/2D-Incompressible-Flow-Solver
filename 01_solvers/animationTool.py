import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re

# Set the color range constant
MIN = -1
MAX = 1

folder_path = "."  # Current folder, change if needed

# Function to extract the numeric part of the filename for natural sorting
def extract_number(filename):
    match = re.search(r'\d+', filename)  # Find the first number in the filename
    return int(match.group()) if match else float('inf')

# Get all CSV files and sort them numerically
csv_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.csv')], key=extract_number)

# Enable interactive mode for real-time plotting
plt.ion()

for filename in csv_files:
    print(f"Processing file: {filename}")  # Debug message

    # Read CSV file and convert it to a numpy array
    data = pd.read_csv(filename, header=None, delimiter=',').astype(float)

    # Create the figure and plot
    fig, ax = plt.subplots()
    cax = ax.imshow(data, cmap='viridis', vmin=MIN, vmax=MAX)
    cbar = plt.colorbar(cax, ax=ax)
    cbar.set_label('Density', rotation=270, labelpad=20)

    plt.pause(0.001)  # Pause to update the plot
    plt.close(fig)  # Close the figure to prevent overlapping
