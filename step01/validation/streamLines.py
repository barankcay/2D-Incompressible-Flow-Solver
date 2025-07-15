import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import os
from matplotlib import gridspec

# Find all 04_*.txt files in subdirectories
vtk_files = []
for root, dirs, files in os.walk('.'):
    for file in files:
        if file.startswith('04_') and file.endswith('.txt'):
            vtk_files.append(os.path.join(root, file))

# Limit to 5 files if more exist
vtk_files = vtk_files[:5]
if not vtk_files:
    raise FileNotFoundError("No 04_*.txt files found in subdirectories")

# Setup figure and custom grid
plt.rcParams.update({'font.size': 18})  # Base font size
fig = plt.figure(figsize=(23, 28))
spec = gridspec.GridSpec(nrows=3, ncols=4, figure=fig, 
                        height_ratios=[1, 1, 1], 
                        hspace=0.3, wspace=0.4)

axes = [
    fig.add_subplot(spec[0, 0:2]),  # Top left
    fig.add_subplot(spec[0, 2:4]),  # Top right
    fig.add_subplot(spec[1, 0:2]),  # Middle left
    fig.add_subplot(spec[1, 2:4]),  # Middle right
    fig.add_subplot(spec[2, 1:3]),  # Bottom center
]

headline = ['(a)', '(b)', '(c)', '(d)', '(e)']

for i, vtk_file in enumerate(vtk_files):
    # Parse VTK file
    points, velocities = [], []
    with open(vtk_file, 'r') as f:
        lines = f.readlines()
    
    # Extract points
    points_start = next(j for j, line in enumerate(lines) 
                       if line.startswith('POINTS')) + 1
    num_points = int(lines[points_start-1].split()[1])
    points = np.array([list(map(float, line.strip().split()[:2])) 
                      for line in lines[points_start:points_start+num_points]])
    
    # Extract velocities
    vectors_start = next(j for j, line in enumerate(lines) 
                        if line.startswith('VECTORS velocity')) + 1
    velocities = np.array([list(map(float, line.strip().split()[:2])) 
                         for line in lines[vectors_start:vectors_start+num_points]])
    
    # Create regular grid and interpolate
    xi = np.linspace(0, 1, 179)  # Explicit 0-1 range
    yi = np.linspace(0, 1, 179)  # Explicit 0-1 range
    xi, yi = np.meshgrid(xi, yi)
    u = griddata(points, velocities[:,0], (xi, yi), method='cubic')
    v = griddata(points, velocities[:,1], (xi, yi), method='cubic')
    
    # Plot streamlines
    ax = axes[i]
    ax.streamplot(xi, yi, u, v, color='black', density=1.7, 
                 linewidth=0.8, arrowstyle='-', arrowsize=1.2,broken_streamlines=False)
    
    # Customize axes
    ax.set_title(f'{headline[i]}', fontsize=22, pad=15)
    ax.set_aspect('equal')
    
    # Set axis ticks from 0 to 1 with custom formatting
    ax.set_xticks(np.linspace(0, 1, 5))
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.tick_params(axis='both', which='major', labelsize=20, length=6, width=1.5)
    
    # Add axis labels with increased font size
    ax.set_xlabel('x', fontsize=22, labelpad=10)
    ax.set_ylabel('y', fontsize=22, labelpad=10)

    # Make spines (axis lines) thicker
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

# Save and show results
plt.tight_layout()
plt.savefig('streamline_comparison.svg', bbox_inches='tight', dpi=300)
plt.show()