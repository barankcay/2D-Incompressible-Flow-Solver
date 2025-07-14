import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def parse_vtk(filename):
    """Parse the VTK file to extract points and velocities"""
    points = []
    velocities = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find POINTS section
    points_start = None
    for i, line in enumerate(lines):
        if line.startswith('POINTS'):
            points_start = i + 1
            num_points = int(line.split()[1])
            break
    
    # Extract points
    for line in lines[points_start:points_start+num_points]:
        coords = list(map(float, line.strip().split()))
        points.append(coords[:2])  # Only take x,y coordinates
    
    # Find VECTORS section
    vectors_start = None
    for i, line in enumerate(lines):
        if line.startswith('VECTORS velocity'):
            vectors_start = i + 1
            break
    
    # Extract velocities
    for line in lines[vectors_start:vectors_start+num_points]:
        vec = list(map(float, line.strip().split()))
        velocities.append(vec[:2])  # Only take u,v components
    
    return np.array(points), np.array(velocities)

def create_streamline_plot(points, velocities, grid_size=(179,179)):
    """Create streamline visualization from parsed data"""
    
    # Create a regular grid
    xi = np.linspace(points[:, 0].min(), points[:, 0].max(), grid_size[0])
    yi = np.linspace(points[:, 1].min(), points[:, 1].max(), grid_size[1])
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate velocities onto the regular grid
    u = griddata(points, velocities[:, 0], (xi, yi), method='cubic')
    v = griddata(points, velocities[:, 1], (xi, yi), method='cubic')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot simple black streamlines
    ax.streamplot(xi, yi, u, v,
                 color='black',
                 density=3,
                 linewidth=0.6,
                 arrowstyle='-',
                 broken_streamlines=False,
                 arrowsize=1.0)
    
    # Set equal aspect ratio
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig('result.svg')
    plt.show()

# Main execution
if __name__ == "__main__":
    # Parse the VTK file
    points, velocities = parse_vtk('04_output.txt')
    
    # Create streamline plot
    create_streamline_plot(points, velocities, grid_size=(179, 179))