
import numpy as np
import matplotlib.pyplot as plt

# Parameters
Nx = 180
Ny = 180

# Load data
filename = "04_VelocityField.txt"
data = np.loadtxt(filename, comments='#')

# Separate u and v
u = data[:, 0::2]  # shape (Ny, Nx)
v = data[:, 1::2]  # shape (Ny, Nx)


# Create coordinate grid (y increases upward)
x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)  # Now strictly increasing (0 to 1)
X, Y = np.meshgrid(x, y)

# Flip u and v to match the new y-direction
u = np.flipud(u)
v = np.flipud(v)



plt.figure(figsize=(8, 8))
plt.streamplot(X, Y, u, v, density=2, color='black',arrowstyle='-', linewidth=0.3,broken_streamlines=False)
plt.title("Streamlines (Velocity Field)")
plt.savefig("heyo.svg")
plt.show()

