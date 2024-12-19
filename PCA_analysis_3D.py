import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function to read XVG file with skip header option
def read_xvg(file_path, skip_header):
    data = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines[skip_header:]:
            if '&' in line:
                break
            else:
                data.append([float(val) for val in line.split()])
    data = np.array(data)
    return data[:, 0], data[:, 1]

# File paths for PC1 and PC2 XVG files
pc1_file = '47-pc1.xvg'
pc2_file = '47-pc2.xvg'

# Read the XVG files
time_pc1, pc1 = read_xvg(pc1_file, skip_header=24)  # Adjust the skip_header value as needed
time_pc2, pc2 = read_xvg(pc2_file, skip_header=24)

time_pc1 = time_pc1/1000
time_pc2 = time_pc2/1000

# Ensure the time arrays match
assert np.array_equal(time_pc1, time_pc2), "Time arrays do not match between PC1 and PC2 files."

# Create a 3D scatter plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Plotting the scatter with time as the Z-axis
scatter = ax.scatter(pc1, pc2, time_pc1, c=time_pc1, cmap='rainbow', edgecolor='k', linewidth = 0.5)

# Add colorbar
cbar = plt.colorbar(scatter)
cbar.set_label('Time (ps)', fontsize=14)
cbar.ax.tick_params(labelsize=12)

# Set labels
ax.set_xlabel('PC1 (nm)', fontsize=14)
ax.set_ylabel('PC2 (nm)', fontsize=14)
ax.set_zlabel('Time (ps)', fontsize=14)

# Set tick label sizes
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='z', labelsize=12)

# Save the plot as an SVG file
plt.savefig('Figures/PCA_time gradient-3D.svg', format='svg')

# Show plot
plt.show()