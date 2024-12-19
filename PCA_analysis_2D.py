import numpy as np
import matplotlib.pyplot as plt

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
pc1_file = 'pc1.xvg'
pc2_file = 'pc2.xvg'

# Read the XVG files
time_pc1, pc1 = read_xvg(pc1_file, skip_header=24)  # Adjust the skip_header value as needed
time_pc2, pc2 = read_xvg(pc2_file, skip_header=24)

time_pc1 = time_pc1/1000
time_pc2 = time_pc2/1000

# Ensure the time arrays match
assert np.array_equal(time_pc1, time_pc2), "Time arrays do not match between PC1 and PC2 files."

# Create a scatter plot with a color gradient based on time
plt.figure(figsize=(10, 8))
scatter = plt.scatter(pc1, pc2, c=time_pc1, cmap='rainbow', edgecolor='k', linewidth = 0.5)

# Add colorbar and set its label size
cbar = plt.colorbar(scatter)
cbar.set_label('Time (ns)', fontsize=14)
cbar.ax.tick_params(labelsize=12)  # Set colorbar label size

# Add labels and title
plt.xlabel('PC1 (nm)', fontsize=16)
plt.ylabel('PC2 (nm)', fontsize=16)

# Set axis tick label sizes
plt.tick_params(axis='x', labelsize=14)  # Set x-axis tick label size
plt.tick_params(axis='y', labelsize=14)  # Set y-axis tick label size

# Show plot
plt.show()