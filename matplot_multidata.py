import matplotlib.pyplot as plt
import numpy as np

# Read the data from XVG files
data1 = np.loadtxt('filename1.xvg', comments=['@', '#'])
data2 = np.loadtxt('filename2.xvg', comments=['@', '#'])
data3 = np.loadtxt('filename3.xvg', comments=['@', '#'])

# Plot the data
plt.plot(data1[:, 0], data1[:, 1], label='label of your data1')
plt.plot(data2[:, 0], data2[:, 1], label='label of your data2')
plt.plot(data3[:, 0], data3[:, 1], label='label of your data3')

# Customize the plot
plt.xlabel('Time (ps)')
plt.ylabel('RMSD')
plt.title('RMSD')
plt.legend()

# Display the plot
plt.show()
