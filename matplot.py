import numpy as np
import matplotlib.pyplot as plt

# Load data from xvg file
data = np.loadtxt('<filename>.xvg', comments=['#', '@'])

# Extract columns from data
x = data[:, 0]
y = data[:, 1]

# Plot data using Matplotlib
plt.plot(x, y)
plt.xlabel('x label')
plt.ylabel('y label')
plt.title('Title of Plot')
plt.show()
