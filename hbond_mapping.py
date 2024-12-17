import re
import matplotlib.pyplot as plt
import numpy as np

# Read the contents of the text file
with open('hbond_9R_300.txt', 'r') as file:
    data = file.read()

# Regular expression pattern to match frame numbers and atom IDs
pattern = r'Current Frame: (\d+)\n<AtomGroup .+ resid (\d+) .+ segid'
# Find all matches using the pattern
matches = np.array(re.findall(pattern, data), dtype=np.int64)
frames1 = matches[:,0]
atoms_ids1 = matches[:,1]

# Define the threshold number
threshold = 400

# Initialize empty lists to store the colors and labels
colors = []
labels = []

# Iterate through the protein_data array
for value in atoms_ids1:
    if value >= threshold:
        colors.append('red')
        labels.append('Rolipram Residues')
    else:
        colors.append('blue')
        labels.append('MMP9 Residues')

plt.scatter(frames1/100, atoms_ids1, color=colors, label=labels[0])

pattern2 = r'Current Frame: (\d+)\n<AtomGroup .+ resid (\d+) .+ segid .+ segid .+ segid'
matches2 = np.array(re.findall(pattern2, data), dtype=np.int64)
frames2 = matches2[:,0]
atoms_ids2 = matches2[:,1]

# Initialize empty lists to store the colors and labels
colors = []
labels = []

# Iterate through the protein_data array
for value in atoms_ids2:
    if value >= threshold:
        colors.append('red')
        labels.append('Rolipram Residues')
    else:
        colors.append('blue')
        labels.append('MMP9 Residues')

plt.scatter(frames2/100, atoms_ids2, color=colors, label=labels[0])
plt.xlabel("Time (ns)")
plt.ylabel('Residue Number')
plt.grid(True, linestyle='--')
plt.legend()
plt.show()
