import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import time
import sys

# Define the file path
Poll_file_path = '/home/sy432/rds/rds-ye_shutong-xcywAxU6Kd0/Pol_model/results/fixedcis/loc_RNAPolII_postPause.txt'
Ser7_file_path = '/home/sy432/rds/rds-ye_shutong-xcywAxU6Kd0/Pol_model/results/fixedcis/500sim_lines.txt'

# Get the directory from the file path to save the PDF in the same location
output_dir = '/home/sy432/rds/rds-ye_shutong-xcywAxU6Kd0/Pol_model/dataprocess/fixedcis/'
output_filename = 'count_ratio_y=1.pdf'
output_path = os.path.join(output_dir, output_filename)

# Initialize lists to store positions for each simulation
Poll_positions = []
Ser7_positions = []

print("Reading Pol II positions...")
sys.stdout.flush()
start_time = time.time()
with open(Poll_file_path) as f:
    for i, line in enumerate(f):
        if i % 1000000 == 0:  # Print progress every million lines
            print(f"Processed {i} lines...")
            sys.stdout.flush()
        if line[0] == '>':
            continue
        # Skip empty lines
        if not line.strip():
            continue
        splitLine = line.strip().split('\t')
        # Get position from splitLine[4]
        pos = int(splitLine[4])  # Position is in the 5th column
        Poll_positions.append(pos)
print(f"Finished reading Pol II positions in {time.time() - start_time:.2f} seconds")
sys.stdout.flush()

print("\nReading Ser7 positions...")
sys.stdout.flush()
start_time = time.time()
with open(Ser7_file_path) as f:
    for i, line in enumerate(f):
        if i % 1000000 == 0:  # Print progress every million lines
            print(f"Processed {i} lines...")
            sys.stdout.flush()
        if line[0] == '>':
            continue
        # Skip empty lines
        if not line.strip():
            continue
        splitLine = line.strip().split('\t')
        pos = int(splitLine[4])  # Position is in the 5th column
        Ser7_positions.append(pos)
print(f"Finished reading Ser7 positions in {time.time() - start_time:.2f} seconds")
sys.stdout.flush()

print("\nConverting to numpy arrays...")
# Convert to numpy array
allPoll_positions = np.array(Poll_positions)
allSer7_positions = np.array(Ser7_positions)

# Save the raw position arrays
print("\nSaving raw position arrays...")
np.save(os.path.join(output_dir, 'Poll_positions.npy'), allPoll_positions)
np.save(os.path.join(output_dir, 'Ser7_positions.npy'), allSer7_positions)
print("Raw position arrays saved successfully!")

print("Creating plots...")
# Create figure with ten subplots (5 original + 5 zoomed)
fig, ((ax1, ax6), (ax2, ax7), (ax3, ax8), (ax4, ax9), (ax5, ax10)) = plt.subplots(5, 2, figsize=(15, 15))

# Plot histogram with actual counts
bin_edges = np.arange(0, 1000, 100)
ax1.hist(allPoll_positions, bins=bin_edges, alpha=0.6, color='blue', label='Count of PolII')
ax1.hist(allSer7_positions, bins=bin_edges, alpha=0.6, color='red', label='Count of Ser7P')
ax1.set_ylabel('Count')
ax1.set_xlabel('Position on Gene (x100 bp)')
ax1.legend()

# Zoomed histogram
ax6.hist(allPoll_positions, bins=bin_edges, alpha=0.6, color='blue', label='Count of PolII')
ax6.hist(allSer7_positions, bins=bin_edges, alpha=0.6, color='red', label='Count of Ser7P')
ax6.set_ylabel('Count')
ax6.set_xlabel('Position on Gene (x100 bp)')
ax6.set_xlim(0, 50)  # Zoom to first 5kb
ax6.legend()

# Calculate and plot moving average for smoother visualization
window_size = 50
x_range = np.arange(0, 1000, 1)

# Calculate histograms
Poll_hist, _ = np.histogram(allPoll_positions, bins=1000, range=(0, 1000))
Ser7_hist, _ = np.histogram(allSer7_positions, bins=1000, range=(0, 1000))

# Calculate moving average
Poll_ma = np.convolve(Poll_hist, np.ones(window_size)/window_size, mode='same')
Ser7_ma = np.convolve(Ser7_hist, np.ones(window_size)/window_size, mode='same')

# Calculate difference and ratios
difference = Ser7_ma - Poll_ma
ratio = Ser7_ma / Poll_ma
log2_ratio = np.log2(ratio)

# Save the processed arrays
print("\nSaving processed arrays...")
np.save(os.path.join(output_dir, 'Poll_ma.npy'), Poll_ma)
np.save(os.path.join(output_dir, 'Ser7_ma.npy'), Ser7_ma)
np.save(os.path.join(output_dir, 'difference.npy'), difference)
np.save(os.path.join(output_dir, 'ratio.npy'), ratio)
np.save(os.path.join(output_dir, 'log2_ratio.npy'), log2_ratio)
print("Arrays saved successfully!")

