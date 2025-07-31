import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.ndimage import gaussian_filter1d

# --- Configuration ---
# Define the directory where the data was saved by iteration.py
# IMPORTANT: Make sure this path is correct for your system
input_dir = 'YOURPATHWAY/TO/TRANS/MODEL/RESULTS/DIRECTORY'

# Smoothing parameters
# Adjust sigma_smoothing_individual for individual PolII/Ser7P density curve smoothness
sigma_smoothing_individual = 5
# Adjust sigma_smoothing_ratio for the smoothness of the ratio and log2 ratio curves
# Note: This sigma is effectively applied to the inputs of the ratio. The ratio itself will be inherently smoothed.
sigma_smoothing_ratio_unused = 5 # This specific sigma is no longer directly used for ratio smoothing

# Small epsilon value to prevent division by zero in ratio calculations
epsilon = 1e-9 # Used for ratios and log2 ratios to handle zero denominators
epsilon_density = 1e-12 # Used for density normalization sums to prevent division by zero

# --- Load Data ---
npz_filename = os.path.join(input_dir, 'cis_polii_ser7p_density_arrays.npz')

try:
    data = np.load(npz_filename)
    positions = data['positions']
    RNApolIIcount_avg = data['polii_density_avg']
    Ser7Pcount_avg = data['ser7p_density_avg']
    gene_length = len(positions)
    print(f"Data loaded successfully from {npz_filename}")
except FileNotFoundError:
    print(f"Error: Data file '{npz_filename}' not found.")
    print("Please ensure 'iteration.py' has been run and generated this file in the specified directory.")
    exit() # Exit the script if the data file is not found

# --- 1. Generate KDE-like curves by smoothing original averaged counts ---
print(f"Applying Gaussian smoothing to individual densities with sigma={sigma_smoothing_individual}...")
smoothed_polii = gaussian_filter1d(RNApolIIcount_avg, sigma=sigma_smoothing_individual)
smoothed_ser7p = gaussian_filter1d(Ser7Pcount_avg, sigma=sigma_smoothing_individual)

# --- Normalize smoothed curves to represent density (sum to 1) ---
# These normalized curves will be used for the first two plots and for the difference plot
normalized_polii = smoothed_polii / (np.sum(smoothed_polii) + epsilon_density)
normalized_ser7p = smoothed_ser7p / (np.sum(smoothed_ser7p) + epsilon_density)
print("Normalized smoothed PolII and Ser7P curves for density plots.")


# --- 2. Calculate Ratio of Normalized Densities (Ser7P Density / PolII Density) ---
# This ratio is inherently smoothed because its components are smoothed and normalized densities.
print("Calculating ratio of Normalized Ser7P Density / Normalized Pol II Density...")
ratio_normalized_ser7p_polii = normalized_ser7p / (normalized_polii + epsilon)

# --- 3. Calculate Log2 Ratio of Normalized Densities ---
print("Calculating Log2 ratio of Normalized Ser7P Density / Normalized Pol II Density...")
log2_ratio_normalized_ser7p_polii = np.log2(ratio_normalized_ser7p_polii + epsilon) # Add epsilon for robustness


# --- 4. Calculate the difference in densities (Ser7P Density - PolII Density) ---
# This difference is already smoothed because its components (normalized_ser7p, normalized_polii) are smoothed.
difference_density_ser7p_polii = normalized_ser7p - normalized_polii
print("Calculated difference in densities (Ser7P Density - PolII Density).")


# --- Plotting ---
plt.figure(figsize=(18, 22)) # Figure size accommodates 5 rows of plots

# Plot 1: Smoothed PolII and Ser7P Density (Full gene length)
plt.subplot(5, 2, 1)
plt.plot(positions, normalized_polii, color='blue', linestyle='-', label='Pol II (density)')
plt.plot(positions, normalized_ser7p, color='red', linestyle='--', label='Ser7P (density)')
plt.xlabel('Position along gene')
plt.ylabel('Density')
plt.title(f'Smoothed Pol II and Ser7P Density (σ={sigma_smoothing_individual})')
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 2: Zoomed-in (first 50 positions) Smoothed PolII and Ser7P Density
plt.subplot(5, 2, 2)
plt.plot(positions[:50], normalized_polii[:50], color='blue', linestyle='-', label='Pol II (density)')
plt.plot(positions[:50], normalized_ser7p[:50], color='red', linestyle='--', label='Ser7P (density)')
plt.xlabel('Position along gene')
plt.ylabel('Density')
plt.title(f'Zoomed-in (Positions 1-50) Smoothed Pol II and Ser7P Density')
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 3: Ratio of Normalized Densities (Full gene length)
plt.subplot(5, 2, 3)
plt.plot(positions, ratio_normalized_ser7p_polii, color='orange')
plt.xlabel('Position along gene')
plt.ylabel('Normalized Ser7P / Normalized Pol II Ratio') # Updated label
plt.title(f'Ratio of Normalized Densities (Smoothed by components, σ={sigma_smoothing_individual})') # Updated title
plt.axhline(y=1, color='gray', linestyle='--', label='Ratio = 1')
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 4: Zoomed-in (first 50 positions) Ratio of Normalized Densities
plt.subplot(5, 2, 4)
plt.plot(positions[:50], ratio_normalized_ser7p_polii[:50], color='orange')
plt.xlabel('Position along gene')
plt.ylabel('Normalized Ser7P / Normalized Pol II Ratio') # Updated label
plt.title(f'Zoomed-in (Positions 1-50) Ratio of Normalized Densities') # Updated title
plt.axhline(y=1, color='gray', linestyle='--', label='Ratio = 1')
plt.legend()
plt.grid(True, alpha=0.3)


# Plot 5: Log2 Ratio of Normalized Densities (Full gene length)
# This replaces the previous smoothed ratio plot with the new ratio
plt.subplot(5, 2, 5)
plt.plot(positions, log2_ratio_normalized_ser7p_polii, color='green')
plt.xlabel('Position along gene')
plt.ylabel('Log2 (Normalized Ser7P / Normalized Pol II Ratio)') # Updated label
plt.title(f'Log2 Ratio of Normalized Densities (Smoothed by components, σ={sigma_smoothing_individual})') # Updated title
plt.axhline(y=0, color='gray', linestyle='--', label='Log2 Ratio = 0 (Ratio = 1)')
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 6: Zoomed-in (first 50 positions) Log2 Ratio of Normalized Densities
plt.subplot(5, 2, 6)
plt.plot(positions[:50], log2_ratio_normalized_ser7p_polii[:50], color='green')
plt.xlabel('Position along gene')
plt.ylabel('Log2 (Normalized Ser7P / Normalized Pol II Ratio)') # Updated label
plt.title(f'Zoomed-in (Positions 1-50) Log2 Ratio of Normalized Densities') # Updated title
plt.axhline(y=0, color='gray', linestyle='--', label='Log2 Ratio = 0 (Ratio = 1)')
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 7 & 8 are now the difference plots
# The original plots 7 & 8 (Log2 of smoothed_ratio_ser7p_polii) are replaced by the new difference plots.
# The previous Plots 9 & 10 (Smoothed Difference Density) are moved to positions 7 & 8.

# Plot 7 (was 9): Difference in Density (Ser7P Density - PolII Density) (Full gene length)
plt.subplot(5, 2, 7)
plt.plot(positions, difference_density_ser7p_polii, color='darkcyan')
plt.xlabel('Position along gene')
plt.ylabel('Difference in Density (Ser7P - Pol II)')
plt.title(f'Difference in Density (Ser7P - Pol II) (σ={sigma_smoothing_individual})')
plt.axhline(y=0, color='gray', linestyle='--', label='Difference = 0')
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 8 (was 10): Zoomed-in (first 50 positions) Difference in Density (Ser7P - PolII)
plt.subplot(5, 2, 8)
plt.plot(positions[:50], difference_density_ser7p_polii[:50], color='darkcyan')
plt.xlabel('Position along gene')
plt.ylabel('Difference in Density (Ser7P - Pol II)')
plt.title(f'Zoomed-in (Positions 1-50) Difference in Density (Ser7P - Pol II)')
plt.axhline(y=0, color='gray', linestyle='--', label='Difference = 0')
plt.legend()
plt.grid(True, alpha=0.3)


plt.tight_layout()

# Save the plot
output_filename_all_plots = os.path.join(input_dir, 'cis_polii_ser7p.pdf')
plt.savefig(output_filename_all_plots, dpi=300, bbox_inches='tight')
print(f"\nAll analysis plots saved to: {output_filename_all_plots}")

plt.show()

# --- Print Summary of First 10 Values ---
print("\n--- Summary of first 10 values ---")
# Updated header for clarity
print(f"{'Position':<10}{'PolII Density':<15}{'Ser7P Density':<15}{'Ratio Densities':<18}{'Log2 Ratio Densities':<23}{'Diff Density':<15}")
for i in range(10):
    print(f"{positions[i]:<10.0f}{normalized_polii[i]:<15.4f}{normalized_ser7p[i]:<15.4f}{ratio_normalized_ser7p_polii[i]:<18.4f}{log2_ratio_normalized_ser7p_polii[i]:<23.4f}{difference_density_ser7p_polii[i]:<15.4f}")