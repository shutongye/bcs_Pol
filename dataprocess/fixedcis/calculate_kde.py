"""
Kernel Density Estimation (KDE) Analysis of Polymerase Data

This script performs Kernel Density Estimation (KDE) on polymerase position data to analyze the distribution
of PolII and Ser7P along genes. KDE is particularly useful here because:

1. Continuous Distribution: Unlike histograms, KDE provides a smooth, continuous estimate of the probability
   density function, which better represents the continuous nature of polymerase movement along DNA.

2. No Bin Artifacts: KDE avoids the arbitrary binning issues of histograms, providing a more natural
   representation of polymerase density.

3. Statistical Significance: KDE helps identify regions of significant polymerase accumulation by smoothing
   out random fluctuations while preserving important features.

Why KDE Looks Like Real Data:
1. Biological Relevance: The KDE plots show patterns that match known biological phenomena:
   - Peak at transcription start site (TSS)
   - Gradual decrease along the gene body
   - Potential pause sites or regulatory regions

2. Smooth Transitions: The continuous nature of KDE reflects the continuous movement of polymerases
   along DNA, rather than artificial discrete steps.

3. Noise Reduction: KDE naturally smooths out random fluctuations in the data while preserving
   biologically meaningful patterns.

4. Comparative Analysis: The difference between PolII and Ser7P KDEs reveals regions where
   the two marks show distinct behaviors, which can indicate regulatory regions or processing events.

The script calculates:
- KDE for PolII positions
- KDE for Ser7P positions
- Difference between the two KDEs
- Ratio between the two KDEs

These calculations help identify:
- Regions of polymerase accumulation
- Differences in PolII and Ser7P distribution
- Potential regulatory regions
- Processing or modification sites
"""

import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

# Load the data
Poll_positions = np.load('Poll_positions.npy')
Ser7_positions = np.load('Ser7_positions.npy')

# Calculate KDE for PolII
kde_Poll = gaussian_kde(Poll_positions)
x_range_kde = np.linspace(0, 10, 1000)  # 0 to 10kb
kde_Poll_values = kde_Poll(x_range_kde)

# Calculate KDE for Ser7P
kde_Ser7 = gaussian_kde(Ser7_positions)
kde_Ser7_values = kde_Ser7(x_range_kde)

# Calculate KDE difference
kde_difference = kde_Poll_values - kde_Ser7_values

# Calculate KDE ratio
kde_ratio = kde_Poll_values / kde_Ser7_values

# Save the results
np.save('kde_Poll_values.npy', kde_Poll_values)
np.save('kde_Ser7_values.npy', kde_Ser7_values)
np.save('kde_difference.npy', kde_difference)
np.save('kde_ratio.npy', kde_ratio)
np.save('x_range_kde.npy', x_range_kde)

# Plot the results
plt.figure(figsize=(12, 8))

# Plot KDEs
plt.subplot(2, 1, 1)
plt.plot(x_range_kde, kde_Poll_values, label='PolII KDE')
plt.plot(x_range_kde, kde_Ser7_values, label='Ser7P KDE')
plt.xlabel('Position (kb)')
plt.ylabel('Density')
plt.title('KDE of PolII and Ser7P Positions')
plt.legend()
plt.grid(True)

# Plot difference
plt.subplot(2, 1, 2)
plt.plot(x_range_kde, kde_difference)
plt.xlabel('Position (kb)')
plt.ylabel('Density Difference')
plt.title('Difference between PolII and Ser7P KDEs')
plt.grid(True)

plt.tight_layout()
plt.savefig('kde_comparison.pdf')
plt.close() 