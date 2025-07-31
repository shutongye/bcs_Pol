import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as p

gene_length = 1000 #length of the gene in 100 bp

input_dir = '/home/sy432/rds/rds-ye_shutong-xcywAxU6Kd0/Pol_model/trans_flagd/d_reset'
f = open(os.path.join(input_dir, '500sim.simulation.bcs'), 'r') #open bcs output

RNApolIIcount_all = np.zeros(gene_length) #counts all PolII from any gene
RNApolIIcount_sim = np.zeros(gene_length)
Ser7Pcount_all = np.zeros(gene_length) #counts all Ser7P from any gene  
Ser7Pcount_sim = np.zeros(gene_length)
sim_iteration = 0
sim_number = 500

for line in f:
    if sim_iteration == sim_number+1:
        break
    
    if line.startswith('>'):
        RNApolIIcount_all += RNApolIIcount_sim #add the current simulation to the total coun
        RNApolIIcount_sim = np.zeros(gene_length) #keeps track of PolII positions on this simulation
        Ser7Pcount_all += Ser7Pcount_sim
        Ser7Pcount_sim = np.zeros(gene_length)
        sim_iteration += 1
        continue

    splitLine = line.rstrip().split('\t')
        # Skip lines that don't have the expected format (like headers or separators)
    if len(splitLine) < 7 or not splitLine[0].replace('.', '').replace('-', '').isdigit():
        continue

    time = float(splitLine[0])
    action = splitLine[1]
    process = splitLine[2]
    i = int(splitLine[4])
    p = int(splitLine[6])

    if action == 'Pol_ii':

        if i > 0:
            if RNApolIIcount_sim[i-1] != 0:
                RNApolIIcount_sim[i-1] -= 1 
            
        if i <= gene_length:
            if RNApolIIcount_sim[i] == 0:
                RNApolIIcount_sim[i] += 1

        if p ==1:
            if i > 0:
                if Ser7Pcount_sim[i-1] != 0:
                    Ser7Pcount_sim[i-1] -= 1
            
            if i <= gene_length:
                if Ser7Pcount_sim[i] == 0:
                    Ser7Pcount_sim[i] +=1


RNApolIIcount_all += RNApolIIcount_sim #add the current simulation to the total coun
Ser7Pcount_all += Ser7Pcount_sim
f.close()

print(RNApolIIcount_all)
print(Ser7Pcount_all)

#kde plot of RNApolIIcount_all
RNApolIIcount_avg = RNApolIIcount_all / sim_number
Ser7Pcount_avg = Ser7Pcount_all / sim_number

positions = np.arange(1, gene_length + 1) 


plt.figure(figsize=(15, 10))

# Plot 1: PolII density (averaged)
plt.subplot(2, 2, 1)
plt.bar(positions, RNApolIIcount_avg, alpha=0.7, color='blue')
plt.xlabel('Position along gene')
plt.ylabel('Average Pol II count per simulation')
plt.title(f'Pol II Density vs Position (Averaged over {sim_number} simulations)')
plt.grid(True, alpha=0.3)

# Plot 2: Ser7P density (averaged)
plt.subplot(2, 2, 2)
plt.bar(positions, Ser7Pcount_avg, alpha=0.7, color='red')
plt.xlabel('Position along gene')
plt.ylabel('Average Ser7P count per simulation')
plt.title(f'Ser7P Density vs Position (p=1 only, averaged over {sim_number} simulations)')
plt.grid(True, alpha=0.3)

# Plot 3: PolII density (total counts)
plt.subplot(2, 2, 3)
plt.bar(positions, RNApolIIcount_all, alpha=0.7, color='darkblue')
plt.xlabel('Position along gene')
plt.ylabel('Total Pol II count across all simulations')
plt.title(f'Pol II Density vs Position (Total counts from {sim_number} simulations)')
plt.grid(True, alpha=0.3)

# Plot 4: Ser7P density (total counts)
plt.subplot(2, 2, 4)
plt.bar(positions, Ser7Pcount_all, alpha=0.7, color='darkred')
plt.xlabel('Position along gene')
plt.ylabel('Total Ser7P count across all simulations')
plt.title(f'Ser7P Density vs Position (Total counts from {sim_number} simulations)')
plt.grid(True, alpha=0.3)

plt.tight_layout()


# Save the plot to the input directory
output_filename = os.path.join(input_dir, 'trans_flagd_500sim.pdf')
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Plot saved to: {output_filename}")

plt.show()

# Save the data arrays to a file (both averaged and total)
data_filename = os.path.join(input_dir, 'trans_polii_ser7p_density_data.txt')
with open(data_filename, 'w') as f:
    f.write("# Position\tAverage_PolII_Count\tAverage_Ser7P_Count\tTotal_PolII_Count\tTotal_Ser7P_Count\n")
    for pos in range(1, gene_length + 1):
        f.write(f"{pos}\t{RNApolIIcount_avg[pos-1]:.6f}\t{Ser7Pcount_avg[pos-1]:.6f}\t{int(RNApolIIcount_all[pos-1])}\t{int(Ser7Pcount_all[pos-1])}\n")

print(f"Data saved to: {data_filename}")

# Also save as numpy arrays for easy loading
npz_filename = os.path.join(input_dir, 'trans_polii_ser7p_density_arrays.npz')
np.savez(npz_filename, 
         positions=np.arange(1, gene_length + 1),
         polii_density_avg=RNApolIIcount_avg,
         ser7p_density_avg=Ser7Pcount_avg,
         polii_density_total=RNApolIIcount_all,
         ser7p_density_total=Ser7Pcount_all)
print(f"NumPy arrays saved to: {npz_filename}")
