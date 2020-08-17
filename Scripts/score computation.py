# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 17:04:01 2020

@author: yann-stanislas.barral@roche.com
"""

import numpy as np
import matplotlib.pyplot as plt


# In[Load compounds ex vivo data]
#Load IC50 data

filename = 'C:/Users/barraly/Documents/PhD/Ion Channel Data/benchmark drugs IC50.csv'
IC50_data = np.loadtxt(filename, delimiter = ';', skiprows = 1, dtype=list)

# Transform IC50_data to numeric values
for row in range(len(IC50_data)):
    for col in range(1, len(IC50_data[row])):
        IC50_data[row, col] = float(IC50_data[row, col])

# Load ex vivo data

filename = 'C:/Users/barraly/Documents/PhD/ex vivo data/20 APs data - 6 compounds/benchmark APD90.csv'
ex_vivo = np.loadtxt(filename, delimiter = ';', skiprows = 1, dtype=list)


# Transform ex_vivo to numeric values
for row in range(len(ex_vivo)):
    for col in range(1, len(ex_vivo[row])):
        ex_vivo[row, col] = float(ex_vivo[row, col])

# In[define the score function]

# Define the function to put weight on the different experimental points
def weight(x):
    f = np.array(x)
    for i in range(0, len(x)):
        if x[i] < 0 :
            f[i] = 3 * np.exp((x[i]) / 5) + 1
        elif x[i] < 20:
            f[i] = 3
        elif x[i] > 20:
            f[i] = 3 * np.exp((20 - x[i]) / 20) + 1
                    
    return f

# Define the function to compute the block for a given IC50 data and conc
def Hill(conc, IC50, hill):
    available = 1 / (1 + np.power(conc/IC50,hill))
    return available

# Define the function to retrieve the indices of the simulated map to 
# Read out the predited APD90 change related to a drug at known concentration
def coordinates(drug_name, conc):
    for drug in range(len(IC50_data)):
        # Look for the drug data for the ex vivo point
        if IC50_data[drug, 0] == drug_name:
            index = drug
    
    # Compute the coordinates
    x = int(Hill(conc,IC50_data[index, 1], IC50_data[index, 2]) * 100)
    y = int(Hill(conc, IC50_data[index, 3], IC50_data[index, 4]) * 100)
    
    # Determine the index of the map to retrieve
    y = 100 - y
    
    return x, y

# Define the function to get the APD90 prolongation for a drug at a given 
# concentration
def APD90_ex_vivo(drug_name, conc):
    for point in range(len(ex_vivo)):
        # Look for the drug data for the ex vivo point
        if ex_vivo[point, 0] == drug_name:
            if ex_vivo[point, 1] == conc:
                APD90 = ex_vivo[point, 2]
    return APD90
    
# Define the function to comnpute the score for a drug
def score_drug(map_sim, drug_name):
    # Read out the Baseline APD90
    # map_sim must be entered as a 101x101 array
    baseline = map_sim[0, 100] 
    
    # Read out the concentrations from the ex vivo data
    concs = []
    for i in range(len(ex_vivo)):
        if ex_vivo[i, 0] == drug_name:
            concs.append(ex_vivo[i, 1])
    
    # Retrieve the simulated APD90 change
    coords = [coordinates(drug_name, concs[0]),
              coordinates(drug_name, concs[1]),
              coordinates(drug_name, concs[2])]
        
    # Summarise that in the APD90_sim variable
    APD90_sim = np.array([map_sim[coords[0][0], coords[0][1]], map_sim[coords[1][0], coords[1][1]], map_sim[coords[2][0], coords[2][1]]]) - baseline
    
    # Read the experimental APD90 change
    APD90_exp = np.array([APD90_ex_vivo(drug_name, concs[0]),
                     APD90_ex_vivo(drug_name, concs[1]),
                     APD90_ex_vivo(drug_name, concs[2])])
    
    # Apply the weighting function
    weight_cloz = weight(APD90_exp)
    
    # Compute the score for the drug
    return np.sum(weight_cloz * np.power(abs(APD90_sim - APD90_exp), 1/2))
    

def score(map_sim):
     
    # Compute the score for the drug
    score_cloz = score_drug(map_sim, 'Clozapine')
    
    # Compute the score for the drug
    score_dof = score_drug(map_sim, 'Dofetilide')
    
    # Compute the score for the drug
    score_nif = score_drug(map_sim, 'Nifedipine')
    
    # Compute the score for the drug
    score_qui = score_drug(map_sim, 'Quinidine')
    
    # Compute the score for the drug
    score_sot = score_drug(map_sim, 'Sotalol')
    
    # Compute the score for the drug
    score_ver = score_drug(map_sim, 'Verapamil')
    
    return np.array([score_cloz, score_dof, score_nif, score_qui, score_sot, score_ver])


# In[Load the data for the models and compute the corresponding scores]
Grd10 = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/Grandi 2010/paper/Grandi.csv', delimiter = ',')
score_Grd10 = score(Grd10)

GrdMann = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/Grandi Mann/Grandi Mann.csv', delimiter = ',')
score_GrdMann = score(GrdMann)

ORd11 = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/ORd 2011/paper/ORd 2011.csv', delimiter = ',')
score_ORd11 = score(ORd11)

ORd_CiPA = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/Ord 2017/Paper/ORd CiPA.csv', delimiter = ',')
score_ORd_CiPA = score(ORd_CiPA)

ORdKrg = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/ORd Krogh/ORd Krogh.csv', delimiter = ',')
score_ORdKrg = score(ORdKrg)

ORdMann = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/ORd Mann/ORd Mann.csv', delimiter = ',')
score_ORdMann = score(ORdMann)

ToRORd = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/ToR ORd/ToR ORd.csv', delimiter = ',')
score_ToRORd = score(ToRORd)

TT = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/Original TT06/paper/TT06.csv', delimiter = ',')
score_TT = score(TT)

TTMann = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/TT06 Mann/TT06 Mann.csv', delimiter = ',')
score_TTMann = score(TTMann)

TTupd = np.loadtxt('C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/TT2.0/TT2.0.csv', delimiter = ',')

score_TTupd = score(TTupd)


# In[Show the results of the benchmark by increasing score]

benchmark = np.array([score_TTupd, score_TT, score_GrdMann, score_TTMann, score_ORdKrg, score_ORdMann, score_Grd10, score_ToRORd, score_ORd_CiPA, score_ORd11])
labels = ['TT 2.0', 'TT06','GrdMann', 'TT Mann', 'ORd Krg', 'ORd Mann', 'Grd10','ToR ORd', 'ORd CiPA', 'ORd11']
width = 0.35       # the width of the bars: can also be len(x) sequence

# Sort by increasing score
benchmark_sorted = np.zeros(np.shape(benchmark))
labels_sorted = labels[:]

# Compute the total score for each model
tot_score = np.sum(benchmark, 1)
for_sorting = np.sum(benchmark, 1)

for i in range(len(tot_score)):
    # Find the index for which the lowest score is achieved
    index = np.where(for_sorting == np.min(for_sorting))[0][0]
    
    # Update the sorted outputs
    benchmark_sorted[i] = benchmark[index]
    labels_sorted[i] = labels[index]
    tot_score[i] = for_sorting[index]
    
    # Delete the minimum to get the next value
    for_sorting = np.delete(for_sorting, index, axis = 0)
    benchmark = np.delete(benchmark, index, axis = 0)
    labels.remove(labels[index])

# Create the figure    
plt.figure(figsize = (15, 10))

# Draw the stacked bar chart for the benchmark
plt.bar(labels_sorted, benchmark_sorted[:, 0], width, label='Clozapine')
plt.bar(labels_sorted, benchmark_sorted[:, 1], width, label='Dofetilide', bottom = benchmark_sorted[:, 0])
plt.bar(labels_sorted, benchmark_sorted[:, 2], width, label='Nifedipine', bottom = np.sum(benchmark_sorted[:, :2], axis = 1))
plt.bar(labels_sorted, benchmark_sorted[:, 3], width, label='Quinidine', bottom = np.sum(benchmark_sorted[:, :3], axis = 1))
plt.bar(labels_sorted, benchmark_sorted[:, 4], width, label='Sotalol', bottom = np.sum(benchmark_sorted[:, :4], axis = 1))
plt.bar(labels_sorted, benchmark_sorted[:, 5], width, label='Verapamil', bottom = np.sum(benchmark_sorted[:, :5], axis = 1))

# Add a textbox to display the score of each model
for i in range(len(labels_sorted)):
    plt.text(-0.25 + i, np.round(tot_score[i]) + 5, str(np.round(tot_score[i])), Fontsize = 15)

plt.legend(fontsize = 12)
plt.xticks(Fontsize = 20, rotation = 45)


