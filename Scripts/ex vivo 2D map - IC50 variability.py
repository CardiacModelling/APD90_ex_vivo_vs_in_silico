# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:30:26 2019

@author: yanral
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import matplotlib.colors as colors


# In[Load compounds ex vivo data]
#Load IC50 data

filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/Variability of IC50/all drugs.csv'
IC50_data = np.loadtxt(filename, delimiter = ';', skiprows = 1, dtype=list)

# Transform IC50_data to numeric values
for row in range(len(IC50_data)):
    for col in range(2, len(IC50_data[row])):
        IC50_data[row, col] = float(IC50_data[row, col])

# Load ex vivo data

filename = 'C:/Users/barraly/Documents/PhD/ex vivo data/20 APs data - 6 compounds/benchmark APD90.csv'
ex_vivo = np.loadtxt(filename, delimiter = ';', skiprows = 1, dtype=list)


# Transform ex_vivo to numeric values
for row in range(len(ex_vivo)):
    for col in range(1, len(ex_vivo[row])):
        ex_vivo[row, col] = float(ex_vivo[row, col])


# In[define the functions for getting the data]
# Define the function to put weight on the different experimental points
def weight(x):
    f = np.array(x)
    for i in range(0, len(x)):
        if x[i] < 0 :
            f[i] = 1
        elif x[i] < 20:
            f[i] = 3
        elif x[i] > 20:
            f[i] = 2.5 * np.exp((20 - x[i]) / 20) + 0.5
    return f

# Define the function to compute the block for a given IC50 data and conc
def Hill(conc, IC50, hill):
    available = 1 / (1 + np.power(conc/IC50,hill))
    return available

# Define the function to get the APD90 prolongation for a drug at a given 
# concentration
def APD90_ex_vivo(drug_name, conc):
    for point in range(len(ex_vivo)):
        # Look for the drug data for the ex vivo point
        if ex_vivo[point, 0] == drug_name:
            if ex_vivo[point, 1] == conc:
                APD90 = ex_vivo[point, 2]
    return APD90

def compute_block(drug_name, conc, channel):
    output = []
    
    # Loop over IC50_data to get all the data
    for row in range(len(IC50_data)):
        # Look for the IC50 data for the drug
        if IC50_data[row, 0] == drug_name and IC50_data[row, 1] == channel:
            # Compute the corresponding rescale of channel conductance
            rescale = Hill(conc, IC50_data[row, 2], IC50_data[row, 3])
            output.append(rescale)
            
    return output


# In[Prepare the points to plot]
exp = []
X = []
Y = []

for point in range(len(ex_vivo)):
    drug = ex_vivo[point, 0]
    conc = ex_vivo[point, 1]
    
    # Loop over the experimental points to retrieve the experimental APD90 change
    exp.append(APD90_ex_vivo(drug, conc))
    
    # Compute the corresponding rescale of hERG and CaV for all of the IC50s
    X.append(compute_block(drug, conc, 'hERG'))
    Y.append(compute_block(drug, conc, 'CaV'))
    

# In[]
cmap = plt.cm.Spectral_r
norm = colors.Normalize(-200, 300)
#norm = colors.Normalize(np.min(data[:,4]), np.max(data[:,4]))


fig, ax = plt.subplots(1,1)
fig.set_figheight(7)
fig.set_figwidth(8)
for point in range(len(exp)):
    # Summarise the IC50 data
    width = np.std(X[point])
    height = np.std(Y[point])
    x = np.mean(X[point])
    y = np.mean(Y[point])
    
    # Plot an ellipse with the summarised IC50 data
    e = patches.Ellipse((x, y), width, height)
    e.set_edgecolor(np.array([0, 0, 0, 1]))
    e.set_facecolor(cmap(norm(exp[point]), alpha=0.8))
    ax.add_artist(e)
'''
    # Add crosses to show where each IC50 data lies
    # For hERG
    for data in range(len(X[point])):
        plt.scatter(X[point][data],
                    y,
                    color = cmap(norm(exp[point])),
                    marker = 'x',
                    LineWidth = 25)

    # For CaV
    for data in range(len(Y[point])):
        plt.scatter(x,
                    Y[point][data],
                    color = cmap(norm(exp[point])),
                    marker = 'x',
                    LineWidth = 25)
'''

plt.ylabel('Available fraction of ICaL', Fontsize = 15)
plt.xlabel('Available fraction of hERG', Fontsize = 15)

plt.xlim([0, 1])
plt.ylim([0, 1])

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
fig.colorbar(sm)


# In[Plot the variability of IC50 versus APD90 change]
plt.figure(figsize = (7, 7))

# Loop over the experimental data points
for point in range(len(exp)):
    # Summarise the IC50 data
    ellipse_width = np.std(X[point])
    ellipse_height = np.std(Y[point])
    variability_IC50 = np.sqrt(ellipse_height **2 + ellipse_width **2)
    excel_1.append(variability_IC50)
    
    # Read out the ex vivo data
    variability_exp = ex_vivo[point, 3]
    
    plt.scatter(variability_IC50, variability_exp, color = 'k')

plt.xlabel('IC50 variability', Fontsize = 15)
plt.ylabel('Experimental APD90 change variability', Fontsize = 15)

