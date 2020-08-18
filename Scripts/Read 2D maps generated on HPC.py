# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 17:30:32 2020


@author: barraly
"""

import csv
import sabs_pkpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clrs
import matplotlib.patches as patches


# In[Load the 2D map of APD90]
# Grd10
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/Grandi 2010/paper/Grandi.csv'

# Grd Mann
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/Grandi Mann/Grandi Mann.csv'

# ORd11
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/ORd 2011/paper/ORd 2011.csv'

# ORd CiPA
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/Ord 2017/Paper/ORd CiPA.csv'

# ORd Mann
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/ORd Mann/ORd Mann.csv'

# ORd Krogh
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/ORd Krogh/ORd Krogh.csv'

# ToR ORd
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/ToR ORd/ToR ORd.csv'

# TT06
filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/Original TT06/paper/TT06.csv'

# TT06 Mann
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/TT06 Mann/TT06 Mann.csv'

# TT2.0
#filename = 'C:/Users/barraly/Documents/PhD/Benchmark of models/2D maps/TT2.0/TT2.0.csv'

APD90_array = np.loadtxt(filename, delimiter=',', skiprows = 0)
num_spaces = 101
percentage_of_IKr_block = np.linspace(0, 1, num_spaces)
percentage_of_ICaL_block = np.linspace(0, 1, num_spaces)

#
#
## In[Restructure the map] 
#index_of_state = 2
#restruct_map = np.zeros((num_spaces, num_spaces))
#for i in range(len(APD90_array)):
#    restruct_map[int(APD90_array[i,0]), int(APD90_array[i, 1])] = APD90_array[i, index_of_state]
#
#
## In[]
#viridis = cm.get_cmap('Spectral_r', 500)
#
#fig, ax = plt.subplots(1,1)
#fig.set_figheight(30)
#fig.set_figwidth(30)
#
#APD90_res = np.array(restruct_map)
#plt.imshow(APD90_res, cmap=viridis)
#
##plt.clim(120, 145)
#plt.colorbar()
#
#plt.title('2D map of APD90 versus ICaL and IKr block,\n for ' + model + ' model and ex vivo data', Fontsize = 45)
#plt.xlabel('Fraction of available IKr', Fontsize=45)
#plt.ylabel('Fraction of available ICaL', Fontsize=45)
#plt.xticks(np.linspace(0,num_spaces-1,11), np.round(np.linspace(np.min(percentage_of_IKr_block),np.max(percentage_of_IKr_block),11),3), 
#           rotation = 90,
#           Fontsize= 45)
#plt.yticks(np.linspace(0,num_spaces-1,11), np.round(1 - np.linspace(np.min(percentage_of_ICaL_block),np.max(percentage_of_ICaL_block),11),3),
#           Fontsize = 45)
#
#
#
#
#
#
#






# In[Load compounds ex vivo data]
#Load IC50 data

filename = 'C:/Users/barraly/Documents/PhD/Ion Channel Data/benchmark drugs IC50.csv'
IC50_data = np.loadtxt(filename, delimiter = ';', skiprows = 1, dtype=list)

# Load ex vivo data

filename = 'C:/Users/barraly/Documents/PhD/ex vivo data/20 APs data - 6 compounds/benchmark APD90.csv'
ex_vivo = np.loadtxt(filename, delimiter = ';', skiprows = 1, dtype=list)

# Define the function to compute the block for a given IC50 data and conc
def Hill(conc, IC50, hill):
    available = 1 / (1 + np.power(conc/IC50,hill))
    return available

# Compute the coordinates of the ex vivo data points

X = []
Y = []

for point in range(len(ex_vivo)):
    for drug in range(len(IC50_data)):
        # Look for the drug data for the ex vivo point
        if IC50_data[drug, 0] == ex_vivo[point, 0]:
            index = drug
    # Add the coordinate for the point to the X list
    X.append(Hill(float(ex_vivo[point, 1]), float(IC50_data[index, 1]), float(IC50_data[index, 2])))
    Y.append(Hill(float(ex_vivo[point, 1]), float(IC50_data[index, 3]), float(IC50_data[index, 4])))
    
    

# In[]
viridis = cm.get_cmap('Spectral_r', 500)
newcolors = viridis(np.linspace(0, 1, 500))

# Put the 0 ms line in black
black = np.array([0, 0, 0, 1])
newcolors[199:201, :] = black
newcmp = clrs.ListedColormap(newcolors)

# Put the +10 ms line in light grey
lightgrey = np.array([0, 0, 0, 0.5])
newcolors[209:211, :] = lightgrey

# Put the -10 ms line in light grey
lightgrey = np.array([0, 0, 0, 0.5])
newcolors[189:191, :] = lightgrey


cmap = newcmp

fig, ax = plt.subplots(1,1)
fig.set_figheight(30)
fig.set_figwidth(30)

APD90_res = np.flipud(APD90_array[:, :])
APD90_res = (APD90_res - APD90_res[-1, -1])
plt.imshow(APD90_res, cmap=newcmp)

plt.clim(-200, 300)
plt.colorbar()

#plt.title('2D map of the APD90 prolongation (in ms) versus ICaL and IKr block,\n for cell model and ex vivo data', Fontsize = 36)
plt.xlabel('Fraction of available IKr', Fontsize=45)
plt.ylabel('Fraction of available ICaL', Fontsize=45)
plt.xlim([0, 100])
plt.ylim([0, 100])
plt.xticks(np.linspace(0,num_spaces-1,11), np.round(np.linspace(np.min(percentage_of_IKr_block),np.max(percentage_of_IKr_block),11),3), 
           rotation = 90,
           Fontsize= 25)
plt.yticks(np.linspace(0,num_spaces-1,11), np.round(np.linspace(np.min(percentage_of_ICaL_block),np.max(percentage_of_ICaL_block),11),3),
           Fontsize = 25)

norm = clrs.Normalize(-200, 300)
cmap = plt.cm.Spectral_r
for k in range(len(ex_vivo)):
    x = X[k] * 100
    y = Y[k] * 100
    radius = float(ex_vivo[k, 3])/10
    
    for i in range(1):
        width = radius*1.25
        height = radius * 1.25
        e = patches.Ellipse((x, y), width, height)
        e.set_edgecolor(color = 'k')
        e.set_facecolor(cmap(norm(np.mean(float(ex_vivo[k, 2])))))
        ax.add_artist(e)    
ax.autoscale()

