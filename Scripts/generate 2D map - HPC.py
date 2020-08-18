# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:06:09 2019

@author: yann-stanislas.barral@roche.com
"""

import sabs_pkpd
import numpy as np
import multiprocessing as mp
import csv

# REPLACE THE PATH BELOW BY THE DESIRED MMT MODEL YOU
# MIGHT WANT TO CHANGE THE DESTINATION CSV FILE FOR THE MAP
# ON LINE 74
mmt = './MMT models/model.mmt'

# Load the model using sabs_pkpd package
sabs_pkpd.constants.s = sabs_pkpd.load_model.load_simulation_from_mmt(mmt)
sabs_pkpd.constants.default_state = sabs_pkpd.constants.s.state()

# Define all of the simulation conditions and parameters
num_spaces = 101
percentage_of_ICaL_block = 1 - np.linspace(0,1,num_spaces)
percentage_of_IKr_block = np.linspace(0,1,num_spaces)
fixed_params_annot = ['drug.ical_rescale', 'drug.ikr_rescale'] # Define the annotation of parameters which value will be changed during simulations
time_max = 1000 # maximal time for readout (in ms)
n_timepoints = 4001 # amount of points sampled for 1 AP
time_samples = np.linspace(0, time_max, n_timepoints) # Time points to sample
pre_run = 1500000 # Duration of pre-pace (in ms)
read_out = 'membrane_potential.V_m' # Define the annotation of parameter to read out from the simulation


def evaluate(i, j):
    # Use quick_simulate function to retrieve the AP with the selected simulation parameters
    AP = sabs_pkpd.run_model.quick_simulate(sabs_pkpd.constants.s, 1000, read_out, 
                                            fixed_params_annot=fixed_params_annot,
                                            fixed_params_values=[percentage_of_ICaL_block[i],percentage_of_IKr_block[j]],
                                            pre_run=pre_run,
                                            time_samples=time_samples)
    
    # Read out the APD90 for the simulated AP with the selected block of hERG/CaV1.2
    output = sabs_pkpd.cardiac.compute_APD(list(AP[0]), time_samples, 50, print_warnings=False)
    q.put([i, j, output])
    
# Initialise the multiprocessing Queue
q = mp.Queue()

# Initialise the array with all the data for output
APD_array = np.zeros((len(percentage_of_ICaL_block), len(percentage_of_IKr_block)))

# Explore the whole block of hERG/CaV1.2 possibilities
for i in range(len(percentage_of_ICaL_block)):
    # Cleaning Processes at each iteration over percentage_of_IKr_block allows to join
    # often enough, not to overload the HPC memory and preserve the performances
    Processes = []
    print(str(i))
    
    # Parallelise the evaluation of the map
    for j in range(len(percentage_of_IKr_block)):
        proc = mp.Process(target=evaluate, args=(i, j,))
        proc.start()
        Processes.append(proc)
        
    # Retrieve the output once all the threads finished running
    for proc in Processes:
        proc.join()
    
    # Save to APD_array the iteration over percentage_of_IKr_block
    while not q.empty():
        temp = q.get()
        APD_array[temp[0], temp[1]] = temp[2]
        
# Write to the CSV file selected below the evaluated map
with open('./model.csv', 'w') as csvfile:
    mywriter = csv.writer(csvfile, delimiter=',')
    for i in range(len(percentage_of_ICaL_block)):
        mywriter.writerow(APD_array[i, :])
