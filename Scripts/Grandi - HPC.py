# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:06:09 2019

@author: yanral
"""

import sabs_pkpd
import numpy as np
import multiprocessing as mp
import csv


mmt = '/pstore/home/barraly/MMT_models/Benchmark/Grd10.mmt'
sabs_pkpd.constants.s = sabs_pkpd.load_model.load_simulation_from_mmt(mmt)
sabs_pkpd.constants.default_state = sabs_pkpd.constants.s.state()

num_spaces = 101
percentage_of_ICaL_block = 1 - np.linspace(0,1,num_spaces)
percentage_of_IKr_block = np.linspace(0,1,num_spaces)
fixed_params_annot = ['drug.ical_rescale', 'drug.ikr_rescale']
time_max = 1000 # ms
n_timepoints = 4001 # index of the AP trigger
time_samples = np.linspace(0, time_max, n_timepoints)
pre_run = 2000000
read_out = 'membrane_potential.V_m'


def evaluate(i, j):
    AP = sabs_pkpd.run_model.quick_simulate(sabs_pkpd.constants.s, 1000, read_out, 
                                            fixed_params_annot=fixed_params_annot,
                                            fixed_params_values=[percentage_of_ICaL_block[i],percentage_of_IKr_block[j]],
                                            pre_run=pre_run,
                                            time_samples=time_samples)
    output = sabs_pkpd.cardiac.compute_APD(list(AP[0]), time_samples, 50, print_warnings=False)
    q.put([i, j, output])
    

q = mp.Queue()

APD_array = np.zeros((len(percentage_of_ICaL_block), len(percentage_of_IKr_block)))

for i in range(len(percentage_of_ICaL_block)):
    Processes = []
    print(str(i))
    
    for j in range(len(percentage_of_IKr_block)):
        proc = mp.Process(target=evaluate, args=(i, j,))
        proc.start()
        Processes.append(proc)
        
    for proc in Processes:
        proc.join()
    
    while not q.empty():
        temp = q.get()
        APD_array[temp[0], temp[1]] = temp[2]
        
with open('/pstore/home/barraly/Codes/Benchmark/Grandi.csv', 'w') as csvfile:
    mywriter = csv.writer(csvfile, delimiter=',')
    for i in range(len(percentage_of_ICaL_block)):
        mywriter.writerow(APD_array[i, :])