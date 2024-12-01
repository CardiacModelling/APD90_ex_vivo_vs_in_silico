import myokit
import numpy as np
import pandas as pd
import sys

# Change the model name as desired
model = '../MMT models/ORd_CiPA.mmt'
print('\nModel for the simulation : ' + model)



############################ LOAD THE MODEL ########################################################################
print('\nLoading the model ...')
m, p, _ = myokit.load(model)

print('Model loaded. Creating the simulation...')
s = myokit.Simulation(m, p)
print('Simulation ready.\n')

initial_state = s.state()
s.set_tolerance(1e-09, 1e-09)



########################### DEFINE THE SIMULATIONS #################################################################
def simulate(s, ikr_rescale, ical_rescale):
    # Rescale the currents
    s.set_constant('drug.ikr_rescale', ikr_rescale)
    s.set_constant('drug.ical_rescale', ical_rescale)
    
    # Replicate the external concentrations in the ex vivo experiments
    s.set_constant(K_o, 4)
    s.set_constant(Na_o, 148.35)
    s.set_constant(Ca_o, 1.8)    
    
    # Reset and run at steady-state
    s.reset()
    s.set_state(initial_state)
    s.pre(1500000)
    output = s.run(2000, log_interval = 0.05)
    APs = np.array(output['membrane.V'])
    
    # Compute the APD for two consecutive APs to avoid alternans
    APD_1 = compute_APD(APs[:20000])
    APD_2 = compute_APD(APs[20000:])
    
    return max(APD_1, APD_2)


def compute_APD(AP):
    # Normalise the AP
    peak = np.percentile(AP, 95)
    RMP = np.mean(AP[-3000:])
    normalised = (np.array(AP) - RMP) / (peak - RMP)
    
    # Deduce the APD as the moment where 90% repolarisation is reached
    APD90 = np.where(normalised[2000:] < 0.1)[0][0] / 20 + 50
    return APD90  

def compute_score(DAPD_sim, DAPD_exp, sigma):
    return np.sum(np.power((DAPD_sim - DAPD_exp) / sigma, 1))
    
    

########################### RUN THE SIMULATIONS AND SAVE FOR THE PHARM DATASET #####################################
# Define the IKr/ICaL rescaling factors
Pharm = True
dataset_pharm = pd.read_csv('../Experimental data/Data for benchmark - Pharm.csv', delimiter = ',')
n_tested = len(dataset_pharm)

print('\nDataset correctly loaded.\n')

# Initialise the output
APD90s = np.zeros(n_tested)

# Run the simulations
for i in range(n_tested):
    # Print progress
    print('Computing condition ' + str(i) + '. IKr rescale : ' + str(1.0 - float(dataset_pharm['IKr block'][i])) + ', ICaL rescale: ' + str(1.0 - float(dataset_pharm['ICaL block'][i])))
    
    # Compute
    APD90s[i] = simulate(s, 1.0 - float(dataset_pharm['IKr block'][i]), 1.0 - float(dataset_pharm['ICaL block'][i]))

# Compute APD90 change from baseline
APD90_base = simulate(s, 1, 1)
DAPD_sim = APD90s - APD90_base

# Compute the score
score = compute_score(DAPD_sim, np.array(dataset_pharm['DAPD']), np.array(dataset_pharm['STD']))
print('The score for the model ' + model + 'is : ' + str(score))

# Save
filename = 'Predictions/DAPD ' + model + ' - Pharm.csv'
np.savetxt(filename, DAPD_sim, delimiter = ',', newline = '\n')



########################### RUN THE SIMULATIONS AND SAVE FOR THE CiPA DATASET #####################################
# Define the IKr/ICaL rescaling factors
Pharm = False
dataset_CiPA = pd.read_csv('../Experimental data/Data for benchmark - CiPA.csv', delimiter = ',')
n_tested = len(dataset_CiPA)

# Initialise the output
APD90s = np.zeros(n_tested)

# Run the simulations
for i in range(n_tested):
    # Print progress
    print('Computing condition ' + str(i) + '. IKr rescale : ' + str(1.0 - float(dataset_CiPA['IKr block'][i])) + ', ICaL rescale: ' + str(1.0 - float(dataset_CiPA['ICaL block'][i])))
    
    # Compute
    APD90s[i] = simulate(s, 1.0 - float(dataset_CiPA['IKr block'][i]), 1.0 - float(dataset_CiPA['ICaL block'][i]))

# Compute APD90 change from baseline
APD90_base = simulate(s, 1, 1)
DAPD_sim = APD90s - APD90_base

# Compute the score
score = compute_score(DAPD_sim, np.array(dataset_CiPA['DAPD']), np.array(dataset_CiPA['STD']))
print('The score for the model ' + model + 'is : ' + str(score))

# Save
filename = 'Predictions/DAPD ' + model + ' - CiPA.csv'
np.savetxt(filename, DAPD_sim, delimiter = ',', newline = '\n')



