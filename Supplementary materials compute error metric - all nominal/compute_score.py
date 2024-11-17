import myokit
import numpy as np
import sys
import pandas as pd

# Read out inputs to the script
model ='../MMT models/Bartolucci_2020.mmt' # Should be replaced by the desired MMT model
print('\nModel for the simulation : ' + model)

# Get the labels of the external concentrations in the MMT model
# These may need to be changed if the model is changed
K_o = 'extracellula.ko'
Na_o = 'extracellula.nao'
Ca_o = 'extracellula.cao'



############################ LOAD THE MODEL ########################################################################
print('\nLoading the model ...')
m, p, _ = myokit.load(model)

print('Model loaded. Creating the simulation...')
s = myokit.Simulation(m, p)
print('Simulation ready.\n')

initial_state = s.state()
s.set_tolerance(1e-10, 1e-10)



########################### DEFINE THE SIMULATIONS #################################################################
def simulate(s, ikr_block, ical_block):
    # Rescale the currents
    s.set_constant('drug.ikr_rescale', 1 - ikr_block)
    s.set_constant('drug.ical_rescale', 1 - ical_block)
    
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
    return np.sum(np.power((DAPD_sim - DAPD_exp) / sigma, 2))
    
    

########################### RUN THE SIMULATIONS AND SAVE FOR THE PHARM DATASET #####################################
# Define the IKr/ICaL rescaling factors
Pharm = True
dataset_pharm = pd.read_csv('../Experimental data/Data for benchmark - all nominal - Pharm.csv', delimiter = ';')
n_tested = len(dataset_pharm)

print('\nDataset correctly loaded.\n')

# Initialise the output
APD90s = np.zeros(n_tested)

# Run the simulations
for i in range(n_tested):
    # Print progress
    print('Computing condition ' + str(i) + '. IKr block : ' + str(float(dataset_pharm['IKr block'][i])) + ', ICaL block: ' + str(float(dataset_pharm['ICaL block'][i])))
    
    # Compute
    APD90s[i] = simulate(s, float(dataset_pharm['IKr block'][i]), float(dataset_pharm['ICaL block'][i]))

# Compute APD90 change from baseline
APD90_base = simulate(s, 0, 0)
DAPD_sim = APD90s - APD90_base

# Compute the score
score = compute_score(DAPD_sim, np.array(dataset_pharm['DAPD']), np.array(dataset_pharm['STD']))
print('The score for the model ' + model + 'is : ' + str(score))

# Save
filename = 'predictions - all nominal/DAPD ' + model + ' - Pharm.csv'
np.savetxt(filename, DAPD_sim, delimiter = ',', newline = '\n')



########################### RUN THE SIMULATIONS AND SAVE FOR THE PHARM DATASET #####################################
# Define the IKr/ICaL rescaling factors
Pharm = False
dataset_CiPA = pd.read_csv('../Experimental data/Data for benchmark - all nominal - CiPA.csv', delimiter = ';')
n_tested = len(dataset_CiPA)

# Initialise the output
APD90s = np.zeros(n_tested)

# Run the simulations
for i in range(n_tested):
    # Print progress
    print('Computing condition ' + str(i) + '. IKr block : ' + str(float(dataset_CiPA['IKr block'][i])) + ', ICaL block: ' + str(float(dataset_CiPA['ICaL block'][i])))
    
    # Compute
    APD90s[i] = simulate(s, float(dataset_CiPA['IKr block'][i]), float(dataset_CiPA['ICaL block'][i]))

# Compute APD90 change from baseline
APD90_base = simulate(s, 0, 0)
DAPD_sim = APD90s - APD90_base

# Compute the score
score = compute_score(DAPD_sim, np.array(dataset_CiPA['DAPD']), np.array(dataset_CiPA['STD']))
print('The score for the model ' + model + 'is : ' + str(score))

# Save
filename = 'predictions - all nominal/DAPD ' + model + ' - CiPA.csv'
np.savetxt(filename, DAPD_sim, delimiter = ',', newline = '\n')



