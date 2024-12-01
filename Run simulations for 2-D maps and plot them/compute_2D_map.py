import myokit
import numpy as np
import sys

"""
To run this script, one needs to pass 3 arguments: 
The input corresponding to the percentage of available IKr --> num==100 <=> 0% IKr block. This can be passed as an array of jobs. This is particularly recommended to parallelise the computation of the 2D maps.
The model name, e.g., "ORd_CiPA". Look for it in the folder "MMT models".
The folder to log out the simulation results, e.g., "ORd_CiPA". This must match one of the sub-folders.
A command line example is" `python compute_2D_map.py 100 "ORd_CiPA" "ORd_CiPA"` for computing the column of the 2D map corresponding to a 0% IKr block.
"""

# Read out inputs to the script
num = float(sys.argv[1]) # Will be parsed for the IKr rescaling factor. 100 = 0% IKr block, 0 = 100% IKr block.
model = '../MMT models/' + str(sys.argv[2]) + '.mmt' # Retrieve the MMT filename
folder = str(sys.argv[3]) # Define the folder where to save the results



############################ LOAD THE MODEL ########################################################################
# Load the model
print('\nLoading the model...')
m, p, _ = myokit.load(model)

print('\nModel loaded. Preparing Simulation...')
s = myokit.Simulation(m, p)
initial_state = s.state()

print('\nSimulation correctly loaded.\n')


########################### DEFINE THE SIMULATIONS ################################################################
def simulate(ikr_rescale, ical_rescale):
    # Reset
    s.reset()
    s.set_state(initial_state)
    
    # Rescale the currents
    s.set_constant('drug.ikr_rescale', ikr_rescale)
    s.set_constant('drug.ical_rescale', ical_rescale)
    
    # Run to steady-state
    s.pre(1500000)
    output = s.run(1000, log_interval = 0.05)
    AP = np.array(output['membrane.V'])
    
    # Compute the APD for two consecutive APs to avoid alternans
    APD = compute_APD(AP)

    return APD


def compute_APD(AP):
    # Normalise the AP
    peak = np.percentile(AP, 95)
    RMP = np.mean(AP[-3000:])
    normalised = (np.array(AP) - RMP) / (peak - RMP)
    
    # Deduce the APD as the moment where 90% repolarisation is reached
    APD90 = np.where(normalised[2000:] < 0.1)[0][0] / 20 + 50
    return APD90 



########################### RUN THE SIMULATIONS AND SAVE ###########################################################
# Define the IKr/ICaL rescaling factors
ikr_rescale = num/100 # From the script input
print('IKr rescale : ' + str(ikr_rescale) + '.')
ical_rescales = np.linspace(0, 1, 101)

# Initialise the output
APD90s = np.zeros(101)

# Run the simulations
for i in range(101):
    print('Computing with ical_rescale: ' + str(ical_rescales[i]) + '...')
    try:
        try:
            s.set_tolerance(1e-10, 1e-10)
            APD90 = simulate(ikr_rescale, ical_rescales[i])
            APD90s[i] = APD90
        except:
            s.set_tolerance(1e-09, 1e-09)
            APD90 = simulate(ikr_rescale, ical_rescales[i])
            APD90s[i] = APD90
    except:
        APD90s[i] = 0

# Save
filename = folder + '/ikr_' + str(ikr_rescale) + '.csv'
np.savetxt(filename, APD90s, delimiter = ',', newline = '\n')
