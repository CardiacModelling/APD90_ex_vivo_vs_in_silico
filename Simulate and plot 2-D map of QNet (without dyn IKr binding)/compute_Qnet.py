import numpy as np
import myokit

# Parse the argument for ICaL rescale
import sys
ical_rescale = float(sys.argv[1])/100



# Load the model
m, p, _ = myokit.load('../MMT models/ORd_CiPA.mmt')

# Set 0.5 Hz pacing
protocol = myokit.Protocol()
protocol.schedule(1, 0, 0.5, 2000, 0)
s = myokit.Simulation(m, protocol)

# Set tolerance
s.set_tolerance(1e-10, 1e-10)

# Save initial state
default_state = s.state()



def compute_AP_and_currents(ikr_rescale, ical_rescale):
    # Set the drug effect
    s.set_constant('drug.ikr_rescale', ikr_rescale)
    s.set_constant('drug.ical_rescale', ical_rescale)
    
    # Pre-run for 1000 pre-paces
    s.reset()
    s.set_state(default_state)
    s.pre(2000000)
    
    # Output every 0.01 ms
    out = s.run(1990, log_interval = 0.01)
    AP = np.array(out['membrane.V'])
    IKr = np.array(out['IKr.IKr'])
    ICaL = np.array(out['ICaL.ICaL'])
    INaL = np.array(out['INaL.INaL'])
    IKs = np.array(out['IKs.IKs'])
    IK1 = np.array(out['IK1.IK1'])
    Ito = np.array(out['Ito.Ito'])
    
    return AP, IKr, ICaL, INaL, IKs, IK1, Ito
    


# Loop over the whole possible IKr rescale
ikr_rescale = np.linspace(0, 1, 101)
QNets = np.zeros(101)
for index in range(101):
    AP, IKr, ICaL, INaL, IKs, IK1, Ito = compute_AP_and_currents(ikr_rescale[index], ical_rescale)
    QNets[index] = np.sum(IKr + ICaL + INaL + IKs + IK1 + Ito) * 0.00001 # *dt for integration



# Save the outputs    
np.savetxt('Simulated QNet/ical_rescale_' + str(ical_rescale) + '.csv', QNets)
