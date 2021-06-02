import numpy as np
def se_integrate():
    samples = 1000
    step = 0.1
    u, r = np.array[samples], np.array[samples]
    u[-1] = r[-1]*np.exp(-r[-1])
    u[-2] = r[-2]*np.exp(-r[-2])
   #integrate inward using Verlet algorithm
    #DO I = SAMPLES-1, 2, -1
     # U(I-1) = 2*U(I) - U(I+1) + STEP**2*(-2*E + 2*V(I) + L*(L+1)/R(I)**2)*U(I)