
from data.constants import *

def dT_dt(t, T):
    dT_dt = ((1/rho_c_Z) * (gamma*(1-alpha_sky)*(1-alpha)*S_0 - transmissivity*SB*(T**4)) 
    + (1/(A_E*rho_c_Z)) * np.pad((-L[0:-1]*k[0:-1]*(T[1:-1] - T[0:-2]) + L[1:]*k[1:]*(T[2:] - T[1:-1])), (1, 1)))

    return dT_dt