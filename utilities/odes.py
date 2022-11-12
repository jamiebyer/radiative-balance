
from data.constants import *
from utilities.volcanicforcing import phi

def dT_dt(t, T, t_max):
    dT_dt = ((1/rho_c_Z) * (gamma*(1-alpha_sky)*(1-alpha)*phi(t,t_max)*S_0 - transmissivity*SB*(T**4)) 
    + (1/(A_E*rho_c_Z)) * np.pad((-L[0:-1]*k[0:-1]*(T[1:-1] - T[0:-2]) + L[1:]*k[1:]*(T[2:] - T[1:-1])), (1, 1)))

    return dT_dt