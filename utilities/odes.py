
from data.constants import *
from utilities.update_functions import *

def dT_dt(t, T):
    dT_dt = ((1/rho_c_Z) * (gamma*(1-alpha_sky)*(1-alpha)*S_0 - transmissivity*SB*(T**4)) 
    + (1/(A_E*rho_c_Z)) * np.pad((-L[0:-1]*k[0:-1]*(T[1:-1] - T[0:-2]) + L[1:]*k[1:]*(T[2:] - T[1:-1])), (1, 1)))

    if t > 500 * 3.154E7 and t < 600 *  3.154E7:
        dT_dt[5] += 0.5E-5

    return dT_dt

def dT_dt_changefrac(t, T):

    ## Update f
    fracs_new = update_f(fractions, T)
    rho_c_Z_new = np.pad(np.dot(fracs_new, LOI_rho_Z_c), (1, 1), constant_values=(1, 1)) # Area averaged product of density, specific heat capacity, thermal scale depth 
    alpha_new = np.pad(np.dot(fracs_new, LOI_alpha), (1, 1)) # Area averaged albedo
    
    ## use these vars in ODE
    dT_dt = ((1/rho_c_Z_new) * (gamma*(1-alpha_sky)*(1-alpha_new)*S_0 - transmissivity*SB*(T**4)) 
    + (1/(A_E*rho_c_Z_new)) * np.pad((-L[0:-1]*k[0:-1]*(T[1:-1] - T[0:-2]) + L[1:]*k[1:]*(T[2:] - T[1:-1])), (1, 1)))

    if t > 500 * 3.154E7 and t < 550 *  3.154E7:
        dT_dt[5] += 0.5E-5

    return dT_dt