
from data.constants import *
from utilities.volcanicforcing import phi, D, gamma, B, beta
from utilities.update_functions import *

def dT_dt(t, T):

    dT_dt = ((1/rho_c_Z) * (gamma*(1-alpha_sky)*(1-alpha)*S_0 - transmissivity*SB*(T**4)) 
    + np.pad((1/(a[1:-1]*A_E*rho_c_Z[1:-1])), (1, 1)) * np.pad((-L[0:-1]*k[0:-1]*(T[1:-1] - T[0:-2]) + L[1:]*k[1:]*(T[2:] - T[1:-1])), (1, 1)))

    return dT_dt

def dT_dt_phi(t,T,t_max):
    # Include radiative forcing in ODE
    dT_dt = ((1/rho_c_Z) * (gamma*(1-alpha_sky)*(1-alpha)*phi(t,t_max)*S_0 - transmissivity*SB*(T**4)) 
    + np.pad((1/(a[1:-1]*A_E*rho_c_Z[1:-1])), (1, 1)) * np.pad((-L[0:-1]*k[0:-1]*(T[1:-1] - T[0:-2]) + L[1:]*k[1:]*(T[2:] - T[1:-1])), (1, 1)))

    return dT_dt

def dT_dt_2d(t, T, lat_res, lon_res, rho_c_Z, alpha, L_lat, L_lon, k, A, gamma):
    # Reshape array since solve_ivp works with 1d arrays.
    T = T.reshape((lat_res+1, lon_res+1))

    # Add last longitudinal temp to start of array, and add first longitudinal temp to end of array.
    # This is for calculating the derivative between the first and last longitudinal zone.
    T[:, 0] = T[:, -2]
    T[:, -1] = T[:, 1]

    dT_dt = ((1/rho_c_Z) * ((gamma*((1-alpha_sky)*(1-alpha)*S_0).T).T - transmissivity*SB*(T**4))) 
    + (1/((A*rho_c_Z.T).T)) * (np.pad((np.multiply(-L_lat[0:-1]*k, (T[1:-1, :] - T[0:-2, :]).T).T + np.multiply(L_lat[1:]*k, (T[2:, :] - T[1:-1, :]).T).T), ((1, 1), (0, 0))) 
    + np.pad(np.multiply(-L_lon[0:-1]*k, (T[:, 1:-1] - T[:, 0:-2])) + np.multiply(L_lon[1:]*k, (T[:, 2:] - T[:, 1:-1])), ((0, 0), (1, 1))))

    return dT_dt.flatten()

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