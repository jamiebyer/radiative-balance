
from data.constants import *

def dT_dt(t, T):
    dT_dt = ((1/rho_c_Z) * (gamma*(1-alpha_sky)*(1-alpha)*S_0 - transmissivity*SB*(T**4)) 
    + np.pad((1/(a[1:-1]*A_E*rho_c_Z[1:-1])), (1, 1)) * np.pad((-L[0:-1]*k[0:-1]*(T[1:-1] - T[0:-2]) + L[1:]*k[1:]*(T[2:] - T[1:-1])), (1, 1)))

    
    #print(1/(a*A_E*rho_c_Z))

    return dT_dt

def dT_dt_2d(t, T, lat_res, lon_res, rho_c_Z, alpha, L_lat, L_lon, k, A, gamma):
    T = T.reshape((lat_res+1, lon_res+1))

    dT_dt = ((1/rho_c_Z) * (gamma*(1-alpha_sky)*(1-alpha)*S_0 - transmissivity*SB*(T**4))) 
    + (1/(A*rho_c_Z)) * (np.pad((np.multiply(-L_lat[0:-1]*k, (T[1:-1, :] - T[0:-2, :]).T).T + np.multiply(L_lat[1:]*k, (T[2:, :] - T[1:-1, :]).T).T), ((1, 1), (0, 0))) 
    + np.pad(np.multiply(-L_lon[0:-1]*k, (T[:, 1:-1] - T[:, 0:-2])) + np.multiply(L_lon[1:]*k, (T[:, 2:] - T[:, 1:-1])), ((0, 0), (1, 1))))

    return dT_dt.flatten()