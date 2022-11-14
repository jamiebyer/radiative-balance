import numpy as np
from data.constants import *

def update_f(f, T):
    '''
    Update the fraction of ice cover in zones 6 and 1 if temperature exceeds melting temperature.
    INPUT
    f - 6 X 3 array, unitless
        array of fractions of land, water, and ice for each zone
    T - float, [K]
        temperature at each step
    OUTPUT
    f - 6 x 3 array, unitless
        Updated array of land-type fractions
    '''
    # change of ice fraction each year (1-decadal_decrease/10)/secs/year
    change_step = 12.6/3.154E7
    fractional_changes_6 = np.array([0.5, 0.5, -1])*change_step
    fractional_changes_1 = np.array([0.9, 0.1, -1])*change_step
    ## Check temp in zone 6, if exceeds melting temperature change ice fraction
    if T[6] >=  273. and f[5, 2] > 0:
        f[5,:] = fractional_changes_6 + f[5,:]
    else:
        pass
    ## Check zone 1, if exceeds melting temperature change ice fraction
    if T[1] >= 273. and f[0, 2] > 0:
        f[0,:] = fractional_changes_1 + f[0,:]
    else:
        pass
    return f

def atmos_emissivity(t):
    """
    Give differing Atmospheric emissivity values to be used to simulate GHG warming. Assume emissivity is 
    constant at pre-industrial levels before 1880, and capture increasing temperature due to GHGs by fitting 
    to global surface temperature anomaly record (GISS Surface Temp Analysis, 2022). 

    INPUT
    t - float or integer, [years]
        Time step for emissivity to be calculated at
    OUTPUT:
    emissivity - float, unitless
        Value of emissivity for given
    """
    T_e = (S_0*(1-0.3)/(4*SB))**(0.25)
    temp_record = temp_anomalies[2:-2, 3] + 287.2
    year = temp_anomalies[2:-2, 0]
    if t > 1880:
        model_to_fit = np.polyfit(year, temp_record, 2)
        fit_to_temp_record = np.poly1d(model_to_fit)
        emissivity = 2*(1-(T_e/fit_to_temp_record(t))**(4))
    else: 
        # Determine unforced temperature value taking mean of temperature from 1880 - 1920
        temp_unforced = np.mean(temp_record[0:39])
        emissivity = 2*(1-(T_e/temp_unforced)**(4))
    return emissivity