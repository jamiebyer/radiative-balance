import numpy as np

def update_f(f, T):
    ## update the fraction of ice cover in zone 6 if temperature exceeds melting temperature
    fractional_changes = np.array([-0.1, 0.05, 0.05])
    ## Check temp in zone 6
    if T[6] >= 235:
       f[5,:] = (1 + fractional_changes) * f[5,:]
    else:
        pass
    return f


def update_f(f, T):
    ## update the fraction of ice cover in zone 6 if temperature exceeds melting temperature
    fractional_changes = np.array([-0.1, 0.05, 0.05])
    ## Check temp in zone 6
    if T[6] >= 235:
       f[5,:] = (1 + fractional_changes) * f[5,:]
    else:
        pass
    return f