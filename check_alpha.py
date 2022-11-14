from utilities.ice_change_funcs import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from utilities.odes import dT_dt_changefrac as dT_dt
from data.constants import *

t_min = 0
t_max = 2000*3.154E7
max_step = 1E6
T_init = [0, 212, 280, 296, 295, 268, 210, 0]

dTdtFunc = lambda t, T: dT_dt(t,T,t_max)
sol = solve_ivp(fun=dTdtFunc, t_span=(t_min, t_max), y0=T_init, method="LSODA", max_step=max_step)
t = sol.t
T = sol.y.T

alpha = np.zeros([len(t), 8])
fracs = np.zeros([len(t), 6, 3])
fracs[0, :, :] = np.genfromtxt('data/model-zonal-fractions.csv', delimiter=',')

for i in range(1, len(t)):
    fracs[i, :, :] = update_f(fracs[i-1, :, :], T[i, :])
    alpha[i, :] =  np.pad(np.dot(fracs[i, :, :], LOI_alpha), (1, 1)) # Area averaged albedo

print(fracs[5, 5, :], fracs[-1, 5, :],'\n', alpha[5, 6], alpha[-1, 6])
t = t[1:]

plt.figure(1, figsize=(10,10))
plt.subplot(311)
plt.plot(t/3.154E7, T[1:, 1:7])
plt.ylabel('temp (k)')
plt.legend(['zone 1', 'zone 2', 'zone 3', 'zone 4', 'zone 5', 'zone 6'])
plt.subplot(312)
plt.plot(t/3.154E7, fracs[1:, 5, 0], label = 'land')
plt.plot(t/3.154E7, fracs[1:, 5, 1], label = 'water')
plt.plot(t/3.154E7, fracs[1:, 5, 2], label = 'ice')
plt.ylabel('fractions')
plt.legend()
plt.subplot(313)
plt.plot(t/3.154E7, alpha[1:, 6])
#plt.legend(['zone 1', 'zone 2', 'zone 3', 'zone 4', 'zone 5', 'zone 6'])
plt.xlabel('time (yr)')
plt.ylabel('alpha')

plt.show()