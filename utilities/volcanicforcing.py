import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

mon2sec = 30*24*3600 # conversion factor between months and seconds
global t_max
t_max = 1.1E11/mon2sec # maximum time [months]
dt_long = 60 # longer time step for extrapolation [months]

## PARAMETERS TO CHANGE ##################################################################################
t_start = 0*mon2sec # time when eruption happens [seconds]

D = 30 # Diffusion rate [(lattitude degrees)^2/month]
Gamma = 0.05  # Precipitation rate [1/month]

a = 15 # Half-width of eruption/initial aerosol column [lattitude degrees]
B = 1 # Maximum relative aerosol concentration
beta = 0.25 # Scale factor when converting relative aerosol concentration to relative radiation blocking

dt = 0.2 # time step [months]
t_end = 12*20 # stop time [months] for shorter computation
###########################################################################################################

N = 50 # Number of spatial intervals
x = np.linspace(-90,90,N+1) # Spatial discretization
h = x[1] - x[0] # Step size [lattitude degrees]

if dt*(2*D/h**2 + Gamma) >= 1: # check that solution will be stable
    print('Solution unstable. dt*(2*D/h^2 + gamma) = ', dt*(2*D/h**2 + Gamma))

def getCDM(): # Central difference matrix
    CDM = np.zeros((N-1,N-1))
    for i in range(N-1):
        CDM[i,i] += -Gamma
        if i != 0:
            CDM[i,i] += -D/h**2
            CDM[i,i-1] = D/h**2
        if i != N-2:
            CDM[i,i] += -D/h**2
            CDM[i,i+1] = D/h**2
    return CDM


def init_dist(x,width,amp): # Define the initial distribution
    A0 = np.zeros(np.shape(x))
    A0temp = np.exp(-width**2/(width**2 - x**2))
    inrange = (x > -width)*(x < width)
    A0[inrange] = A0temp[inrange]
    return amp*np.exp(1)*A0

def plot_forcing(t,A,phi_k):
    fig,ax = plt.subplots(2,figsize=(10,10))
    # plot spatial distribution
    sampling_time_indices = [0,5,10,25,50,200]
    zone_xticks = [-90,-60,-30,0,30,60,90]
    legend = []
    for i in sampling_time_indices:
        ax[0].plot(x,A[:,i])
        legend.append('t = {} months'.format(t[i]))
    ax[0].legend(legend)
    ax[0].set_title("Aerosol distributions, D = {}, $\gamma$ = {}, B = {}".format(D,Gamma,B))
    ax[0].set_xlabel('Latitude degrees')
    ax[0].set_ylabel('Relative aerosol concentration')
    ax[0].set_xticks(zone_xticks)
    # plot temporal distribution
    for xx in range(3):
        ax[1].plot(t/12,phi_k[xx,:])
    ax[1].legend(['Zones 1 & 6', 'Zones 2 & 5', 'Zones 3 & 4'])
    ax[1].set_title(r"Zonal radiation reductions, $\beta$ = {}".format(beta))
    ax[1].set_xlabel('Time (yr)')
    ax[1].set_ylabel('Relative incoming radiation')
    plt.show()
    figname = "figures/forcing/forcing_{}_{}_{}_{}.png".format(D,Gamma,B,beta)
    fig.savefig(figname)


def create_phi_funcs(tt_max):
    if tt_max > t_max*mon2sec:
        t_end_long = tt_max/mon2sec
    else:
        t_end_long = t_max
    t = np.arange(0,t_end,dt) # time discretization
    
    CDM = getCDM()
    A = np.zeros((np.size(x), np.size(t)))
    A[:,0] = init_dist(x,a,B) # Setting the initial condition
    for j in range(1,np.size(t)): # Solve the pde
        dAdt = np.zeros(np.shape(x))
        dAdt[1:N] = np.dot(CDM,A[1:N,j-1])
        dAdt[0] = dAdt[1]
        dAdt[N] = dAdt[N-1]
        A[:,j] = A[:,j-1] + dAdt*dt # Step through time
    
    phi_continuous = 1 - beta*A # Calculate relative radiation reduction
    phi_k = np.zeros((6,np.size(t)))
    for i in range(np.size(t)):
        for j in range(-3,3):
            xrange = (x>=j*30)*(x<=(j+1)*30)
            phi_k[j+3,i] = np.mean(phi_continuous[xrange,i]) # Calculate zonal average
    
    # plot_forcing(t,A,phi_k)

    t_span_long = np.arange(t_end,t_end_long,dt_long)
    phi_extension = np.array(list(map(lambda k: 1-(1-phi_k[k,-1])*np.exp(Gamma*t[-1])*np.exp(-Gamma*t_span_long),range(6))))
    phi_k = np.append(phi_k,phi_extension,axis=1)
    t = np.append(t,t_span_long)

    phi_funcs = []
    for i in range(6):
        phi_funcs.append(interp1d(t*mon2sec,phi_k[i,:]))
    return phi_funcs

global phi_k
phi_k = create_phi_funcs(t_max)

def phi(tt, tt_max):
    global t_max
    global phi_k
    if tt_max > t_max*mon2sec:
        phi_k = create_phi_funcs(tt_max)
        t_max = tt_max/mon2sec + 1

    phi_funcs = phi_k
    phi_k_t = np.ones(8)
    if tt < t_start:
        return phi_k_t
    t = tt - t_start
    for k in range(1,7):
        phi_k_t[k] = phi_funcs[k-1](t)
    return phi_k_t
    
    



