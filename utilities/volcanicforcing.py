import numpy as np
from scipy.interpolate import interp1d

mon2sec = 30*24*3600 # conversion factor between months and seconds

## PARAMETERS TO CHANGE ##################################################################################
t_start = 0*mon2sec # time when eruption happens [seconds]

D = 30 # Diffusion rate [(lattitude degrees)^2/month]
gamma = 0.01  # Precipitation rate [1/month]

a = 15 # Half-width of eruption/initial aerosol column [lattitude degrees]
B = 0.5 # Maximum relative aerosol concentration
beta = 0.75 # Scale factor when converting relative aerosol concentration to relative radiation blocking

dt = 0.2 # time step [months]
t_end = 12*10 # stop time [months] for shorter computation

dt_long = 60 # longer time step for extrapolation [months]
t_max = 2E11/mon2sec # maximum time [months]
###########################################################################################################

N = 50 # Number of spatial intervals
x = np.linspace(-90,90,N+1) # Spatial discretization
h = x[1] - x[0] # Step size [lattitude degrees]

if dt*(2*D/h**2 + gamma) >= 1: # check that solution will be stable
    print('Solution unstable. dt*(2*D/h^2 + gamma) = ', dt*(2*D/h**2 + gamma))

def getCDM(): # Central different matrix
    CDM = np.zeros((N-1,N-1))
    for i in range(N-1):
        CDM[i,i] += -gamma
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
    
    t_span_long = np.arange(t_end,t_end_long,dt_long)
    phi_extension = np.array(list(map(lambda k: 1-(1-phi_k[k,-1])*np.exp(gamma*t[-1])*np.exp(-gamma*t_span_long),range(6))))
    phi_k = np.append(phi_k,phi_extension,axis=1)
    t = np.append(t,t_span_long)

    phi_funcs = []
    for i in range(6):
        phi_funcs.append(interp1d(t*mon2sec,phi_k[i,:]))
    return phi_funcs

phi_k = create_phi_funcs(t_max)

def phi(tt, tt_max):
    if tt_max > t_max*mon2sec:
        phi_funcs = create_phi_funcs(tt_max)
    else:
        phi_funcs = phi_k
    phi_k_t = np.ones(8)
    if tt < t_start:
        return phi_k_t
    t = tt - t_start
    for k in range(1,7):
        phi_k_t[k] = phi_funcs[k-1](t)
    return phi_k_t
    
    


