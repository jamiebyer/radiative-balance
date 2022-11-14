import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import xarray as xr
from utilities.odes import dT_dt, dT_dt_2d, dT_dt_phi, D, Gamma, B, beta
from scipy.integrate import solve_ivp
# from utilities.write_from_models import get_2d_constants
#from mpl_toolkits.basemap import Basemap
import plotly.graph_objects as go
import plotly.express as px
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_model():
    # Time interval for integration (in seconds)
    t_min = 0
    # t_max = 1E11
    t_max = 3E8
    max_step = 1E6

    # Initial temperatures
    T_init = [0, 212, 280, 296, 295, 268, 210, 0]

    sol = solve_ivp(fun=dT_dt, t_span=(t_min, t_max), y0=T_init, method="LSODA", max_step=max_step)
    t = sol.t
    T = sol.y.T

    plt.plot(t / 3.154E7, T[:, 1:-1])
    
    plt.title("radiative balance")
    plt.xlabel("time (yr)")
    plt.ylabel("temperature (K)")
    plt.legend(["zone 1", "zone 2", "zone 3", "zone 4", "zone 5", "zone 6"])
    
    plt.tight_layout()
    plt.show()

def plot_model_withforcing():
    # Time interval for integration (in seconds)
    mon2sec = 30*24*3600
    t_min = 0
    t_max = 20*12*mon2sec
    max_step = 1E6

    # Initial temperatures
    T_init = [0, 212, 280, 296, 295, 268, 210, 0]

    dTdtFunc = lambda t, T: dT_dt_phi(t,T,t_max)
    sol = solve_ivp(fun=dTdtFunc, t_span=(t_min, t_max), y0=T_init, method="LSODA", max_step=max_step)
    t = sol.t
    T = sol.y.T

    plt.plot(t / 3.154E7, T[:, 1:-1])
    
    plt.title("Radiative balance with forcing")
    plt.xlabel("Time (yr)")
    plt.ylabel("Temperature (K)")
    plt.legend(["zone 1", "zone 2", "zone 3", "zone 4", "zone 5", "zone 6"])
    
    plt.tight_layout()
    figname = "figures/radiative-balance-forcing/radiative-balance-F{}_{}_{}_{}.png".format(D, Gamma, B, beta)
    plt.savefig(figname)
    plt.show()

'''
def plot_2d_model():
    # Resolution
    lat_res = 50
    lon_res = 50

    # Time interval for integration (in seconds)
    t_min = 0
    t_max = 1E9
    max_step = 1E6

    # Initial temperatures
    T_init = np.pad(np.full((lat_res-1, lon_res-1), 250), ((1, 1), (1, 1))).flatten()

    # Call solver
    sol = solve_ivp(fun=dT_dt_2d, t_span=(t_min, t_max), y0=T_init, method="LSODA", max_step=max_step, args=get_2d_constants(lat_res, lon_res))
    t = sol.t
    T = sol.y.T
    
    T = T.reshape((len(t), lat_res+1, lon_res+1))

    # Time vs. Temperature plot
    # Loop over each zone to plot
    for i in range(lat_res-1) :
        for j in range(lon_res-1):
            plt.plot(t / 3.154E7, T[:, i+1, j+1])
    
    plt.title("radiative balance")
    plt.xlabel("time (yr)")
    plt.ylabel("temperature (K)")
    
    plt.tight_layout()
    plt.show()

    # Map plot. Lat, lon axes.
    lats = np.linspace(-90, 90, lat_res)
    lons = np.linspace(0, 360, lon_res)

    plot_lats = lats[:-1] + (lats[1:]-lats[:-1]) / 2
    plot_lons = lons[:-1] + (lons[1:]-lons[:-1]) / 2
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    div = make_axes_locatable(ax)
    cax = div.append_axes('right', '5%', '5%')

    vmin, vmax = np.min(T[-1, 1:-1, 1:-1]), np.max(T[-1, 1:-1, 1:-1])

    cf = ax.pcolormesh(plot_lons, plot_lats, T[-1, 1:-1, 1:-1], vmin=vmin, vmax=vmax)
    cb = fig.colorbar(cf, cax=cax)
    tx = ax.set_title('t = ' + str(t[-1] / 3.154E7) + " yr")
    ax.set_ylabel("lat")
    ax.set_xlabel("lon")

    plt.show()

def animate_2d_model(save=False):
    # Resolution
    lat_res = 50
    lon_res = 50

    # Time interval for integration (in seconds)
    t_min = 0
    t_max = 1E9
    max_step = 1E6

    # Initial temperatures
    T_init = np.pad(np.full((lat_res-1, lon_res-1), 250), ((1, 1), (1, 1))).flatten()

    sol = solve_ivp(
        fun=dT_dt_2d, 
        t_span=(t_min, t_max), 
        y0=T_init, 
        method="LSODA", 
        max_step=max_step, 
        args=get_2d_constants(lat_res, lon_res)
        )
    t = sol.t
    T = sol.y.T

    lats = np.linspace(-90, 90, lat_res)
    lons = np.linspace(0, 360, lon_res)

    T = T.reshape((len(t), lat_res+1, lon_res+1))

    plot_lats = lats[:-1] + (lats[1:]-lats[:-1]) / 2
    plot_lons = lons[:-1] + (lons[1:]-lons[:-1]) / 2

    X, Y = np.meshgrid(plot_lats, plot_lons)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    div = make_axes_locatable(ax)
    cax = div.append_axes('right', '5%', '5%')

    vmin, vmax = np.min(T[-1, 1:-1, 1:-1]), np.max(T[-1, 1:-1, 1:-1])

    cf = ax.pcolormesh(plot_lons, plot_lats, T[0, 1:-1, 1:-1], vmin=vmin, vmax=vmax)
    cb = fig.colorbar(cf, cax=cax)
    tx = ax.set_title('t = ' + str(t[0]) + "s")
    ax.set_ylabel("lat")
    ax.set_xlabel("lon")

    def animate(i):
        cf = ax.pcolormesh(plot_lons, plot_lats, T[i, 1:-1, 1:-1], vmin=vmin, vmax=vmax)
        cax.cla()
        fig.colorbar(cf, cax=cax)
        tx.set_text('t = {0}s'.format(t[i]))

    anim = animation.FuncAnimation(fig, animate, frames=len(t), interval=100, repeat=False)

    if save:
        f = r"./animation.gif" 
        writergif = animation.PillowWriter(fps=30) 
        anim.save(f, writer=writergif)

    plt.show()

def plot_albedo():
    lat_res = 50
    lon_res = 50
    _, _, _, alpha, _, _, _, _, _ = get_2d_constants(lat_res, lon_res)
    lats = np.linspace(-90, 90, lat_res)
    lons = np.linspace(0, 360, lon_res)
    plot_lats = lats[:-1] + (lats[1:]-lats[:-1]) / 2
    plot_lons = lons[:-1] + (lons[1:]-lons[:-1]) / 2
    plt.pcolormesh(plot_lons, plot_lats, alpha[1:-1, 1:-1])
    cbar = plt.colorbar()
    cbar.set_label('albedo', rotation=90)
    plt.title("map of albedo")
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.tight_layout()
    plt.show()

def plot_maps():
    land = "sftlf"
    land_ice = "sftgif"
    sea_ice = "siconca"

    i = 0
    for var in [land, land_ice, sea_ice]:
        ds = xr.open_zarr("data/models/"+var+".zarr")
        if var == sea_ice:
            ds = ds.isel(time=0)
        X, Y = np.meshgrid(ds[var]["lon"], ds[var]["lat"])

        i += 1
        plt.subplot(3, 1, i)
        plt.contourf(X, Y, ds[var].values)
        plt.title(var)
        plt.tight_layout()

    plt.show()
    '''