import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import xarray as xr
from utilities.odes import dT_dt, dT_dt_2d
from scipy.integrate import solve_ivp
from utilities.write_from_models import get_2d_constants
from mpl_toolkits.basemap import Basemap
import plotly.graph_objects as go
import plotly.express as px
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_model():
    # Time interval for integration (in seconds)
    t_min = 0
    t_max = 1E11
    max_step = 1E6

    # Initial temperatures
    T_init = [0, 212, 280, 296, 295, 268, 210, 0]

    sol = solve_ivp(fun=dT_dt, t_span=(t_min, t_max), y0=T_init, method="LSODA", max_step=max_step)
    t = sol.t
    T = sol.y.T

    plt.plot(t / 3.154E7, T[:, 1:-1])
    #plt.plot(t / 3.154E7, T[:, 1])
    
    plt.title("radiative balance")
    plt.xlabel("time (yr)")
    plt.ylabel("temperature (K)")
    #plt.legend(["zone 1", "zone 2", "zone 3", "zone 4", "zone 5", "zone 6"])
    
    plt.tight_layout()
    plt.show()

def plot_2d_model():
    # Resolution
    lat_res = 7
    lon_res = 7

    # Time interval for integration (in seconds)
    t_min = 0
    #t_max = 1E11
    t_max = 1E-5
    max_step = 1E5# = 0.1#1E6

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

    plot_lats = lats[1:] + (lats[1:]-lats[:-1]) / 2
    plot_lons = lons[1:] + (lons[1:]-lons[:-1]) / 2

    X, Y = np.meshgrid(plot_lats, plot_lons)
    plt.contourf(X, Y, T[-1, 1:-1, 1:-1])
    
    plt.title("radiative balance")
    plt.xlabel("lat")
    plt.ylabel("lon")
    
    plt.tight_layout()
    plt.show()

def animate_2d_model():
    # Resolution
    lat_res = 50
    lon_res = 50

    # Time interval for integration (in seconds)
    t_min = 0
    t_max = 1E5 #1E11
    max_step = 1E3 #1E6

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

    vmin, vmax = 249.8, 250.2

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

        #plt.title("radiative balance")
        #ax.set_ylabel("lon")
        #ax.set_xlabel("lat")    
        #ax.tight_layout()    

    ani = FuncAnimation(fig, animate, frames=len(t), interval=100, repeat=False)

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