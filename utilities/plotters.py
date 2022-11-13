import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from utilities.odes import *
from scipy.integrate import solve_ivp


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
    
    plt.title("radiative balance")
    plt.xlabel("time (yr)")
    plt.ylabel("temperature (K)")
    plt.legend(["zone 1", "zone 2", "zone 3", "zone 4", "zone 5", "zone 6"])
    
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

def plot_model_fracchanges():
    # Time interval for integration (in seconds)
    t_min = 0
    t_max = 1E11
    max_step = 1E6

    # Initial temperatures
    T_init = [0, 212, 280, 296, 295, 268, 210, 0]

    sol = solve_ivp(fun=dT_dt_changefrac, t_span=(t_min, t_max), y0=T_init, method="LSODA", max_step=max_step)
    t = sol.t
    T = sol.y.T

    plt.plot(t / 3.154E7, T[:, 1:-1])
    
    plt.title("radiative balance")
    plt.xlabel("time (yr)")
    plt.ylabel("temperature (K)")
    plt.legend(["zone 1", "zone 2", "zone 3", "zone 4", "zone 5", "zone 6"])
    
    plt.tight_layout()
    plt.show()