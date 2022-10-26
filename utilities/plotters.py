import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

def plot_models():
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