import xarray as xr
import numpy as np
import csv

land = "sftlf"
land_ice = "sftgif"
sea_ice = "siconca"

# Open data for land, land ice, sea ice
ds_land = xr.open_zarr("data/models/"+land+".zarr")
ds_land_ice = xr.open_zarr("data/models/"+land_ice+".zarr")
ds_sea_ice = xr.open_zarr("data/models/"+sea_ice+".zarr").isel(time=0)

# Array to store zonal fractions
f = [[],[],[]]

for i in np.arange(-90, 60+1, 30):
    land_frac = ds_land.sel(lat=slice(i, i+30))[land].values.mean()
    land_ice_frac = ds_land_ice.sel(lat=slice(i, i+30))[land_ice].values.mean()
    sea_ice_frac = ds_sea_ice.sel(lat=slice(i, i+30))[sea_ice].values.mean()

    f[0].append(land_frac - land_ice_frac) # land
    f[1].append(100 - land_frac - sea_ice_frac) # ocean
    f[2].append(land_ice_frac + sea_ice_frac) # ice

# Write zonal fractions to csv
with open("./data/model-zonalfractions.csv","w+") as out_csv:
    csvWriter = csv.writer(out_csv, delimiter=',')
    csvWriter.writerows(np.array(f).T)
