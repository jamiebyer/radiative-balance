import xarray as xr
import numpy as np
import csv

land = "sftlf"
land_ice = "sftgif"
sea_ice = "siconca"
grid_area = "areacella"

# Open data for land, land ice, sea ice
ds_land = xr.open_zarr("data/models/"+land+".zarr")
ds_land_ice = xr.open_zarr("data/models/"+land_ice+".zarr")
ds_sea_ice = xr.open_zarr("data/models/"+sea_ice+".zarr").isel(time=slice(0, 12)).mean(dim="time")
ds_grid_area = xr.open_zarr("data/models/"+grid_area+".zarr")

# Array to store zonal fractions
f = [[],[],[]]

for i in np.arange(-90, 60+1, 30):
    zone_area = ds_grid_area.sel(lat=slice(i, i+30))[grid_area].values.mean()
    zone_weighted_area = zone_area / zone_area.mean()

    land_frac = (ds_land*zone_weighted_area).sel(lat=slice(i, i+30))[land].values.mean()
    land_ice_frac = (ds_land_ice*zone_weighted_area).sel(lat=slice(i, i+30))[land_ice].values.mean()
    sea_ice_frac = (ds_sea_ice*zone_weighted_area).sel(lat=slice(i, i+30))[sea_ice].values.mean()

    f[0].append((land_frac - land_ice_frac)/100) # land
    f[1].append((100 - land_frac - sea_ice_frac)/100) # ocean
    f[2].append((land_ice_frac + sea_ice_frac)/100) # ice

# Write zonal fractions to csv
with open("./data/model-zonal-fractions.csv","w+") as out_csv:
    csvWriter = csv.writer(out_csv, delimiter=',')
    csvWriter.writerows(np.array(f).T)
