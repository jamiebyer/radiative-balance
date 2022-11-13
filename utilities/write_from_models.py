import xarray as xr
import numpy as np
import csv
from data.constants import *

# Variable names for datasets.
land = "sftlf"
land_ice = "sftgif"
sea_ice = "siconca"
grid_area = "areacella"

# Open data for land, land ice, sea ice. Average sea ice over one year.
ds_land = xr.open_zarr("data/models/"+land+".zarr")
ds_land_ice = xr.open_zarr("data/models/"+land_ice+".zarr")
ds_sea_ice = xr.open_zarr("data/models/"+sea_ice+".zarr").isel(time=slice(0, 12)).mean(dim="time")

# Area of each cell for the other datasets.
ds_grid_area = xr.open_zarr("data/models/"+grid_area+".zarr")

# Interpolate sea ice to have same lat lon as other datasets
ds_sea_ice = ds_sea_ice.interp(lat=ds_grid_area["lat"], lon=ds_grid_area["lon"])

def write_1d_fractions():
    # Write 1d zonal fractions to a csv.

    # Array to store zonal fractions
    f = [[],[],[]]

    # Loop over each zone latitude to add fractions to f
    for i in np.arange(-90, 60+1, 30):
        # Weigh fractions by area, since area varies over a zone
        zone_area = ds_grid_area.sel(lat=slice(i, i+30))[grid_area].values
        zone_weighted_area = zone_area / zone_area.mean()

        # Get mean land and ice fractions from datasets
        land_frac = (zone_weighted_area*ds_land.sel(lat=slice(i, i+30))[land].values).mean()
        land_ice_frac = (zone_weighted_area*ds_land_ice.sel(lat=slice(i, i+30))[land_ice].values).mean()
        sea_ice_frac = (zone_weighted_area*ds_sea_ice.sel(lat=slice(i, i+30))[sea_ice].values).mean()

        f[0].append((land_frac - land_ice_frac)/100) # land
        f[1].append((100 - land_frac - sea_ice_frac)/100) # ocean
        f[2].append((land_ice_frac + sea_ice_frac)/100) # ice

    # Write zonal fractions to csv
    with open("./data/model-zonal-fractions.csv","w+") as out_csv:
        csvWriter = csv.writer(out_csv, delimiter=',')
        csvWriter.writerows(np.array(f).T)

def get_2d_constants(lat_res, lon_res):
    # Return constants needed for 2d model

    # Initialize rho_c_Z, alpha arrays.
    rho_c_Z = np.ones((lat_res-1, lon_res-1))
    alpha = np.zeros((lat_res-1, lon_res-1))

    lats = np.linspace(-90, 90, lat_res)
    lons = np.linspace(0, 360, lon_res)

    # Loop over each zone.
    for i in range(len(lats)-1):
        for j in range(len(lons)-1):
            # Get weighting factor based on area across zone.
            zone_area = ds_grid_area.sel(lat=slice(lats[i], lats[i+1]), lon=slice(lons[j], lons[j+1]))[grid_area].values
            zone_weighted_area = zone_area / zone_area.mean()

            # Get mean land, ice fractions from datasets
            land_frac = (zone_weighted_area*ds_land.sel(lat=slice(lats[i], lats[i+1]), lon=slice(lons[j], lons[j+1]))[land].values).mean()
            land_ice_frac = (zone_weighted_area*ds_land_ice.sel(lat=slice(lats[i], lats[i+1]), lon=slice(lons[j], lons[j+1]))[land_ice].values).mean()
            sea_ice_frac = (zone_weighted_area*ds_sea_ice.sel(lat=slice(lats[i], lats[i+1]), lon=slice(lons[j], lons[j+1]))[sea_ice].values).mean()

            f_land = (land_frac - land_ice_frac)/100
            f_ocean = (100 - land_frac - sea_ice_frac)/100
            f_ice = (land_ice_frac + sea_ice_frac)/100

            f = np.array([f_land, f_ocean, f_ice])
            f = f.clip(min=0) # Remove negative values

            if np.sum(f) > 1:
                # If the total fraction is more than 1, scale it down
                f = (1/(np.sum(f)))*f

            # Update rho_c_Z, alpha arrays from f values
            rho_c_Z[i, j] = np.dot(f, LOI_rho_Z_c)
            alpha[i, j] = np.dot(f, LOI_alpha)


    rho_c_Z = np.pad(rho_c_Z, ((1, 1), (1, 1)), constant_values=(1, 1))
    alpha = np.pad(alpha, ((1, 1), (1, 1)))

    # Constant value for k
    k = 1E7 # W m-1 K-1
    
    # Get length of lat, lon boundaries.
    L_lat = np.zeros(lat_res)
    L_lon = np.zeros(lon_res)

    L_lat[1:] = 2*np.pi*R_E*np.cos(lats[1:])
    L_lon[1:] = 2*np.pi*R_E*np.cos(lons[1:])

    # Area of each spherical latitudinal zone
    A_sphere_zones = A_E * (1/2)*(np.sin(np.deg2rad(lats[1:])) - np.sin(np.deg2rad(lats[0:-1])))
    A = A_sphere_zones / lon_res # Area of each zone, split by lat and lon
    A = np.pad(A, (1, 1))

    # Get area of each circular latitudinal zone
    theta = np.deg2rad(2*(lats + 90))
    A_segments = (R_E**2 / 2) * (theta - np.sin(theta))

    A_disk = A_segments[1:] - A_segments[:-1]

    # Calculate gamma for each latitudinal zone
    gamma = np.zeros(lat_res + 1)
    gamma[1:-1] = A_disk / A_sphere_zones

    return lat_res, lon_res, rho_c_Z, alpha, L_lat, L_lon, k, A, gamma

#write_1d_fractions()