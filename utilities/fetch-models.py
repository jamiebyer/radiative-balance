import xarray as xr
import pooch
import pandas as pd
import fsspec
from pathlib import Path
import json

#get esm datastore
odie = pooch.create(
    # Use the default cache folder for the operating system
    path="./.cache",
    base_url="https://storage.googleapis.com/cmip6/",
    # The registry specifies the files that can be fetched
    registry={
        "pangeo-cmip6.csv": "a19a343d57f3137d11ce89e07e4c66586ed02e3d396361bf209eb7a324328b89",
    },
)

file_path = odie.fetch("pangeo-cmip6.csv")
df_og = pd.read_csv(file_path)

fs = fsspec.filesystem("filecache", target_protocol='gs', target_options={'anon': True}, cache_storage='./.cache/files/')


def save_model(var_id, mod_id, exp_id):
    model_path = Path('data/models/'+var_id+'.zarr')
    
    query = "variable_id=='"+var_id+"' & experiment_id=='"+exp_id+"' & source_id=='"+mod_id+"'"
    df = df_og.query(query)
    #print(df.drop_duplicates(['source_id'])['source_id'])
    zstore_url = df["zstore"].values[0]
    the_mapper=fs.get_mapper(zstore_url)
    ds = xr.open_zarr(the_mapper, consolidated=True)

    # write data to zarr
    ds.to_zarr(model_path, mode='w')


land = "sftlf" # Percentage of the grid cell occupied by land (including lakes) [%]
land_ice = "sftgif" # Land Ice Area Percentage [%]
sea_ice_atm = "siconca" # Sea-ice Area Percentage (Atmospheric Grid) [%]
grid_area = "areacella" # Grid-Cell Area for Atmospheric Grid Variables [m2]

mod_id = "CanESM5" # for land, land_ice, grid area
# mod_id = "CESM2" #for sea_ice_atm
exp_id = "piControl"

#save_model(grid_area, mod_id, exp_id)
