# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Convert configs files to netcdf
#
# Saves lots of I/O overhead.

# %%
import os

import numpy as np
import pandas as pd
import xarray as xr

from fair import FAIR
from fair.io import read_properties

# %%
calibration_versions = ['1.5.0']

for calver in calibration_versions:
    df_configs = pd.read_csv(f'../data/fair-parameters/fair-calibrate-{calver}/calibrated_constrained_parameters.csv', index_col=0)
    df_defaults = pd.read_csv(f'../data/fair-parameters/fair-calibrate-{calver}/species_configs_properties.csv', index_col=0)

    # I think all is required here is to create a dummy instance of fair and read all in from CSV
    f = FAIR(ch4_method="Thornhill2021")
    f.define_time(2023, 2023, 1)
    f.define_scenarios(['dummy'])
    f.define_configs(df_configs.index)
    species, properties = read_properties(f'../data/fair-parameters/fair-calibrate-{calver}/species_configs_properties.csv')
    f.define_species(species, properties)
    f.allocate()
    f.fill_species_configs(f'../data/fair-parameters/fair-calibrate-{calver}/species_configs_properties.csv')
    f.override_defaults(f'../data/fair-parameters/fair-calibrate-{calver}/calibrated_constrained_parameters.csv')

    outdir = f"../output/calibration_binaries/fair-calibrate-{calver}"
    os.makedirs(outdir, exist_ok = True)
    f.climate_configs.to_netcdf(os.path.join(outdir, 'energy_balance_parameters.nc'))
    f.species_configs.to_netcdf(os.path.join(outdir, 'species_configs.nc'))

# %%
