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

# %%
import os

import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import xarray as xr

from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties

# %%
scenarios = ['counterfactual']

# %%
df_configs = pd.read_csv('../data/fair-parameters/fair-calibrate-1.5.0/calibrated_constrained_parameters.csv', index_col=0)

# %%
df_solar = pd.read_csv('../data/forcing/solar_erf_timebounds_2022zero.csv', index_col=0)
df_volcanic = pd.read_csv('../data/forcing/volcanic_erf_timebounds_2020_constant.csv', index_col=0)

# %%
# this could be parallel - let's see how painful in serial
# define problem

temp_out = np.ones((276, len(scenarios), 841)) * np.nan
faer_out = np.ones((276, len(scenarios), 841)) * np.nan
labels = []
weights = np.ones(52)
weights[0] = 0.5
weights[-1] = 0.5

for iscen, scenario in enumerate(tqdm(scenarios)):
    f = FAIR(ch4_method="Thornhill2021")
    f.define_time(1750, 2025, 1)
    f.define_scenarios([scenario])
    f.define_configs(df_configs.index)
    species, properties = read_properties('../data/fair-parameters/fair-calibrate-1.5.0/species_configs_properties.csv')
    f.define_species(species, properties)
    f.allocate()

    f.fill_from_csv('../data/emissions/historical_emissions_1750-2024.csv')

    # fill external forcings
    fill(f.forcing, df_solar.values[:, None] * df_configs['forcing_scale[Solar]'].values[None, :], specie='Solar')
    fill(f.forcing, df_volcanic.values[:, None] * df_configs['forcing_scale[Volcanic]'].values[None, :], specie='Volcanic')

    # override_defaults is very slow so read in xarray directly
    ds_sc = xr.load_dataset("../output/calibration_binaries/fair-calibrate-1.5.0/species_configs.nc")
    ds_ebm = xr.load_dataset("../output/calibration_binaries/fair-calibrate-1.5.0/energy_balance_parameters.nc")
    f.species_configs = ds_sc
    f.climate_configs = ds_ebm

    initialise(f.concentration, f.species_configs["baseline_concentration"])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)
    initialise(f.ocean_heat_content_change, 0)
    
    f.run(progress=False)
    temp_out[:, iscen, :] = f.temperature[:, 0, :, 0]# - np.average(f.temperature[100:152, 0, :, 0], weights=weights, axis=0)
    faer_out[:, iscen, :] = f.forcing[:, 0, :, 54] + f.forcing[:, 0, :, 55]
    labels.append(scenario)

# %%
pl.plot(np.median(temp_out[:, 0, :], axis=1));

# %%
np.median(temp_out[:, 0, :], axis=1)

# %%
# os.makedirs('../output/results', exist_ok=True)

# %%
ds = xr.Dataset(
    data_vars = dict(
        temperature_anomaly_rel_1750 = (['timebound', 'scenario', 'config'], temp_out),
        aerosol_forcing = (['timebound', 'scenario', 'config'], faer_out),
    ),
    coords = dict(
        timebound = np.arange(1750, 2102),
        scenario = labels,
        config = df_configs.index,
    ),
)
ds.to_netcdf('../output/results/scenariomip-cmip7.nc')

# %%
