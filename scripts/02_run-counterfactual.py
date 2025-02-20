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
# temperature rebasing: only used for reporting here
weights = np.ones(52)
weights[0] = 0.5
weights[-1] = 0.5

# %%
# output arrays
temp_out = np.ones((276, 2, 841)) * np.nan
faer_out = np.ones((276, 2, 841)) * np.nan
fnat_out = np.ones((276, 2, 841)) * np.nan
fant_out = np.ones((276, 2, 841)) * np.nan

# %%
# full scenario names
labels = ['counterfactual', 'counterfactual_internal_variability']

# %%
df_configs = pd.read_csv('../data/fair-parameters/fair-calibrate-1.5.0/calibrated_constrained_parameters.csv', index_col=0)

# %%
df_solar = pd.read_csv('../data/forcing/solar_erf_timebounds_1750-2025.csv', index_col=0)
df_volcanic = pd.read_csv('../data/forcing/volcanic_erf_timebounds_1750-2025.csv', index_col=0)

# %%
df_solar['counterfactual'].values

# %% [markdown]
# ## Counterfactual: no internal variability

# %%
scenarios = ['counterfactual']

# %%
# this could be parallel - let's see how painful in serial
# define problem

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
    fill(f.forcing, df_solar['counterfactual'].values[:, None, None] * df_configs['forcing_scale[Solar]'].values[None, None, :], specie='Solar')
    fill(f.forcing, df_volcanic['counterfactual'].values[:, None, None] * df_configs['forcing_scale[Volcanic]'].values[None, None :], specie='Volcanic')

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
    fnat_out[:, iscen, :] = f.forcing[:, 0, :, 52] + f.forcing[:, 0, :, 53]
    fant_out[:, iscen, :] = f.forcing_sum[:, 0, :] - fnat_out[:, iscen, :]

# %%
pl.plot(np.median(temp_out[:, 0, :] - np.average(f.temperature[100:152, 0, :, 0], weights=weights, axis=0), axis=1));

# %%
os.makedirs('../output/results', exist_ok=True)

# %%
np.median(temp_out[250:, 0, :] - np.average(f.temperature[100:152, 0, :, 0], weights=weights, axis=0), axis=1)

# %% [markdown]
# ## Counterfactual, with internal variability

# %%
scenarios = ['counterfactual']

# %%
# this could be parallel - let's see how painful in serial
# define problem

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
    fill(f.forcing, df_solar['counterfactual'].values[:, None, None] * df_configs['forcing_scale[Solar]'].values[None, None, :], specie='Solar')
    fill(f.forcing, df_volcanic['counterfactual'].values[:, None, None] * df_configs['forcing_scale[Volcanic]'].values[None, None, :], specie='Volcanic')

    # override_defaults is very slow so read in xarray directly
    ds_sc = xr.load_dataset("../output/calibration_binaries/fair-calibrate-1.5.0/species_configs.nc")
    ds_ebm = xr.load_dataset("../output/calibration_binaries/fair-calibrate-1.5.0/energy_balance_parameters.nc")
    f.species_configs = ds_sc
    f.climate_configs = ds_ebm

    # turn on internal variability this time
    fill(f.climate_configs['stochastic_run'], True)
    fill(f.climate_configs['use_seed'], True)

    initialise(f.concentration, f.species_configs["baseline_concentration"])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)
    initialise(f.ocean_heat_content_change, 0)
    
    f.run(progress=False)
    temp_out[:, 1, :] = f.temperature[:, 0, :, 0]
    faer_out[:, 1, :] = f.forcing[:, 0, :, 54] + f.forcing[:, 0, :, 55]
    fnat_out[:, 1, :] = f.forcing[:, 0, :, 52] + f.forcing[:, 0, :, 53]
    fant_out[:, 1, :] = f.forcing_sum[:, 0, :] - fnat_out[:, 1, :]

# %%
pl.plot(np.median(temp_out[:, 1, :] - np.average(f.temperature[100:152, 0, :, 0], weights=weights, axis=0), axis=1));

# %%
np.median(temp_out[250:, 1, :] - np.average(f.temperature[100:152, 0, :, 0], weights=weights, axis=0), axis=1)

# %%
ds = xr.Dataset(
    data_vars = dict(
        temperature_anomaly_rel_1750 = (['timebound', 'scenario', 'config'], temp_out),
        aerosol_forcing = (['timebound', 'scenario', 'config'], faer_out),
        natural_forcing = (['timebound', 'scenario', 'config'], fnat_out),
        anthropogenic_forcing = (['timebound', 'scenario', 'config'], fant_out),
    ),
    coords = dict(
        timebound = np.arange(1750, 2026),
        scenario = labels,
        config = df_configs.index,
    ),
)
ds.to_netcdf('../output/results/counterfactual.nc')

# %%
