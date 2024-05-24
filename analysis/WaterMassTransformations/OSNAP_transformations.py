import doralite
import gfdl_utils.core as gu
from CM4Xutils import *

import xbudget
import regionate
import xwmb
import warnings

def save_wmb(name, lons, lats, grid, model, save_maps=False):

    region = regionate.GriddedRegion(
        name,
        lons,
        lats,
        grid
    )
    
    budgets_dict = xbudget.load_preset_budget(model="MOM6_3Donly")
    xbudget.collect_budgets(grid, budgets_dict)
    
    lam = "sigma2"
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)

        print("Computing WMT integrals")
        # Spatially-integrated water mass budget terms
        kwargs = {"greater_than":True, "default_bins":True}
        wmb = xwmb.WaterMassBudget(
            grid,
            budgets_dict,
            region
        )
        wmb.mass_budget(lam, integrate=True, along_section=False, **kwargs)

        print("Computing WMT transports")
        # Convergent transports, resolved along the boundary
        wmb_along = xwmb.WaterMassBudget(
            grid,
            budgets_dict,
            region
        )
        wmb_along.mass_budget(lam, integrate=True, along_section=True, **kwargs)
        wmb.wmt["convergent_mass_transport_along"] = (
                wmb_along.grid._ds['convergent_mass_transport_along'].compute()
        )
        wmb.wmt.to_zarr(f"../../data/wmb_{model}_{region.name}_2010-2024.zarr", mode="w")

        if save_maps:
            # Spatial maps of the time-mean water mass budget terms
            print("Computing WMT maps")
            wmb_maps = xwmb.WaterMassBudget(
                grid,
                budgets_dict,
                region
            )
            wmb_maps.mass_budget(lam, integrate=False, along_section=False, **kwargs)
            wmt_timemean_maps = wmb_maps.wmt.sel(xh=slice(-70, 10), yh=slice(53, 70)).drop_dims("time_bounds").mean("time")
            wmt_timemean_maps.to_zarr(f"../../data/wmb_maps_{model}_{region.name}_2010-2024.zarr", mode="w")

grid_dict = {}
models = exp_dict.keys()
for model in models:
    ds = concat_scenarios([
        load_wmt_ds(model, interval="2010", dmget=True).isel(time_bounds=slice(None, -1)),
        load_wmt_ds(model, interval="2015", dmget=True).isel(time_bounds=slice(None, -1)),
        load_wmt_ds(model, interval="2020", dmget=True)
    ])

    grid_dict[model] = make_grid(ds)

# OSNAP-EAST
name = "IrmingerIceland"
lons = np.array([-39.,  -20.,  -14,  -7.0, -4.5,  -4.5, -14.7, -28.00, -30.54, -44.90])
lats = np.array([ 68,     65,   65, 62.25,   58,  56.5,  58.0,  58.05,  58.86,  60.30])
for model, grid in grid_dict.items():
    save_wmb(name, lons, lats, grid, model)

# OSNAP-WEST
name = "Labrador"
lons = np.array([-56.8775, -52.0956, -49.8604, -47.6107, -44.8000, -50,  -56.0, -61.5,  -63.5, -63.5])
lats = np.array([52.0166,   52.6648,  53.5577,  58.8944,  60.4000,  66.5, 66.5,  66.0,   66.0,  57.5])
for model, grid in grid_dict.items():
    save_wmb(name, lons, lats, grid, model)
