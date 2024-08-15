import doralite
import gfdl_utils.core as gu
from CM4Xutils import *

import xbudget
import regionate
import xwmb
import warnings

def save_wmb_maps(name, lons, lats, grid, model, sigma2_target):

    region = regionate.GriddedRegion(
        name,
        lons,
        lats,
        grid
    )

    years = [grid._ds["time"][0].dt.year.values, grid._ds["time"][-1].dt.year.values]
    
    budgets_dict = xbudget.load_preset_budget(model="MOM6_3Donly")
    xbudget.collect_budgets(grid, budgets_dict)
    
    lam = "sigma2"
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)
        kwargs = {"greater_than":True, "default_bins":True}
        wmb_maps = xwmb.WaterMassBudget(
            grid,
            budgets_dict,
            region
        )
        wmb_maps.mass_budget(lam, integrate=False, along_section=False, **kwargs)
        wmt_timemean_maps = wmb_maps.wmt[["boundary_fluxes"]].sel(xh=slice(-70, 10), yh=slice(53, 70)).sel(sigma2_l_target=sigma2_target).mean("time").compute()
        wmt_timemean_maps.to_zarr(f"../../data/wmb_boundary_flux_maps_{model}_{region.name}_{years[0]}-{years[1]}.zarr", mode="w")

grid_dict = {}
models = ["CM4Xp125"]#exp_dict.keys()
for model in models:
    for interval in ["2010", "2015", "2020"]:
        ds = load_wmt_ds(model, interval=interval, dmget=True)
        grid_dict[model] = make_grid(ds)
    
        # OSNAP-EAST
        name = "IrmingerIceland"
        lons = np.array([-39.,  -20.,  -14,  -7.0, -4.5,  -4.5, -14.7, -28.00, -30.54, -44.90])
        lats = np.array([ 68,     65,   65, 62.25,   58,  56.5,  58.0,  58.05,  58.86,  60.30])
        for model, grid in grid_dict.items():
            save_wmb_maps(name, lons, lats, grid, model, sigma2_target=36.4)
        
        # OSNAP-WEST
        name = "Labrador"
        lons = np.array([-56.8775, -52.0956, -49.8604, -47.6107, -44.8000, -50,  -56.0, -61.5,  -63.5, -63.5])
        lats = np.array([52.0166,   52.6648,  53.5577,  58.8944,  60.4000,  66.5, 66.5,  66.0,   66.0,  57.5])
        for model, grid in grid_dict.items():
            save_wmb_maps(name, lons, lats, grid, model, sigma2_target=36.6)
