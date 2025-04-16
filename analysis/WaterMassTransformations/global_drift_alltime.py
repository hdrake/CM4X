import doralite
import gfdl_utils.core as gu
from CM4Xutils import *

import xbudget
import regionate
import xwmb
import warnings

def save_wmb(grid, model, interval):

    budgets_dict = xbudget.load_preset_budget(model="MOM6_drift")
    xbudget.collect_budgets(grid, budgets_dict)
    
    lam = "sigma2"
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)

        print("Computing WMT integrals")
        # Spatially-integrated water mass budget terms
        kwargs = {"greater_than":True, "default_bins":True}
        wmb = xwmb.WaterMassBudget(
            grid,
            budgets_dict
        )
        wmb.mass_budget(lam, integrate=True, along_section=False, **kwargs)
        wmb.wmt.chunk({"time":1}).to_zarr(f"../../data/wmb_{model}_global_drift_{interval}.zarr", mode="w")

grid_dict = {}
models = exp_dict.keys()
for model in models:
    for interval in np.arange(1750, 2350, 5):
        grid = load_wmt_grid(model, interval=str(interval), dmget=True)
        save_wmb(grid, model, interval)
