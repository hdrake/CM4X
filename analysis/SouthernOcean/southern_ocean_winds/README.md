Notebooks for analysis of southern ocean winds in GFDL models as contribution to the high-res model development team.

Authors: Katherine Turner - University of Arizona / GFDL

CONTENTS:

-- create_zarrs.ipynb -- preprocessing for 5-year netcdf files into appropriately-chunked zarr files

-- westerlies_piC.ipynb -- comparison of average preindustrial control westerly fields (odiv-1, odiv-210 (?), odiv-209)

-- westerlies_hist.ipynb -- comparison of average satellite-era westerly fields (odiv-2, odiv-231, odiv-255) against ERA5 wind fields

    ! Note ! odiv-255 has been run up to 1960 as of 2 August 2023 -- need to update files as more years become available

-- seasonal_hist.ipynb -- as in westerlies_hist.ipynb, but now taking seasonal averages

-- moments.ipynb -- using daily fields for zonal winds to show standard deviations (and potentially higher moments?) as some sort of proxy for storminess?
