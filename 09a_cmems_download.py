"""
09a_cmems_download.py

Download CMEMS North-West Shelf model currents for the ST90 survey box
(NWSHELF_ANALYSISFORECAST_PHY_004_013, 1.5 km, hourly instantaneous 3-D uo/vo,
plus the daily-mean 3-D field for the tide-free residual current).

Needs the Copernicus Marine toolbox and a saved login:
    python -m pip install --user copernicusmarine
    copernicusmarine login          # once; stores credentials in ~/.copernicusmarine
Run from the project root:
    python 09a_cmems_download.py
Files go to output_ST90/cmems/ (git-ignored).
"""
import os
import copernicusmarine as cm

out = os.path.join("output_ST90", "cmems")
os.makedirs(out, exist_ok=True)

box = dict(minimum_longitude=-1.6, maximum_longitude=1.1,
           minimum_latitude=58.8, maximum_latitude=60.15,
           minimum_depth=0, maximum_depth=80)
# daytime schools: 28 Jun 03:00 - 5 Jul 21:00 (sonar clock); pad 3 h for clock offsets
t0, t1 = "2026-06-27T22:00:00", "2026-07-06T02:00:00"

jobs = [("cmems_mod_nws_phy-cur_anfc_1.5km-3D_PT1H-i", "nws_cur_hourly.nc", ["uo", "vo"]),
        ("cmems_mod_nws_phy-cur_anfc_1.5km-3D_P1D-m", "nws_cur_daily.nc", ["uo", "vo"])]
for ds, fn, var in jobs:
    if os.path.exists(os.path.join(out, fn)):
        print("exists:", fn)
        continue
    cm.subset(dataset_id=ds, variables=var, start_datetime=t0, end_datetime=t1,
              output_directory=out, output_filename=fn, disable_progress_bar=True, **box)
    print("downloaded:", fn)
