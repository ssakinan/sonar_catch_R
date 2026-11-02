"""
09a_cmems_download.py

Download CMEMS North-West Shelf model currents for a survey box
(NWSHELF_ANALYSISFORECAST_PHY_004_013, 1.5 km, hourly instantaneous 3-D uo/vo,
plus the daily-mean 3-D field for the tide-free residual current).

Needs the Copernicus Marine toolbox and a saved login:
    python -m pip install --user copernicusmarine
    copernicusmarine login          # once; stores credentials in ~/.copernicusmarine
Run from the project root, with no arguments for the ST90 June-July 2026 box:
    python 09a_cmems_download.py
or for another dataset (output folder, start, end, lon min/max, lat min/max, max depth):
    python 09a_cmems_download.py output 2025-10-22T22:00:00 2025-10-29T12:00:00 -0.2 1.3 59.8 60.85 150
Files go to <output folder>/cmems/ (git-ignored).
"""
import os
import sys
import copernicusmarine as cm

if len(sys.argv) > 1:
    base, t0, t1 = sys.argv[1:4]
    lon0, lon1, lat0, lat1, zmax = map(float, sys.argv[4:9])
else:
    # ST90 daytime schools: 28 Jun 03:00 - 5 Jul 21:00 (sonar clock); pad 3 h for clock offsets
    base, t0, t1 = "output_ST90", "2026-06-27T22:00:00", "2026-07-06T02:00:00"
    lon0, lon1, lat0, lat1, zmax = -1.6, 1.1, 58.8, 60.15, 80

out = os.path.join(base, "cmems")
os.makedirs(out, exist_ok=True)
box = dict(minimum_longitude=lon0, maximum_longitude=lon1, minimum_latitude=lat0, maximum_latitude=lat1,
           minimum_depth=0, maximum_depth=zmax)

jobs = [("cmems_mod_nws_phy-cur_anfc_1.5km-3D_PT1H-i", "nws_cur_hourly.nc", ["uo", "vo"]),
        ("cmems_mod_nws_phy-cur_anfc_1.5km-3D_P1D-m", "nws_cur_daily.nc", ["uo", "vo"])]
for ds, fn, var in jobs:
    if os.path.exists(os.path.join(out, fn)):
        print("exists:", fn)
        continue
    cm.subset(dataset_id=ds, variables=var, start_datetime=t0, end_datetime=t1,
              output_directory=out, output_filename=fn, disable_progress_bar=True, **box)
    print("downloaded:", fn)
