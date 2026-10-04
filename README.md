# Sonar (LSSS Profos) x catch – R workflow

Run in order (paths/settings in `cfg` at the top of `00_setup_functions.R`):

| Script | What it answers | Main outputs |
|---|---|---|
| `00_setup_functions.R` | readers, effort, smoothing, circular stats | – |
| `01_overview_maps.R` | where/when did we fish and record sonar? | 01a trip map, 01b haul vs sonar timeline, 01c sonar zoom map |
| `02_sonar_catch_link.R` | what did the sonar see before/during each tow? | 02a timeline, 02b depth-vs-range, 02c haul zooms, 02d catch vs indices, `haul_sonar_indices.csv` |
| `03_school_track_smoothing.R` | denoise school tracks, speed + significance | 03a–03e, `school_tracks.rds`, `school_movement.csv`, `school_intervals.csv` |
| `04_school_movement_patterns.R` | population-level movement patterns | 04a vectors map, 04b roses, 04c speed covariates, 04d gridded field, 04e Figure-6 style |

To scale up: export Profos pp/ps/ppe for the whole trip into `cfg$profos_dir`
(subfolders are read recursively). Script 02 then gives one row per haul.

Caveats found in the October 2025 sample:
- `Center.dep` = range x sin(|tilt|) exactly (R² = 1). It is the beam-centre depth, not a measured school depth.
- Raw ping-to-ping centroid noise is ~4–5 m per axis (it grows with range). Raw speeds are meaningless; smoothed speeds stabilise for tau ≥ 20–30 s.
- Catch-file times are UTC. Profos `Ship.speed` is in m/s.
