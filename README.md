# Sonar (LSSS Profos) x catch – R workflow

Two sonar datasets are supported as **profiles** (settings in `cfg` and `profiles` at the top of `00_setup_functions.R`):

| Profile | Sonar | Data | Trip | Outputs |
|---|---|---|---|---|
| `SU90` (default) | Simrad SU90 | `C:/Drive/Sonar/October2025/profos` (23 Oct 2025) | 2025K, mackerel | `output/` |
| `ST90` | Simrad ST90 | `D:/SONAR/Reports/profos/ST90June2026` (28 Jun – 5 Jul 2026) | 2026E, herring | `output_ST90/` |

Choose the profile before sourcing a script, e.g. in the RStudio console:

```r
sonar_profile <- "ST90"
source("00b_school_qc.R")   # once, before the others (ST90 only)
source("01_overview_maps.R")
```

Without `sonar_profile` (or the `SONAR_PROFILE` environment variable) the SU90 profile is used.

| Script | What it answers | Main outputs |
|---|---|---|
| `00_setup_functions.R` | settings, profiles, readers, effort, QC, smoothing, circular stats, animation helpers | – |
| `00b_school_qc.R` | which schools are unrealistic or outliers? (profiles with `cfg$qc`, i.e. ST90) | `school_qc.csv`, qc_a–qc_d |
| `01_overview_maps.R` | where/when did we fish and record sonar? | 01a trip map, 01b haul vs sonar timeline, 01c sonar zoom map |
| `02_sonar_catch_link.R` | what did the sonar see before/during each tow? | 02a timeline (+ one per day for multi-day data), 02b depth-vs-range, 02c haul zooms, 02d catch vs indices, `haul_sonar_indices.csv` |
| `03_school_track_smoothing.R` | denoise school tracks, speed + significance | 03a–03e, `school_tracks.rds`, `school_movement.csv`, `school_intervals.csv` |
| `04_school_movement_patterns.R` | population-level movement patterns | 04a vectors map, 04b roses, 04c speed covariates, 04d gridded field, 04e Figure-6 style |
| `05_report.R` | Word report, SU90 | `output/SCH302_sonar_catch_report.docx` |
| `05b_report_ST90.R` | Word report, ST90 (incl. QC) | `output_ST90/SCH302_ST90_sonar_catch_report.docx` |
| `06_school_animation.R` | all smoothed tracks from a common start time | `animation/schools_*.mp4` |
| `07_vessel_school_animation.R` | vessel and schools in real time | `animation/vessel_schools_*.mp4` (one per day for multi-day data) |
| `08_daytime_migration.R` | is there a migration signal in daytime movement? (net drift after removing tide and vessel effect) | `migration/` figures 08a–08f, `08_migration_summary.png`, `08_daytime_flow.mp4`, `SCH302_ST90_daytime_migration.docx`, estimate tables |
| `09a_cmems_download.py` | download CMEMS North-West Shelf model currents (hourly + daily mean, 0–75 m) for the survey box | `output_ST90/cmems/*.nc` |
| `09_cmems_comparison.R` | is the daytime drift current or swimming? (model current at each school; swimming = ground − current; clock check against the model tide) | `migration/` figures 09a–09f, `09_cmems_summary.png`, `cmems_estimates.csv`, `SCH302_ST90_cmems_comparison.docx` |

The CMEMS download needs the Copernicus Marine toolbox and a saved login (`python -m pip install --user copernicusmarine`, then `copernicusmarine login` once).

The report scripts re-run the analysis scripts without plot titles, so the captions go below each figure (`SKIP_FIGURES=1` reuses existing figures).

## School quality control (ST90)

`00b_school_qc.R` drops pings with Sv < −70 dB or tilt shallower than −3°, then excludes schools that are:
faster than 2.5 m/s (smoothed), moving with the ship (heading within 20° of the ship's, > 0.5 × ship speed while steaming > 3 m/s),
showing tracker jumps (> 30 m/s between cleaned pings), closer than 150 m on average, weaker than −60 dB, or extreme in size.
Thresholds are in `profiles$ST90$qc`; the flags per school are in `output_ST90/school_qc.csv`.
For the June–July 2026 data this excluded 642 of 6,000 schools; the speed distribution had a second mode at the
vessel's steaming speed (6–7 m/s), i.e. echoes tied to the ship.

## Animations

Settings can be overridden without editing the scripts, e.g.

```r
anim_settings <- list(views = "follow", skip_no_school_m = 1000, clock_step_s = 30)
source("07_vessel_school_animation.R")
```

- `06`: every school starts at t = 0; views `origin` (common start point) and `map` (observed positions, movement exaggerated). Length = 99th percentile of track duration.
- `07`: the vessel sails its track and schools appear at their real detection time and position; school tracks play 10× slower than the clock and movement is exaggerated (10× map, 5× follow). `skip_no_school_m` jumps over stretches without schools near the vessel; multi-day data are split into one video per day.
- Both need ffmpeg on the PATH and render frames in parallel.

## Caveats found in the data

- `Center.dep` = range × sin(|tilt|) exactly (R² = 1) for both sonars: it is the beam-centre depth, not a measured school depth.
- Raw ping-to-ping centroid noise is ~3–5 m per axis; raw speeds are meaningless, smoothed speeds stabilise for tau ≥ 20–30 s.
- Profos `Ship.speed` is in m/s.
- Catch-file times: equal to the sonar clock for trip 2025K; **1 h behind** the sonar clock for trip 2026E (`catch_offset_h = 1` in the ST90 profile).
- The ST90 export has no ppe file (track/effort only from pings with schools; indices per hour instead of per km), spaces in the File column (handled by the reader), a ~590 m range limit and only the forward ±90° sector.
- Haul ids such as `2026E-41` are read by `fread` as numbers unless `colClasses = list(character = "haul_id")` is given.
