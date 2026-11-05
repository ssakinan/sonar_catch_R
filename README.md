# Sonar (LSSS Profos) x catch – R workflow

Two sonar datasets are supported as **profiles** (settings in `cfg` and `profiles` at the top of `00_setup_functions.R`):

| Profile | Sonar | Data | Trip | Outputs |
|---|---|---|---|---|
| `SU90` (default) | Simrad SU90 | `D:/SONAR/Reports/2025/profos` (23–29 Oct 2025) | 2025K, mackerel | `output/` |
| `ST90` | Simrad ST90 | `D:/SONAR/Reports/profos/ST90June2026` (28 Jun – 5 Jul 2026) | 2026E, herring | `output_ST90/` |
| `SU90_2026` | Simrad SU90 | `D:/SONAR/Reports/profos/SU90June2026` (same period) | 2026E, herring | `output_SU90_2026/` |

Choose the profile before sourcing a script, e.g. in the RStudio console:

```r
sonar_profile <- "ST90"
source("00b_school_qc.R")   # once, before the others
source("01_overview_maps.R")
```

Without `sonar_profile` (or the `SONAR_PROFILE` environment variable) the SU90 profile is used.

| Script | What it answers | Main outputs |
|---|---|---|
| `00_setup_functions.R` | settings, profiles, readers, effort, QC, smoothing, circular stats, animation helpers | – |
| `00b_school_qc.R` | which schools are unrealistic or outliers? (same rules for all profiles) | `school_qc.csv`, qc_a–qc_d |
| `01_overview_maps.R` | where/when did we fish and record sonar? | 01a trip map, 01b haul vs sonar timeline, 01c sonar zoom map |
| `02_sonar_catch_link.R` | what did the sonar see before/during each tow? | 02a timeline (+ one per day for multi-day data), 02b depth-vs-range, 02c haul zooms, 02d catch vs indices, `haul_sonar_indices.csv` |
| `03_school_track_smoothing.R` | denoise school tracks, speed + significance | 03a–03e, `school_tracks.rds`, `school_movement.csv`, `school_intervals.csv` |
| `04_school_movement_patterns.R` | population-level movement patterns | 04a vectors map, 04b roses, 04c speed covariates, 04d gridded field, 04e Figure-6 style |
| `05_report.R` | Word report for the chosen profile (QC, sonar vs catch, tracks, movement) | `output/SCH302_sonar_catch_report.docx` (SU90), `output_ST90/SCH302_ST90_sonar_catch_report.docx` (ST90) |
| `06_school_animation.R` | all smoothed tracks from a common start time | `animation/schools_*.mp4` |
| `07_vessel_school_animation.R` | vessel and schools in real time | `animation/vessel_schools_*.mp4` (one per day for multi-day data) |
| `08_daytime_migration.R` | is there a migration signal in daytime movement? (net drift after removing tide and vessel effect) | `migration/` figures 08a–08f, `08_migration_summary.png`, `08_daytime_flow.mp4`, `SCH302_ST90_daytime_migration.docx`, estimate tables |
| `09a_cmems_download.py` | download CMEMS North-West Shelf model currents (hourly + daily mean, 0–75 m) for the survey box | `output_ST90/cmems/*.nc` |
| `09_cmems_comparison.R` | is the daytime drift current or swimming? (model current at each school; swimming = ground − current; clock check against the model tide) | `migration/` figures 09a–09f, `09_cmems_summary.png`, `cmems_estimates.csv`, `SCH302_ST90_cmems_comparison.docx` |

| `10b_matched_schools.R` | schools tracked by both ST90 and SU90: overlap, beam geometry and depth bounds, Sv/size differences, velocity error per sonar | `output_compare/` figures 10g–10k, `10b_matched_summary.png`, `SCH302_matched_schools.docx` |
| `12_report_figures.R` | figures for chapter 4 of the pelAcoustics 2026 report | `output_compare/report_chapter4/` |
| `11_project_summary.R` | two-page summary of all analyses (run after 00b–10) | `SCH302_sonar_summary.docx` |
| `10_compare_sonars.R` | ST90 vs SU90 on the same period: detections, school properties, schools tracked by both, drift/swimming, catch (run 00b, 01–04, 08, 09 for both profiles first; for SU90_2026 run 09 with `CMEMS_LAG_H=1` to use the ST90 clock offset) | `output_compare/` figures 10a–10f, `10_sonar_comparison_summary.png`, `SCH302_ST90_vs_SU90.docx` |

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
- Catch-file times: equal to the sonar clock for trip 2025K; **1 h behind** the ST90 clock for trip 2026E (`catch_offset_h = 1` in the ST90 profile).
- On trip 2026E the SU90 clock runs **1 h 03 min 38 s behind** the ST90 (ship positions in the two exports align at this shift, stable over the trip); the SU90_2026 profile moves its times onto the ST90 clock (`clock_shift_s = 3818`). The SU90 and catch-log clocks are therefore probably near UTC and the ST90 about 1 h ahead.
- The ST90 export has no ppe file (track/effort only from pings with schools; indices per hour instead of per km), spaces in the File column (handled by the reader), a ~590 m range limit and only the forward ±90° sector.
- Haul ids such as `2026E-41` are read by `fread` as numbers unless `colClasses = list(character = "haul_id")` is given.
