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

Extra scripts:
- `05_report.R` rebuilds the Word report (`output/SCH302_sonar_catch_report.docx`). It re-runs 01–04 without plot titles, so the captions go below each figure.
- `06_school_animation.R` makes MP4 animations (`output/animation/`) of all smoothed school tracks, with every school starting at t = 0. It writes two views: a common start point, and the observed positions with movement exaggerated. It needs ffmpeg on the PATH and takes about 4 minutes per view.
- `07_vessel_school_animation.R` makes real-time MP4s (`output/animation/vessel_schools_*.mp4`): the vessel sails its track and schools appear at their real detection time and position. The `map` view shows the whole area; the `follow` view is centred on the vessel. To keep schools visible, their tracks play 10× slower than the clock and their movement is exaggerated (10× on the map, 5× when following the vessel). Rendering runs in parallel and takes about 3 minutes.
  To jump over stretches without schools near the vessel, e.g. for the follow view only:
  `anim_settings <- list(views = "follow", skip_no_school_m = 1000); source("07_vessel_school_animation.R")`
  This writes `vessel_schools_follow_skip1000m.mp4`.
