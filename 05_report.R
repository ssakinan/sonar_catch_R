################################################################################
# 05_report.R
#
# Builds a Word report (output/SCH302_sonar_catch_report.docx) summarising the
# analyses of scripts 01-04, with every figure followed by a caption paragraph.
#
#   1. Re-runs 01-04 in "report mode" (SONAR_REPORT=1): identical figures but
#      without titles/subtitles, written to output/report_figs/.
#   2. Assembles the document with officer + flextable.
#
# Numbers in the text are computed from the outputs. The interpretation was
# written for the 23 Oct 2025 sample (hauls 7-8) -- revise the wording when
# Profos exports for more days are added.
################################################################################

rerun_figures <- !identical(Sys.getenv("SKIP_FIGURES"), "1")   # set SKIP_FIGURES=1 to reuse figures

if (rerun_figures) {
  Sys.setenv(SONAR_REPORT = "1")
  for (s in c("01_overview_maps.R", "02_sonar_catch_link.R",
              "03_school_track_smoothing.R", "04_school_movement_patterns.R")) {
    message("Running ", s, " in report mode ...")
    source(s, local = new.env())
  }
}
Sys.setenv(SONAR_REPORT = "1")
source("00_setup_functions.R")      # sets cfg$out_dir to output/report_figs
library(officer)
library(flextable)

fig_dir  <- cfg$out_dir
out_docx <- file.path(dirname(fig_dir), "SCH302_sonar_catch_report.docx")

# ---- numbers for the text ---------------------------------------------------------
hauls <- read_catch(cfg$catch_file, trips = cfg$trip)
ship  <- ship_track(read_profos(cfg$profos_dir, "ppe"))
pp    <- read_profos(cfg$profos_dir, "pp")
ind   <- fread(file.path(fig_dir, "haul_sonar_indices.csv"))
mov   <- fread(file.path(fig_dir, "school_movement.csv"))
cst   <- fread(file.path(fig_dir, "movement_circular_stats.csv"))
mm    <- mov[moving == TRUE]

f1 <- function(x) formatC(x, format = "f", digits = 1, big.mark = ",")
f2 <- function(x) formatC(x, format = "f", digits = 2)
f0 <- function(x) formatC(x, format = "d", big.mark = ",")
pv <- function(p) if (p < 0.001) "p < 0.001" else paste("p =", formatC(p, digits = 2, format = "g"))
cs <- function(g) cst[group == g]

n_sch_all  <- uniqueN(pp[!is.na(Id)]$sid)
rec_h      <- sum(ship$dt, na.rm = TRUE) / 3600
dist_km    <- sum(ship$dist_m) / 1000
sp_by_sun  <- mm[, .(med = median(speed_lm), n = .N), by = sun]
w_day_night <- wilcox.test(speed_lm ~ sun, mm[sun %in% c("day", "night")], exact = FALSE)$p.value
rho_area   <- suppressWarnings(cor.test(log(mm$Area), mm$speed_lm, method = "spearman"))
rho_range  <- suppressWarnings(cor.test(mm$range_m, mm$speed_lm, method = "spearman"))
mov[, v_along := vx * sin(ship_bearing * pi / 180) + vy * cos(ship_bearing * pi / 180)]
fit_art    <- summary(lm(v_along ~ ship_speed, mov[moving == TRUE]))$coefficients
raw_step   <- {
  m <- merge_pp(pp)
  m[, v := c(NA, sqrt(diff(to_xy(lon, lat)[, 1])^2 + diff(to_xy(lon, lat)[, 2])^2) /
               diff(as.numeric(t))), by = sid]
  median(m$v[is.finite(m$v)])
}
noise_m <- {   # same estimator as script 03
  m <- merge_pp(pp)[sid %in% mov$sid]
  nz <- m[, {xy <- to_xy(lon, lat); .(sx = noise_sd(xy[, 1]), sy = noise_sd(xy[, 2]))}, by = sid]
  median(c(nz$sx, nz$sy), na.rm = TRUE)
}
fits_tilt <- merge_pp(pp)[, .(r2 = summary(lm(depth ~ range_m))$r.squared,
                              max_depth = max(depth)), by = tilt]
hz <- hauls[haul_id %in% ind$haul_id]
spd <- function(s) f2(sp_by_sun[sun == s, med])

# ---- document helpers ---------------------------------------------------------------
fp_cap_lab <- fp_text(bold = TRUE, font.size = 9.5)
fp_cap     <- fp_text(font.size = 9.5)
fp_bold    <- fp_text(bold = TRUE)
fig_no <- 0; tab_no <- 0

doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(orient = "portrait", width = 21 / 2.54, height = 29.7 / 2.54),
  page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1)))

h1 <- function(txt) doc <<- body_add_par(doc, txt, style = "heading 1")
h2 <- function(txt) doc <<- body_add_par(doc, txt, style = "heading 2")
para <- function(...) doc <<- body_add_par(doc, paste0(...), style = "Normal")
lead <- function(head, ...) doc <<- body_add_fpar(doc, fpar(ftext(paste0(head, " "), fp_bold),
                                                            ftext(paste0(...))), style = "Normal")

add_fig <- function(file, caption, max_w = 6.25, max_h = 8.3) {
  fig_no <<- fig_no + 1
  d <- dim(png::readPNG(file.path(fig_dir, file)))
  w <- max_w; h <- w * d[1] / d[2]
  if (h > max_h) { h <- max_h; w <- h * d[2] / d[1] }
  doc <<- body_add_fpar(doc, fpar(external_img(file.path(fig_dir, file), width = w, height = h),
                                  fp_p = fp_par(text.align = "center", keep_with_next = TRUE,
                                                padding.top = 6)))
  doc <<- body_add_fpar(doc, fpar(ftext(sprintf("Figure %d. ", fig_no), fp_cap_lab),
                                  ftext(caption, fp_cap),
                                  fp_p = fp_par(padding.bottom = 12)),
                        style = "Image Caption")
  invisible(fig_no)
}
add_table_caption <- function(caption) {
  tab_no <<- tab_no + 1
  doc <<- body_add_fpar(doc, fpar(ftext(sprintf("Table %d. ", tab_no), fp_cap_lab),
                                  ftext(caption, fp_cap),
                                  fp_p = fp_par(keep_with_next = TRUE, padding.top = 6)),
                        style = "Table Caption")
}

# ---- title -------------------------------------------------------------------------------
doc <- body_add_fpar(doc, fpar(ftext("Omnidirectional sonar and midwater trawl catches on SCH302",
                                     fp_text(bold = TRUE, font.size = 20))))
doc <- body_add_fpar(doc, fpar(ftext(sprintf(
  "Exploratory analysis of LSSS Profos school detections, trip %s (sonar sample %s)",
  cfg$trip, format(min(ship$t), "%d %B %Y")), fp_text(font.size = 13, color = "grey30"))))
doc <- body_add_fpar(doc, fpar(ftext(paste("Report generated", format(Sys.Date(), "%d %B %Y"),
                                           "with the sonar_catch_R scripts"),
                                     fp_text(font.size = 9.5, color = "grey40"))))

# ---- 1 summary ----------------------------------------------------------------------------
h1("1. Summary")
para("This report links school detections from the vessel's omnidirectional sonar (processed with ",
     "the LSSS Profos school tracker) to the catch of the midwater trawl hauls, and describes how the ",
     "detected schools moved. It uses the Profos exports currently available for ",
     format(min(ship$t), "%d %B %Y"), ", which cover hauls ", paste(hz$haul, collapse = " and "), ". ",
     "The workflow is written for the whole trip: adding further Profos exports updates every ",
     "figure and table.")
lead("Profos depth is not a measured depth.",
     "The school depth reported by Profos equals range x sin(tilt) exactly (R² = ",
     paste(formatC(fits_tilt$r2, format = "f", digits = 4), collapse = " and "),
     " for the two tilts), i.e. the depth of the beam centre at the school's range. ",
     "At the tilts used (", paste(sort(unique(fits_tilt$tilt)), collapse = "° and "),
     "°) the sonar cannot observe schools deeper than about ",
     paste(round(fits_tilt$max_depth), collapse = "–"), " m, while the hauls fished from ",
     paste(range(hz$headline_depth), collapse = "–"), " m. Sonar school depth can therefore not be ",
     "compared with fishing depth; the echosounder is needed for that.")
lead("Raw school positions are noisy; smoothing is essential.",
     "The centroid of a tracked school jitters by about ", f1(noise_m), " m per axis from ping to ping ",
     "(more at long range), which makes raw ping-to-ping speeds meaningless (median ", f1(raw_step),
     " m/s). Smoothing over 20–30 s or more gives stable speeds; the median school speed is ",
     f2(median(mov$speed_lm)), " m/s, consistent with the LSSS estimate (r = ",
     f2(cor(mov$lsss_speed, mov$speed_lm, use = "complete")), ").")
lead("Most schools were really moving.",
     f0(nrow(mm)), " of ", f0(nrow(mov)), " tracked schools (", round(100 * nrow(mm) / nrow(mov)),
     "%) had a speed significantly above zero. Schools were faster at night (median ", spd("night"),
     " m/s) than by day (", spd("day"), " m/s; Wilcoxon ", pv(w_day_night), ") and larger schools ",
     "swam faster (Spearman ρ = ", f2(rho_area$estimate), ", ", pv(rho_area$p.value), ").")
lead("Directions were local, not general.",
     "Across all schools there was no common swimming direction (R̄ = ", f2(cs("all")$Rbar), ", ",
     pv(cs("all")$rayleigh_p), "), but schools near each haul shared a direction: ",
     paste(sprintf("haul %s %.0f° (R̄ = %s, %s)", sub(".*-", "", sub("haul: ", "", cst[grepl("^haul", group)]$group)),
                   cst[grepl("^haul", group)]$mean_dir, f2(cst[grepl("^haul", group)]$Rbar),
                   sapply(cst[grepl("^haul", group)]$rayleigh_p, pv)), collapse = "; "), ".")
lead("Sonar and catch.",
     "Both hauls were pure mackerel (", paste(round(hz$catch), collapse = " and "), " t). ",
     "School encounter rates were high while towing in both hauls (",
     paste(f1(ind[phase == "tow", schools_per_10km]), collapse = " and "),
     " schools per 10 km sailed). With two hauls no relationship between sonar indices and catch can ",
     "be tested; the per-haul indices (Table 1) are produced so that this can be done once the whole ",
     "trip is exported.")

# ---- 2 data -------------------------------------------------------------------------------
h1("2. Data")
para("Catch data: ", nrow(hauls), " hauls of trip ", cfg$trip, " (", format(min(hauls$shoot_time), "%d %b"),
     " – ", format(max(hauls$haul_time), "%d %b %Y"), "), ", f0(round(sum(hauls$catch))), " t in total. ",
     "The dominant species was mackerel in ", sum(hauls$dom_species == "mac"), " hauls, blue whiting in ",
     sum(hauls$dom_species == "whb"), " and herring in ", sum(hauls$dom_species == "her"), ". ",
     "Shoot and haul times are in UTC: the sonar ship track passes the logged shoot and haul ",
     "positions with no time offset, within the one-minute rounding of the positions.")
para("Sonar data: Profos exports for ", format(min(ship$t), "%d %B %Y"), " (", f1(rec_h), " h of recording, ",
     f0(round(dist_km)), " km sailed, ", f0(nrow(ship)), " pings, ", f0(n_sch_all), " schools with ",
     f0(nrow(pp[!is.na(Id)])), " detections). Three file types are used: pp (one row per school per ping), ",
     "ps (one row per school, the LSSS summary) and ppe (one row per ping including pings without ",
     "schools, used as the ship track and sonar effort). Ship speed in the Profos files is in m/s.")

# ---- 3 methods ----------------------------------------------------------------------------
h1("3. Methods")
h2("3.1 Linking sonar to hauls")
para("Each ping and school detection was assigned to a haul phase: search (the ", cfg$pre_shoot_h,
     " h before shoot) or tow (shoot to haul). Schools were included if they were within ",
     cfg$buffer_km, " km of the ship track during the tow. Sonar effort was taken from the ppe file as the ",
     "number of pings, the recording time (gaps longer than ", cfg$gap_s, " s count as sonar off) and ",
     "the distance sailed while recording. Per haul and phase the following indices were computed: ",
     "number of distinct schools, schools per 10 km sailed, a school backscatter index (10 log10 of the ",
     "summed linear Sv x school area per ping; relative, uncalibrated), the backscatter-weighted school ",
     "depth, and the fraction of school backscatter inside the trawl layer (headline depth to headline ",
     "+ ", cfg$default_vert_opening, " m; the vertical opening is not recorded in the catch file).")
h2("3.2 Cleaning and smoothing school tracks")
para("The approach of Figure_6_R_code.R was extended from one school to all schools. Schools with fewer than ",
     cfg$min_pings, " pings or a track shorter than ", cfg$min_duration, " s were excluded. When Profos ",
     "reported the same school twice in one ping (a split school), the detections were merged into one ",
     "area-weighted position instead of adding random jitter to the time stamps. Sudden tracker jumps were ",
     "removed with a two-dimensional Hampel filter (running median over ", 2 * cfg$hampel_k + 1,
     " pings; points more than ", cfg$hampel_nsig, " robust SD and at least ", cfg$hampel_floor,
     " m from it were removed). The positioning noise was estimated from second differences of the ",
     "positions, which cancel constant swimming velocity.")
para("Two smoothers were compared: a smoothing spline with about one degree of freedom per ", cfg$tau_s,
     " s (a physical time scale, rather than the fixed lambda of Figure_6, whose effect depends on track ",
     "length), and a constant-velocity Kalman filter with Rauch–Tung–Striebel smoother (process noise ",
     cfg$sigma_acc, " m/s², measurement noise as estimated), which also gives the uncertainty of the velocity.")
h2("3.3 Movement metrics")
para("For each school the velocity was estimated by linear regression of the cleaned positions on time. ",
     "Its standard error was inflated for lag-1 autocorrelation of the residuals, and a school was classed as ",
     "moving when its speed exceeded twice this standard error. Directions are compass bearings (0° = north, ",
     "clockwise). Directions were summarised with circular statistics: mean direction, mean resultant length ",
     "R̄ (0 = no common direction, 1 = identical directions) and the Rayleigh test of uniformity. Light ",
     "conditions were derived from solar elevation (day > 0°, twilight 0 to −6°, night < −6°). ",
     "Velocities over ", cfg$interval_s, "-s intervals of the smoothed tracks reproduce the Figure 6 summary for all schools.")

# ---- 4 results ----------------------------------------------------------------------------
doc <- body_add_break(doc)
h1("4. Results")
h2("4.1 Fishing activity and sonar coverage")
para("The trip fished mackerel east of Shetland and east of Orkney in late October and blue whiting ",
     "west of the Hebrides in early November (Figure 1). Profos data have so far been exported for one day, ",
     "which overlaps hauls 7 and 8 (Figure 2). On that day the vessel searched and towed in an area of about ",
     "20 x 25 km where schools were detected throughout (Figure 3).")
add_fig("01a_trip_map.png", paste(
  "Fishing positions of trip", cfg$trip, "(SCH302). Lines join the shoot and haul positions; circles at",
  "the haul position are sized by catch (t) and coloured by dominant species (mac = mackerel, whb = blue",
  "whiting, her = herring); numbers are haul numbers. The short black line east of Shetland is the ship",
  "track for the period with exported sonar data."))
add_fig("01b_trip_timeline.png", paste(
  "Hauls of trip", cfg$trip, "over time. Each bar spans shoot to haul time and is placed at the haul's",
  "catch (t); colour = dominant species, numbers = haul numbers. Yellow bands mark periods with exported",
  "Profos sonar data."))
add_fig("01c_sonar_zoom_map.png", paste0(
  "Ship track and detected schools on ", format(min(ship$t), "%d %B %Y"), ". Track colour: orange = search ",
  "(", cfg$pre_shoot_h, " h before a shoot), red = tow (shoot to haul), grey = other. Circles are the mean ",
  "position of each school track, sized by mean school area (m²) and coloured by Profos depth (see ",
  "Section 4.2 for what this depth represents). Red arrows join the logged shoot and haul positions; ",
  "labels give haul number, dominant species and catch. Map cropped at ", cfg$map_lat_min, "°N."))

h2("4.2 What the sonar saw before and during the tows")
para("School detections occurred in bursts that follow the ship's track through aggregations rather than ",
     "a steady rate (Figure 4). For haul 7 schools were detected both while searching (",
     f1(ind[haul == 7 & phase == "search", schools_per_10km]), " per 10 km) and while towing (",
     f1(ind[haul == 7 & phase == "tow", schools_per_10km]), " per 10 km). Before haul 8 few schools were ",
     "seen within the buffer of the later tow track (", f1(ind[haul == 8 & phase == "search", schools_per_10km]),
     " per 10 km), whereas during the tow the rate was the highest of the day (",
     f1(ind[haul == 8 & phase == "tow", schools_per_10km]), " per 10 km; Table 1, Figure 6).")
para("Figure 5 shows why the Profos depths must be interpreted with care: every detection lies on the ",
     "beam-centre line of its tilt. The deeper schools seen around haul 8 are therefore schools at longer ",
     "range with the −7° tilt rather than evidence of schools moving deeper. Both hauls were fished below ",
     "the depth the sonar covers at these tilts, so the fraction of school backscatter in the trawl layer ",
     "is zero by construction. The same applies, more strongly, to the blue whiting hauls (headline depths ",
     paste(range(hauls[dom_species == "whb", headline_depth]), collapse = "–"), " m).")
add_fig("02a_timeline_schools_vs_tows.png", paste0(
  "School detections against time. Upper panel: each school detection (one per ping) by time and Profos ",
  "depth, coloured by mean volume backscattering strength Sv (dB) and sized by school area (m²). Shading: ",
  "orange = ", cfg$pre_shoot_h, " h search period before shoot, red = tow. Dashed boxes mark the trawl layer ",
  "from headline depth to headline + ", cfg$default_vert_opening, " m (assumed vertical opening). Labels give ",
  "haul number, catch and the share of the dominant species. Lower panel: distinct schools per hour of ",
  "sonar recording in 10-min bins (bars) and ship speed (line, right axis)."))
add_fig("02b_depth_vs_range_tilt.png", paste0(
  "Profos school depth against horizontal range from the ship, by transducer tilt. The detections fall on ",
  "straight lines (R² = 1.000): Profos depth = range x sin(|tilt|). Solid lines = beam-centre depth with a ",
  cfg$transducer_depth, " m transducer depth; shaded bands = ±", cfg$vert_beam_half, "° vertical beam ",
  "(approximate); dashed lines = headline depths of the hauls."))
add_fig("02c_haul_zoom_maps.png", paste0(
  "Schools within ", cfg$buffer_km, " km of the tow track for each haul with sonar data. Orange = ship track ",
  "in the ", cfg$pre_shoot_h, " h before shoot, red = ship track while towing (gaps = no sonar recording), ",
  "dashed = straight line between logged shoot and haul positions. Circles = school mean positions, sized ",
  "by school area and coloured by Profos depth."))

add_table_caption(paste0(
  "Sonar indices per haul and phase (schools within ", cfg$buffer_km, " km of the tow track). Rec. h = ",
  "hours of sonar recording; km = distance sailed while recording; SBI = school backscatter index ",
  "(relative, dB); depth = backscatter-weighted Profos depth (beam-centre depth, see Figure 5)."))
tab <- ind[, .(Haul = haul, Phase = phase, `Catch (t)` = round(catch), Species = dom_species,
               `Rec. h` = round(rec_h, 1), km = round(dist_km, 1), Schools = n_schools,
               `Schools / 10 km` = round(schools_per_10km, 1), `SBI (dB)` = round(sbi_db, 1),
               `Depth (m)` = round(depth_w), `Ship speed (m/s)` = round(ship_speed, 1))]
ft <- flextable(tab)
ft <- theme_booktabs(ft)
ft <- fontsize(ft, size = 9, part = "all")
ft <- bold(ft, part = "header")
ft <- align(ft, align = "center", part = "all")
ft <- width(ft, width = c(0.45, 0.6, 0.6, 0.6, 0.5, 0.5, 0.6, 0.7, 0.6, 0.6, 0.75))
doc <- body_add_flextable(doc, ft, align = "center")
para("")

h2("4.3 Cleaning and smoothing school tracks")
para(f0(nrow(mov)), " of ", f0(n_sch_all), " schools passed the quality criteria. The raw centroids scatter ",
     "around a smooth path (Figure 7); the ", f0(sum(mov$n_outliers)), " outlying detections removed by the ",
     "Hampel filter are mostly tracker jumps at the start or end of a track. The spline and Kalman ",
     "smoothers give nearly identical paths. Raw ping-to-ping speeds fluctuate between 0 and more than ",
     "5 m/s, while the smoothed speeds vary slowly (Figure 8).")
para("Positioning noise increases with range (Figure 9, left). Mean speed falls steeply as the smoothing ",
     "time scale grows from 5 to 20 s and is nearly constant beyond 30 s (Figure 9, right): this plateau is ",
     "where noise has been removed without smoothing away real movement, and supports the ",
     cfg$tau_s, "-s default. Short tracks cannot resolve speed reliably (Figure 10, left); our speeds agree ",
     "with the LSSS speeds in the ps file (Figure 10, right).")
para("Relative to the ship, schools tended to move ahead and to starboard (Figure 11). The component of ",
     "school velocity along the ship's heading did not depend significantly on ship speed (slope ",
     f2(fit_art[2, 1]), ", ", pv(fit_art[2, 4]), "), so there is no clear sign of a positioning bias, but ",
     "avoidance of the vessel cannot be excluded.")
add_fig("03a_raw_vs_smoothed_tracks.png", paste(
  "Raw and smoothed tracks of the six longest-tracked schools (local coordinates, m; panel headers give",
  "the Profos run and school Id). Dark grey points = raw centroids after merging same-ping duplicates; red",
  "points = outliers removed; blue line = smoothing spline; red line = Kalman smoother; triangle = first",
  "detection."))
add_fig("03b_speed_raw_vs_smoothed.png", paste(
  "Speed over time for the same six schools. Grey = raw ping-to-ping speed; blue = spline speed; red =",
  "Kalman speed with ±2 SD band. The y-axis is truncated at 6 m/s."))
add_fig("03c_noise_and_smoothing_sensitivity.png", paste(
  "Left: centroid positioning noise per school (from second differences) against mean range from the",
  "ship, with loess fit and 95% confidence band. Right: mean school speed against the spline smoothing",
  "time scale τ (median and interquartile range over schools; log x-axis); dashed line = median",
  "straight-line (regression) speed."))
add_fig("03d_movement_significance.png", paste(
  "Left: straight-line speed of each school against track duration (log scale) with ±2 SE error bars;",
  "black = speed significantly above zero. Right: speed from this analysis against the LSSS Profos",
  "school speed (ps file); dashed = 1:1 line."))
add_fig("03e_ship_motion_artefact_check.png", paste(
  "Check for ship-related artefacts, moving schools only. Left: school velocity component along the",
  "ship's heading against ship speed, with linear fit and 95% confidence band. Right: swimming direction",
  "relative to the ship's heading."))

h2("4.4 School movement patterns")
para("Moving schools swam in all directions over the area as a whole (Figure 12; Figure 13: R̄ = ",
     f2(cs("all")$Rbar), ", ", pv(cs("all")$rayleigh_p), "). Directions were more consistent during tows (",
     round(cs("phase: tow")$mean_dir), "°, R̄ = ", f2(cs("phase: tow")$Rbar), ", ", pv(cs("phase: tow")$rayleigh_p),
     ") and around twilight (", round(cs("sun: twilight")$mean_dir), "°, R̄ = ", f2(cs("sun: twilight")$Rbar),
     ", n = ", cs("sun: twilight")$n, "), and within the haul areas (Section 1). The gridded field (Figure 15) ",
     "shows that neighbouring schools often moved the same way while different parts of the area showed ",
     "different directions, i.e. local movement rather than a common migration.")
para("Speed was highest before dawn and lowest in the middle of the day (Figure 14, bottom; night median ",
     spd("night"), " m/s, day ", spd("day"), " m/s). Larger schools were faster, whereas speed was not ",
     "related to range from the ship (Spearman ρ = ", f2(rho_range$estimate), ", ", pv(rho_range$p.value),
     "), which suggests that the speed estimates are not driven by range-dependent noise. Light condition ",
     "and fishing phase are confounded on a single day, so the day–night difference needs confirming ",
     "with more days.")
add_fig("04a_movement_vector_map.png", paste0(
  "Movement of schools with significant speed (n = ", nrow(mm), "). Arrows start at each school's mean ",
  "position and show the distance swum in 10 min at its estimated velocity, coloured by speed; dark lines = ",
  "Kalman-smoothed tracks; grey = ship track; red dashed = logged tow lines. Cropped at ", cfg$map_lat_min, "°N."))
add_fig("04b_direction_roses.png", paste(
  "Swimming directions (compass bearings, 22.5° bins) of moving schools by ship activity (top: search,",
  "tow, other) and light condition (bottom: day, twilight, night). Panel headers give the number of",
  "schools, mean direction, mean resultant length R̄ (0 = no common direction, 1 = identical directions)",
  "and the Rayleigh test p-value."))
add_fig("04c_speed_vs_covariates.png", paste(
  "Speed of moving schools against range from the ship (top left), school area (top right, log scale) and",
  "time of day (bottom; red = tows), coloured by light condition; black lines = loess fits with 95%",
  "confidence bands."))
add_fig("04d_gridded_movement_field.png", paste0(
  "Mean school movement per ", 1, " km grid cell. Fill = number of moving schools in the cell; arrow = mean ",
  "velocity (relative length, longest arrow = 40% of a cell); arrow thickness = R̄, i.e. how consistently the ",
  "schools in a cell moved in the same direction (thin = one school or no agreement); grey = ship track. ",
  "Cropped at ", cfg$map_lat_min, "°N."))
add_fig("04e_figure6_style_all_schools.png", paste0(
  "Figure 6-style summary for all moving schools, using ", cfg$interval_s, "-s intervals of the smoothed ",
  "tracks. Left: distribution of interval speed and direction. Right: interval direction against speed ",
  "(number of intervals per bin, log colour scale). A linear box plot is not suitable for directions ",
  "(0° and 360° are the same direction); Figure 13 gives the circular summary."))

# ---- 5 caveats + next steps ------------------------------------------------------------------
h1("5. Limitations and next steps")
lead("One day of sonar data.", "All sonar results come from ", format(min(ship$t), "%d %B %Y"),
     " and two mackerel hauls. Export Profos for the rest of the trip (all *-pp, -ps and -ppe files in ",
     "one folder) and re-run the scripts; Figure 2 shows which hauls then gain sonar coverage.")
lead("Depth.", "Profos gives no independent school depth. Use the echosounder (EK80) for vertical ",
     "distribution, or periods with steeper tilt, before comparing schools with fishing depth. For the ",
     "deep blue whiting hauls the omnidirectional sonar mainly informs on the surface layer.")
lead("Relative indices.", "Sv, school area and the backscatter index are uncalibrated and relative; they ",
     "are suitable for comparing hauls on the same vessel and settings, not for biomass.")
lead("Trawl geometry.", "Vertical net opening is not recorded in the catch file; an opening of ",
     cfg$default_vert_opening, " m was assumed. Net sonar or trawl-sensor data would sharpen the depth comparison.")
lead("Tracking.", "Profos can split or merge schools and restart Ids between runs (handled by run-specific ",
     "school keys). Smoothed tracks are very straight (median straightness ",
     f2(median(mov$straightness, na.rm = TRUE)), ") partly because the Kalman model assumes nearly ",
     "constant velocity over a track of about a minute.")
lead("Statistics with more hauls.", "With the whole trip, model catch (or catch rate) against sonar ",
     "indices from the search and tow phases, e.g. with a GAM including species and time of day, and test ",
     "the day–night speed difference and vessel-avoidance hypothesis with day and haul as random effects.")

# ---- appendix --------------------------------------------------------------------------------
h1("Appendix: scripts")
para("All analyses are in the RStudio project sonar_catch_R (GitHub: ssakinan/sonar_catch_R). Settings are ",
     "in cfg at the top of 00_setup_functions.R. Run the scripts in order; 05_report.R re-runs 01–04 ",
     "without plot titles and rebuilds this document.")
apx <- data.table(
  Script = c("00_setup_functions.R", "01_overview_maps.R", "02_sonar_catch_link.R",
             "03_school_track_smoothing.R", "04_school_movement_patterns.R", "05_report.R"),
  Content = c("Settings, file readers, sonar effort, smoothing, circular statistics",
              "Trip map, haul and sonar timeline, sonar zoom map (Figures 1–3)",
              "School timeline, depth geometry, haul zooms, per-haul indices (Figures 4–6, Table 1)",
              "Outlier removal, noise estimate, spline and Kalman smoothing, speed significance (Figures 7–11)",
              "Movement vectors, direction roses, speed covariates, gridded field (Figures 12–16)",
              "This report"))
ft2 <- flextable(apx)
ft2 <- theme_booktabs(ft2)
ft2 <- fontsize(ft2, size = 9, part = "all")
ft2 <- bold(ft2, part = "header")
ft2 <- width(ft2, width = c(2.2, 4.0))
doc <- body_add_flextable(doc, ft2, align = "left")

print(doc, target = out_docx)
Sys.unsetenv("SONAR_REPORT")
message("Report written: ", normalizePath(out_docx))
