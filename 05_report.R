################################################################################
# 05_report.R
#
# Word report for one sonar profile: school quality control, sonar vs catch,
# track smoothing and movement patterns.
#   sonar_profile <- "SU90"; source("05_report.R")   # SU90, 23-29 Oct 2025 (mackerel)
#   sonar_profile <- "ST90"; source("05_report.R")   # ST90, 28 Jun - 5 Jul 2026 (herring)
# Output: output/SCH302_sonar_catch_report.docx (SU90) or
#         output_ST90/SCH302_ST90_sonar_catch_report.docx (ST90)
#
#   1. Re-runs 00b and 01-04 in report mode (figures without titles; captions
#      go below each figure in the document). Set SKIP_FIGURES=1 to reuse them.
#   2. Builds the document with officer + flextable.
#
# Numbers and most statements are computed from the outputs; check the wording
# if the data change.
################################################################################

if (!exists("sonar_profile")) sonar_profile <- Sys.getenv("SONAR_PROFILE", "SU90")
rerun_figures <- !identical(Sys.getenv("SKIP_FIGURES"), "1")
if (rerun_figures) {
  Sys.setenv(SONAR_REPORT = "1")
  for (s in c("00b_school_qc.R", "01_overview_maps.R", "02_sonar_catch_link.R",
              "03_school_track_smoothing.R", "04_school_movement_patterns.R")) {
    message("Running ", s, " in report mode ...")
    source(s, local = new.env())
  }
}
Sys.setenv(SONAR_REPORT = "1")
source("00_setup_functions.R")      # cfg$out_dir -> output_ST90/report_figs
library(officer)
library(flextable)

fig_dir  <- cfg$out_dir
base_dir <- dirname(fig_dir)
out_docx <- file.path(base_dir, if (cfg$label == "ST90") "SCH302_ST90_sonar_catch_report.docx" else "SCH302_sonar_catch_report.docx")

# ---- numbers for the text ---------------------------------------------------------------
hauls <- read_catch(cfg$catch_file, trips = cfg$trip)
ship  <- read_ship()
raw   <- read_profos(cfg$profos_dir, "pp")
qc    <- fread(file.path(base_dir, "school_qc.csv"))
pp    <- apply_school_qc(raw)
ind   <- fread(file.path(fig_dir, "haul_sonar_indices.csv"), colClasses = list(character = "haul_id"))  # "2026E-41" would parse as a number
mov   <- fread(file.path(fig_dir, "school_movement.csv"))
cst   <- fread(file.path(fig_dir, "movement_circular_stats.csv"))
mm    <- mov[moving == TRUE]
q     <- cfg$qc

f0 <- function(x) formatC(round(x), format = "d", big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1, big.mark = ",")
f2 <- function(x) formatC(x, format = "f", digits = 2)
pv <- function(p) if (p < 0.001) "p < 0.001" else paste("p =", formatC(p, digits = 2, format = "g"))
cs <- function(g) cst[group == g]

n_det_raw  <- nrow(raw[!is.na(Id)])
n_sch_raw  <- uniqueN(raw[!is.na(Id)]$sid)
n_ping_drop <- n_det_raw - nrow(qc_pings(raw)[!is.na(Id)])
n_excl     <- sum(qc$excluded); n_kept <- uniqueN(pp[!is.na(Id)]$sid)
rule_n     <- qc[, lapply(.SD, sum), .SDcols = patterns("^r_")]
qt         <- qc[tracked == TRUE]
days       <- sort(unique(as.Date(ship$t)))
cover_h    <- ship[, diff(range(as.numeric(t))), by = seg][, sum(V1)] / 3600
m_all      <- merge_pp(pp)
fits_tilt  <- m_all[, .(n = .N, r2 = summary(lm(depth ~ range_m))$r.squared, max_depth = max(depth)), by = tilt][order(tilt)]
rel_sector <- m_all[, mean(rel_bearing > 90 & rel_bearing < 270)]
noise_m    <- { nz <- m_all[sid %in% mov$sid, { xy <- to_xy(lon, lat); .(sx = noise_sd(xy[, 1]), sy = noise_sd(xy[, 2])) }, by = sid]
                median(c(nz$sx, nz$sy), na.rm = TRUE) }
mov[, v_along := vx * sin(ship_bearing * pi / 180) + vy * cos(ship_bearing * pi / 180)]
fit_art    <- summary(lm(v_along ~ ship_speed, mov[moving == TRUE]))$coefficients
rho_area   <- suppressWarnings(cor.test(log(mm$Area), mm$speed_lm, method = "spearman"))
rho_range  <- suppressWarnings(cor.test(mm$range_m, mm$speed_lm, method = "spearman"))
mm[, hour := as.numeric(format(as.POSIXct(t_start, tz = "UTC"), "%H")) +
             as.numeric(format(as.POSIXct(t_start, tz = "UTC"), "%M")) / 60]
g_hour     <- mgcv::gam(speed_lm ~ s(hour, bs = "cc", k = 12), data = mm, knots = list(hour = c(0, 24)))
hh         <- data.frame(hour = seq(0, 23.9, 0.1)); hh$fit <- predict(g_hour, hh)
peak_h     <- hh$hour[which.max(hh$fit)]; low_h <- hh$hour[which.min(hh$fit)]
haul_st    <- cst[grepl("^haul", group)]
hz         <- hauls[haul_id %in% ind$haul_id]
ind[, tow_h := as.numeric(difftime(haul_time, shoot_time, units = "hours"))]
ind[, rate := catch / tow_h]
cors <- rbindlist(lapply(c("search", "tow"), function(ph) rbindlist(lapply(c("schools_per_h", "sbi_h_db", "depth_w"), function(v) {
  rbindlist(lapply(c("catch", "rate"), function(y) {
    d <- ind[phase == ph & is.finite(get(v))]
    ct <- suppressWarnings(cor.test(d[[y]], d[[v]], method = "spearman"))
    data.table(phase = ph, index = v, response = y, n = nrow(d), rho = ct$estimate, p = ct$p.value)
  }))
}))))
max_abs_rho <- max(abs(cors$rho)); min_p <- min(cors$p)
busy_day <- names(which.max(table(as.Date(pp$t))))
n_below  <- sum(hz$headline_depth > max(m_all$depth))   # tows with headline deeper than any sonar detection
shallow  <- hz[headline_depth <= max(m_all$depth) - 10]                              # well inside the sonar's depth range
edge     <- hz[headline_depth > max(m_all$depth) - 10 & headline_depth <= max(m_all$depth)]  # only the deepest detections
species  <- if (cfg$trip == "2025K") "mackerel" else "herring"
sp_tab   <- hz[, .N, by = dom_species][order(-N)]
sp_name  <- c(mac = "mackerel", her = "herring", hom = "horse mackerel", whb = "blue whiting")
sp_txt   <- paste0(sp_tab$N, " ", ifelse(is.na(sp_name[sp_tab$dom_species]), sp_tab$dom_species, sp_name[sp_tab$dom_species]), collapse = ", ")
# second speed mode at steaming speed? (as many tracks at 5-8 m/s as at 2.5-4 m/s: ST90 0.93, SU90 0.22-0.25)
ship_mode <- qt[, sum(speed_lm > 5 & speed_lm < 8) / sum(speed_lm > 2.5 & speed_lm < 4)] > 0.5
ping_int  <- median(diff(sort(unique(as.numeric(raw$t)))))
has_space <- any(grepl(" ", raw$File))
rel_dir   <- circ_summary((mm$bearing_lm - mm$ship_bearing) %% 360)
day_share <- mov[, .(day = mean(sun == "day"), night = mean(sun == "night"))]
catch_iqr <- round(quantile(hz$catch, c(0.25, 0.75)), -1)

# ---- document helpers ---------------------------------------------------------------------
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
                                  fp_p = fp_par(text.align = "center", keep_with_next = TRUE, padding.top = 6)))
  doc <<- body_add_fpar(doc, fpar(ftext(sprintf("Figure %d. ", fig_no), fp_cap_lab), ftext(caption, fp_cap),
                                  fp_p = fp_par(padding.bottom = 12)), style = "Image Caption")
}
add_table <- function(tab, caption, widths) {
  tab_no <<- tab_no + 1
  doc <<- body_add_fpar(doc, fpar(ftext(sprintf("Table %d. ", tab_no), fp_cap_lab), ftext(caption, fp_cap),
                                  fp_p = fp_par(keep_with_next = TRUE, padding.top = 6)), style = "Table Caption")
  ft <- flextable(tab); ft <- theme_booktabs(ft); ft <- fontsize(ft, size = 8.5, part = "all")
  ft <- bold(ft, part = "header"); ft <- align(ft, align = "center", part = "all"); ft <- width(ft, width = widths)
  doc <<- body_add_flextable(doc, ft, align = "center")
  doc <<- body_add_par(doc, "", style = "Normal")
}
fig_no_next <- function(k = 1) fig_no + k

# ---- title ----------------------------------------------------------------------------------
doc <- body_add_fpar(doc, fpar(ftext(sprintf("%s sonar and %s catches on SCH302", cfg$label, species), fp_text(bold = TRUE, font.size = 20))))
doc <- body_add_fpar(doc, fpar(ftext(sprintf(
  "School quality control, sonar versus catch and school movement, trip %s, %s – %s",
  cfg$trip, format(min(days), "%d %B"), format(max(days), "%d %B %Y")), fp_text(font.size = 13, color = "grey30"))))
doc <- body_add_fpar(doc, fpar(ftext(paste("Report generated", format(Sys.Date(), "%d %B %Y"),
                                           paste0("with the sonar_catch_R scripts (profile ", sonar_profile, ")")),
                                     fp_text(font.size = 9.5, color = "grey40"))))

# ---- 1 summary ------------------------------------------------------------------------------
h1("1. Summary")
para("This report analyses school detections from the Simrad ", cfg$label, " omnidirectional sonar on SCH302, ",
     "processed with the LSSS Profos school tracker, for ", length(days), " days (",
     format(min(days), "%d %B"), " – ", format(max(days), "%d %B %Y"), ") during which the vessel fished mainly ",
     species, ". ", if (cfg$label == "ST90") paste("The ST90 works at a lower frequency and longer range than the SU90, and its",
     "targets here are mainly herring schools. ") else "",
     "Unrealistic or outlier schools are removed before analysis, with the same rules for all datasets.")
lead("Unrealistic schools were removed.",
     f0(n_excl), " of ", f0(n_sch_raw), " schools (", round(100 * n_excl / n_sch_raw), "%) were excluded. ",
     "The largest group appeared to swim faster than ", q$speed_max, " m/s (", rule_n$r_speed, " schools)",
     if (ship_mode) paste(", with a second mode at 6–7 m/s that matches the vessel's steaming speed: echoes that move",
                          "with the ship rather than fish") else ", forming a tail of implausibly fast tracks",
     ". After exclusion the median school speed is ", f2(median(mov$speed_lm)), " m/s and 95% of schools swim ",
     "slower than ", f1(quantile(mov$speed_lm, 0.95)), " m/s, which is realistic for ", species, ".")
lead("Profos depth is a beam-centre depth.",
     "Profos school depth equals range x sin(tilt) exactly (R² > 0.9999 ",
     "for every tilt). With the export limited to about ",
     round(max(m_all$range_m), -1), " m range, no school deeper than ", round(max(m_all$depth)), " m could be seen; ",
     if (n_below > 0) paste0(n_below, " of ", nrow(hz), " hauls fished with the headline deeper than that")
     else paste0("all ", nrow(hz), " hauls fished with the headline within that depth range"),
     " (median headline depth ", round(median(hz$headline_depth)), " m).")
lead("Sonar abundance did not predict catch.",
     "Across ", nrow(hz), " hauls, neither the number of schools nor their backscatter, during the search ",
     "before shooting or during the tow, was related to catch or catch per hour of towing (all |ρ| ≤ ",
     f2(max_abs_rho), ", all p ≥ ", formatC(min_p, digits = 2, format = "f"), ").")
sig_h <- haul_st[rayleigh_p < 0.05]
lead(if (nrow(sig_h) >= nrow(haul_st) / 2) "Schools moved locally in the same direction." else "Little common direction.",
     "Over the whole area there was little common swimming direction (R̄ = ", f2(cs("all")$Rbar), "); within ",
     nrow(sig_h), " of ", nrow(haul_st), " haul areas schools moved together (Rayleigh p < 0.05",
     if (nrow(sig_h)) paste0(", R̄ = ", paste(unique(f2(range(sig_h$Rbar))), collapse = "–")) else "", "). ",
     "In daylight the mean direction was ", round(cs("sun: day")$mean_dir), "° (R̄ = ", f2(cs("sun: day")$Rbar), ", ",
     pv(cs("sun: day")$rayleigh_p), "), at night ", round(cs("sun: night")$mean_dir), "° (R̄ = ", f2(cs("sun: night")$Rbar),
     ", ", pv(cs("sun: night")$rayleigh_p), ").")
lead("Speed varied through the day.",
     "Pooled over all days, school speed peaked around ", sprintf("%02.0f:00", peak_h), " (sonar clock) and was lowest ",
     "around ", sprintf("%02.0f:00", low_h), ". ",
     if (rho_area$p.value < 0.05) paste0(if (rho_area$estimate > 0) "Larger" else "Smaller", " schools swam slightly faster")
     else "Speed was not related to school size", " (Spearman ρ = ", f2(rho_area$estimate), ", ", pv(rho_area$p.value), ").")
if (cfg$catch_offset_h != 0) lead("Clock offset.",
     "Haul times in the catch file are one hour behind the sonar clock for this trip: the vessel track passes the ",
     "logged shoot and haul positions only after adding 1 h (median distance 0.6 and 0.9 km, against 4.3 and 6.9 km ",
     "without correction). All haul times were shifted by +1 h.") else lead("Clock.",
     "Haul times in the catch file match the sonar clock for this trip (the ship track passes the logged shoot positions ",
     "without a time shift).")

# ---- 2 data ---------------------------------------------------------------------------------
h1("2. Data")
para("Sonar: one Profos run (", paste(unique(raw$run), collapse = ", "), ") with a pp file (", f0(n_det_raw), " school detections of ",
     f0(n_sch_raw), " schools) and a ps file (one row per school). There is no ppe export, so pings without schools are ",
     "missing: the vessel track and sonar effort can only be reconstructed from pings that detected a school (",
     length(unique(ship$seg)), " track sections covering ", f0(cover_h), " h). Effort-based indices are therefore expressed ",
     "per hour of the search or tow window, assuming the sonar ran throughout. ",
     if (has_space) "The File column of this export contains spaces, which required a more robust file reader. " else "",
     "Ping interval is about ", f1(ping_int), " s and tilts were mostly ",
     paste(fits_tilt[n == max(n)]$tilt), "° (", f0(max(fits_tilt$n)), " detections). All detections lie within ±90° of the bow (",
     round(100 * (1 - rel_sector)), "%), so the sector behind the vessel is not represented in this export.")
para("Catch: ", nrow(hz), " hauls of trip ", cfg$trip, " overlap the sonar data (hauls ", min(hz$haul), "–", max(hz$haul),
     "); dominant species: ", sp_txt, " (", f0(sum(hz$catch)), " t in total, ", f0(min(hz$catch)), "–", f0(max(hz$catch)),
     " t per haul).")

# ---- 3 QC -----------------------------------------------------------------------------------
h1("3. Quality control: unrealistic and outlier schools")
para("Before analysis, ", f0(n_ping_drop), " detections were removed at ping level: Sv below ", q$ping_sv_min,
     " dB (far below ", species, " school values) or tilt shallower than ", q$ping_tilt_max, "° (near-horizontal beam, ",
     "dominated by surface echoes). Schools were then tested against six rules (Table 1). Speed, ship-alignment and ",
     "jump rules need a track and were applied to the ", f0(nrow(qt)), " schools with at least ", cfg$min_pings,
     " pings over ", cfg$min_duration, " s; the other rules were applied to all schools.")
para(if (ship_mode) paste0("The speed distribution of tracked schools is bimodal (Figure ", fig_no_next(), "): a main mode below 1 m/s and a ",
     "second mode at 6–7 m/s. The second mode coincides with the ship's steaming speed, and these 'schools' head almost ",
     "exactly where the ship heads (Figure ", fig_no_next(2), "), i.e. they are echoes tied to the vessel. The ", q$speed_max,
     " m/s threshold separates the two modes. ") else paste0("Most tracked schools swim slower than 1 m/s, with a tail of faster ",
     "tracks (Figure ", fig_no_next(), "); some of these move in the ship's direction at a large fraction of its speed (Figure ",
     fig_no_next(2), "), i.e. they are echoes tied to the vessel. "), "Of the excluded tracked schools the median speed was ",
     f1(median(qt[excluded == TRUE]$speed_lm)), " m/s, against ", f2(median(qt[excluded == FALSE]$speed_lm)), " m/s for ",
     "those kept. Excluded schools are spread over the whole area (Figure ", fig_no_next(4), "), so the exclusion does not ",
     "remove a particular region.")
qtab <- data.table(
  Rule = c("Speed", "Moves with ship", "Track jump", "Near field", "Weak echo", "Size outlier"),
  Criterion = c(sprintf("Smoothed speed > %s m/s", q$speed_max),
                sprintf("Heading within %s° of ship, speed > %s x ship speed, ship > %s m/s", q$ship_hd_max, q$ship_frac, q$ship_speed_min),
                sprintf("Step between cleaned pings > %s m/s", q$jump_max),
                sprintf("Mean range < %s m", q$range_min),
                sprintf("Mean Sv < %s dB", q$school_sv_min),
                sprintf("|robust z| of log area > %s", q$area_z_max)),
  Reason = c(paste("Faster than", species, "sustain"), "Echo fixed to the vessel", "Tracker merged or switched schools",
             "Wake, bubbles, surface echoes", paste("Too weak for a", species, "school"), "Implausible size"),
  Schools = c(rule_n$r_speed, rule_n$r_ship, rule_n$r_jump, rule_n$r_near, rule_n$r_weak, rule_n$r_area))
add_table(qtab, sprintf("Exclusion rules and number of schools flagged by each (rules overlap; %s schools excluded in total).",
                        f0(n_excl)), c(1.0, 2.6, 1.9, 0.7))
add_fig("qc_a_speed.png", sprintf(paste(
  "Left: smoothed speed of the %s tracked schools (log scale), kept (grey) and excluded (red); dashed line = %s m/s",
  "limit. Right: school speed against ship speed during the track; dotted = 1:1 line, dashed = speed limit."),
  f0(nrow(qt)), q$speed_max))
add_fig("qc_b_ship_alignment.png", sprintf(paste(
  "School speed relative to ship speed (log scale) against the difference between school and ship heading.",
  "The dashed box marks the 'moves with ship' rule (heading within %s°, more than %s x ship speed; the rule also",
  "requires the ship to steam faster than %s m/s)."), q$ship_hd_max, q$ship_frac, q$ship_speed_min))
add_fig("qc_d_reasons.png", "Excluded schools by the combination of rules that flagged them.")
add_fig("qc_c_map.png", "Mean positions of kept (grey) and excluded (red) schools.")

# ---- 4 methods (short) -----------------------------------------------------------------------
h1("4. Methods")
para("Each school detection was assigned to a haul phase: search (the ",
     cfg$pre_shoot_h, " h before shoot) or tow (shoot to haul), and schools within ", cfg$buffer_km, " km of the ship ",
     "track during the tow were linked to the haul. Indices per haul and phase: schools per hour of the window, school ",
     "backscatter per hour (10 log10 of summed linear Sv x area; relative), backscatter-weighted Profos depth and the ",
     "fraction of backscatter in the trawl layer (headline to headline + ", cfg$default_vert_opening, " m, assumed). ",
     "Tracks were cleaned with a Hampel filter and smoothed with a spline (about one degree of freedom per ", cfg$tau_s,
     " s) and a constant-velocity Kalman smoother; a school counted as moving when its speed exceeded twice its ",
     "standard error. Directions were summarised with circular statistics (mean direction, mean resultant length R̄, ",
     "Rayleigh test). Movement maps use ", cfg$grid_km, "-km grid cells.")

# ---- 5 results --------------------------------------------------------------------------------
doc <- body_add_break(doc)
h1("5. Results")
h2("5.1 Fishing and sonar coverage")
para(if (cfg$trip == "2026E") "The trip first fished horse mackerel west of Ireland and then herring east of Shetland and Orkney, where the "
     else "The trip fished mackerel east of Shetland and Orkney in late October and blue whiting west of the Hebrides in early November; the ",
     cfg$label, " data were recorded in the first part (Figures ", fig_no_next(), "–", fig_no_next(3), ").")
add_fig("01a_trip_map.png", paste("Fishing positions of trip", cfg$trip, "(SCH302). Lines join shoot and haul positions;",
  "circles at the haul position are sized by catch and coloured by dominant species (her = herring, hom = horse",
  "mackerel, mac = mackerel, whb = blue whiting); numbers are haul numbers. The black line is the", cfg$label, "ship track."))
add_fig("01b_trip_timeline.png", paste("Hauls of trip", cfg$trip, "over time (bars = shoot to haul, placed at the catch;",
  "colour = dominant species) and periods with", cfg$label, "data (yellow).",
  if (cfg$catch_offset_h != 0) sprintf("Haul times shifted by %+g h to the sonar clock.", cfg$catch_offset_h) else ""))
add_fig("01c_sonar_zoom_map.png", paste0("Ship track and mean positions of the ", f0(n_kept), " schools kept after quality ",
  "control, sized by school area and coloured by Profos depth. Track colour: orange = search (", cfg$pre_shoot_h,
  " h before shoot), red = tow, grey = other; red arrows = logged tow lines with haul number, species and catch."))

h2("5.2 Sonar versus catch")
para("School detections and tows over the whole period are shown in Figure ", fig_no_next(), ", and for the busiest day (",
     format(as.Date(busy_day), "%d %B"), ") in Figure ", fig_no_next(2), "; timelines for every day are in the ",
     basename(base_dir), " folder. ", if (n_below > 0) paste0(n_below, " of ", nrow(hz), " tows fished with the headline deeper than any school the ",
     "sonar could detect at these tilts (Figure ", fig_no_next(3), "), so for those hauls the fraction of school ",
     "backscatter in the trawl layer is zero. ") else paste0("All tows fished with the headline within the depths the sonar could see ",
     "at these tilts (Figure ", fig_no_next(3), "). "),
     if (nrow(shallow)) paste0(nrow(shallow), " tows (", paste0("haul ", shallow$haul, ": ", shallow$headline_depth, " m", collapse = ", "),
                               ") fished well within the depths seen by the sonar") else "No tow fished well within the depths seen by the sonar",
     if (nrow(edge)) paste0("; the tows at ", paste(unique(range(edge$headline_depth)), collapse = "–"),
                            " m reach only the deepest detections, at the edge of the range") else "", ".")
para("Schools were found near every tow (Figure ", fig_no_next(4), ", Table 2). The median rate was ",
     f1(median(ind[phase == "search"]$schools_per_h)), " schools per hour while searching and ",
     f1(median(ind[phase == "tow"]$schools_per_h)), " per hour while towing. Neither catch nor catch per hour of towing ",
     "was related to the school rate, the school backscatter or the Profos depth of either phase (Figure ", fig_no_next(5),
     "; Spearman |ρ| ≤ ", f2(max_abs_rho), ", n = ", nrow(hz), "). Likely reasons: catch per haul is limited by the ",
     "skipper's decision and vessel capacity (half of the hauls were ", catch_iqr[1], "–", catch_iqr[2], " t); the sonar does not see the depth layer that ",
     "was fished; and counts within ", cfg$buffer_km, " km describe the area rather than the fish that entered the net.")
add_fig("02a_timeline_schools_vs_tows.png", paste0("School detections against time for the whole period. Upper panel: ",
  "each detection by time and Profos depth, coloured by Sv and sized by school area; red shading = tows, orange = ",
  cfg$pre_shoot_h, " h search before shoot, dashed boxes = trawl layer; labels = haul numbers. Lower panel: distinct ",
  "schools per hour of recording (10-min bins) and ship speed (line, right axis)."))
add_fig(sprintf("02a_timeline_%s.png", busy_day), paste0("As the previous figure, for ", format(as.Date(busy_day), "%d %B %Y"),
  "; labels give haul number, catch and share of the dominant species."))
add_fig("02b_depth_vs_range_tilt.png", paste0("Profos school depth against range by tilt. Detections fall exactly on the ",
  "beam-centre lines (depth = range x sin(|tilt|)); bands = ±", cfg$vert_beam_half, "° vertical beam (approximate); ",
  "dashed lines = headline depths of the hauls."))
add_fig("02c_haul_zoom_maps.png", paste0("Schools within ", cfg$buffer_km, " km of each tow (panel title = haul number and ",
  "catch). Orange = ship track while searching, red = while towing (shown only where schools were detected, as there is no ",
  "ppe export), dashed arrow = logged tow line; circles = school mean positions coloured by Profos depth. Each panel has ",
  "equal east–west and north–south extent."))
ind[!is.finite(sbi_h_db), sbi_h_db := NA]
ind[, `:=`(catch_t = round(catch), tow_hr = round(tow_h, 1))]
itab <- dcast(ind, haul + catch_t + tow_hr ~ phase, value.var = c("schools_per_h", "sbi_h_db"))
setnames(itab, c("Haul", "Catch (t)", "Tow (h)", "Schools/h search", "Schools/h tow", "SBI/h search (dB)", "SBI/h tow (dB)"))
for (j in 4:7) set(itab, j = j, value = round(itab[[j]], 1))
add_table(itab, paste0("Sonar indices per haul: schools per hour and school backscatter per hour (SBI; relative) ",
  "within ", cfg$buffer_km, " km of the tow track, in the ", cfg$pre_shoot_h, " h search before shoot and during the tow. ",
  "Blank = no schools in that phase."), c(0.5, 0.75, 0.6, 1.0, 0.9, 1.1, 1.0))
add_fig("02d_catch_vs_sonar_indices.png", "Catch per haul against the sonar indices of the search (circles) and tow (triangles) phases, with linear fits and 95% confidence bands per phase.")

h2("5.3 Cleaning and smoothing school tracks")
para(f0(nrow(mov)), " schools passed the track criteria. Centroid positioning noise was about ", f1(noise_m), " m per axis; ",
     "the smoothed tracks follow the raw positions closely (Figures ", fig_no_next(), "–", fig_no_next(2), "), and mean speed ",
     "becomes stable for smoothing time scales of about 20–30 s or more (Figure ", fig_no_next(3), "). ",
     round(100 * mean(mov$moving)), "% of schools had a speed significantly above zero. Our speeds agree with the LSSS ",
     "estimates for most schools (Figure ", fig_no_next(4), ").")
para("Even after exclusion, school velocity along the ship's heading increased slightly with ship speed (slope ",
     f2(fit_art[2, 1]), ", ", pv(fit_art[2, 4]), "; Figure ", fig_no_next(5), "). The effect is small (about ",
     f1(fit_art[2, 1] * 5), " m/s at 5 m/s ship speed) but indicates either remaining vessel-related echoes or avoidance; ",
     "relative to the ship's heading, moving schools swam on average towards ", round(rel_dir$mean_dir), "° (R̄ = ", f2(rel_dir$Rbar), ").")
add_fig("03a_raw_vs_smoothed_tracks.png", "Raw and smoothed tracks of the six longest-tracked schools (local coordinates, m). Dark points = raw centroids; red = outliers removed; blue = spline; red line = Kalman smoother; triangle = first detection.")
add_fig("03b_speed_raw_vs_smoothed.png", "Speed over time for the same schools: raw ping-to-ping (grey), spline (blue) and Kalman (red, ±2 SD band). Y-axis truncated at 6 m/s.")
add_fig("03c_noise_and_smoothing_sensitivity.png", "Left: positioning noise per school against range (loess fit). Right: mean school speed against the spline smoothing time scale (median and interquartile range; log axis); dashed = median regression speed.")
add_fig("03d_movement_significance.png", "Left: school speed against track duration with ±2 SE bars; black = speed significantly above zero. Right: speed from this analysis against the LSSS Profos school speed; dashed = 1:1.")
add_fig("03e_ship_motion_artefact_check.png", "Left: school velocity along the ship's heading against ship speed, with linear fit. Right: swimming direction relative to the ship's heading (moving schools).")

h2("5.4 School movement patterns")
para("Moving schools (n = ", f0(nrow(mm)), ") swam in all directions over the area as a whole (Figures ", fig_no_next(), "–",
     fig_no_next(2), "); the gridded field (Figure ", fig_no_next(4), ") shows where neighbouring schools shared a direction. R̄ within haul areas ",
     "ranged from ", f2(min(haul_st$Rbar)), " to ", f2(max(haul_st$Rbar)), ". Speed ",
     if (rho_range$p.value < 0.05) "varied" else "did not vary", " with range from the ship (ρ = ", f2(rho_range$estimate),
     ") and peaked around ", sprintf("%02.0f:00", peak_h), " (Figure ", fig_no_next(3), "). ",
     round(100 * day_share$day), "% of tracked schools were seen in daylight and ", round(100 * day_share$night), "% at night.")
add_fig("04a_movement_vector_map.png", sprintf("Movement of schools with significant speed (n = %s). Arrows start at each school's mean position and show the distance swum in 10 min, coloured by speed; dark lines = smoothed tracks; grey = ship track; red dashed = logged tow lines.", f0(nrow(mm))))
add_fig("04b_direction_roses.png", "Swimming directions (22.5° bins) by ship activity (top) and light condition (bottom). Panel headers give n, mean direction, R̄ and the Rayleigh test p-value.")
add_fig("04c_speed_vs_covariates.png", "School speed against range (top left), school area (top right, log scale) and hour of day pooled over all days (bottom, sonar clock; cyclic GAM with 95% band), coloured by light condition.")
add_fig("04d_gridded_movement_field.png", sprintf("Mean school movement per %s km cell. Fill = number of moving schools; arrow = mean velocity (relative length); arrow thickness = R̄ (agreement in direction); grey = ship track.", cfg$grid_km))
add_fig("04e_figure6_style_all_schools.png", sprintf("Figure 6-style summary using %s-s intervals of the smoothed tracks. Left: interval speed and direction. Right: direction against speed (intervals per bin, log colour scale). Figure %d gives the circular summary of directions.", cfg$interval_s, fig_no - 3))   # caption is evaluated after add_fig increments fig_no

# ---- 6 limitations ----------------------------------------------------------------------------
h1("6. Limitations and next steps")
lead("No ppe export.", "Effort and the ship track come only from pings with schools. Exporting the ppe file from LSSS would ",
     "give true recording time and distance sailed, and a complete ship track.")
lead("Range and sector.", sprintf("Detections stop at about %s m and only the forward ±90° sector is present, although the sonar can reach further. Re-processing with a longer range and full sector would exploit the instrument and show schools nearer the fishing depth.", round(max(m_all$range_m), -1)))
lead("Depth.", "Profos depth is the beam-centre depth; the echosounder is needed to compare schools with the trawl depth.")
if (cfg$catch_offset_h != 0) lead("Clock.", "The catch log and the sonar clock differ by 1 h on this trip; which one is UTC should be confirmed, as it shifts the hour-of-day results by 1 h.")
lead("Thresholds.", paste0("The exclusion thresholds were set from the ST90 herring data and are used unchanged for all datasets. They are stored in cfg$qc (00_setup_functions.R) and the flags per school in ", basename(base_dir), "/school_qc.csv, so they can be reviewed and changed."))
lead("Catch.", "Catch per haul is a coarse response. Linking sonar to catch may work better with net-sonar or trawl-sensor data, or with sonar backscatter scaled to the trawl path rather than a 5-km buffer.")

h1("Appendix: scripts and animations")
para("Run with sonar_profile <- \"", sonar_profile, "\": 00b_school_qc.R (exclusions), 01–04 (figures), 05_report.R (this ",
     "report), 06_school_animation.R and 07_vessel_school_animation.R (animations in ", basename(base_dir), "/animation; the vessel ",
     "animations are split into one video per day), 08 (daytime migration) and 09 (CMEMS currents).")

print(doc, target = out_docx)
Sys.unsetenv("SONAR_REPORT")
message("Report written: ", normalizePath(out_docx))
