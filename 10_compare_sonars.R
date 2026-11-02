################################################################################
# 10_compare_sonars.R
#
# ST90 versus SU90 on the same trip and period (28 Jun - 5 Jul 2026).
# Run first, for both profiles "ST90" and "SU90_2026": 00b, 01-04, 08 and 09.
#   source("10_compare_sonars.R")
#
#   1. detections: schools per hour over time, QC exclusions
#   2. school properties: range, depth, Sv, area, speed
#   3. schools seen by both sonars: tracks matched in time and space,
#      position offset and agreement of the smoothed velocities
#   4. movement: daytime direction, net drift (08) and swimming (09)
#   5. catch: sonar indices per haul from both sonars against catch
#
# Outputs: output_compare/ figures 10a-10f, 10_sonar_comparison_summary.png,
# compare_*.csv and a Word summary.
################################################################################

sonar_profile <- "ST90"
source("00_setup_functions.R")   # run from the RStudio project root
library(officer); library(flextable)

match_dist_m <- 100   # max distance between the two tracks at the common mid time
match_min_s  <- 30    # min time overlap of the two tracks
sonars <- c(ST90 = "output_ST90", SU90 = "output_SU90_2026")
cols_s <- c(ST90 = "#1b7f79", SU90 = "#d9822b")
out <- "output_compare"
dir.create(out, showWarnings = FALSE)

rd <- function(file, ...) rbindlist(lapply(names(sonars), function(s)
  fread(file.path(sonars[[s]], file), ...)[, sonar := s]), fill = TRUE)
qc  <- rd("school_qc.csv")
mov <- rd("school_movement.csv")
mov[, `:=`(t_start = as.POSIXct(t_start, tz = "UTC"), t_end = as.POSIXct(t_end, tz = "UTC"))]
ind <- rd("haul_sonar_indices.csv", colClasses = list(character = "haul_id"))
mig <- rd("migration/migration_estimates.csv")
cme <- rd("migration/cmems_estimates.csv")
trk <- rbindlist(lapply(names(sonars), function(s)
  readRDS(file.path(sonars[[s]], "school_tracks.rds"))[keep == TRUE, .(sid, t, lon_kf, lat_kf)][, sonar := s]))

# common period: hours in which at least one sonar recorded schools
t_rng <- range(mov$t_start)

# ---- 1. detections over time -------------------------------------------------------------------------
mov[, hr := as.POSIXct(floor(as.numeric(t_start) / 3600) * 3600, origin = "1970-01-01", tz = "UTC")]
hours <- CJ(hr = seq(min(mov$hr), max(mov$hr), by = 3600), sonar = names(sonars))
ph <- mov[, .N, by = .(hr, sonar)][hours, on = .(hr, sonar)][is.na(N), N := 0L]
phw <- dcast(ph, hr ~ sonar, value.var = "N")
both_on <- phw[ST90 > 0 | SU90 > 0]
rho_hour <- cor.test(both_on$ST90, both_on$SU90, method = "spearman", exact = FALSE)

# ---- 2. school properties -------------------------------------------------------------------------------
summ <- mov[, .(schools = .N, moving = mean(moving), day = mean(sun == "day"),
                range_med = median(range_m), range_95 = quantile(range_m, 0.95),
                depth_med = median(depth), Sv_med = median(Sv), area_med = median(Area),
                speed_med = median(speed_lm), speed_95 = quantile(speed_lm, 0.95),
                dur_med = median(duration_s)), by = sonar]
qcs <- qc[, .(detected = .N, excluded = sum(excluded)), by = sonar]
summ <- qcs[summ, on = "sonar"]
fwrite(summ, file.path(out, "compare_summary.csv"))
print(summ)

# ---- 3. schools seen by both sonars -----------------------------------------------------------------------
a <- mov[sonar == "ST90", .(sid_st = sid, s0 = t_start, s1 = t_end)]
b <- mov[sonar == "SU90", .(sid_su = sid, s0 = t_start, s1 = t_end)]
setkey(a, s0, s1)
cand <- foverlaps(b, a, by.x = c("s0", "s1"), by.y = c("s0", "s1"), nomatch = NULL)
cand[, `:=`(o0 = pmax(s0, i.s0), o1 = pmin(s1, i.s1))]
cand <- cand[as.numeric(o1 - o0, units = "secs") >= match_min_s]
cand[, tm := as.numeric(o0) + as.numeric(o1 - o0, units = "secs") / 2]
pos_at <- function(tr, ids, tm) {
  tr <- tr[sid %in% ids]; sp <- split(tr, tr$sid)
  t(mapply(function(id, tq) { x <- sp[[id]]; tn <- as.numeric(x$t)
    c(approx(tn, x$lon_kf, tq, rule = 2)$y, approx(tn, x$lat_kf, tq, rule = 2)$y) }, ids, tq = tm, USE.NAMES = FALSE))
}
pst <- pos_at(trk[sonar == "ST90"], cand$sid_st, cand$tm)
psu <- pos_at(trk[sonar == "SU90"], cand$sid_su, cand$tm)
cand[, `:=`(dx = (psu[, 1] - pst[, 1]) * 111320 * cos(pst[, 2] * pi / 180), dy = (psu[, 2] - pst[, 2]) * 111320)]
cand[, dist := sqrt(dx^2 + dy^2)]
# one-to-one: closest pairs first
cand <- cand[dist <= match_dist_m][order(dist)]
used_st <- used_su <- character(0); keep <- logical(nrow(cand))
for (i in seq_len(nrow(cand))) if (!(cand$sid_st[i] %in% used_st) && !(cand$sid_su[i] %in% used_su)) {
  keep[i] <- TRUE; used_st <- c(used_st, cand$sid_st[i]); used_su <- c(used_su, cand$sid_su[i]) }
pairs <- cand[keep]
vcols <- c("sid", "vx", "vy", "speed_lm", "bearing_lm", "depth", "Sv", "Area", "range_m", "moving", "sun", "duration_s")
pairs <- merge(pairs, setnames(mov[sonar == "ST90", ..vcols], vcols, paste0(vcols, "_st")), by = "sid_st")
pairs <- merge(pairs, setnames(mov[sonar == "SU90", ..vcols], vcols, paste0(vcols, "_su")), by = "sid_su")
pairs[, ddir := ((bearing_lm_su - bearing_lm_st + 180) %% 360) - 180]
fwrite(pairs, file.path(out, "compare_matched_schools.csv"))
n_pairs <- nrow(pairs)
cor_v <- c(x = cor(pairs$vx_st, pairs$vx_su), y = cor(pairs$vy_st, pairs$vy_su))
mv_pairs <- pairs[moving_st & moving_su]
ddir_mv <- circ_summary(mv_pairs$ddir %% 360)
med_abs_ddir <- median(abs(mv_pairs$ddir))
off <- pairs[, .(dx = median(dx), dy = median(dy), dist = median(dist))]
# overlap between sonars: share of ST90 / SU90 tracked schools that have a partner
share_st <- n_pairs / mov[sonar == "ST90", .N]; share_su <- n_pairs / mov[sonar == "SU90", .N]
# simultaneous recording: share of ST90 schools at times when the SU90 also recorded
message(sprintf("Matched schools: %d (%.0f%% of ST90, %.0f%% of SU90); velocity r: east %.2f, north %.2f",
                n_pairs, 100 * share_st, 100 * share_su, cor_v["x"], cor_v["y"]))

# ---- 4. movement -----------------------------------------------------------------------------------------------
keyrows <- function(x, pat) x[grepl(pat, model)]
drift <- rbind(keyrows(mig, "\\(main\\)")[, what := "Net drift over ground (08)"],
               keyrows(cme, "^Swimming = ")[, what := "Swimming relative to CMEMS (09)"],
               keyrows(cme, "^CMEMS current at schools")[, what := "CMEMS current at the schools (09)"], fill = TRUE)
drift[, what := factor(what, unique(what))]
fwrite(drift, file.path(out, "compare_drift.csv"))
day_dirs <- mov[sun == "day" & moving == TRUE, .(n = .N, mean_dir = circ_summary(bearing_lm)$mean_dir,
                                                 Rbar = circ_summary(bearing_lm)$Rbar), by = sonar]

# ---- 5. catch ------------------------------------------------------------------------------------------------
ind[, tow_h := as.numeric(difftime(as.POSIXct(haul_time, tz = "UTC"), as.POSIXct(shoot_time, tz = "UTC"), units = "hours"))]
ind[, rate := catch / tow_h]
cors <- ind[, rbindlist(lapply(c("schools_per_h", "sbi_h_db"), function(v) rbindlist(lapply(c("catch", "rate"), function(y) {
  x <- .SD[is.finite(get(v))]
  ct <- suppressWarnings(cor.test(x[[y]], x[[v]], method = "spearman"))
  data.table(index = v, response = y, n = nrow(x), rho = ct$estimate, p = ct$p.value)
})))), by = .(sonar, phase)]
fwrite(cors, file.path(out, "compare_catch_correlations.csv"))
iw <- dcast(ind, haul_id + phase + catch ~ sonar, value.var = c("schools_per_h", "sbi_h_db"))
rho_haul <- iw[, .(rho_n = cor(schools_per_h_ST90, schools_per_h_SU90, method = "spearman", use = "complete.obs"),
                   rho_sbi = cor(sbi_h_db_ST90, sbi_h_db_SU90, method = "spearman", use = "complete.obs")), by = phase]

# ================================ FIGURES ================================================================
lab_s <- scale_colour_manual(values = cols_s, name = NULL)
fill_s <- scale_fill_manual(values = cols_s, name = NULL)

# 10a: detections per hour ----------------------------------------------------------------------------------
p10a1 <- ggplot(ph, aes(hr, N, colour = sonar)) + geom_line(alpha = 0.8) + lab_s +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%d %b") +
  labs(x = NULL, y = "Tracked schools per hour", title = "Detections over time") + theme(legend.position = "bottom")
p10a2 <- ggplot(both_on, aes(ST90, SU90)) + geom_abline(colour = "grey70", linetype = 2) +
  geom_point(alpha = 0.4, size = 1.2) + coord_equal() +
  labs(x = "ST90 schools per hour", y = "SU90 schools per hour", title = "Same hours",
       subtitle = sprintf("Spearman rho = %.2f, %d hours", rho_hour$estimate, nrow(both_on)))
ggsave(file.path(out, "10a_detections.png"), p10a1 + p10a2 + plot_layout(widths = c(2.4, 1)), width = 13, height = 4.8, dpi = 200)

# 10b: school properties -----------------------------------------------------------------------------------
pl <- melt(mov[, .(sonar, `Range (m)` = range_m, `Beam-centre depth (m)` = depth, `Sv (dB)` = Sv,
                   `log10 area (m2)` = log10(Area), `Speed (m/s)` = speed_lm, `Track duration (min)` = duration_s / 60)],
           id.vars = "sonar")
p10b <- ggplot(pl, aes(value, colour = sonar, fill = sonar)) + geom_density(alpha = 0.15) +
  facet_wrap(~ variable, scales = "free", ncol = 3) + lab_s + fill_s +
  labs(x = NULL, y = "Density", title = "School properties (tracked schools after QC)") + theme(legend.position = "bottom")
ggsave(file.path(out, "10b_school_properties.png"), p10b, width = 11, height = 6.5, dpi = 200)

# 10c: schools seen by both sonars -------------------------------------------------------------------------
vl <- rbind(pairs[, .(comp = "East velocity (m/s)", st = vx_st, su = vx_su)], pairs[, .(comp = "North velocity (m/s)", st = vy_st, su = vy_su)])
lim_v <- quantile(abs(c(vl$st, vl$su)), 0.995)
p10c1 <- ggplot(vl, aes(st, su)) + geom_abline(colour = "grey60", linetype = 2) +
  geom_point(alpha = 0.3, size = 0.8) + facet_wrap(~ comp) + coord_equal(xlim = c(-lim_v, lim_v), ylim = c(-lim_v, lim_v)) +
  labs(x = "ST90", y = "SU90", title = "Same school, two sonars: smoothed velocity",
       subtitle = sprintf("%d matched schools; r = %.2f (east), %.2f (north)", n_pairs, cor_v["x"], cor_v["y"]))
p10c2 <- ggplot(mv_pairs, aes(ddir)) + geom_histogram(breaks = seq(-180, 180, 15), fill = "grey40", colour = "white") +
  scale_x_continuous(breaks = seq(-180, 180, 90)) +
  labs(x = "SU90 - ST90 direction (deg)", y = "Schools", title = "Direction difference",
       subtitle = sprintf("Both moving (n = %d); median |diff| %.0f deg", nrow(mv_pairs), med_abs_ddir))
p10c3 <- ggplot(pairs, aes(dx, dy)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(alpha = 0.25, size = 0.8) + annotate("point", x = off$dx, y = off$dy, colour = "firebrick", size = 3) + coord_equal() +
  labs(x = "East offset (m)", y = "North offset (m)", title = "SU90 - ST90 position",
       subtitle = sprintf("Median offset %.0f m E, %.0f m N", off$dx, off$dy))
ggsave(file.path(out, "10c_matched_schools.png"), p10c1 + p10c2 + p10c3 + plot_layout(widths = c(2, 1, 1)),
       width = 15, height = 5, dpi = 200)

# 10d: daytime directions -----------------------------------------------------------------------------------
dd <- mov[sun == "day" & moving == TRUE]
p10d <- ggplot(dd, aes(bearing_lm, fill = sonar)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 360, 15), colour = "white") +
  geom_vline(data = day_dirs, aes(xintercept = mean_dir), colour = "black", linewidth = 0.8) +
  coord_polar(start = 0) + facet_wrap(~ sonar) + fill_s +
  scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270), labels = c("N", "E", "S", "W")) +
  geom_text(data = day_dirs, aes(x = 45, y = Inf, label = sprintf("n = %d\nmean %.0f deg, Rbar %.2f", n, mean_dir, Rbar)),
            inherit.aes = FALSE, vjust = 1.2, size = 3) +
  labs(x = NULL, y = NULL, title = "Direction of daytime moving schools (observed, tide not removed)") +
  theme(legend.position = "none", axis.text.y = element_blank())
ggsave(file.path(out, "10d_daytime_directions.png"), p10d, width = 10, height = 5.5, dpi = 200)

# 10e: net drift and swimming ---------------------------------------------------------------------------------
lim_d <- max(abs(drift[, c(vx, vy)])) * 86.4 * 1.5
p10e <- ggplot(drift) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_segment(aes(0, 0, xend = vx * 86.4, yend = vy * 86.4, colour = sonar), arrow = arrow(length = unit(2.5, "mm")),
               linewidth = 1.2) +
  ggrepel::geom_text_repel(aes(vx * 86.4, vy * 86.4, colour = sonar,
                label = sprintf("%.1f km/d, %.0f deg\n(%.0f-%.0f)", km_day, dir, dir_lo, dir_hi)),
            size = 2.8, min.segment.length = 0, box.padding = 0.6, seed = 1, show.legend = FALSE) +
  facet_wrap(~ what) + coord_equal(xlim = c(-lim_d, lim_d), ylim = c(-lim_d, lim_d)) + lab_s +
  labs(x = "East (km/day)", y = "North (km/day)", title = "Daytime net movement estimated from each sonar",
       subtitle = "Tide and vessel effect removed; numbers = speed, direction and 95% bootstrap interval of direction") +
  theme(legend.position = "bottom")
ggsave(file.path(out, "10e_drift_comparison.png"), p10e, width = 13, height = 5.5, dpi = 200)

# 10f: catch ---------------------------------------------------------------------------------------------------
p10f1 <- ggplot(iw, aes(schools_per_h_ST90, schools_per_h_SU90)) + geom_abline(colour = "grey70", linetype = 2) +
  geom_point(aes(size = catch), alpha = 0.6) + facet_wrap(~ phase) + scale_size_area(max_size = 5, name = "Catch (t)") +
  labs(x = "ST90 schools per hour", y = "SU90 schools per hour", title = "Per haul: the two sonars")
cl <- copy(cors)[index == "schools_per_h"][, lab := paste(phase, response)]
p10f2 <- ggplot(cors, aes(rho, interaction(response, index, sep = " ~ "), colour = sonar)) +
  geom_vline(xintercept = 0, colour = "grey70") + geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ phase) + lab_s + xlim(-1, 1) +
  labs(x = "Spearman rho", y = NULL, title = "Catch vs sonar index") + theme(legend.position = "bottom")
ggsave(file.path(out, "10f_catch.png"), p10f1 + p10f2 + plot_layout(widths = c(1.3, 1)), width = 14, height = 5.5, dpi = 200)

# summary figure ------------------------------------------------------------------------------------------------
p_sum <- (p10a1 + labs(title = "A  Tracked schools per hour")) /
  ((p10c1 + labs(title = "B  Same school, two sonars", subtitle = NULL)) + (p10c2 + labs(title = "C  Direction difference", subtitle = NULL)) +
     plot_layout(widths = c(2, 1))) /
  (p10e + labs(title = "D  Daytime net drift, swimming and current", subtitle = NULL)) +
  plot_layout(heights = c(1, 1.2, 1.3)) +
  plot_annotation(title = sprintf("ST90 vs SU90, SCH302 %s - %s: %d matched schools, velocity r = %.2f / %.2f (east / north)",
                                  format(t_rng[1], "%d %b"), format(t_rng[2], "%d %b %Y"), n_pairs, cor_v["x"], cor_v["y"]))
ggsave(file.path(out, "10_sonar_comparison_summary.png"), p_sum, width = 13, height = 15, dpi = 200)
message("Figures written to ", out)

# ================================ WORD SUMMARY ==========================================================
f0 <- function(x) formatC(round(x), format = "d", big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
S <- function(s, v) summ[sonar == s][[v]]
D <- function(s, w) drift[sonar == s & what == w]
ang <- function(a, b) abs(((a - b + 180) %% 360) - 180)
wg <- levels(drift$what)
fp_cap_lab <- fp_text(bold = TRUE, font.size = 9.5); fp_cap <- fp_text(font.size = 9.5); fp_bold <- fp_text(bold = TRUE)
fig_no <- 0
doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(orient = "portrait", width = 21 / 2.54, height = 29.7 / 2.54),
  page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1)))
h1 <- function(x) doc <<- body_add_par(doc, x, style = "heading 1")
para <- function(...) doc <<- body_add_par(doc, paste0(...), style = "Normal")
lead <- function(h, ...) doc <<- body_add_fpar(doc, fpar(ftext(paste0(h, " "), fp_bold), ftext(paste0(...))), style = "Normal")
add_fig <- function(file, caption, max_w = 6.25, max_h = 8.3) {
  fig_no <<- fig_no + 1
  dd <- dim(png::readPNG(file.path(out, file))); wd <- max_w; ht <- wd * dd[1] / dd[2]
  if (ht > max_h) { ht <- max_h; wd <- ht * dd[2] / dd[1] }
  doc <<- body_add_fpar(doc, fpar(external_img(file.path(out, file), width = wd, height = ht),
                                  fp_p = fp_par(text.align = "center", keep_with_next = TRUE, padding.top = 6)))
  doc <<- body_add_fpar(doc, fpar(ftext(sprintf("Figure %d. ", fig_no), fp_cap_lab), ftext(caption, fp_cap),
                                  fp_p = fp_par(padding.bottom = 12)), style = "Image Caption")
}

doc <- body_add_fpar(doc, fpar(ftext("ST90 versus SU90 on the same herring trip", fp_text(bold = TRUE, font.size = 18))))
doc <- body_add_fpar(doc, fpar(ftext(sprintf("SCH302, trip 2026E, %s – %s. Generated %s with 10_compare_sonars.R",
  format(t_rng[1], "%d %B"), format(t_rng[2], "%d %B %Y"), format(Sys.Date(), "%d %B %Y")), fp_text(font.size = 10, color = "grey40"))))

h1("Main results")
lead("Clocks.", sprintf(paste(
  "The two exports do not share a clock: the SU90 time stamps run 1 h %d min %02d s behind the ST90 (found by aligning the ship",
  "positions written into the two exports; they then agree to about 15 m, with the same offset on every day). The SU90",
  "times were moved onto the ST90 clock before the analysis (clock_shift_s in the SU90_2026 profile). As the catch log is",
  "about 1 h behind the ST90, the SU90 clock and the catch log are probably both close to UTC, and the ST90 clock about one",
  "hour ahead."), profiles$SU90_2026$clock_shift_s %/% 60 - 60, profiles$SU90_2026$clock_shift_s %% 60))
lead("Detections.", sprintf(paste(
  "Both exports were processed identically (same quality control, smoothing and tracking). The ST90 gave %s tracked schools",
  "(%s detected, %s excluded by QC) and the SU90 %s (%s detected, %s excluded). Hourly numbers of tracked schools from the",
  "two sonars are correlated (Spearman ρ = %s over %d hours with data)."),
  f0(S("ST90", "schools")), f0(S("ST90", "detected")), f0(S("ST90", "excluded")),
  f0(S("SU90", "schools")), f0(S("SU90", "detected")), f0(S("SU90", "excluded")), f2(rho_hour$estimate), nrow(both_on)))
lead("School properties.", sprintf(paste(
  "Median range %s m (ST90) and %s m (SU90); median beam-centre depth %s and %s m (this reflects the tilt settings, not the",
  "fish); median Sv %s and %s dB; median area %s and %s m²; median smoothed speed %s and %s m/s."),
  f0(S("ST90", "range_med")), f0(S("SU90", "range_med")), f0(S("ST90", "depth_med")), f0(S("SU90", "depth_med")),
  f1(S("ST90", "Sv_med")), f1(S("SU90", "Sv_med")), f0(S("ST90", "area_med")), f0(S("SU90", "area_med")),
  f2(S("ST90", "speed_med")), f2(S("SU90", "speed_med"))))
lead("Same schools.", sprintf(paste(
  "%s schools were tracked by both sonars at the same time within %d m (%.0f%% of ST90 and %.0f%% of SU90 tracks).",
  "Their smoothed velocities agree with r = %s (east) and %s (north); for pairs where both tracks were significantly moving,",
  "the median direction difference is %.0f°. The median position offset SU90 − ST90 is %.0f m east and %.0f m north.",
  "For the same schools the beam-centre depth is %.0f m (ST90) and %.0f m (SU90) in the median, confirming that Profos depth",
  "reflects the tilt of each sonar, not the depth of the fish."),
  f0(n_pairs), match_dist_m, 100 * share_st, 100 * share_su, f2(cor_v["x"]), f2(cor_v["y"]), med_abs_ddir, off$dx, off$dy,
  median(pairs$depth_st), median(pairs$depth_su)))
lead("Migration signal.", sprintf(paste(
  "Daytime net drift over ground: ST90 %s km/day towards %.0f°, SU90 %s km/day towards %.0f° (%.0f° apart).",
  "Swimming relative to the CMEMS current: ST90 %s km/day towards %.0f°, SU90 %s km/day towards %.0f°."),
  f1(D("ST90", wg[1])$km_day), D("ST90", wg[1])$dir, f1(D("SU90", wg[1])$km_day), D("SU90", wg[1])$dir,
  ang(D("ST90", wg[1])$dir, D("SU90", wg[1])$dir),
  f1(D("ST90", wg[2])$km_day), D("ST90", wg[2])$dir, f1(D("SU90", wg[2])$km_day), D("SU90", wg[2])$dir))
lead("Catch.", sprintf(paste(
  "Per haul, the two sonars' school rates are correlated (Spearman ρ = %s in the search phase, %s during the tow). Neither",
  "sonar index predicts catch well: the largest |ρ| with catch or catch rate is %s for the ST90 and %s for the SU90."),
  f2(rho_haul[phase == "search"]$rho_n), f2(rho_haul[phase == "tow"]$rho_n),
  f2(max(abs(cors[sonar == "ST90"]$rho))), f2(max(abs(cors[sonar == "SU90"]$rho)))))
add_fig("10_sonar_comparison_summary.png", "Summary. A: tracked schools per hour from each sonar. B: smoothed velocity of the same school tracked by both sonars. C: direction difference for pairs where both tracks are significantly moving. D: daytime net drift over ground, swimming relative to the CMEMS current and the model current at the schools, estimated separately from each sonar.")

h1("Details")
add_fig("10a_detections.png", "Tracked schools per hour over the trip (left) and in the same hours (right; dashed = 1:1).")
add_fig("10b_school_properties.png", "Distributions of school properties for the tracked schools of each sonar after the same quality control.")
add_fig("10c_matched_schools.png", sprintf("Schools tracked by both sonars (overlap ≥ %d s, positions within %d m at the common mid time; one-to-one, closest first). Left: smoothed velocity; middle: direction difference for pairs where both tracks move significantly; right: position offset (red = median).", match_min_s, match_dist_m))
add_fig("10d_daytime_directions.png", "Direction of significantly moving daytime schools (observed velocity, tide not removed); black line = mean direction.")
add_fig("10e_drift_comparison.png", "Daytime net drift over ground (script 08), swimming relative to CMEMS (09) and the CMEMS current at the schools (09), estimated separately from each sonar.")
add_fig("10f_catch.png", "Left: schools per hour from the two sonars per haul (point size = catch). Right: Spearman correlation of catch and catch rate with schools per hour and backscatter per hour.")
tab <- summ[, .(Sonar = sonar, Detected = f0(detected), Excluded = f0(excluded), Tracked = f0(schools),
                `Moving %` = f0(100 * moving), `Range m` = f0(range_med), `Depth m` = f0(depth_med),
                `Sv dB` = f1(Sv_med), `Speed m/s` = f2(speed_med))]
doc <- body_add_fpar(doc, fpar(ftext("Table 1. ", fp_cap_lab), ftext("School statistics per sonar (medians).", fp_cap),
                               fp_p = fp_par(keep_with_next = TRUE, padding.top = 6)), style = "Table Caption")
ft <- flextable(tab); ft <- theme_booktabs(ft); ft <- fontsize(ft, size = 9, part = "all"); ft <- bold(ft, part = "header")
doc <- body_add_flextable(doc, autofit(ft), align = "center")

h1("Caveats")
lead("Exports, not instruments.", "Both datasets are LSSS Profos exports: range, sector, thresholds and tilt settings used during recording and processing affect detections as much as the instruments do. Differences in numbers and depths should be read as differences between these exports.")
lead("Same thresholds.", "The quality control thresholds were tuned on the ST90 data and applied unchanged to the SU90, so that both are treated the same; they are not necessarily optimal for the SU90.")
lead("Matching.", "Schools are matched on position and time only. Two nearby schools could be swapped, and a large school may be split differently by the two sonars.")
print(doc, target = file.path(out, "SCH302_ST90_vs_SU90.docx"))
message("10: done -> ", normalizePath(out))
