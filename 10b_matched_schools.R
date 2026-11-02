################################################################################
# 10b_matched_schools.R
#
# Schools tracked by BOTH the ST90 and the SU90 at the same place and time
# (June-July 2026). Each pair is the same fish seen by two instruments, so
# differences within a pair are instrument and processing differences.
# Run 10_compare_sonars.R first (it builds the pairs).
#   source("10b_matched_schools.R")
#
#   1. which schools are seen by both: share matched by range, Sv and size
#   2. geometry: the two sonars' tilts, beam overlap and what it implies for depth
#   3. echo level, school size and range of the same school
#   4. velocity: agreement, and the velocity error of one sonar estimated from
#      the paired differences (error variance = var(difference) / 2)
#   5. timing: which sonar picks a school up first
#
# Outputs: output_compare/ figures 10g-10k, 10b_matched_summary.png,
# compare_matched_stats.csv and SCH302_matched_schools.docx.
################################################################################

sonar_profile <- "ST90"
source("00_setup_functions.R")   # run from the RStudio project root
library(officer); library(flextable)

half_beam <- cfg$vert_beam_half   # deg, assumed half vertical beamwidth for both sonars
out  <- "output_compare"
cols_s <- c(ST90 = "#1b7f79", SU90 = "#d9822b")

p <- fread(file.path(out, "compare_matched_schools.csv"))
mov <- rbind(fread("output_ST90/school_movement.csv")[, sonar := "ST90"],
             fread("output_SU90_2026/school_movement.csv")[, sonar := "SU90"], fill = TRUE)
mov[, matched := sid %in% c(p$sid_st, p$sid_su)]
trk <- rbind(readRDS("output_ST90/school_tracks.rds")[, .(sid, tilt)],
             readRDS("output_SU90_2026/school_tracks.rds")[, .(sid, tilt)])
tl <- trk[, .(tilt = median(tilt)), by = sid]
mov <- tl[mov, on = "sid"]
p <- merge(p, setnames(copy(tl), c("sid_st", "tilt_st")), by = "sid_st")
p <- merge(p, setnames(copy(tl), c("sid_su", "tilt_su")), by = "sid_su")

# ---- 1. which schools are seen by both (hours in which both sonars recorded) -----------------------------
mov[, hr := floor(as.numeric(as.POSIXct(t_start, tz = "UTC")) / 3600)]
both_hr <- mov[, uniqueN(sonar), by = hr][V1 == 2, hr]
mb <- mov[hr %in% both_hr]
share <- mb[, .(n = .N, matched = mean(matched)), by = sonar]
bin_share <- function(x, by, breaks, lab) {
  x[, .(n = .N, share = mean(matched)), by = .(sonar, bin = cut(get(by), breaks))][!is.na(bin)][, var := lab][]
}
shares <- rbind(bin_share(mb, "range_m", seq(150, 650, 100), "Range (m)"),
                bin_share(mb, "Sv", c(-70, -50, -45, -40, -35, -30, 0), "Sv (dB)"),
                bin_share(mb[, lA := log10(Area)], "lA", c(0, 3, 3.25, 3.5, 3.75, 4, 6), "log10 area (m2)"))
tilt_mode <- mov[, .N, by = .(sonar, tilt)][order(-N)][, .SD[1], by = sonar]
shares[, var := factor(var, c("Range (m)", "Sv (dB)", "log10 area (m2)"))]

# ---- 2. geometry: beam overlap ------------------------------------------------------------------------------
# beam of each sonar covers tilt +- half_beam; a school seen by both must reach into both beams
p[, `:=`(st_top = -(tilt_st + half_beam), st_bot = -(tilt_st - half_beam),     # depression angles (deg, positive down)
         su_top = -(tilt_su + half_beam), su_bot = -(tilt_su - half_beam))]
p[, `:=`(ov_top = pmax(st_top, su_top), ov_bot = pmin(st_bot, su_bot))]
p[, overlap := ov_bot > ov_top]
p[, rng := (range_m_st + range_m_su) / 2]
p[, `:=`(dep_lo = rng * sin(pmax(0, ov_top) * pi / 180), dep_hi = rng * sin(ov_bot * pi / 180))]
p[overlap == FALSE, `:=`(dep_lo = rng * sin(ov_bot * pi / 180),     # gap between the beams:
                         dep_hi = rng * sin(ov_top * pi / 180))]     # the school spans at least this
p[, dep_mid := (dep_lo + dep_hi) / 2]
share_overlap <- mean(p$overlap)

# ---- 3. echo level, size, range -------------------------------------------------------------------------------
p[, `:=`(dSv = Sv_su - Sv_st, lA = log10(Area_su / Area_st), dR = range_m_su - range_m_st)]
dsv_fit <- lm(dSv ~ range_m_st, p)
cor_sv <- cor(p$Sv_st, p$Sv_su); cor_a <- cor(log(p$Area_st), log(p$Area_su))

# ---- 4. velocity error from the pairs -----------------------------------------------------------------------------
p[, dur := pmin(duration_s_st, duration_s_su)]
p[, dur_bin := cut(dur, c(0, 60, 120, Inf), labels = c("< 1 min", "1-2 min", "> 2 min"))]
err <- function(d) c(sd = sd(d) / sqrt(2), robust = mad(d) / sqrt(2))   # per-sonar error, per component
dv <- c(p$vx_su - p$vx_st, p$vy_su - p$vy_st)
err_all <- err(dv)
err_dur <- p[, .(n = .N, err_robust = mad(c(vx_su - vx_st, vy_su - vy_st)) / sqrt(2),
                 r = cor(c(vx_st, vy_st), c(vx_su, vy_su))), by = dur_bin][order(dur_bin)]
bias <- c(x = mean(p$vx_su - p$vx_st), y = mean(p$vy_su - p$vy_st))
mv <- p[moving_st & moving_su]

# ---- 5. timing --------------------------------------------------------------------------------------------------
p[, dt_onset := as.numeric(as.POSIXct(s0, tz = "UTC")) - as.numeric(as.POSIXct(i.s0, tz = "UTC"))]   # ST90 start - SU90 start (s)
first_st <- mean(p$dt_onset < 0)

stats <- data.table(
  item = c("pairs", "ST90 share matched", "SU90 share matched", "ST90 modal tilt", "SU90 modal tilt", "beams overlap (share of pairs)",
           "median depth bound low (m)", "median depth bound high (m)", "median Sv SU90 - ST90 (dB)", "Sv difference per 100 m range (dB)",
           "cor Sv", "median log10 area ratio SU90/ST90", "cor log area", "median range SU90 - ST90 (m)", "median position distance (m)",
           "velocity error per sonar, robust (m/s)", "velocity error per sonar, SD (m/s)", "mean vx diff", "mean vy diff",
           "median onset ST90 - SU90 (s)", "share ST90 first"),
  value = c(nrow(p), share[sonar == "ST90"]$matched, share[sonar == "SU90"]$matched, tilt_mode[sonar == "ST90"]$tilt,
            tilt_mode[sonar == "SU90"]$tilt, share_overlap, median(p$dep_lo), median(p$dep_hi), median(p$dSv),
            coef(dsv_fit)[[2]] * 100, cor_sv, median(p$lA), cor_a, median(p$dR), median(p$dist), err_all[["robust"]], err_all[["sd"]],
            bias[["x"]], bias[["y"]], median(p$dt_onset), first_st))
fwrite(stats, file.path(out, "compare_matched_stats.csv"))
print(stats[, .(item, value = signif(value, 3))])

# ================================ FIGURES ================================================================
lab_s <- scale_colour_manual(values = cols_s, name = NULL)

# 10g: share of schools seen by both
p10g <- ggplot(shares[n >= 20], aes(bin, share, colour = sonar, group = sonar)) +
  geom_line() + geom_point(aes(size = n)) + facet_wrap(~ var, scales = "free_x") + lab_s +
  scale_size_area(max_size = 3, name = "Schools") + scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Share also tracked by the other sonar", title = "Which schools are seen by both sonars",
       subtitle = "Tracked schools in hours when both sonars recorded; bins with at least 20 schools") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(out, "10g_matched_share.png"), p10g, width = 12, height = 5, dpi = 200)

# 10h: geometry -- beams at the modal tilts and the matched schools
rr <- seq(0, 600, 10)
beam <- rbindlist(lapply(names(cols_s), function(s) { tlt <- tilt_mode[sonar == s]$tilt
  data.table(sonar = s, r = rr, centre = rr * sin(-tlt * pi / 180), lo = rr * sin(-(tlt + half_beam) * pi / 180),
             hi = rr * sin(-(tlt - half_beam) * pi / 180)) }))
pg <- rbind(p[, .(id = .I, r = range_m_st, depth = depth_st, sonar = "ST90")], p[, .(id = .I, r = range_m_su, depth = depth_su, sonar = "SU90")])
p10h <- ggplot() +
  geom_ribbon(data = beam, aes(r, ymin = pmax(lo, 0), ymax = hi, fill = sonar), alpha = 0.15) +
  geom_line(data = beam, aes(r, centre, colour = sonar), linetype = 2) +
  geom_line(data = pg, aes(r, depth, group = id), colour = "grey70", linewidth = 0.2) +
  geom_point(data = pg, aes(r, depth, colour = sonar), size = 0.8, alpha = 0.7) +
  scale_y_reverse() + lab_s + scale_fill_manual(values = cols_s, guide = "none") +
  labs(x = "Range (m)", y = "Depth (m)", title = "Same school, two beams",
       subtitle = sprintf("Bands = beams at the modal tilts (ST90 %d°, SU90 %d°, ±%g° assumed); points = Profos beam-centre depth of each matched school, lines join the pair",
                          tilt_mode[sonar == "ST90"]$tilt, tilt_mode[sonar == "SU90"]$tilt, half_beam)) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
ggsave(file.path(out, "10h_beam_geometry.png"), p10h, width = 10, height = 5.5, dpi = 200)

# 10i: echo level and size
p10i1 <- ggplot(p, aes(Sv_st, Sv_su)) + geom_abline(colour = "grey60", linetype = 2) +
  geom_abline(intercept = median(p$dSv), slope = 1, colour = "firebrick") + geom_point(alpha = 0.35, size = 0.9) + coord_equal() +
  labs(x = "ST90 Sv (dB)", y = "SU90 Sv (dB)", title = "School Sv",
       subtitle = sprintf("r = %.2f; red = 1:1 + %.1f dB (median difference)", cor_sv, median(p$dSv)))
p10i2 <- ggplot(p, aes(range_m_st, dSv)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(alpha = 0.35, size = 0.9) + geom_smooth(method = "lm", formula = y ~ x, colour = "firebrick") +
  labs(x = "Range (m, ST90)", y = "SU90 - ST90 Sv (dB)", title = "Sv difference against range",
       subtitle = sprintf("%+.1f dB per 100 m", coef(dsv_fit)[[2]] * 100))
p10i3 <- ggplot(p, aes(Area_st, Area_su)) + geom_abline(colour = "grey60", linetype = 2) + geom_point(alpha = 0.35, size = 0.9) +
  scale_x_log10() + scale_y_log10() + coord_equal() +
  labs(x = "ST90 area (m2)", y = "SU90 area (m2)", title = "School area",
       subtitle = sprintf("r(log) = %.2f; SU90 median %.0f%% of ST90", cor_a, 100 * 10^median(p$lA)))
ggsave(file.path(out, "10i_echo_level_size.png"), p10i1 + p10i2 + p10i3, width = 15, height = 5, dpi = 200)

# 10j: velocity agreement and error
ve <- rbind(p[, .(dur_bin, comp = "East", st = vx_st, su = vx_su)], p[, .(dur_bin, comp = "North", st = vy_st, su = vy_su)])
p10j1 <- ggplot(ve, aes(abs(su - st), colour = dur_bin)) + stat_ecdf() +
  coord_cartesian(xlim = c(0, 1.5)) + scale_colour_viridis_d(end = 0.85, name = "Shorter track") +
  labs(x = "|SU90 - ST90| velocity component (m/s)", y = "Cumulative share", title = "Velocity difference by track length")
p10j2 <- ggplot(err_dur, aes(dur_bin, err_robust)) + geom_col(fill = "grey40", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f m/s\nn = %d", err_robust, n)), vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Shorter of the two tracks", y = "Velocity error of one sonar (m/s)", title = "Estimated velocity error",
       subtitle = "Robust SD of paired difference / sqrt(2), per component")
ggsave(file.path(out, "10j_velocity_error.png"), p10j1 + p10j2 + plot_layout(widths = c(1.4, 1)), width = 12, height = 5, dpi = 200)

# 10k: onset timing and position
p10k1 <- ggplot(p, aes(pmax(pmin(dt_onset, 120), -120))) + geom_histogram(binwidth = 10, fill = "grey40", colour = "white") +
  geom_vline(xintercept = 0, colour = "firebrick") +
  labs(x = "ST90 start - SU90 start (s; clipped at ±120)", y = "Pairs", title = "Who sees the school first",
       subtitle = sprintf("Median %+.0f s; ST90 first in %.0f%% of pairs", median(p$dt_onset), 100 * first_st))
p10k2 <- ggplot(p, aes(dR)) + geom_histogram(binwidth = 20, fill = "grey40", colour = "white") + geom_vline(xintercept = 0, colour = "firebrick") +
  labs(x = "SU90 - ST90 mean range (m)", y = "Pairs", title = "Range of the same school",
       subtitle = sprintf("Median %+.0f m", median(p$dR)))
ggsave(file.path(out, "10k_timing_range.png"), p10k1 + p10k2, width = 11, height = 4.5, dpi = 200)

p_sum <- (p10g + labs(title = "A  Share of schools seen by both sonars", subtitle = NULL) + theme(legend.position = "right")) /
  ((p10h + labs(title = "B  Two beams, one school", subtitle = NULL)) + (p10i1 + labs(title = "C  School Sv", subtitle = NULL)) +
     plot_layout(widths = c(1.6, 1))) /
  ((p10i2 + labs(title = "D  Sv difference vs range")) + (p10j2 + labs(title = "E  Velocity error of one sonar", subtitle = NULL))) +
  plot_annotation(title = sprintf("%d schools tracked by both ST90 and SU90 (within %d m, same time), 28 Jun - 5 Jul 2026", nrow(p), 100))
ggsave(file.path(out, "10b_matched_summary.png"), p_sum, width = 14, height = 15, dpi = 200)
message("Figures written to ", out)

# ================================ WORD ==========================================================================
f0 <- function(x) formatC(round(x), format = "d", big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
S <- function(i) stats[item == i]$value
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

doc <- body_add_fpar(doc, fpar(ftext("The same schools seen by two sonars", fp_text(bold = TRUE, font.size = 18))))
doc <- body_add_fpar(doc, fpar(ftext(sprintf("SCH302, ST90 and SU90, 28 June – 5 July 2026. Generated %s with 10b_matched_schools.R",
  format(Sys.Date(), "%d %B %Y")), fp_text(font.size = 10, color = "grey40"))))
para(sprintf(paste(
  "The ST90 and SU90 recorded at the same time on this trip (after aligning the SU90 clock by 1 h 03 min 38 s). %s schools were",
  "tracked by both: smoothed tracks overlapping for at least 30 s and centred within 100 m of each other at the common mid time",
  "(median %.0f m apart). Because each pair is the same fish at the same moment, differences within a pair measure how the two",
  "instruments and their Profos processing differ."), f0(nrow(p)), S("median position distance (m)")))

h1("Main findings")
lead("Few schools are seen by both, mainly because the beams point at different depths.", sprintf(paste(
  "In hours when both sonars recorded, %.0f%% of ST90 and %.0f%% of SU90 tracked schools were also tracked by the other sonar.",
  "The ST90 was mostly tilted at %d° and the SU90 at %d°: at 400 m their beam centres lie at about %.0f m and %.0f m depth.",
  "With a vertical beam of about ±%g°, the beams overlap only between about %g° and %g° below the horizontal, so most schools fall into one beam only.",
  "The share seen by both rises with school size and echo strength: large, dense schools reach into both beams."),
  100 * S("ST90 share matched"), 100 * S("SU90 share matched"), S("ST90 modal tilt"), S("SU90 modal tilt"),
  400 * sin(-S("ST90 modal tilt") * pi / 180), 400 * sin(-S("SU90 modal tilt") * pi / 180), half_beam,
  -S("SU90 modal tilt") - half_beam, -S("ST90 modal tilt") + half_beam))
lead("Depth.", sprintf(paste(
  "For the same school the Profos depth is %.0f m (ST90) and %.0f m (SU90) in the median: a pure tilt effect, as Profos",
  "reports the beam-centre depth. A school seen by both sonars must reach into the depth band where the two beams overlap",
  "(beams overlap for %.0f%% of pairs with ±%g° beams): at the pairs' ranges this band lies between about %.0f and %.0f m",
  "(medians). A compact school lies within it; a taller one at least extends into it. This is a better depth indication than",
  "either sonar's beam-centre depth alone, but it depends on the assumed beamwidth."),
  median(p$depth_st), median(p$depth_su), 100 * share_overlap, half_beam, S("median depth bound low (m)"), S("median depth bound high (m)")))
lead("Echo level.", sprintf(paste(
  "The SU90 reports a %.1f dB higher school Sv than the ST90 for the same school (median), and the difference grows by",
  "%.1f dB per 100 m of range. Sv from the two sonars is only moderately correlated (r = %.2f). The constant part is",
  "consistent with calibration and frequency differences; the range trend points to different TVG/absorption settings or",
  "to the school sitting off-axis in the shallower ST90 beam at longer range. Sv from the two sonars is therefore not",
  "interchangeable without an inter-calibration."), S("median Sv SU90 - ST90 (dB)"), S("Sv difference per 100 m range (dB)"), cor_sv))
lead("Size and range.", sprintf(paste(
  "The SU90 school area is %.0f%% of the ST90 area in the median (r = %.2f on log scale), and its school centre is %.0f m",
  "closer to the ship: each beam cuts a different slice of the school."),
  100 * 10^S("median log10 area ratio SU90/ST90"), cor_a, -S("median range SU90 - ST90 (m)")))
lead("Velocity.", sprintf(paste(
  "The two sonars agree on average (mean differences %+.2f and %+.2f m/s east and north), so there is no bias between them.",
  "The scatter between them gives the velocity error of a single sonar: about %.2f m/s per component (robust estimate),",
  "falling from %.2f m/s for tracks shorter than 1 min to %.2f m/s for tracks longer than 2 min. This is about half of the",
  "typical school speed, so individual school directions are uncertain, while averages over many schools (as in the",
  "migration analysis) are reliable. For pairs where both tracks were significantly moving the median direction",
  "difference is %.0f°."),
  bias[["x"]], bias[["y"]], S("velocity error per sonar, robust (m/s)"), err_dur[dur_bin == "< 1 min"]$err_robust,
  err_dur[dur_bin == "> 2 min"]$err_robust, median(abs(mv$ddir))))
lead("Timing.", sprintf(paste(
  "Both sonars start tracking a shared school at almost the same moment (median %+.0f s, ST90 first in %.0f%% of pairs)."),
  S("median onset ST90 - SU90 (s)"), 100 * S("share ST90 first")))
lead("Practical consequences.", paste(
  "(1) Use matching tilts if the two sonars are meant to count the same schools, or deliberately different tilts to cover",
  "two depth layers. (2) Paired detections give school depth bounds that neither sonar gives alone. (3) Inter-calibrate",
  "before comparing Sv or school abundance between the sonars. (4) Treat single-school velocities as ±0.2–0.35 m/s;",
  "rely on population averages. (5) Both work in overlapping nominal frequency bands (SU90 20–30 kHz, ST90 14–24 kHz); pinging",
  "together may cause cross-talk, which could explain part of the SU90's larger tracker jumps and is worth checking."))
add_fig("10b_matched_summary.png", "Summary. A: share of tracked schools also tracked by the other sonar, by range, Sv and size. B: beams at the modal tilts and the beam-centre depths of matched schools (lines join the pair). C: school Sv from both sonars. D: Sv difference against range. E: velocity error of one sonar by track length.")

h1("Details")
add_fig("10g_matched_share.png", "Share of tracked schools also tracked by the other sonar (hours when both recorded), by range, Sv and school area; point size = number of schools.")
add_fig("10h_beam_geometry.png", sprintf("Beams of the two sonars at their modal tilts (bands, ±%g° assumed; dashed = beam centre) and the Profos beam-centre depth of each matched school from each sonar, joined by a line.", half_beam))
add_fig("10i_echo_level_size.png", "School Sv from both sonars (dashed = 1:1, red = 1:1 plus the median difference), Sv difference against range (linear fit) and school area.")
add_fig("10j_velocity_error.png", "Left: cumulative distribution of the absolute difference in each velocity component, by the shorter track's duration. Right: implied velocity error of one sonar (robust SD of the paired difference divided by √2).")
add_fig("10k_timing_range.png", "Left: difference in start time of the two tracks. Right: difference in mean range of the same school.")
tab <- stats[, .(Quantity = item, Value = ifelse(abs(value) < 1 & !grepl("tilt", item), f2(value), f1(value)))]
doc <- body_add_fpar(doc, fpar(ftext("Table 1. ", fp_cap_lab), ftext("Statistics of the matched schools (also in compare_matched_stats.csv).", fp_cap),
                               fp_p = fp_par(keep_with_next = TRUE, padding.top = 6)), style = "Table Caption")
ft <- flextable(tab); ft <- theme_booktabs(ft); ft <- fontsize(ft, size = 9, part = "all"); ft <- bold(ft, part = "header")
doc <- body_add_flextable(doc, width(ft, width = c(3.8, 1)), align = "center")

h1("Caveats")
lead("Beam width.", sprintf("The depth bounds assume a ±%g° vertical beam for both sonars; the real beamwidths depend on frequency and settings and should be taken from the installation.", half_beam))
lead("Matching.", "Pairs are matched on position and time only (within 100 m, overlap ≥ 30 s, closest first). In dense aggregations two neighbouring schools can be paired wrongly, which inflates the differences.")
lead("Profos settings.", "Thresholds and school detection settings in Profos may differ between the two exports; they affect area and Sv as much as the instruments do.")
print(doc, target = file.path(out, "SCH302_matched_schools.docx"))
message("10b: done -> ", normalizePath(out))
