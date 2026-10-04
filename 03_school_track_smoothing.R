################################################################################
# 03_school_track_smoothing.R
#
# Taming the erratic school detections. Generalises Figure_6_R_code.R from
# one school (Id 1721) to ALL schools, with these changes:
#
#  1. Same-ping duplicates are MERGED (area-weighted) instead of jittering
#     timestamps with runif() -- no random numbers, reproducible.
#  2. Tracker jumps are removed with a 2-D Hampel filter (running median).
#  3. The centroid positioning noise (m) is ESTIMATED from the data
#     (2nd differences), so you know how much smoothing is justified.
#  4. Two smoothers, compared side by side:
#       spline : smooth.spline with ~1 df per tau_s seconds (physical scale,
#                instead of a fixed lambda whose effect varies per track)
#       kalman : constant-velocity Kalman filter + RTS smoother, using the
#                estimated noise; gives velocity AND its uncertainty.
#  5. Velocity significance: is the school really moving, or is the
#     "speed" just noise? (speed > 2 SE, SE corrected for autocorrelation)
#  6. Artefact check: apparent school velocity vs ship velocity.
#
# Outputs: output/school_tracks.rds (point level), output/school_movement.csv
#          (one row per school), output/school_intervals.csv (Figure_6 style),
#          figures 03a-03e.
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root

pp <- merge_pp(read_profos(cfg$profos_dir, "pp"))
ps <- read_profos(cfg$profos_dir, "ps")

# ---- 1. QC: keep schools with enough pings / duration ----------------------------
ss <- school_summary(pp)
good <- ss[n_pings >= cfg$min_pings & duration_s >= cfg$min_duration, sid]
message(length(good), " of ", nrow(ss), " schools pass QC (>= ", cfg$min_pings,
        " pings, >= ", cfg$min_duration, " s)")
pp <- pp[sid %in% good]

# ---- 2. estimate positioning noise -------------------------------------------------
noise <- pp[, {xy <- to_xy(lon, lat); .(sx = noise_sd(xy[, 1]), sy = noise_sd(xy[, 2]),
                                        range_m = mean(range_m), n = .N)}, by = sid]
sigma_meas <- median(c(noise$sx, noise$sy), na.rm = TRUE)
message(sprintf("Estimated centroid positioning noise: %.1f m per axis (median over schools)",
                sigma_meas))

# ---- 3. smooth every school ---------------------------------------------------------
tracks <- rbindlist(lapply(split(pp, by = "sid"), smooth_school,
                           sigma_meas = sigma_meas), fill = TRUE)
saveRDS(tracks, file.path(cfg$out_dir, "school_tracks.rds"))

# ---- 4. track-level metrics -----------------------------------------------------------
mov <- tracks[keep == TRUE, track_metrics(.SD), by = sid]
mov[, n_outliers := tracks[keep == FALSE, .N, by = sid][mov, on = "sid", x.N]]
mov[is.na(n_outliers), n_outliers := 0L]
mov <- merge(ss, mov, by = "sid", suffixes = c("", "_trk"))
# add LSSS's own speed/heading for comparison
mov <- merge(mov, ps[, .(sid, lsss_speed = Speed, lsss_heading = Heading)],
             by = "sid", all.x = TRUE)
sp <- sun_phase(mov$t_start, mov$lon, mov$lat)
mov[, `:=`(sun_elev = sp$sun_elev, sun = sp$sun)]
fwrite(mov, file.path(cfg$out_dir, "school_movement.csv"))

ints <- tracks[, interval_velocity(.SD), by = sid]
fwrite(ints, file.path(cfg$out_dir, "school_intervals.csv"))

cat("\nSpeeds (m/s), schools passing QC:\n")
print(summary(mov[, .(speed_lm, speed_kf, lsss_speed, straightness)]))
cat("Share of schools with significant movement:", round(mean(mov$moving), 2), "\n")

# ---- 03a  raw vs smoothed for example schools --------------------------------------------
ex <- mov[order(-n_pings)][1:min(6, .N), sid]
te <- tracks[sid %in% ex]
# equal x/y span per panel (coord_equal is not allowed with free facet scales)
eq <- te[, {h <- max(diff(range(x)), diff(range(y))) / 2
            .(x = mean(range(x)) + c(-h, h), y = mean(range(y)) + c(-h, h))}, by = sid]
p3a <- ggplot(te, aes(x, y)) +
  geom_blank(data = eq) +
  geom_path(colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = keep), size = 1) +
  geom_path(aes(x_spl, y_spl, linetype = "spline"), colour = "steelblue", linewidth = 0.8,
            data = te[keep == TRUE]) +
  geom_path(aes(x_kf, y_kf, linetype = "kalman"), colour = "tomato", linewidth = 0.8,
            data = te[keep == TRUE]) +
  geom_point(data = te[keep == TRUE, .SD[1], by = sid], aes(x_kf, y_kf), shape = 17, size = 3) +
  facet_wrap(~ sid, scales = "free") +
  theme(aspect.ratio = 1) +
  scale_colour_manual(values = c(`TRUE` = "grey30", `FALSE` = "red"),
                      labels = c(`TRUE` = "kept", `FALSE` = "outlier"), name = "Ping") +
  scale_linetype_manual(values = c(spline = 1, kalman = 1), name = "Smoother",
                        guide = guide_legend(override.aes = list(colour = c("tomato", "steelblue")))) +
  labs(x = "East (m)", y = "North (m)",
       title = "Raw school centroids vs smoothed tracks (6 longest tracks)",
       subtitle = sprintf("Triangle = start. Positioning noise ~ %.1f m; spline tau = %s s; Kalman sigma_acc = %s m/s2",
                          sigma_meas, cfg$tau_s, cfg$sigma_acc))
ggsave(file.path(cfg$out_dir, "03a_raw_vs_smoothed_tracks.png"), p3a, width = 12, height = 8, dpi = 200)

# ---- 03b  speed over time for the same schools: raw vs smoothed --------------------------
te[, raw_speed := c(NA, sqrt(diff(x)^2 + diff(y)^2) / diff(tt)), by = sid]
p3b <- ggplot(te[keep == TRUE], aes(tt)) +
  geom_line(aes(y = raw_speed, colour = "raw ping-to-ping"), alpha = 0.6) +
  geom_line(aes(y = sqrt(vx_spl^2 + vy_spl^2), colour = "spline")) +
  geom_ribbon(aes(ymin = pmax(0, sqrt(vx_kf^2 + vy_kf^2) - 2 * v_sd_kf),
                  ymax = sqrt(vx_kf^2 + vy_kf^2) + 2 * v_sd_kf), fill = "tomato", alpha = 0.2) +
  geom_line(aes(y = sqrt(vx_kf^2 + vy_kf^2), colour = "kalman (+/- 2 SD)")) +
  facet_wrap(~ sid, scales = "free_x") +
  coord_cartesian(ylim = c(0, 6)) +
  scale_colour_manual(values = c(`raw ping-to-ping` = "grey50", spline = "steelblue",
                                 `kalman (+/- 2 SD)` = "tomato"), name = NULL) +
  labs(x = "Seconds since first detection", y = "Speed (m/s)",
       title = "Raw ping-to-ping speed is dominated by noise; smoothed speeds are plausible")
ggsave(file.path(cfg$out_dir, "03b_speed_raw_vs_smoothed.png"), p3b, width = 12, height = 7, dpi = 200)

# ---- 03c  noise vs range + effect of the smoothing scale -----------------------------------
p3c1 <- ggplot(noise, aes(range_m, (sx + sy) / 2)) +
  geom_point(alpha = 0.6) + geom_smooth(method = "loess", formula = y ~ x, colour = "tomato") +
  labs(x = "Mean range from ship (m)", y = "Positioning noise (m)",
       title = "Centroid noise vs range")
# sensitivity: median speed for different spline time scales
sens <- rbindlist(lapply(c(5, 10, 20, 30, 60, 120), function(tau) {
  v <- tracks[keep == TRUE, {
    sx <- smooth_spline_track(tt, x, tau); sy <- smooth_spline_track(tt, y, tau)
    .(speed = mean(sqrt(sx$vel^2 + sy$vel^2)))
  }, by = sid]
  data.table(tau_s = tau, median_speed = median(v$speed), q25 = quantile(v$speed, .25),
             q75 = quantile(v$speed, .75))
}))
p3c2 <- ggplot(sens, aes(tau_s, median_speed)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.2) +
  geom_line(colour = "steelblue") + geom_point(colour = "steelblue") +
  geom_hline(yintercept = median(mov$speed_lm), linetype = 2) +
  scale_x_log10() +
  labs(x = "Spline smoothing scale tau (s, log)", y = "Mean school speed (m/s)",
       title = "Speed depends on smoothing",
       subtitle = "Plateau = noise removed; dashed = straight-line (regression) speed")
p3c <- p3c1 + p3c2
ggsave(file.path(cfg$out_dir, "03c_noise_and_smoothing_sensitivity.png"), p3c, width = 11, height = 4.5, dpi = 200)

# ---- 03d  is the movement real? ---------------------------------------------------------------
p3d1 <- ggplot(mov, aes(duration_s, speed_lm, colour = moving)) +
  geom_errorbar(aes(ymin = pmax(0, speed_lm - 2 * speed_se), ymax = speed_lm + 2 * speed_se),
                alpha = 0.4, width = 0) +
  geom_point() + scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"),
                      name = "Speed > 2 SE") +
  labs(x = "Track duration (s, log)", y = "Speed (m/s)",
       title = "Short tracks cannot resolve speed")
p3d2 <- ggplot(mov, aes(lsss_speed, speed_lm, colour = moving)) +
  geom_abline(linetype = 2) + geom_point() +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60"), guide = "none") +
  labs(x = "LSSS Profos speed (ps file, m/s)", y = "This script (m/s)",
       title = "Comparison with LSSS school speed")
ggsave(file.path(cfg$out_dir, "03d_movement_significance.png"), p3d1 + p3d2,
       width = 11, height = 4.5, dpi = 200)

# ---- 03e  artefact check: does apparent school motion follow the ship? -----------------------
mov[, `:=`(v_along = vx * sin(ship_bearing * pi / 180) + vy * cos(ship_bearing * pi / 180),
           v_across = vx * cos(ship_bearing * pi / 180) - vy * sin(ship_bearing * pi / 180),
           rel_dir = (bearing_lm - ship_bearing) %% 360)]
fit_art <- lm(v_along ~ ship_speed, mov[moving == TRUE])
print(summary(fit_art)$coefficients)
p3e1 <- ggplot(mov[moving == TRUE], aes(ship_speed, v_along)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point() + geom_smooth(method = "lm", formula = y ~ x) +
  labs(x = "Ship speed during track (m/s)", y = "School velocity along ship heading (m/s)",
       title = "Artefact check",
       subtitle = sprintf("Slope = %.2f (p = %.2g). A clear +/- slope suggests positioning bias",
                          coef(fit_art)[2], summary(fit_art)$coefficients[2, 4]))
p3e2 <- ggplot(mov[moving == TRUE], aes(rel_dir)) +
  geom_histogram(breaks = seq(0, 360, 20), fill = "steelblue", colour = "white") +
  coord_polar(start = 0) +
  scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270),
                     labels = c("ahead", "stbd", "astern", "port")) +
  labs(x = NULL, y = NULL, title = "School direction relative to ship heading")
ggsave(file.path(cfg$out_dir, "03e_ship_motion_artefact_check.png"), p3e1 + p3e2,
       width = 11, height = 5, dpi = 200)

message("03: tracks, movement tables and figures written to ", cfg$out_dir)
