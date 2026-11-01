################################################################################
# 00b_school_qc.R
#
# Flags unrealistic / outlier schools so that scripts 01-07 can exclude them.
# Run once per profile before the other scripts, e.g.
#   sonar_profile <- "ST90"; source("00b_school_qc.R")
#
# Ping rules (dropped before anything else):
#   Sv < ping_sv_min, tilt shallower than ping_tilt_max (near-horizontal beam)
# School rules (any one excludes the school):
#   weak     mean Sv < school_sv_min                 -> not a herring school
#   near     mean range < range_min                  -> near field / wake / surface
#   area     |robust z of log area| > area_z_max     -> implausible size
#   speed    smoothed speed > speed_max              -> faster than herring sustain
#   ship     heading within ship_hd_max of the ship's heading, speed > ship_frac x
#            ship speed, ship steaming > ship_speed_min -> echo moving with the ship
#   jump     step between consecutive cleaned pings > jump_max m/s -> tracker jump
# speed, ship and jump need a track (>= cfg$min_pings, >= cfg$min_duration);
# shorter schools are only tested with weak, near and area.
#
# Outputs (cfg$out_dir): school_qc.csv, qc_a..qc_d figures.
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root
q <- cfg$qc
if (is.null(q)) stop("No QC thresholds (cfg$qc) for profile ", cfg$label)

raw <- read_profos(cfg$profos_dir, "pp")
n_ping0 <- nrow(raw[!is.na(Id)])
pp <- qc_pings(raw)
message(sprintf("Ping rules: %d of %d detections dropped (Sv < %s dB or tilt > %s deg)",
                n_ping0 - nrow(pp[!is.na(Id)]), n_ping0, q$ping_sv_min, q$ping_tilt_max))
m  <- merge_pp(pp)
ss <- school_summary(m)

# ---- tracks for the movement-based rules --------------------------------------------
trk_ids <- ss[n_pings >= cfg$min_pings & duration_s >= cfg$min_duration, sid]
mt <- m[sid %in% trk_ids]
nz <- mt[, { xy <- to_xy(lon, lat); .(sx = noise_sd(xy[, 1]), sy = noise_sd(xy[, 2])) }, by = sid]
sigma_meas <- median(c(nz$sx, nz$sy), na.rm = TRUE)
message(sprintf("Smoothing %d tracks (noise %.1f m) ...", length(trk_ids), sigma_meas))
tr <- rbindlist(lapply(split(mt, by = "sid"), smooth_school, sigma_meas = sigma_meas), fill = TRUE)
mv <- tr[keep == TRUE, track_metrics(.SD), by = sid]
jumps <- tr[keep == TRUE][order(sid, t), .(max_step_speed = max(sqrt(diff(x)^2 + diff(y)^2) / diff(tt))), by = sid]
mv <- merge(mv, jumps, by = "sid")
mv[, hd_diff := abs(((bearing_lm - ship_bearing) + 180) %% 360 - 180)]

# ---- flags -----------------------------------------------------------------------------------
qc <- merge(ss, mv[, .(sid, speed_lm, bearing_lm, ship_speed_trk = ship_speed, hd_diff, max_step_speed)],
            by = "sid", all.x = TRUE)
la <- log(qc$Area)
qc[, `:=`(
  tracked  = sid %in% trk_ids,
  r_weak   = Sv < q$school_sv_min,
  r_near   = range_m < q$range_min,
  r_area   = abs((la - median(la)) / mad(la)) > q$area_z_max,
  r_speed  = !is.na(speed_lm) & speed_lm > q$speed_max,
  r_ship   = !is.na(hd_diff) & hd_diff < q$ship_hd_max & speed_lm > q$ship_frac * ship_speed_trk &
             ship_speed_trk > q$ship_speed_min,
  r_jump   = !is.na(max_step_speed) & max_step_speed > q$jump_max)]
rules <- c("r_weak", "r_near", "r_area", "r_speed", "r_ship", "r_jump")
qc[, excluded := Reduce(`|`, .SD), .SDcols = rules]
qc[, reason := apply(.SD, 1, function(r) paste(sub("r_", "", rules[r]), collapse = "+")), .SDcols = rules]
fwrite(qc, file.path(cfg$out_dir, "school_qc.csv"))

cnt <- qc[, lapply(.SD, sum), .SDcols = rules]
print(cnt)
message(sprintf("Excluded %d of %d schools (%d of %d tracked)", sum(qc$excluded), nrow(qc),
                sum(qc[tracked == TRUE]$excluded), sum(qc$tracked)))
print(qc[tracked == TRUE, .(n = .N, median_speed = median(speed_lm), q95_speed = quantile(speed_lm, 0.95)),
         by = excluded])

# ---- figures ---------------------------------------------------------------------------------
qt <- qc[tracked == TRUE]
qt[, status := factor(ifelse(excluded, "excluded", "kept"), c("kept", "excluded"))]
st_cols <- c(kept = "grey35", excluded = "tomato")

pa1 <- ggplot(qt, aes(speed_lm, fill = status)) +
  geom_histogram(bins = 60, boundary = 0, colour = NA) +
  geom_vline(xintercept = q$speed_max, linetype = 2) +
  scale_x_log10() + scale_fill_manual(values = st_cols, name = NULL) +
  labs(x = "Smoothed school speed (m/s, log)", y = "Schools", title = "School speed",
       subtitle = sprintf("Dashed = %s m/s limit", q$speed_max))
pa2 <- ggplot(qt, aes(ship_speed_trk, speed_lm, colour = status)) +
  geom_point(size = 0.8, alpha = 0.6) + geom_abline(linetype = 3) +
  geom_hline(yintercept = q$speed_max, linetype = 2) +
  scale_colour_manual(values = st_cols, name = NULL, guide = "none") +
  labs(x = "Ship speed during track (m/s)", y = "School speed (m/s)", title = "School vs ship speed",
       subtitle = "Dotted = 1:1 line")
ggsave(file.path(cfg$out_dir, "qc_a_speed.png"), pa1 + pa2 + plot_layout(guides = "collect"),
       width = 11, height = 4.5, dpi = 200)

pb <- ggplot(qt, aes(hd_diff, speed_lm / pmax(ship_speed_trk, 0.1), colour = status)) +
  geom_point(size = 0.8, alpha = 0.6) +
  annotate("rect", xmin = 0, xmax = q$ship_hd_max, ymin = q$ship_frac, ymax = Inf,
           fill = NA, colour = "tomato", linetype = 2) +
  scale_y_log10() + scale_colour_manual(values = st_cols, name = NULL) +
  labs(x = "Difference between school and ship heading (deg)", y = "School speed / ship speed (log)",
       title = "Echoes moving with the ship",
       subtitle = sprintf("Box = 'ship' rule (also requires ship speed > %s m/s)", q$ship_speed_min))
ggsave(file.path(cfg$out_dir, "qc_b_ship_alignment.png"), pb, width = 7, height = 5, dpi = 200)

qc[, status := factor(ifelse(excluded, "excluded", "kept"), c("kept", "excluded"))]
pc <- ggplot(qc[order(status)], aes(lon, lat, colour = status)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = st_cols, name = NULL) + coord_quickmap() +
  labs(x = "Longitude", y = "Latitude", title = "Kept and excluded schools (mean positions)")
ggsave(file.path(cfg$out_dir, "qc_c_map.png"), pc, width = 8, height = 7, dpi = 200)

reasons <- qc[excluded == TRUE, .N, by = reason][order(-N)]
pd <- ggplot(reasons, aes(reorder(reason, N), N)) + geom_col(fill = "tomato") +
  geom_text(aes(label = N), hjust = -0.2, size = 3) + coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Rule(s) triggered", y = "Excluded schools", title = "Why schools were excluded")
ggsave(file.path(cfg$out_dir, "qc_d_reasons.png"), pd, width = 7, height = 4.5, dpi = 200)

message("00b: school_qc.csv and qc figures written to ", cfg$out_dir)
