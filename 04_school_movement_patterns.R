################################################################################
# 04_school_movement_patterns.R
#
# School movement patterns across all tracked schools (run 03 first).
#   Fig 04a  Map of movement vectors (smoothed tracks + arrow = 2 min of
#            swimming), ship track underneath, tow lines.
#   Fig 04b  Rose diagrams of swimming direction, by fishing phase and by
#            sun phase, with circular stats (mean direction, Rbar, Rayleigh p).
#   Fig 04c  Speed vs depth / school size / time of day.
#   Fig 04d  Gridded mean movement field (cells of grid_km): population-
#            level flow, robust to individual noisy tracks.
#   Fig 04e  Figure_6-style panel for ALL schools (interval speeds and
#            directions), using smoothed positions.
#   Table    output/movement_circular_stats.csv
#
# Only schools whose movement is significant (moving == TRUE) are used for
# directions: a "direction" of a school that is not resolvably moving is noise.
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root

grid_km <- 1   # cell size for the gridded movement field

hauls  <- read_catch(cfg$catch_file, trips = cfg$trip)
mov    <- fread(file.path(cfg$out_dir, "school_movement.csv"))
tracks <- readRDS(file.path(cfg$out_dir, "school_tracks.rds"))
ints   <- fread(file.path(cfg$out_dir, "school_intervals.csv"))
ship   <- ship_track(read_profos(cfg$profos_dir, "ppe"))

mov[, t_start := as.POSIXct(t_start, tz = "UTC")]
mov <- cbind(mov, label_phase(mov$t_start, hauls))
mov[, sun := factor(sun, c("day", "twilight", "night"))]
mov[, phase := factor(phase, c("search", "tow", "other"))]
mm <- mov[moving == TRUE]
hz <- hauls[haul_time >= min(ship$t) & shoot_time <= max(ship$t)]

# ---- 04a  movement vector map ----------------------------------------------------------
arrow_s <- 600
mm[, `:=`(lon_end = lon + vx * arrow_s / (111320 * cos(lat * pi / 180)),
          lat_end = lat + vy * arrow_s / 111320)]
p4a <- ggplot() +
  geom_path(data = ship[seq(1, .N, 5)], aes(lon, lat, group = seg), colour = "grey75") +
  geom_segment(data = hz, aes(shootlon, shootlat, xend = haullon, yend = haullat),
               colour = "red", linetype = 2) +
  geom_path(data = tracks[keep == TRUE], aes(lon_kf, lat_kf, group = sid),
            colour = "grey30", linewidth = 0.3) +
  geom_segment(data = mm, aes(lon, lat, xend = lon_end, yend = lat_end, colour = speed_lm),
               arrow = arrow(length = unit(1.5, "mm")), linewidth = 0.7) +
  scale_colour_viridis_c(option = "plasma", name = "Speed (m/s)") +
  coord_quickmap() +
  labs(x = "Longitude", y = "Latitude",
       title = "School movement vectors (arrow = distance swum in 10 min)",
       subtitle = "Grey = ship track; dark lines = smoothed school tracks; red dashed = tow lines")
ggsave(file.path(cfg$out_dir, "04a_movement_vector_map.png"), p4a, width = 9, height = 8, dpi = 200)

# ---- 04b  rose diagrams + circular stats -------------------------------------------------
stats <- rbind(
  cbind(group = "all", mm[, circ_summary(bearing_lm)]),
  mm[, circ_summary(bearing_lm), by = .(group = paste("phase:", phase))],
  mm[, circ_summary(bearing_lm), by = .(group = paste("sun:", sun))],
  mm[, circ_summary(bearing_lm), by = .(group = paste("haul:", haul_id))][!grepl("NA", group)]
)
print(stats)
fwrite(stats, file.path(cfg$out_dir, "movement_circular_stats.csv"))

rose <- function(d, facet, st) {
  lab <- st[, setNames(sprintf("%s
n=%d, mean=%.0f deg, Rbar=%.2f, p=%.2g",
                               sub(".*: ", "", group), n, mean_dir, Rbar, rayleigh_p),
                       sub(".*: ", "", group))]
  ggplot(d, aes(bearing_lm)) +
    geom_histogram(breaks = seq(0, 360, 22.5), fill = "steelblue", colour = "white") +
    coord_polar(start = 0) +
    scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270),
                       labels = c("N", "E", "S", "W")) +
    facet_wrap(as.formula(paste("~", facet)), nrow = 1, labeller = as_labeller(lab)) +
    labs(x = NULL, y = "Schools")
}
p4b <- rose(mm, "phase", stats[grepl("^phase", group)]) /
  rose(mm, "sun", stats[grepl("^sun", group)]) +
  plot_annotation(title = "Swimming direction of schools with significant movement",
                  subtitle = "Rbar: 0 = random, 1 = all same direction; p = Rayleigh test of uniformity")
ggsave(file.path(cfg$out_dir, "04b_direction_roses.png"), p4b, width = 10, height = 8, dpi = 200)

# ---- 04c  speed vs covariates ------------------------------------------------------------
p_d <- ggplot(mm, aes(range_m, speed_lm)) + geom_point(aes(colour = sun), alpha = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "black") +
  labs(x = "Range from ship (m)", y = "Speed (m/s)", title = "Speed vs range",
       subtitle = "(range, not depth: Profos depth = range x tilt)")
p_a <- ggplot(mm, aes(Area, speed_lm)) + geom_point(aes(colour = sun), alpha = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "black") + scale_x_log10() +
  labs(x = "School area (m2, log)", y = "Speed (m/s)", title = "Speed vs school size")
p_t <- ggplot(mm, aes(t_start, speed_lm)) +
  geom_rect(data = hz, aes(xmin = shoot_time, xmax = haul_time, ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 0.1, inherit.aes = FALSE) +
  geom_point(aes(colour = sun), alpha = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "black", span = 0.5) +
  labs(x = "Time (UTC)", y = "Speed (m/s)", title = "Speed through the day (red = tows)")
p4c <- (p_d + p_a) / p_t + plot_layout(guides = "collect") &
  scale_colour_manual(values = c(day = "gold3", twilight = "orchid", night = "navy"),
                      drop = FALSE)
ggsave(file.path(cfg$out_dir, "04c_speed_vs_covariates.png"), p4c, width = 11, height = 8, dpi = 200)

# ---- 04d  gridded movement field ------------------------------------------------------------
lon0 <- mean(mm$lon); lat0 <- mean(mm$lat)
xy <- to_xy(mm$lon, mm$lat, lon0, lat0)
mm[, `:=`(gx = floor(xy[, 1] / (grid_km * 1000)), gy = floor(xy[, 2] / (grid_km * 1000)))]
grid <- mm[, .(n = .N, vx = mean(vx), vy = mean(vy),
               Rbar = circ_summary(bearing_lm)$Rbar), by = .(gx, gy)]
cc <- from_xy((grid$gx + 0.5) * grid_km * 1000, (grid$gy + 0.5) * grid_km * 1000, lon0, lat0)
grid[, `:=`(lon = cc[, 1], lat = cc[, 2], speed = sqrt(vx^2 + vy^2))]
grid[is.na(Rbar), Rbar := 0]   # single-school cells: direction agreement undefined
sc <- 0.4 * grid_km * 1000 / max(grid$speed)   # longest arrow = 40% of a cell
grid[, `:=`(lon_end = lon + vx * sc / (111320 * cos(lat * pi / 180)),
            lat_end = lat + vy * sc / 111320)]
p4d <- ggplot(grid) +
  geom_path(data = ship[seq(1, .N, 5)], aes(lon, lat, group = seg), colour = "grey80") +
  geom_tile(aes(lon, lat, fill = n), width = grid_km * 1000 / (111320 * cos(lat0 * pi / 180)),
            height = grid_km * 1000 / 111320, alpha = 0.5) +
  geom_segment(aes(lon, lat, xend = lon_end, yend = lat_end, linewidth = Rbar),
               arrow = arrow(length = unit(1.5, "mm"))) +
  scale_fill_viridis_c(name = "Schools", option = "mako", direction = -1) +
  scale_linewidth(range = c(0.2, 1.2), limits = c(0, 1), name = "Rbar") +
  coord_quickmap() +
  labs(x = "Longitude", y = "Latitude",
       title = sprintf("Mean school movement per %s km cell", grid_km),
       subtitle = "Arrow = mean velocity (relative length); thick = schools agree on direction")
ggsave(file.path(cfg$out_dir, "04d_gridded_movement_field.png"), p4d, width = 9, height = 8, dpi = 200)

# ---- 04e  Figure_6 style for all schools (interval level) --------------------------------------
ints <- ints[sid %in% mm$sid]
pA <- ggplot(melt(ints, measure.vars = c("speed", "bearing")), aes(variable, value, fill = variable)) +
  geom_boxplot(width = 0.4, outlier.size = 0.4) +
  facet_wrap(~ variable, scales = "free",
             labeller = as_labeller(c(speed = "Speed (m/s)", bearing = "Direction (deg)"))) +
  scale_fill_manual(values = c("steelblue", "tomato"), guide = "none") +
  labs(x = NULL, y = NULL, title = "A") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
pB <- ggplot(ints, aes(bearing, speed)) +
  geom_bin2d(bins = c(36, 30)) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", name = "Intervals") +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  labs(x = "Direction (deg)", y = "Speed (m/s)", title = "B")
ggsave(file.path(cfg$out_dir, "04e_figure6_style_all_schools.png"), pA + pB,
       width = 12, height = 5.5, dpi = 200)

message("04: movement figures + movement_circular_stats.csv written to ", cfg$out_dir)
