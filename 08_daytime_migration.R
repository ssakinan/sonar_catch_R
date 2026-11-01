################################################################################
# 08_daytime_migration.R
#
# Is there a migration signal in the DAYTIME school movements?
# (run 00b and 03 first; intended for the ST90 profile)
#   sonar_profile <- "ST90"; source("08_daytime_migration.R")
#
# The ground velocity of a tracked school (vx east, vy north; m/s) is modelled as
#
#   v = net drift + semidiurnal tide + vessel effect + noise
#   vx = a_x + c_x cos(wt) + s_x sin(wt) + b_x1 ship_vx + b_x2 ship_vy
#   vy = a_y + c_y cos(wt) + s_y sin(wt) + b_y1 ship_vx + b_y2 ship_vy
#
#   w = 2 pi / 12.42 h (M2 tide; with 8 days of data M2 and S2 cannot be separated)
#   ship terms absorb the small apparent dependence of school velocity on ship
#   velocity found in script 03 (vessel-related echoes or avoidance).
#   (a_x, a_y) = net residual drift = the migration candidate.
#
# Uncertainty: block bootstrap over 1-hour blocks (schools seen in the same hour
# are not independent). Robustness: model variants, leave-one-day-out, daily and
# spatial (grid) residual vectors, and a progressive vector diagram.
#
# Outputs (cfg$out_dir/migration): 08a-08f figures, 08_migration_summary.png,
# 08_daytime_flow.mp4, migration_estimates.csv, and a Word summary.
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root
library(officer); library(flextable)

B          <- 1000                 # bootstrap replicates
w          <- 2 * pi / (12.4206 * 3600)
grid_km_m  <- 5                    # cell size for the spatial residual map
min_cell_n <- 10                   # schools per cell to draw an arrow
make_anim  <- !identical(Sys.getenv("SKIP_ANIM"), "1")   # SKIP_ANIM=1: figures + document only
set.seed(1)

out <- file.path(sub("/report_figs$", "", cfg$out_dir), "migration")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- data: daytime tracked schools -----------------------------------------------------------
m <- fread(file.path(sub("/report_figs$", "", cfg$out_dir), "school_movement.csv"))
m[, t := as.POSIXct(t_start, tz = "UTC") + duration_s / 2]      # track mid time
d <- m[sun == "day"]
d[, `:=`(tt = as.numeric(t), hb = floor(as.numeric(t) / 3600), day = as.Date(t))]
d[, `:=`(cw = cos(w * tt), sw = sin(w * tt))]
message(sprintf("Daytime tracked schools: %d (%d with significant movement), %d hour blocks, %d days",
                nrow(d), sum(d$moving), uniqueN(d$hb), uniqueN(d$day)))

# ---- model fitting + block bootstrap ----------------------------------------------------------
terms_full <- c("cw", "sw", "ship_vx", "ship_vy")
fit_flow <- function(x, terms) {
  X <- if (length(terms)) cbind(1, as.matrix(x[, ..terms])) else matrix(1, nrow(x), 1)
  cx <- lm.fit(X, x$vx)$coefficients; cy <- lm.fit(X, x$vy)$coefficients
  names(cx) <- names(cy) <- c("a", terms)
  list(x = cx, y = cy)
}
boot_flow <- function(x, terms, B) {
  blocks <- split(seq_len(nrow(x)), x$hb)
  t(replicate(B, {
    idx <- unlist(blocks[sample(length(blocks), replace = TRUE)], use.names = FALSE)
    f <- fit_flow(x[idx], terms); c(ax = f$x[["a"]], ay = f$y[["a"]])
  }))
}
summ_vec <- function(ax, ay, bt, label, n) {
  sp <- sqrt(ax^2 + ay^2); dir <- bearing_of(ax, ay)
  bsp <- sqrt(bt[, "ax"]^2 + bt[, "ay"]^2); bdir <- bearing_of(bt[, "ax"], bt[, "ay"])
  ddev <- ((bdir - dir + 180) %% 360) - 180
  data.table(model = label, n = n, vx = ax, vy = ay, speed = sp, km_day = sp * 86.4, dir = dir,
             speed_lo = quantile(bsp, 0.025), speed_hi = quantile(bsp, 0.975),
             dir_lo = (dir + quantile(ddev, 0.025)) %% 360, dir_hi = (dir + quantile(ddev, 0.975)) %% 360,
             p_zero = mean(bt[, "ax"] * ax + bt[, "ay"] * ay <= 0))   # share of replicates pointing away
}
variant <- function(x, terms, label) {
  f <- fit_flow(x, terms); bt <- boot_flow(x, terms, B)
  list(fit = f, boot = bt, sum = summ_vec(f$x[["a"]], f$y[["a"]], bt, label, nrow(x)))
}
message("Fitting models and bootstrapping ...")
v_raw   <- variant(d, character(0), "Raw mean velocity")
v_tide  <- variant(d, c("cw", "sw"), "Tide removed")
v_main  <- variant(d, terms_full, "Tide + vessel removed (main)")
v_slow  <- variant(d[ship_speed_trk < 2.5], terms_full, "Main, ship speed < 2.5 m/s")
v_mov   <- variant(d[moving == TRUE], terms_full, "Main, moving schools only")
loo <- rbindlist(lapply(sort(unique(d$day)), function(dd) {
  x <- d[day != dd]; f <- fit_flow(x, terms_full); bt <- boot_flow(x, terms_full, 200)
  summ_vec(f$x[["a"]], f$y[["a"]], bt, paste("Without", format(dd, "%d %b")), nrow(x))
}))
est <- rbind(v_raw$sum, v_tide$sum, v_main$sum, v_slow$sum, v_mov$sum, loo)
fwrite(est, file.path(out, "migration_estimates.csv"))
print(est[, .(model, n, speed = round(speed, 3), km_day = round(km_day, 1), dir = round(dir),
              dir_lo = round(dir_lo), dir_hi = round(dir_hi), p_zero)])

fm <- v_main$fit
m2_amp <- c(x = sqrt(fm$x[["cw"]]^2 + fm$x[["sw"]]^2), y = sqrt(fm$y[["cw"]]^2 + fm$y[["sw"]]^2))
r2 <- function(comp) { X <- cbind(1, as.matrix(d[, ..terms_full])); f <- lm.fit(X, d[[comp]])
  1 - sum(f$residuals^2) / sum((d[[comp]] - mean(d[[comp]]))^2) }
tide_p <- c(x = anova(lm(vx ~ ship_vx + ship_vy, d), lm(vx ~ cw + sw + ship_vx + ship_vy, d))[2, "Pr(>F)"],
            y = anova(lm(vy ~ ship_vx + ship_vy, d), lm(vy ~ cw + sw + ship_vx + ship_vy, d))[2, "Pr(>F)"])

# de-tided, vessel-corrected velocity of every school (= net drift + its own residual)
pred_nonnet <- function(f, x, comp) as.matrix(x[, ..terms_full]) %*% f[[comp]][terms_full]
d[, `:=`(rx = vx - pred_nonnet(fm, .SD, "x")[, 1], ry = vy - pred_nonnet(fm, .SD, "y")[, 1])]
d[, `:=`(tide_x = fm$x[["cw"]] * cw + fm$x[["sw"]] * sw, tide_y = fm$y[["cw"]] * cw + fm$y[["sw"]] * sw)]

# ---- daily residual vectors -------------------------------------------------------------------------
daily <- d[, {
  bl <- split(seq_len(.N), hb)
  bt <- t(replicate(500, { i <- unlist(bl[sample(length(bl), replace = TRUE)]); c(mean(rx[i]), mean(ry[i])) }))
  .(n = .N, hours = length(bl), rx = mean(rx), ry = mean(ry), vx = mean(vx), vy = mean(vy),
    rx_lo = quantile(bt[, 1], .025), rx_hi = quantile(bt[, 1], .975),
    ry_lo = quantile(bt[, 2], .025), ry_hi = quantile(bt[, 2], .975))
}, by = day][order(day)]
daily[, `:=`(dir = bearing_of(rx, ry), km_day = sqrt(rx^2 + ry^2) * 86.4)]
fwrite(daily, file.path(out, "migration_daily.csv"))
daily_dirs <- circ_summary(daily[n >= 50]$dir)

# ---- spatial residual field ----------------------------------------------------------------------------
lon0 <- mean(d$lon); lat0 <- mean(d$lat)
xy <- to_xy(d$lon, d$lat, lon0, lat0)
d[, `:=`(gx = floor(xy[, 1] / (grid_km_m * 1000)), gy = floor(xy[, 2] / (grid_km_m * 1000)))]
cells <- d[, .(n = .N, rx = mean(rx), ry = mean(ry)), by = .(gx, gy)][n >= min_cell_n]
cc <- from_xy((cells$gx + 0.5) * grid_km_m * 1000, (cells$gy + 0.5) * grid_km_m * 1000, lon0, lat0)
cells[, `:=`(lon = cc[, 1], lat = cc[, 2], dir = bearing_of(rx, ry), spd = sqrt(rx^2 + ry^2))]
cell_dirs <- circ_summary(cells$dir, w = cells$n)
share_w <- cells[, sum(n[dir > 180 & dir < 360]) / sum(n)]   # share of schools in cells drifting westward

# ---- progressive vector diagram (hourly means, observed hours only) ------------------------------------
hourly <- d[, .(t = mean(t), n = .N, vx = mean(vx), vy = mean(vy), rx = mean(rx), ry = mean(ry),
                tide_x = mean(tide_x), tide_y = mean(tide_y), lon = mean(lon), lat = mean(lat)), by = hb][order(hb)]
hourly[, `:=`(px_raw = cumsum(vx) * 3.6, py_raw = cumsum(vy) * 3.6,      # km (1 h per block)
              px_res = cumsum(rx) * 3.6, py_res = cumsum(ry) * 3.6)]
obs_h <- nrow(hourly)

# ================================ FIGURES ================================================================
sch_arrow_s <- 7200   # school arrows: distance covered in 2 h
k_lon <- function(lat) 111320 * cos(lat * pi / 180)
day_cols <- setNames(scales::hue_pal()(nrow(daily)), format(daily$day, "%d %b"))
ship_bg <- read_ship()[seq(1, .N, 10)]

# 08a: flow map -----------------------------------------------------------------------------------------
cell_scale <- 0.8 * grid_km_m * 1000 / max(cells$spd)
cells[, `:=`(lon_e = lon + rx * cell_scale / k_lon(lat), lat_e = lat + ry * cell_scale / 111320)]
d[, `:=`(lon_e = lon + vx * sch_arrow_s / k_lon(lat), lat_e = lat + vy * sch_arrow_s / 111320)]
p8a <- ggplot() +
  geom_path(data = ship_bg, aes(lon, lat, group = seg), colour = "grey85") +
  geom_segment(data = d, aes(lon, lat, xend = lon_e, yend = lat_e), colour = "grey55", alpha = 0.35,
               linewidth = 0.25) +
  geom_segment(data = cells, aes(lon, lat, xend = lon_e, yend = lat_e, colour = dir > 180),
               arrow = arrow(length = unit(1.6, "mm")), linewidth = 0.9) +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue"),
                      labels = c(`TRUE` = "westward", `FALSE` = "eastward"), name = "Net drift") +
  coord_quickmap() +
  labs(x = "Longitude", y = "Latitude", title = "Daytime school movement",
       subtitle = sprintf("Grey = each school (distance in 2 h, observed velocity); arrows = mean net drift per %d-km cell (tide and vessel effect removed, n >= %d)",
                          grid_km_m, min_cell_n))
ggsave(file.path(out, "08a_flow_map.png"), p8a, width = 9, height = 8, dpi = 200)

# 08b: tidal signal -----------------------------------------------------------------------------------------
hl <- melt(hourly[, .(t, n, east = vx, north = vy, tide_east = tide_x + fm$x[["a"]], tide_north = tide_y + fm$y[["a"]])],
           id.vars = c("t", "n"))
hl[, comp := ifelse(grepl("east", variable), "East velocity (m/s)", "North velocity (m/s)")]
hl[, kind := ifelse(grepl("tide", variable), "fit", "obs")]
p8b1 <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(data = hl[kind == "obs"], aes(t, value, size = n), colour = "grey30", alpha = 0.6) +
  geom_line(data = hl[kind == "fit"], aes(t, value, group = interaction(comp, as.Date(t))), colour = "tomato") +
  facet_wrap(~ comp, ncol = 1) + scale_size_area(max_size = 2.5, name = "Schools\nper hour") +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%d %b") +
  labs(x = NULL, y = NULL, title = "Hourly mean school velocity and fitted tide",
       subtitle = "Red = net drift + semidiurnal (12.42 h) tidal component")
ph <- seq(0, 2 * pi, length = 100)
ell <- data.table(x = fm$x[["cw"]] * cos(ph) + fm$x[["sw"]] * sin(ph), y = fm$y[["cw"]] * cos(ph) + fm$y[["sw"]] * sin(ph))
p8b2 <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(data = hourly, aes(vx, vy), colour = "grey50", alpha = 0.5, size = 1) +
  geom_path(data = ell, aes(x + fm$x[["a"]], y + fm$y[["a"]]), colour = "tomato", linewidth = 0.8) +
  annotate("segment", x = 0, y = 0, xend = fm$x[["a"]], yend = fm$y[["a"]], colour = "black",
           arrow = arrow(length = unit(2, "mm")), linewidth = 1) +
  coord_equal() +
  labs(x = "East (m/s)", y = "North (m/s)", title = "Tidal ellipse",
       subtitle = "Grey = hourly means\nred = fitted tidal cycle\narrow = net drift")
ggsave(file.path(out, "08b_tidal_signal.png"), p8b1 + p8b2 + plot_layout(widths = c(2.2, 1)),
       width = 13, height = 5.5, dpi = 200)

# 08c: daily net drift vectors -------------------------------------------------------------------------------
daily[, lab := format(day, "%d %b")]
lim <- max(abs(unlist(daily[n >= 20, .(rx_lo, rx_hi, ry_lo, ry_hi)]))) * 86.4 * 1.1
p8c <- ggplot(daily[n >= 20]) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_errorbar(aes(x = rx * 86.4, ymin = ry_lo * 86.4, ymax = ry_hi * 86.4, colour = lab), width = 0, alpha = 0.6) +
  geom_errorbarh(aes(y = ry * 86.4, xmin = rx_lo * 86.4, xmax = rx_hi * 86.4, colour = lab), height = 0, alpha = 0.6) +
  geom_segment(aes(x = 0, y = 0, xend = rx * 86.4, yend = ry * 86.4, colour = lab),
               arrow = arrow(length = unit(2, "mm")), linewidth = 0.9) +
  ggrepel::geom_text_repel(aes(rx * 86.4, ry * 86.4, label = sprintf("%s (n=%d)", lab, n), colour = lab),
            min.segment.length = 0, seed = 1,
            vjust = -0.8, size = 3, show.legend = FALSE) +
  annotate("segment", x = 0, y = 0, xend = fm$x[["a"]] * 86.4, yend = fm$y[["a"]] * 86.4,
           arrow = arrow(length = unit(3, "mm")), linewidth = 1.6) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  scale_colour_manual(values = day_cols, guide = "none") +
  labs(x = "East (km/day)", y = "North (km/day)", title = "Net drift per day",
       subtitle = "Coloured = daily mean (tide and vessel effect removed)\nwith 95% bootstrap intervals; black = all days")
ggsave(file.path(out, "08c_daily_vectors.png"), p8c, width = 7.5, height = 7.5, dpi = 200)

# 08d: estimates across models -------------------------------------------------------------------------------
est[, model := factor(model, rev(model))]
est[, grp := ifelse(grepl("^Without", model), "Leave one day out", "Model variants")]
p8d1 <- ggplot(est, aes(km_day, model, colour = grp)) +
  geom_errorbarh(aes(xmin = speed_lo * 86.4, xmax = speed_hi * 86.4), height = 0.3) + geom_point(size = 2) +
  scale_colour_manual(values = c("Model variants" = "black", "Leave one day out" = "grey50"), name = NULL) +
  labs(x = "Net drift (km/day)", y = NULL, title = "Net drift speed")
unwrap <- function(x, ref) ref + ((x - ref + 180) %% 360) - 180
est[, `:=`(dir_u = unwrap(dir, v_main$sum$dir), lo_u = unwrap(dir_lo, v_main$sum$dir), hi_u = unwrap(dir_hi, v_main$sum$dir))]
p8d2 <- ggplot(est, aes(dir_u, model, colour = grp)) +
  geom_errorbarh(aes(xmin = lo_u, xmax = hi_u), height = 0.3) + geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0, 540, 15), labels = function(b) b %% 360) +
  scale_colour_manual(values = c("Model variants" = "black", "Leave one day out" = "grey50"), guide = "none") +
  labs(x = "Direction of net drift (deg, compass)", y = NULL, title = "Net drift direction") +
  theme(axis.text.y = element_blank())
ggsave(file.path(out, "08d_net_drift_estimates.png"), p8d1 + p8d2 + plot_layout(guides = "collect"),
       width = 12, height = 5, dpi = 200)

# 08e: progressive vector diagram --------------------------------------------------------------------------------
pvd <- rbind(hourly[, .(t, x = px_raw, y = py_raw, what = "Observed velocity")],
             hourly[, .(t, x = px_res, y = py_res, what = "Tide and vessel effect removed")])
pvd <- rbind(data.table(t = min(hourly$t), x = 0, y = 0, what = unique(pvd$what)), pvd)
p8e <- ggplot(pvd, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_path(aes(colour = what, linetype = what), linewidth = 0.8) +
  geom_point(data = pvd[, .SD[.N], by = what], aes(colour = what), size = 2.5) +
  scale_colour_manual(values = c("Observed velocity" = "grey50", "Tide and vessel effect removed" = "black"), name = NULL) +
  scale_linetype_manual(values = c("Observed velocity" = 2, "Tide and vessel effect removed" = 1), name = NULL) +
  coord_equal() +
  labs(x = "Cumulative east displacement (km)", y = "Cumulative north displacement (km)",
       title = "Progressive vector diagram",
       subtitle = sprintf("Hourly mean school velocity integrated over the %d daytime hours with data", obs_h)) +
  theme(legend.position = "bottom")
ggsave(file.path(out, "08e_progressive_vector.png"), p8e, width = 7.5, height = 7.5, dpi = 200)

# 08f: residual direction rose ----------------------------------------------------------------------------------------
d[, rdir := bearing_of(rx, ry)]
rc <- circ_summary(d$rdir)
p8f <- ggplot(d, aes(rdir)) +
  geom_histogram(breaks = seq(0, 360, 15), fill = "steelblue", colour = "white") +
  annotate("segment", x = rc$mean_dir, xend = rc$mean_dir, y = 0, yend = Inf, colour = "firebrick", linewidth = 1) +
  coord_polar(start = 0) +
  scale_x_continuous(limits = c(0, 360), breaks = c(0, 90, 180, 270), labels = c("N", "E", "S", "W")) +
  labs(x = NULL, y = "Schools", title = "Direction of individual schools (tide and vessel effect removed)",
       subtitle = sprintf("n = %d, mean %.0f deg (red), Rbar = %.2f", nrow(d), rc$mean_dir, rc$Rbar))
ggsave(file.path(out, "08f_residual_direction_rose.png"), p8f, width = 6.5, height = 6.5, dpi = 200)

# summary figure ----------------------------------------------------------------------------------------------------
ms <- v_main$sum
p_sum <- (p8a + labs(title = "A  Net drift per 5-km cell", subtitle = NULL) + theme(legend.position = "bottom")) +
  (p8c + labs(title = "B  Net drift per day", subtitle = NULL)) +
  (p8e + labs(title = "C  Progressive vector diagram", subtitle = NULL)) +
  (p8b2 + labs(title = "D  Tidal ellipse and net drift", subtitle = NULL)) +
  plot_layout(ncol = 2) +
  plot_annotation(title = sprintf("Daytime herring schools: net drift %.1f km/day towards %.0f deg (95%% CI %.0f-%.0f deg) after removing the tide",
                                  ms$km_day, ms$dir, ms$dir_lo, ms$dir_hi),
                  subtitle = sprintf("%d schools, %d hour blocks, %s - %s", nrow(d), uniqueN(d$hb),
                                     format(min(d$day), "%d %b"), format(max(d$day), "%d %b %Y")))
ggsave(file.path(out, "08_migration_summary.png"), p_sum, width = 14, height = 13, dpi = 200)
message("Figures written to ", out)

# ================================ ANIMATION ============================================================
if (make_anim) {
  hours <- hourly$hb
  map_lim <- list(x = range(c(d$lon, d$lon_e)), y = range(c(d$lat, d$lat_e)))
  pvd_lim <- range(c(pvd$x, pvd$y))
  cmp_lim <- 0.6                  # m/s, fixed scale of the compass panel
  draw_hour <- function(k) {
    hk <- hours[k]
    cur <- d[hb == hk]; prev <- d[hb < hk & hb >= hk - 3]
    hm <- hourly[hb == hk]
    pm <- ggplot() +
      geom_path(data = ship_bg, aes(lon, lat, group = seg), colour = "grey90") +
      geom_segment(data = prev, aes(lon, lat, xend = lon_e, yend = lat_e), colour = "grey70", alpha = 0.4, linewidth = 0.3) +
      geom_segment(data = cur, aes(lon, lat, xend = lon_e, yend = lat_e, colour = speed_lm),
                   arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.6) +
      scale_colour_viridis_c(option = "plasma", limits = c(0, 2), oob = scales::squish, guide = "none") +
      coord_quickmap(xlim = map_lim$x, ylim = map_lim$y) +
      labs(x = NULL, y = NULL, title = format(hm$t, "%d %b %H:00 UTC (sonar clock)"),
           subtitle = sprintf("%d schools this hour (arrows = 2 h of movement; grey = previous 3 h)", hm$n))
    vec <- data.table(what = factor(c("Observed mean", "Tide", "Net drift"), c("Observed mean", "Tide", "Net drift")),
                      x = c(hm$vx, hm$tide_x, hm$rx), y = c(hm$vy, hm$tide_y, hm$ry))
    vec[, `:=`(x = pmax(pmin(x, cmp_lim), -cmp_lim), y = pmax(pmin(y, cmp_lim), -cmp_lim))]
    rings <- rbindlist(lapply(c(0.2, 0.4, 0.6), function(r) data.table(r = r, a = seq(0, 2 * pi, length = 90))))
    pc <- ggplot() +
      geom_path(data = rings, aes(r * sin(a), r * cos(a), group = r), colour = "grey88") +
      geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
      geom_segment(data = vec, aes(0, 0, xend = x, yend = y, colour = what), linewidth = 1.3,
                   arrow = arrow(length = unit(2.5, "mm"))) +
      annotate("segment", x = 0, y = 0, xend = fm$x[["a"]], yend = fm$y[["a"]], colour = "black", linetype = 3) +
      scale_colour_manual(values = c("Observed mean" = "grey55", "Tide" = "tomato", "Net drift" = "black"), name = NULL) +
      coord_equal(xlim = c(-cmp_lim, cmp_lim), ylim = c(-cmp_lim, cmp_lim)) +
      labs(x = "East (m/s)", y = "North (m/s)", title = "This hour",
           subtitle = "Rings 0.2, 0.4, 0.6 m/s; dotted = overall net drift") +
      theme(legend.position = "bottom")
    pp_ <- ggplot() +
      geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
      geom_path(data = pvd[t <= hm$t & what == "Observed velocity"], aes(x, y), colour = "grey60", linetype = 2) +
      geom_path(data = pvd[t <= hm$t & what != "Observed velocity"], aes(x, y), colour = "black", linewidth = 0.8) +
      geom_point(data = pvd[t <= hm$t & what != "Observed velocity"][.N], aes(x, y), colour = "firebrick", size = 3) +
      coord_equal(xlim = pvd_lim, ylim = pvd_lim) +
      labs(x = "East (km)", y = "North (km)", title = "Cumulative net drift",
           subtitle = "black = tide and vessel removed; grey dashed = observed")
    pm + pc + pp_ + plot_layout(widths = c(1.3, 0.9, 1))
  }
  fdir <- file.path(tempdir(), "frames_migration")
  message(sprintf("Rendering %d hourly frames ...", length(hours)))
  render_frames(draw_hour, length(hours), fdir, width_px = 2200, height_px = 900, res = 130)
  encode_mp4(fdir, file.path(out, "08_daytime_flow.mp4"), fps = 3)
}

# ================================ WORD SUMMARY ==========================================================
f0 <- function(x) formatC(round(x), format = "d", big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
compass <- function(b) c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")[round(b / 22.5) %% 16 + 1]
loo_dir <- range(unwrap(loo$dir, ms$dir)) %% 360
west_days <- daily[n >= 50, sum(dir > 180 & dir < 360)]; n_days50 <- daily[, sum(n >= 50)]
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

doc <- body_add_fpar(doc, fpar(ftext("Daytime herring school movement: is there a migration signal?", fp_text(bold = TRUE, font.size = 18))))
doc <- body_add_fpar(doc, fpar(ftext(sprintf("SCH302, %s sonar, trip %s, %s – %s. Generated %s with 08_daytime_migration.R",
  cfg$label, cfg$trip, format(min(d$day), "%d %B"), format(max(d$day), "%d %B %Y"), format(Sys.Date(), "%d %B %Y")),
  fp_text(font.size = 10, color = "grey40"))))

h1("1. Main result")
lead("Net drift.", sprintf(paste(
  "After removing the semidiurnal tide and the small vessel-related effect, daytime schools show a net drift of",
  "%s km/day (%s m/s) towards %.0f° (%s; 95%% bootstrap interval %.0f–%.0f°, speed %s–%s km/day)."),
  f1(ms$km_day), f2(ms$speed), ms$dir, compass(ms$dir), ms$dir_lo, ms$dir_hi, f1(ms$speed_lo * 86.4), f1(ms$speed_hi * 86.4)),
  if (ms$p_zero == 0) sprintf(" None of the %d bootstrap replicates pointed away from this direction.", B)
  else sprintf(" In %.1f%% of %d bootstrap replicates the drift pointed away from this direction.", 100 * ms$p_zero, B))
lead("Consistency.", sprintf(paste(
  "The direction is stable when any single day is left out (%.0f–%.0f°), when only slow ship speeds or only clearly moving",
  "schools are used, and on %d of %d days with at least 50 schools the daily drift pointed westward.",
  "Across %d grid cells of %d km, %.0f%% of schools were in cells with a westward net drift (cell directions R̄ = %s)."),
  loo_dir[1], loo_dir[2], west_days, n_days50, nrow(cells), grid_km_m, 100 * share_w, f2(cell_dirs$Rbar)))
lead("Tide.", sprintf(paste(
  "School velocities over the ground oscillate with the semidiurnal tide (amplitude %s m/s east–west and %s m/s",
  "north–south; F-test p < 10^-%d). The tidal amplitude is %s–%s times the net drift, so within a few hours the tide",
  "dominates the movement of schools over ground, but over a tidal cycle it largely cancels."), f2(m2_amp["x"]), f2(m2_amp["y"]),
  floor(-log10(max(max(tide_p), 1e-300))), f1(min(m2_amp) / ms$speed), f1(max(m2_amp) / ms$speed)))
lead("Individual schools.", sprintf(paste(
  "Individual schools still swim in all directions (R̄ = %s for de-tided directions): the drift is a weak population-level bias",
  "on top of a lot of local movement, not a coherent stream of schools."), f2(rc$Rbar)))
lead("Interpretation.", paste(
  "The data are consistent with a slow westward to south-westward displacement of the herring schools of the order of",
  sprintf("%s km/day during daylight in late June – early July.", f1(ms$km_day)),
  "This is suggestive of a migration signal but not proof of one: see Section 4 for the caveats, in particular the unknown",
  "residual (non-tidal) current, which the sonar cannot separate from swimming.",
  "In direction, a south-westward drift from this area points towards the Orkney–Shetland spawning grounds of North Sea",
  "autumn-spawning herring, which spawn there in August–September; a pre-spawning movement is therefore one plausible",
  "explanation, but these data alone cannot establish it."))
add_fig("08_migration_summary.png", sprintf(paste(
  "Summary. A: mean net drift per %d-km cell (tide and vessel effect removed; red = westward, blue = eastward), over the",
  "individual school vectors (grey). B: net drift per day with 95%% bootstrap intervals (black = all days). C: progressive",
  "vector diagram of hourly mean velocity (dashed = observed, solid = tide and vessel removed). D: hourly mean velocities",
  "(grey), fitted tidal ellipse (red) and net drift (arrow)."), grid_km_m))

h1("2. Data and approach")
para(sprintf(paste(
  "Data: the %s daytime schools (sun above the horizon) among the quality-controlled, smoothed school tracks of script 03",
  "(%s with speed significantly above zero), %d hour blocks on %d days. Each school contributes one velocity vector: the",
  "straight-line velocity of its smoothed track over ground (east and north components)."), f0(nrow(d)), f0(sum(d$moving)),
  uniqueN(d$hb), uniqueN(d$day)))
para(paste(
  "Model: for each component, velocity = net drift + semidiurnal tidal harmonic (period 12.42 h, free amplitude and phase) +",
  "a linear effect of ship velocity (absorbs the small apparent dependence of school velocity on ship motion found earlier) +",
  "noise. The intercept is the net drift: the mean ground velocity that remains when the tide is averaged out and the ship",
  "is not moving. Schools seen in the same hour are not independent, so confidence intervals come from a block bootstrap",
  sprintf("over 1-hour blocks (%d replicates).", B),
  "Robustness was checked with model variants, by leaving out one day at a time, per day and per grid cell, and with a",
  "progressive vector diagram that integrates the hourly mean velocities over time."))

h1("3. Results")
add_fig("08b_tidal_signal.png", "Left: hourly mean school velocity (points, sized by number of schools) and the fitted net drift + semidiurnal tide (red). Right: hourly means in velocity space with the fitted tidal ellipse (red) and the net drift (arrow).")
add_fig("08c_daily_vectors.png", "Net drift per day (tide and vessel effect removed) with 95% bootstrap intervals over hour blocks; black arrow = all days. Days with fewer than 20 schools are not shown.")
add_fig("08d_net_drift_estimates.png", "Net drift speed (left) and direction (right) with 95% bootstrap intervals for the model variants (black) and leaving out one day at a time (grey).")
add_fig("08a_flow_map.png", sprintf("Daytime school vectors (grey: distance covered in 2 h at the observed velocity) and mean net drift per %d-km cell with at least %d schools (red = westward, blue = eastward).", grid_km_m, min_cell_n))
add_fig("08e_progressive_vector.png", sprintf("Progressive vector diagram: hourly mean velocities integrated over the %d daytime hours with data. Dashed = observed velocity; solid = tide and vessel effect removed. Gaps between hours (night, no schools) are not included, so the end point is the displacement accumulated during observed hours only.", obs_h))
add_fig("08f_residual_direction_rose.png", "Directions of individual daytime schools after removing the tide and vessel effect; red line = mean direction.")
est_tab <- est[, .(Estimate = as.character(model), n = n, `km/day` = f1(km_day), `95% CI km/day` = paste0(f1(speed_lo * 86.4), "–", f1(speed_hi * 86.4)),
                   `Direction (°)` = round(dir), `95% CI (°)` = paste0(round(dir_lo), "–", round(dir_hi)))]
doc <- body_add_fpar(doc, fpar(ftext("Table 1. ", fp_cap_lab), ftext("Net drift estimates (block bootstrap, 1-h blocks).", fp_cap),
                               fp_p = fp_par(keep_with_next = TRUE, padding.top = 6)), style = "Table Caption")
ft <- flextable(est_tab); ft <- theme_booktabs(ft); ft <- fontsize(ft, size = 8.5, part = "all"); ft <- bold(ft, part = "header")
ft <- width(ft, width = c(2.2, 0.5, 0.6, 1.0, 0.8, 0.9))
doc <- body_add_flextable(doc, ft, align = "center")

h1("4. Caveats")
lead("Current versus swimming.", "The sonar measures movement over ground. The model removes the oscillating tide, but a residual (non-tidal) current would appear as net drift. Comparing the result with a tidal and residual current model for the area (e.g. CMEMS North West Shelf) is the key next step before calling the drift migration.")
lead("Vessel sampling.", "Schools are only seen within about 600 m of the vessel and in the forward sector, and the vessel goes where the fish are. The ship-velocity terms remove the direct effect of ship motion on apparent school velocity, but where and when schools were sampled is set by fishing.")
lead("Daytime only.", "With daylight from about 03 to 21 UTC, the semidiurnal tide is sampled unevenly and diurnal tides cannot be separated from the daily cycle of school behaviour. Eight days are not enough to separate M2 and S2.")
lead("Clock.", "Times are the sonar clock; the catch log differs by 1 h on this trip. This does not affect the drift estimate, but it shifts the tidal phase and the hour labels by 1 h if the sonar clock is not UTC.")
lead("Scale.", "A drift of a few km/day accumulates to tens of km over the period, comparable to the size of the fishing area. Longer series (more trips, other years) are needed to see whether this is a recurring seasonal movement.")

h1("5. Outputs")
para("Folder ", normalizePath(out, winslash = "/"), ": figures 08a–08f and 08_migration_summary.png; 08_daytime_flow.mp4 (one frame per",
     " daytime hour with schools: school vectors on the map, that hour's observed mean, tidal and net-drift velocity in the",
     " middle, cumulative net drift on the right);",
     " migration_estimates.csv and migration_daily.csv.")
print(doc, target = file.path(out, "SCH302_ST90_daytime_migration.docx"))
message("08: done -> ", normalizePath(out))
