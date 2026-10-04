################################################################################
# 00_setup_functions.R
#
# Shared settings + helper functions for combining LSSS Profos school exports
# (omnidirectional sonar, SCH302) with commercial midwater-trawl catch data.
#
# Source this file at the top of every analysis script:
#   source("00_setup_functions.R")   # run from the RStudio project root
#
# Profos export types (one set per LSSS "run", e.g. 2025-1023-0303-*):
#   *-pp.txt  : one row per school per ping   (school track points)
#   *-ps.txt  : one row per school            (LSSS school summary)
#   *-ppe.txt : one row per ping, incl. pings WITHOUT schools (Id = N/A)
#               -> this is the sonar EFFORT / ship track file.
#
# Facts established from the October 2025 sample (23 Oct 2025):
#   * Ship.speed in Profos is in m/s (matches speed derived from GPS track).
#   * Catch-file shoot/haul times are UTC (ship track passes the logged
#     shoot/haul positions at offset 0 h, within the 1-minute position rounding).
#   * A school Id can have 2 detections in the same ping (split school):
#     these are merged (area-weighted) instead of jittering timestamps.
#   * Raw ping-to-ping school centroid steps imply a median "speed" of
#     ~3 m/s -> mostly positioning noise; smoothing is essential.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(geosphere)
  library(scales)
})

# ------------------------------------------------------------------------------
# Settings -- edit these
# ------------------------------------------------------------------------------
cfg <- list(
  profos_dir   = "C:/Drive/Sonar/October2025/profos",  # all *-pp/ps/ppe.txt here
  catch_file   = "C:/Drive/Sonar/SCH302_2025K_2026E.csv",
  trip         = "2025K",
  out_dir      = "output",   # relative to the project root

  # sonar geometry (for depth-coverage diagnostics)
  transducer_depth = 6,      # m
  vert_beam_half   = 4,      # deg, half the vertical beamwidth (approx.)

  # school quality control
  min_pings     = 8,         # pings per school track
  min_duration  = 20,        # s
  gap_s         = 30,        # ship-track gap (s) treated as "sonar off"

  # smoothing
  hampel_k      = 3,         # half-window (pings) for outlier filter
  hampel_nsig   = 3,
  hampel_floor  = 10,        # m, never flag deviations smaller than this
  tau_s         = 30,        # spline: ~1 degree of freedom per tau_s seconds
  sigma_acc     = 0.03,      # Kalman: process noise, m s^-2
  interval_s    = 10,        # Figure_6-style interval length (s)

  # sonar <-> haul matching
  pre_shoot_h   = 2,         # hours of searching before shoot to include
  buffer_km     = 5,         # max distance of school from the tow track
  default_vert_opening = 30, # m, used when vert_opening is NA

  # map display (plots only -- does not remove data from any analysis)
  map_lat_min   = 60.12      # southern limit of sonar maps; NULL = no limit
)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

species_cols <- c(her = "#1f77b4", mac = "#2ca02c", whb = "#9467bd",
                  hom = "#ff7f0e", boc = "#8c564b", had = "#e377c2",
                  pok = "#7f7f7f")

# Map coordinates for the sonar zoom maps, cropped at cfg$map_lat_min.
# Pass the lon/lat of everything drawn so the view fits the remaining area.
map_coord <- function(lon, lat, lat_min = cfg$map_lat_min) {
  if (is.null(lat_min)) return(coord_quickmap())
  keep <- !is.na(lat) & lat >= lat_min
  coord_quickmap(xlim = range(lon[keep], na.rm = TRUE),
                 ylim = c(lat_min, max(lat[keep], na.rm = TRUE)))
}

theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey92")))

# ------------------------------------------------------------------------------
# Readers
# ------------------------------------------------------------------------------

# Read and stack all Profos exports of one type in a directory.
# Adds: run (file prefix), sid (unique school key across runs), t (POSIXct UTC)
read_profos <- function(dir, type = c("pp", "ps", "ppe")) {
  type  <- match.arg(type)
  files <- list.files(dir, pattern = paste0("-", type, "\\.txt$"),
                      full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop("No *-", type, ".txt files in ", dir)

  dt <- rbindlist(lapply(files, function(f) {
    x <- fread(f, na.strings = "N/A", fill = TRUE)
    x[, run := sub(paste0("-", type, "\\.txt$"), "", basename(f))]
    x
  }), fill = TRUE)

  if (type == "ps") {
    dt[, t_start := as.POSIXct(paste(StartDate, StartTime), tz = "UTC",
                               format = "%Y-%m-%d %H:%M:%OS")]
    dt[, t_stop  := as.POSIXct(paste(StopDate, StopTime), tz = "UTC",
                               format = "%Y-%m-%d %H:%M:%OS")]
  } else {
    dt[, t := as.POSIXct(paste(Date, Time), tz = "UTC",
                         format = "%Y-%m-%d %H:%M:%OS")]
  }
  dt[!is.na(Id), sid := paste(run, Id, sep = "_")]
  dt[]
}

# Haul-level table from the catch file: one row per haul, with catch per
# species in wide format and the dominant species.
read_catch <- function(file, trips = NULL) {
  d <- fread(file)
  if (!is.null(trips)) d <- d[trip %in% trips]

  hauls <- unique(d[, .(trip, haul, shoot_time, haul_time,
                        shootlon, shootlat, haullon, haullat,
                        headline_depth, vert_opening, water_depth,
                        duration_fishing, speed_fishing, catch)],
                  by = c("trip", "haul"))

  sp <- unique(d[, .(trip, haul, species, catch_sp)])
  sp[is.na(catch_sp), catch_sp := 0]
  dom <- sp[order(-catch_sp),
            .(dom_species = species[1],
              dom_frac    = catch_sp[1] / sum(catch_sp)),
            by = .(trip, haul)]
  wide <- dcast(sp, trip + haul ~ species, value.var = "catch_sp",
                fill = 0, fun.aggregate = sum)

  hauls <- merge(hauls, dom, by = c("trip", "haul"))
  hauls <- merge(hauls, wide, by = c("trip", "haul"))
  hauls[, `:=`(vert_opening = as.numeric(vert_opening),
               headline_depth = as.numeric(headline_depth))]
  hauls[is.na(vert_opening), vert_opening := cfg$default_vert_opening]
  hauls[, haul_id := paste(trip, haul, sep = "-")]
  setorder(hauls, shoot_time)
  hauls[]
}

# ------------------------------------------------------------------------------
# Ship track / sonar effort (from ppe)
# ------------------------------------------------------------------------------

# One row per ping. seg = continuous recording segment (new seg after gap_s).
# dist_m = distance sailed since previous ping (0 across gaps).
ship_track <- function(ppe, gap_s = cfg$gap_s) {
  sh <- unique(ppe[, .(t, lon = Ship.lon, lat = Ship.lat,
                       speed = Ship.speed, heading = Ship.heading,
                       tilt = Trans.tilt)], by = "t")
  setorder(sh, t)
  sh[, dt := c(NA_real_, diff(as.numeric(t)))]
  sh[, seg := cumsum(is.na(dt) | dt > gap_s)]
  sh[, dist_m := c(0, distHaversine(cbind(lon[-.N], lat[-.N]),
                                    cbind(lon[-1], lat[-1])))]
  sh[is.na(dt) | dt > gap_s, `:=`(dist_m = 0, dt = NA_real_)]
  sh[]
}

# Label each ping (or any time vector) with the haul it belongs to and phase:
# "search" = up to pre_h hours before shoot, "tow" = shoot..haul_time.
label_phase <- function(t, hauls, pre_h = cfg$pre_shoot_h) {
  out <- data.table(haul_id = NA_character_, phase = "other")[rep(1, length(t))]
  for (i in seq_len(nrow(hauls))) {
    h <- hauls[i]
    tow <- t >= h$shoot_time & t <= h$haul_time
    pre <- t >= h$shoot_time - pre_h * 3600 & t < h$shoot_time &
      out$phase == "other"
    out[tow, `:=`(haul_id = h$haul_id, phase = "tow")]
    out[pre, `:=`(haul_id = h$haul_id, phase = "search")]
  }
  out
}

# ------------------------------------------------------------------------------
# Geometry
# ------------------------------------------------------------------------------

# Local metric coordinates (m) around (lon0, lat0). Equirectangular is
# indistinguishable from the aeqd projection used in Figure_6 at the
# few-km scale of a school track, and needs no sf transform per school.
to_xy <- function(lon, lat, lon0 = mean(lon), lat0 = mean(lat)) {
  R <- 6371008.8
  cbind(x = (lon - lon0) * pi / 180 * R * cos(lat0 * pi / 180),
        y = (lat - lat0) * pi / 180 * R)
}
from_xy <- function(x, y, lon0, lat0) {
  R <- 6371008.8
  cbind(lon = lon0 + x / (R * cos(lat0 * pi / 180)) * 180 / pi,
        lat = lat0 + y / R * 180 / pi)
}

# Compass bearing (0 = N, clockwise) of a velocity vector (vx east, vy north).
# Equivalent to atan2 + convert.heading.angle() in Figure_6_R_code.R.
bearing_of <- function(vx, vy) (atan2(vx, vy) * 180 / pi) %% 360

# Circular statistics (degrees). Rbar = mean resultant length (0 = random
# directions, 1 = all identical). Rayleigh p-value via Zar's approximation.
circ_summary <- function(b, w = NULL) {
  b <- b[!is.na(b)]; n <- length(b)
  if (n < 2) return(data.table(n = n, mean_dir = NA_real_, Rbar = NA_real_,
                               rayleigh_p = NA_real_))
  if (is.null(w)) w <- rep(1, n)
  r  <- b * pi / 180
  C  <- sum(w * cos(r)) / sum(w); S <- sum(w * sin(r)) / sum(w)
  Rb <- sqrt(C^2 + S^2)
  Z  <- n * Rb^2
  p  <- exp(sqrt(1 + 4 * n + 4 * (n^2 - (n * Rb)^2)) - (1 + 2 * n))
  data.table(n = n, mean_dir = (atan2(S, C) * 180 / pi) %% 360,
             Rbar = Rb, rayleigh_p = min(1, max(0, p)))
}

# Solar elevation (deg) and day/twilight/night label.
sun_phase <- function(t, lon, lat) {
  el <- suncalc::getSunlightPosition(
    data = data.frame(date = t, lat = lat, lon = lon))$altitude * 180 / pi
  ph <- fifelse(el > 0, "day", fifelse(el < -6, "night", "twilight"))
  list(sun_elev = el, sun = factor(ph, c("day", "twilight", "night")))
}

# ------------------------------------------------------------------------------
# School data cleaning
# ------------------------------------------------------------------------------

# Merge multiple detections of the same school in the same ping
# (area-weighted centroid; Sv averaged in the linear domain).
# Adds range (m) and bearing relative to ship heading.
merge_pp <- function(pp) {
  pp <- pp[!is.na(Id) & !is.na(Center.lon)]
  pp[, sv_lin := 10^(Sv.mean / 10)]
  m <- pp[, .(lon   = weighted.mean(Center.lon, Area),
              lat   = weighted.mean(Center.lat, Area),
              depth = weighted.mean(Center.dep, Area),
              Area  = sum(Area),
              sv_lin = weighted.mean(sv_lin, Area),
              AlongBeamSize = sum(AlongBeamSize),
              AlongRingSize = max(AlongRingSize),
              ship_lon = Ship.lon[1], ship_lat = Ship.lat[1],
              ship_speed = Ship.speed[1], ship_heading = Ship.heading[1],
              tilt = Trans.tilt[1], n_parts = .N),
          by = .(run, Id, sid, t)]
  m[, Sv := 10 * log10(sv_lin)]
  m[, backscatter := sv_lin * Area]   # school "area backscatter" (relative)
  xy <- to_xy(m$lon, m$lat, m$ship_lon, m$ship_lat)
  m[, range_m := sqrt(xy[, 1]^2 + xy[, 2]^2)]
  m[, rel_bearing := (bearing_of(xy[, 1], xy[, 2]) - ship_heading) %% 360]
  setorder(m, sid, t)
  m[]
}

# Per-school summary from merged pp (used for matching with hauls/maps).
school_summary <- function(m) {
  m[, .(run = run[1], Id = Id[1],
        t_start = min(t), t_end = max(t),
        duration_s = as.numeric(max(t)) - as.numeric(min(t)),
        n_pings = .N,
        lon = mean(lon), lat = mean(lat),
        depth = weighted.mean(depth, backscatter),
        Area = mean(Area),
        Sv = 10 * log10(mean(sv_lin)),
        backscatter = mean(backscatter),
        range_m = mean(range_m),
        ship_speed = mean(ship_speed)),
    by = sid]
}

# ------------------------------------------------------------------------------
# Track smoothing
# ------------------------------------------------------------------------------

# Hampel-type outlier flag on a 2-D track: a point is an outlier if its
# distance to the running-median position exceeds nsig robust SDs.
hampel_xy <- function(x, y, k = cfg$hampel_k, nsig = cfg$hampel_nsig,
                      floor_m = cfg$hampel_floor) {
  n <- length(x)
  if (n < 2 * k + 1) return(rep(TRUE, n))
  mx <- runmed(x, 2 * k + 1, endrule = "median")
  my <- runmed(y, 2 * k + 1, endrule = "median")
  d  <- sqrt((x - mx)^2 + (y - my)^2)
  d <= max(floor_m, nsig * 1.4826 * median(d))
}

# Robust estimate of centroid measurement noise (m, per axis) from second
# differences: for white noise var(2nd diff) = 6 sigma^2, and a constant
# swimming velocity cancels out.
noise_sd <- function(x) {
  d2 <- diff(x, differences = 2)
  if (length(d2) < 3) return(NA_real_)
  1.4826 * median(abs(d2 - median(d2))) / sqrt(6)
}

# Smoothing spline with a PHYSICAL smoothing scale: ~1 df per tau_s seconds.
# (Figure_6 used a fixed lambda, whose effect depends on track length and
# time scaling, so the same lambda smooths short and long tracks differently.)
# Returns smoothed position and velocity (analytic derivative).
smooth_spline_track <- function(tt, z, tau_s = cfg$tau_s) {
  n  <- length(z)
  df <- min(n - 1, max(2, 1 + diff(range(tt)) / tau_s))
  f  <- smooth.spline(tt, z, df = df)
  list(pos = predict(f, tt)$y, vel = predict(f, tt, deriv = 1)$y)
}

# Constant-velocity Kalman filter + Rauch-Tung-Striebel smoother (one axis).
# Handles irregular ping intervals; gives velocity and its uncertainty.
#   sigma_meas : centroid positioning noise (m)
#   sigma_acc  : how much the school may change velocity (m s^-2)
kalman_cv_1d <- function(tt, z, sigma_meas, sigma_acc = cfg$sigma_acc,
                         v0_sd = 1) {
  n  <- length(z)
  xp <- xf <- matrix(0, n, 2)
  Pp <- Pf <- array(0, c(2, 2, n))
  Fl <- vector("list", n)
  R  <- sigma_meas^2
  x  <- c(z[1], 0)
  P  <- diag(c(R, v0_sd^2))
  for (i in seq_len(n)) {
    if (i > 1) {
      dt <- tt[i] - tt[i - 1]
      Fm <- matrix(c(1, 0, dt, 1), 2)
      Q  <- sigma_acc^2 * matrix(c(dt^3 / 3, dt^2 / 2, dt^2 / 2, dt), 2)
      x  <- as.vector(Fm %*% x)
      P  <- Fm %*% P %*% t(Fm) + Q
      Fl[[i]] <- Fm
    }
    xp[i, ] <- x; Pp[, , i] <- P
    K <- P[, 1] / (P[1, 1] + R)
    x <- x + K * (z[i] - x[1])
    P <- P - K %*% t(P[1, ])
    xf[i, ] <- x; Pf[, , i] <- P
  }
  xs <- xf; Ps <- Pf
  for (i in (n - 1):1) {
    C <- Pf[, , i] %*% t(Fl[[i + 1]]) %*% solve(Pp[, , i + 1])
    xs[i, ] <- xf[i, ] + as.vector(C %*% (xs[i + 1, ] - xp[i + 1, ]))
    Ps[, , i] <- Pf[, , i] + C %*% (Ps[, , i + 1] - Pp[, , i + 1]) %*% t(C)
  }
  list(pos = xs[, 1], vel = xs[, 2], vel_sd = sqrt(pmax(Ps[2, 2, ], 0)))
}

# Full per-school pipeline: QC -> outlier removal -> spline + Kalman smoothing.
# Returns the point-level track (raw, cleaned, smoothed) for one school.
smooth_school <- function(s, sigma_meas) {
  s   <- copy(s)
  lon0 <- mean(s$lon); lat0 <- mean(s$lat)
  xy  <- to_xy(s$lon, s$lat, lon0, lat0)
  s[, `:=`(x = xy[, 1], y = xy[, 2], tt = as.numeric(t) - as.numeric(t[1]))]
  s[, keep := hampel_xy(x, y)]
  k <- s[keep == TRUE]
  if (nrow(k) < 4) return(NULL)

  sx <- smooth_spline_track(k$tt, k$x); sy <- smooth_spline_track(k$tt, k$y)
  kx <- kalman_cv_1d(k$tt, k$x, sigma_meas)
  ky <- kalman_cv_1d(k$tt, k$y, sigma_meas)

  k[, `:=`(x_spl = sx$pos, y_spl = sy$pos, vx_spl = sx$vel, vy_spl = sy$vel,
           x_kf = kx$pos, y_kf = ky$pos, vx_kf = kx$vel, vy_kf = ky$vel,
           v_sd_kf = sqrt((kx$vel_sd^2 + ky$vel_sd^2) / 2))]
  ll <- from_xy(k$x_kf, k$y_kf, lon0, lat0)
  k[, `:=`(lon_kf = ll[, 1], lat_kf = ll[, 2])]
  ll <- from_xy(k$x_spl, k$y_spl, lon0, lat0)
  k[, `:=`(lon_spl = ll[, 1], lat_spl = ll[, 2])]
  # keep outliers in the output (flagged) for plotting
  rbind(k, s[keep == FALSE], fill = TRUE)[order(t)]
}

# Track-level movement metrics for one smoothed school track.
#  speed_lm / bearing_lm : straight-line regression of cleaned positions on time
#  speed_se              : SE inflated for lag-1 autocorrelation of residuals
#  moving                : speed significantly > 0 (speed > 2 SE)
#  straightness          : net displacement / smoothed path length (0..1)
track_metrics <- function(k) {
  k <- k[keep == TRUE]
  fx <- lm(x ~ tt, k); fy <- lm(y ~ tt, k)
  vx <- coef(fx)[2]; vy <- coef(fy)[2]
  infl <- function(f) {
    r <- resid(f); rho <- if (length(r) > 3) acf(r, 1, plot = FALSE)$acf[2] else 0
    rho <- max(0, min(rho, 0.95)); sqrt((1 + rho) / (1 - rho))
  }
  sex <- summary(fx)$coef[2, 2] * infl(fx); sey <- summary(fy)$coef[2, 2] * infl(fy)
  sp  <- sqrt(vx^2 + vy^2)
  se  <- sqrt((vx * sex)^2 + (vy * sey)^2) / max(sp, 1e-9)
  path <- sum(sqrt(diff(k$x_kf)^2 + diff(k$y_kf)^2))
  net  <- sqrt((last(k$x_kf) - k$x_kf[1])^2 + (last(k$y_kf) - k$y_kf[1])^2)
  ship_vx <- mean(k$ship_speed * sin(k$ship_heading * pi / 180))
  ship_vy <- mean(k$ship_speed * cos(k$ship_heading * pi / 180))
  data.table(vx = vx, vy = vy, speed_lm = sp, speed_se = se,
             bearing_lm = bearing_of(vx, vy),
             moving = sp > 2 * se,
             speed_kf = mean(sqrt(k$vx_kf^2 + k$vy_kf^2)),
             path_m = path, net_m = net,
             straightness = ifelse(path > 0, net / path, NA_real_),
             ship_vx = ship_vx, ship_vy = ship_vy,
             ship_speed = sqrt(ship_vx^2 + ship_vy^2),
             ship_bearing = bearing_of(ship_vx, ship_vy),
             rmse_raw_vs_kf = sqrt(mean((k$x - k$x_kf)^2 + (k$y - k$y_kf)^2)),
             n_outliers = NA_integer_)
}

# Figure_6-style interval velocities (regression within interval_s windows),
# but applied to the SMOOTHED positions of any school.
interval_velocity <- function(k, interval_s = cfg$interval_s) {
  k <- k[keep == TRUE]
  k[, int := floor(tt / interval_s)]
  k[, if (.N > 1) {
    fx <- coef(lm(x_kf ~ tt))[2]; fy <- coef(lm(y_kf ~ tt))[2]
    .(t_mid = mean(t), vx = fx, vy = fy, speed = sqrt(fx^2 + fy^2),
      bearing = bearing_of(fx, fy), lon = mean(lon_kf), lat = mean(lat_kf),
      depth = mean(depth))
  }, by = int]
}
