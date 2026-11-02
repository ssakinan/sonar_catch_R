################################################################################
# 09_cmems_comparison.R
#
# Is the daytime net drift of script 08 swimming or current?
# Compares school ground velocities with CMEMS North-West Shelf model currents
# (NWSHELF_ANALYSISFORECAST_PHY_004_013, 1.5 km, hourly, tides included).
# (run 00b, 03 and 09a_cmems_download.py first; intended for the ST90 profile)
#   sonar_profile <- "ST90"; source("09_cmems_comparison.R")
#
# For every daytime school the model current (u, v) is interpolated to the
# school's position, track mid time and depth (bilinear in space, linear in
# depth and time). Then
#   swimming velocity = school ground velocity - model current
# and the net drift model of script 08 (intercept + M2 harmonic + ship velocity)
# is fitted to the observed velocity, to the model current at the schools and to
# the swimming velocity, with the same 1-hour block bootstrap.
#
# The sonar clock is checked against the model tide: the offset between sonar
# clock and UTC is the lag that best aligns the school velocities with the
# model current.
#
# Outputs (cfg$out_dir/migration): 09a-09f figures, 09_cmems_summary.png,
# cmems_estimates.csv, and a Word summary.
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root
library(officer); library(flextable); library(ncdf4)

B          <- 1000                 # bootstrap replicates
w          <- 2 * pi / (12.4206 * 3600)
grid_km_m  <- 5                    # cell size for the spatial map
min_cell_n <- 10                   # schools per cell to draw an arrow
lags_h     <- seq(-3, 3, 0.25)     # candidate sonar clock - UTC offsets
map_depth  <- 25                   # depth (m) of the model residual field on the map
set.seed(1)

base <- sub("/report_figs$", "", cfg$out_dir)
out  <- file.path(base, "migration")
nc_dir <- file.path(base, "cmems")
if (!dir.exists(nc_dir)) nc_dir <- "output_ST90/cmems"   # SU90_2026: same box and period as ST90
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(file.path(nc_dir, "nws_cur_hourly.nc")))
  stop("CMEMS files missing: run  python 09a_cmems_download.py  from the project root first")

# ---- CMEMS fields --------------------------------------------------------------------------------
read_nc <- function(file) {
  nc <- nc_open(file); on.exit(nc_close(nc))
  list(lon = ncvar_get(nc, "longitude"), lat = ncvar_get(nc, "latitude"), dep = ncvar_get(nc, "depth"),
       time = as.numeric(ncvar_get(nc, "time")), u = ncvar_get(nc, "uo"), v = ncvar_get(nc, "vo"))   # [lon, lat, depth, time]
}
hr <- read_nc(file.path(nc_dir, "nws_cur_hourly.nc"))
dy <- read_nc(file.path(nc_dir, "nws_cur_daily.nc"))
dy$time <- dy$time + 43200        # daily means are stamped at 00:00, centre them at noon

# quadrilinear interpolation, ignoring land/sea-bed (NA) corners
interp4 <- function(f, var, lon, lat, dep, t) {
  pos <- function(g, x) { i <- findInterval(x, g, all.inside = TRUE)
    list(i = i, f = pmin(pmax((x - g[i]) / (g[i + 1] - g[i]), 0), 1)) }
  px <- pos(f$lon, lon); py <- pos(f$lat, lat); pz <- pos(f$dep, pmin(pmax(dep, min(f$dep)), max(f$dep)))
  pt <- pos(f$time, t)
  num <- den <- numeric(length(lon))
  for (a in 0:1) for (b in 0:1) for (c in 0:1) for (e in 0:1) {
    wt <- (if (a) px$f else 1 - px$f) * (if (b) py$f else 1 - py$f) *
          (if (c) pz$f else 1 - pz$f) * (if (e) pt$f else 1 - pt$f)
    val <- f[[var]][cbind(px$i + a, py$i + b, pz$i + c, pt$i + e)]
    ok <- !is.na(val); num[ok] <- num[ok] + wt[ok] * val[ok]; den[ok] <- den[ok] + wt[ok]
  }
  ifelse(den > 0, num / den, NA_real_)
}
# depth-mean current over the model levels (0-75 m) at the school position and time
interp_dmean <- function(f, var, lon, lat, t) {
  z <- sapply(f$dep, function(dd) interp4(f, var, lon, lat, rep(dd, length(lon)), t))
  rowMeans(z, na.rm = TRUE)
}

# ---- data: daytime tracked schools (as in 08) ----------------------------------------------------
m <- fread(file.path(base, "school_movement.csv"))
m[, t := as.POSIXct(t_start, tz = "UTC") + duration_s / 2]      # track mid time, sonar clock
d <- m[sun == "day"]
d[, `:=`(tt = as.numeric(t), hb = floor(as.numeric(t) / 3600), day = as.Date(t))]
d[, `:=`(cw = cos(w * tt), sw = sin(w * tt))]
model_at <- function(lag_h, var, depth = d$depth, f = hr) interp4(f, var, d$lon, d$lat, depth, d$tt - lag_h * 3600)

# ---- clock check: which sonar clock - UTC offset aligns schools with the model tide? -------------
message("Clock check ...")
lag_scan <- rbindlist(lapply(lags_h, function(L) {
  x <- d[, .(hb, vx, vy, u = model_at(L, "u"), v = model_at(L, "v"))][, lapply(.SD, mean), by = hb]
  # vector correlation of hourly means (components centred separately)
  ob <- c(x$vx - mean(x$vx), x$vy - mean(x$vy)); mo <- c(x$u - mean(x$u), x$v - mean(x$v))
  data.table(lag_h = L, r = cor(ob, mo), slope = sum(ob * mo) / sum(mo^2))
}))
best_lag <- lag_scan[which.max(r), lag_h]
lag_use  <- round(best_lag)                   # clock offsets are whole hours
if (nzchar(Sys.getenv("CMEMS_LAG_H"))) lag_use <- as.numeric(Sys.getenv("CMEMS_LAG_H"))   # fixed offset, e.g. from another sonar
message(sprintf("Best sonar clock - UTC offset: %+.2f h (r = %.2f); using %+d h", best_lag, max(lag_scan$r), lag_use))

# ---- model current at each school ------------------------------------------------------------------
d[, `:=`(u = model_at(lag_use, "u"), v = model_at(lag_use, "v"))]                        # school depth
d[, `:=`(u0 = model_at(lag_use, "u", 0), v0 = model_at(lag_use, "v", 0))]               # surface
d[, `:=`(ud = interp_dmean(hr, "u", lon, lat, tt - lag_use * 3600),
         vd = interp_dmean(hr, "v", lon, lat, tt - lag_use * 3600))]                     # depth mean 0-75 m
d[, `:=`(ur = model_at(lag_use, "u", f = dy), vr = model_at(lag_use, "v", f = dy))]     # daily mean (residual)
d[, `:=`(ux_lag0 = model_at(0, "u"), vy_lag0 = model_at(0, "v"))]
n_na <- d[, sum(is.na(u) | is.na(ud) | is.na(ur))]
if (n_na) message(n_na, " schools outside the model sea area dropped")
d <- d[!is.na(u) & !is.na(ud) & !is.na(ur)]
message(sprintf("Daytime schools with model current: %d, %d hour blocks", nrow(d), uniqueN(d$hb)))

# ---- fits (same model and bootstrap as 08) ----------------------------------------------------------
terms_full <- c("cw", "sw", "ship_vx", "ship_vy")
fit_flow <- function(x, terms, cx = "vx", cy = "vy") {
  X <- if (length(terms)) cbind(1, as.matrix(x[, ..terms])) else matrix(1, nrow(x), 1)
  ax <- lm.fit(X, x[[cx]])$coefficients; ay <- lm.fit(X, x[[cy]])$coefficients
  names(ax) <- names(ay) <- c("a", terms)
  list(x = ax, y = ay)
}
boot_flow <- function(x, terms, B, cx = "vx", cy = "vy") {
  blocks <- split(seq_len(nrow(x)), x$hb)
  t(replicate(B, {
    idx <- unlist(blocks[sample(length(blocks), replace = TRUE)], use.names = FALSE)
    f <- fit_flow(x[idx], terms, cx, cy); c(ax = f$x[["a"]], ay = f$y[["a"]])
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
variant <- function(x, label, cx = "vx", cy = "vy", terms = terms_full) {
  f <- fit_flow(x, terms, cx, cy); bt <- boot_flow(x, terms, B, cx, cy)
  list(fit = f, boot = bt, sum = summ_vec(f$x[["a"]], f$y[["a"]], bt, label, nrow(x)))
}
d[, `:=`(sx = vx - u, sy = vy - v, sx0 = vx - u0, sy0 = vy - v0, sxd = vx - ud, syd = vy - vd,
         sx_lag0 = vx - ux_lag0, sy_lag0 = vy - vy_lag0)]

message("Fitting and bootstrapping ...")
v_obs   <- variant(d, "Schools: ground velocity (as 08)")
v_mod   <- variant(d, "CMEMS current at schools (hourly, same fit)", "u", "v")
v_res   <- variant(d, "CMEMS daily-mean current at schools", "ur", "vr", character(0))
v_swim  <- variant(d, "Swimming = ground - CMEMS (main)", "sx", "sy")
v_swimd <- variant(d, "Swimming, depth-mean current 0-75 m", "sxd", "syd")
v_swim0 <- variant(d, "Swimming, surface current", "sx0", "sy0")
v_swim_lag0 <- variant(d, "Swimming, sonar clock taken as UTC", "sx_lag0", "sy_lag0")
v_swim_slow <- variant(d[ship_speed_trk < 2.5], "Swimming, ship speed < 2.5 m/s", "sx", "sy")
v_swim_mov  <- variant(d[moving == TRUE], "Swimming, moving schools only", "sx", "sy")
loo <- rbindlist(lapply(sort(unique(d$day)), function(dd) {
  x <- d[day != dd]; f <- fit_flow(x, terms_full, "sx", "sy"); bt <- boot_flow(x, terms_full, 200, "sx", "sy")
  summ_vec(f$x[["a"]], f$y[["a"]], bt, paste("Swimming, without", format(dd, "%d %b")), nrow(x))
}))
est <- rbind(v_obs$sum, v_mod$sum, v_res$sum, v_swim$sum, v_swimd$sum, v_swim0$sum, v_swim_lag0$sum,
             v_swim_slow$sum, v_swim_mov$sum, loo)
fwrite(est, file.path(out, "cmems_estimates.csv"))
print(est[, .(model, n, speed = round(speed, 3), km_day = round(km_day, 1), dir = round(dir),
              dir_lo = round(dir_lo), dir_hi = round(dir_hi), p_zero)])

# tidal ellipses: observed vs model at the same schools
ellipse <- function(f, label) { ph <- seq(0, 2 * pi, length = 100)
  data.table(what = label, x = f$x[["a"]] + f$x[["cw"]] * cos(ph) + f$x[["sw"]] * sin(ph),
             y = f$y[["a"]] + f$y[["cw"]] * cos(ph) + f$y[["sw"]] * sin(ph)) }
amp <- function(f) c(x = sqrt(f$x[["cw"]]^2 + f$x[["sw"]]^2), y = sqrt(f$y[["cw"]]^2 + f$y[["sw"]]^2))
pha <- function(f) c(x = atan2(f$x[["sw"]], f$x[["cw"]]), y = atan2(f$y[["sw"]], f$y[["cw"]])) * 12.4206 / (2 * pi)  # h
amp_obs <- amp(v_obs$fit); amp_mod <- amp(v_mod$fit); amp_swim <- amp(v_swim$fit)
pha_dif <- ((pha(v_obs$fit) - pha(v_mod$fit) + 6.21) %% 12.42) - 6.21       # h, obs later = positive

# hourly means for time series, regression and progressive vectors
hourly <- d[, .(t = mean(t), n = .N, vx = mean(vx), vy = mean(vy), u = mean(u), v = mean(v),
                sx = mean(sx), sy = mean(sy)), by = hb][order(hb)]
adv <- c(x = coef(lm(vx ~ u, hourly))[["u"]], y = coef(lm(vy ~ v, hourly))[["v"]])  # slope obs on model
r2_tide <- 1 - c(x = sum(hourly[, (vx - u - mean(vx - u))^2]) / sum(hourly[, (vx - mean(vx))^2]),
                 y = sum(hourly[, (vy - v - mean(vy - v))^2]) / sum(hourly[, (vy - mean(vy))^2]))
hourly[, `:=`(ox = cumsum(vx) * 3.6, oy = cumsum(vy) * 3.6, mx = cumsum(u) * 3.6, my = cumsum(v) * 3.6,
              wx = cumsum(sx) * 3.6, wy = cumsum(sy) * 3.6)]                  # km (1 h per block)
obs_h <- nrow(hourly)

# per-day swimming vectors (tide residual + vessel effect removed from swimming velocity)
fs <- v_swim$fit
pred_nonnet <- function(f, x) cbind(as.matrix(x[, ..terms_full]) %*% f$x[terms_full],
                                    as.matrix(x[, ..terms_full]) %*% f$y[terms_full])
pn <- pred_nonnet(fs, d); d[, `:=`(qx = sx - pn[, 1], qy = sy - pn[, 2])]
pn <- pred_nonnet(v_obs$fit, d); d[, `:=`(rx = vx - pn[, 1], ry = vy - pn[, 2])]
daily <- d[, .(n = .N, rx = mean(rx), ry = mean(ry), qx = mean(qx), qy = mean(qy), ur = mean(ur), vr = mean(vr)),
           by = day][order(day)]
daily[, `:=`(dir_obs = bearing_of(rx, ry), dir_swim = bearing_of(qx, qy), dir_cur = bearing_of(ur, vr))]
near_dir <- function(b, ref) abs(((b - ref + 180) %% 360) - 180) < 90           # within 90 deg of the overall direction
swim_days_w <- daily[n >= 50, sum(near_dir(dir_swim, v_swim$sum$dir))]; n_days50 <- daily[, sum(n >= 50)]

# ---- spatial: model residual field and per-cell drift ------------------------------------------------------------
lon0 <- mean(d$lon); lat0 <- mean(d$lat)
xy <- to_xy(d$lon, d$lat, lon0, lat0)
d[, `:=`(gx = floor(xy[, 1] / (grid_km_m * 1000)), gy = floor(xy[, 2] / (grid_km_m * 1000)))]
cells <- d[, .(n = .N, rx = mean(rx), ry = mean(ry), qx = mean(qx), qy = mean(qy), ur = mean(ur), vr = mean(vr)),
           by = .(gx, gy)][n >= min_cell_n]
cc <- from_xy((cells$gx + 0.5) * grid_km_m * 1000, (cells$gy + 0.5) * grid_km_m * 1000, lon0, lat0)
cells[, `:=`(lon = cc[, 1], lat = cc[, 2])]
cell_dirs_swim <- circ_summary(bearing_of(cells$qx, cells$qy), w = cells$n)
share_w_swim <- cells[, sum(n[near_dir(bearing_of(qx, qy), v_swim$sum$dir)]) / sum(n)]

# period-mean model residual at map_depth (mean of the daily means over the school days)
iz <- which.min(abs(dy$dep - map_depth))
it <- which(dy$time >= min(d$tt) - 43200 & dy$time <= max(d$tt) + 43200)
fld <- CJ(i = seq(1, length(dy$lon), 3), j = seq(1, length(dy$lat), 3))
fld[, `:=`(lon = dy$lon[i], lat = dy$lat[j],
           ur = apply(dy$u[cbind(rep(i, each = length(it)), rep(j, each = length(it)), iz, it)] |>
                      matrix(nrow = length(it)), 2, mean),
           vr = apply(dy$v[cbind(rep(i, each = length(it)), rep(j, each = length(it)), iz, it)] |>
                      matrix(nrow = length(it)), 2, mean))]
fld <- fld[!is.na(ur) & !is.na(vr)]

# ================================ FIGURES ================================================================
k_lon <- function(lat) 111320 * cos(lat * pi / 180)
ship_bg <- read_ship()[seq(1, .N, 10)]
cols3 <- c("Schools (ground)" = "black", "CMEMS current" = "steelblue", "Swimming" = "firebrick")

# 09a: clock check ----------------------------------------------------------------------------------------
p9a <- ggplot(lag_scan, aes(lag_h, r)) +
  geom_vline(xintercept = lag_use, colour = "firebrick", linetype = 2) +
  geom_line() + geom_point(size = 1.5) +
  scale_x_continuous(breaks = -3:3) +
  labs(x = "Assumed sonar clock - UTC (h)", y = "Vector correlation, hourly means",
       title = "Clock check against the model tide",
       subtitle = sprintf("Best alignment %+.2f h (r = %.2f); dashed = offset used (%+d h)", best_lag, max(lag_scan$r), lag_use))
ggsave(file.path(out, "09a_clock_check.png"), p9a, width = 7, height = 4.5, dpi = 200)

# 09b: hourly time series --------------------------------------------------------------------------------
hl <- melt(hourly[, .(t, n, `East|Schools (ground)` = vx, `North|Schools (ground)` = vy,
                      `East|CMEMS current` = u, `North|CMEMS current` = v)], id.vars = c("t", "n"))
hl[, c("comp", "what") := tstrsplit(variable, "|", fixed = TRUE)]
hl[, comp := paste(comp, "velocity (m/s)")]
p9b <- ggplot(hl, aes(t, value, colour = what)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(aes(group = interaction(what, as.Date(t))), alpha = 0.7) +
  geom_point(data = hl[what == "Schools (ground)"], aes(size = n), alpha = 0.5) +
  facet_wrap(~ comp, ncol = 1) +
  scale_colour_manual(values = cols3, name = NULL) + scale_size_area(max_size = 2.5, guide = "none") +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%d %b") +
  labs(x = "Sonar clock", y = NULL, title = "Hourly mean school velocity and model current at the same schools",
       subtitle = sprintf("Model at school position, depth and time (clock offset %+d h); slope schools on model: east %.2f, north %.2f",
                          lag_use, adv["x"], adv["y"])) +
  theme(legend.position = "bottom")
ggsave(file.path(out, "09b_tide_timeseries.png"), p9b, width = 13, height = 6, dpi = 200)

# 09c: vector comparison ---------------------------------------------------------------------------------
vec <- rbind(v_obs$sum[, what := "Schools (ground)"], v_mod$sum[, what := "CMEMS current"],
             v_swim$sum[, what := "Swimming"])
cloud <- rbind(data.table(v_obs$boot, what = "Schools (ground)"), data.table(v_mod$boot, what = "CMEMS current"),
               data.table(v_swim$boot, what = "Swimming"))
lim <- max(abs(unlist(cloud[, .(ax, ay)]))) * 86.4 * 1.1
p9c <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(data = cloud, aes(ax * 86.4, ay * 86.4, colour = what), alpha = 0.08, size = 0.6) +
  geom_segment(data = vec, aes(0, 0, xend = vx * 86.4, yend = vy * 86.4, colour = what),
               arrow = arrow(length = unit(2.5, "mm")), linewidth = 1.3) +
  geom_point(data = v_res$sum, aes(vx * 86.4, vy * 86.4), shape = 4, size = 3, stroke = 1.2, colour = "steelblue") +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  scale_colour_manual(values = cols3, name = NULL) +
  labs(x = "East (km/day)", y = "North (km/day)", title = "Net drift: schools, current and swimming",
       subtitle = "Arrows = estimate, dots = bootstrap replicates;\ncross = CMEMS daily-mean current at the schools") +
  theme(legend.position = "bottom")
ggsave(file.path(out, "09c_vector_comparison.png"), p9c, width = 7.5, height = 8, dpi = 200)

# 09d: map ------------------------------------------------------------------------------------------------
fld_scale <- 0.9 * 3 * 1500 / quantile(sqrt(fld$ur^2 + fld$vr^2), 0.95)
fld[, `:=`(lon_e = lon + ur * fld_scale / k_lon(lat), lat_e = lat + vr * fld_scale / 111320)]
cell_scale <- 1.6 * grid_km_m * 1000 / max(sqrt(c(cells$rx, cells$qx)^2 + c(cells$ry, cells$qy)^2))
cl <- rbind(cells[, .(lon, lat, lon_e = lon + rx * cell_scale / k_lon(lat), lat_e = lat + ry * cell_scale / 111320, what = "Schools (ground)")],
            cells[, .(lon, lat, lon_e = lon + qx * cell_scale / k_lon(lat), lat_e = lat + qy * cell_scale / 111320, what = "Swimming")])
p9d <- ggplot() +
  geom_segment(data = fld, aes(lon, lat, xend = lon_e, yend = lat_e), colour = "steelblue", alpha = 0.6,
               arrow = arrow(length = unit(1, "mm")), linewidth = 0.35) +
  geom_path(data = ship_bg, aes(lon, lat, group = seg), colour = "grey80") +
  geom_segment(data = cl, aes(lon, lat, xend = lon_e, yend = lat_e, colour = what),
               arrow = arrow(length = unit(1.6, "mm")), linewidth = 0.8) +
  scale_colour_manual(values = cols3, name = NULL) +
  coord_quickmap(xlim = range(d$lon) + c(-0.1, 0.1), ylim = range(d$lat) + c(-0.05, 0.05)) +
  labs(x = "Longitude", y = "Latitude", title = "Model residual current and school net drift",
       subtitle = sprintf("Blue = CMEMS mean current at %d m, %s - %s (daily means, tide removed)\nper %d-km cell: black = school net drift, red = swimming (arrows not to the same scale)",
                          round(dy$dep[iz]), format(min(d$day), "%d %b"), format(max(d$day), "%d %b"), grid_km_m)) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
ggsave(file.path(out, "09d_current_map.png"), p9d, width = 9, height = 8, dpi = 200)

# 09e: progressive vector diagram -------------------------------------------------------------------------
pvd <- rbind(hourly[, .(hb, x = ox, y = oy, what = "Schools (ground)")], hourly[, .(hb, x = mx, y = my, what = "CMEMS current")],
             hourly[, .(hb, x = wx, y = wy, what = "Swimming")])
pvd <- rbind(data.table(hb = min(hourly$hb) - 1, x = 0, y = 0, what = names(cols3)), pvd)
p9e <- ggplot(pvd, aes(x, y, colour = what)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_path(linewidth = 0.8) + geom_point(data = pvd[, .SD[.N], by = what], size = 2.5) +
  scale_colour_manual(values = cols3, name = NULL) + coord_equal() +
  labs(x = "Cumulative east displacement (km)", y = "Cumulative north displacement (km)",
       title = "Progressive vector diagram",
       subtitle = sprintf("Hourly means over the %d daytime hours with schools; swimming = schools - current", obs_h)) +
  theme(legend.position = "bottom")
ggsave(file.path(out, "09e_progressive_vector.png"), p9e, width = 7.5, height = 7.5, dpi = 200)

# 09f: estimates across variants --------------------------------------------------------------------------
unwrap <- function(x, ref) ref + ((x - ref + 180) %% 360) - 180
ef <- copy(est)[, model := factor(model, rev(model))]
ef[, grp := fcase(grepl("without", model), "Leave one day out", grepl("^Swimming", model), "Swimming",
                  default = "Schools / current")]
ref_dir <- v_swim$sum$dir
ef[, `:=`(dir_u = unwrap(dir, ref_dir), lo_u = unwrap(dir_lo, ref_dir), hi_u = unwrap(dir_hi, ref_dir))]
gcol <- c("Schools / current" = "steelblue", "Swimming" = "firebrick", "Leave one day out" = "grey55")
p9f1 <- ggplot(ef, aes(km_day, model, colour = grp)) +
  geom_errorbarh(aes(xmin = speed_lo * 86.4, xmax = speed_hi * 86.4), height = 0.3) + geom_point(size = 2) +
  scale_colour_manual(values = gcol, name = NULL) + labs(x = "Net speed (km/day)", y = NULL, title = "Speed")
p9f2 <- ggplot(ef, aes(dir_u, model, colour = grp)) +
  geom_errorbarh(aes(xmin = lo_u, xmax = hi_u), height = 0.3) + geom_point(size = 2) +
  scale_x_continuous(breaks = seq(-360, 720, 45), labels = function(b) b %% 360) +
  scale_colour_manual(values = gcol, guide = "none") +
  labs(x = "Direction (deg, compass)", y = NULL, title = "Direction") + theme(axis.text.y = element_blank())
ggsave(file.path(out, "09f_estimates.png"), p9f1 + p9f2 + plot_layout(guides = "collect"), width = 13, height = 6, dpi = 200)

# tidal ellipse comparison (summary panel D)
ell <- rbind(ellipse(v_obs$fit, "Schools (ground)"), ellipse(v_mod$fit, "CMEMS current"))
p9g <- ggplot(ell, aes(x, y, colour = what)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_path(linewidth = 0.9) + coord_equal() + scale_colour_manual(values = cols3, name = NULL) +
  labs(x = "East (m/s)", y = "North (m/s)", title = "Fitted tidal ellipses") + theme(legend.position = "bottom")

# summary figure --------------------------------------------------------------------------------------------
sw <- v_swim$sum; mo <- v_mod$sum; ob <- v_obs$sum
p_sum <- (p9d + labs(title = "A  Model residual current and cell net drift", subtitle = NULL)) +
  (p9c + labs(title = "B  Net drift: schools, current, swimming", subtitle = NULL)) +
  (p9e + labs(title = "C  Progressive vector diagram", subtitle = NULL)) +
  (p9g + labs(title = "D  Tidal ellipses (12.42 h)")) +
  plot_layout(ncol = 2) +
  plot_annotation(title = sprintf("Schools %.1f km/day to %.0f deg; CMEMS current at the schools %.1f km/day to %.0f deg; swimming %.1f km/day to %.0f deg (95%% CI %.0f-%.0f deg)",
                                  ob$km_day, ob$dir, mo$km_day, mo$dir, sw$km_day, sw$dir, sw$dir_lo, sw$dir_hi),
                  subtitle = sprintf("%d daytime schools, %d hour blocks, %s - %s; model sampled at school position, depth and time (sonar clock %+d h from UTC)",
                                     nrow(d), uniqueN(d$hb), format(min(d$day), "%d %b"), format(max(d$day), "%d %b %Y"), lag_use))
ggsave(file.path(out, "09_cmems_summary.png"), p_sum, width = 15, height = 13.5, dpi = 200)
message("Figures written to ", out)

# ================================ WORD SUMMARY ==========================================================
f0 <- function(x) formatC(round(x), format = "d", big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
compass <- function(b) c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")[round(b / 22.5) %% 16 + 1]
vtxt <- function(s) sprintf("%s km/day towards %.0f° (%s; 95%% CI %s–%s km/day, %.0f–%.0f°)", f1(s$km_day), s$dir, compass(s$dir),
                            f1(s$speed_lo * 86.4), f1(s$speed_hi * 86.4), s$dir_lo, s$dir_hi)
loo_dir <- range(unwrap(loo$dir, sw$dir)) %% 360
ang <- function(a, b) abs(((a - b + 180) %% 360) - 180)
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

species <- if (cfg$trip == "2025K") "mackerel" else "herring"
doc <- body_add_fpar(doc, fpar(ftext(sprintf("Daytime %s drift versus CMEMS model currents", species), fp_text(bold = TRUE, font.size = 18))))
doc <- body_add_fpar(doc, fpar(ftext(sprintf("SCH302, %s sonar, trip %s, %s – %s. Generated %s with 09_cmems_comparison.R",
  cfg$label, cfg$trip, format(min(d$day), "%d %B"), format(max(d$day), "%d %B %Y"), format(Sys.Date(), "%d %B %Y")),
  fp_text(font.size = 10, color = "grey40"))))
h1("Main result")
explains <- ang(mo$dir, ob$dir) < 45 && mo$speed > 0.5 * ob$speed      # current carries most of the drift
swim_clear <- sw$p_zero < 0.025 && ang(sw$dir_lo, sw$dir_hi) < 90 && sw$dir_lo != sw$dir_hi
lead(if (explains) "The current explains most of the drift." else "The current does not explain the drift.", sprintf(paste(
  "At the same places, depths and times as the schools, the CMEMS model current has a net component of %s, %.0f° from",
  "the school drift over ground (%s). The daily-mean (tide-free) model current at the schools gives %s km/day towards %.0f°.",
  "Over the survey area the model residual current at %.0f m has a median speed of %s m/s."),
  vtxt(mo), ang(mo$dir, ob$dir), vtxt(ob), f1(v_res$sum$km_day), v_res$sum$dir, dy$dep[iz], f2(median(sqrt(fld$ur^2 + fld$vr^2)))))
lead("Swimming.", sprintf(paste(
  "Subtracting the model current from the school velocities gives a net swimming velocity of %s, %.0f° from the ground drift.",
  "In %.1f%% of the %d bootstrap replicates it pointed the other way.%s"), vtxt(sw), ang(sw$dir, ob$dir), 100 * sw$p_zero, B,
  if (swim_clear) "" else " The swimming direction is therefore not well determined."))
lead("Robustness.", sprintf(paste(
  "The swimming direction ranges from %.0f° to %.0f° for the depth-mean (0–75 m) and surface currents, with the sonar",
  "clock taken as UTC, for slow ship speeds and for clearly moving schools, and %.0f–%.0f° when any single day is left out.",
  "On %d of %d days with at least 50 schools the daily swimming vector was within 90° of the overall direction, and %.0f%% of",
  "schools were in %d-km cells where it was (cell directions R̄ = %s)."),
  min(unwrap(est[grepl("^Swimming", model) & !grepl("without", model)]$dir, sw$dir)) %% 360,
  max(unwrap(est[grepl("^Swimming", model) & !grepl("without", model)]$dir, sw$dir)) %% 360,
  loo_dir[1], loo_dir[2], swim_days_w, n_days50, 100 * share_w_swim, grid_km_m, f2(cell_dirs_swim$Rbar)))
lead("Tide.", sprintf(paste(
  "The model tide at the schools has a 12.42-h amplitude of %s m/s east and %s m/s north; the tidal signal in the school",
  "velocities is %s and %s m/s. Hourly mean school velocity rises by %s (east) and %s (north) m/s per m/s of model current,",
  "so schools %s. The phase agrees within %s h."),
  f2(amp_mod["x"]), f2(amp_mod["y"]), f2(amp_obs["x"]), f2(amp_obs["y"]), f2(adv["x"]), f2(adv["y"]),
  if (mean(adv) > 0.85) "move largely with the tidal stream"
  else "follow the tidal stream only partly (either they compensate part of it by swimming, or the sonar tracks underestimate it)",
  f1(max(abs(pha_dif)))))
lead("Clock.", sprintf(paste(
  "School velocities line up best with the model tide when the sonar clock is taken as %+.2f h from UTC (r = %s, versus %s",
  "at 0 h); the analysis used %+d h%s."),
  best_lag, f2(max(lag_scan$r)), f2(lag_scan[lag_h == 0, r]), lag_use,
  if (nzchar(Sys.getenv("CMEMS_LAG_H"))) " (set from another sonar on the same trip)" else ""))
lead("Interpretation.", if (explains) paste(
  sprintf("The daytime net displacement of the %s schools is largely what the model current would carry them;", species),
  sprintf("the swimming that remains (%s km/day) is small compared with the drift and the uncertainty of the model's", f1(sw$km_day)),
  "residual current, so these data show no clear directed movement relative to the water.")
  else paste(
  sprintf("Within the limits of the model, the daytime displacement of the %s schools is not passive drift with the", species),
  sprintf("residual current; it is consistent with directed swimming of about %s km/day towards %.0f° relative to the water.", f1(sw$km_day), sw$dir),
  if (species == "herring") "The direction is towards the Orkney–Shetland spawning grounds, which fits a pre-spawning movement, but a population-level bias on schools that individually swim in all directions is not proof of migration." else ""))

h1("Data and approach")
para(sprintf(paste(
  "Model: CMEMS North-West European Shelf analysis and forecast (NWSHELF_ANALYSISFORECAST_PHY_004_013, version %s),",
  "1.5 km grid, hourly instantaneous currents with tides, levels 0–75 m; and the daily-mean currents of the same product",
  "(tide largely averaged out). Downloaded with 09a_cmems_download.py."), "202511"))
para(sprintf(paste(
  "Sampling: for each of the %s daytime schools of script 08 the model current was interpolated to the school's position,",
  "track mid time and beam-centre depth (median %.0f m). Swimming velocity = school ground velocity − model current.",
  "The net drift model of script 08 (intercept + 12.42-h tidal harmonic + ship velocity terms, 1-h block bootstrap with %d",
  "replicates) was fitted to the school ground velocity, to the model current at the schools and to the swimming velocity,",
  "so all three use the same samples in time and space."), f0(nrow(d)), median(d$depth), B))
para(sprintf(paste(
  "Clock: the offset between the sonar clock and UTC was estimated as the time shift of the model that best correlates the",
  "hourly mean school velocities with the model current at the schools (shifts of %+g to %+g h in steps of 0.25 h)."),
  min(lags_h), max(lags_h)))

h1("Results")
add_fig("09_cmems_summary.png", sprintf(paste(
  "Summary. A: CMEMS mean residual current at %.0f m (blue) and, per %d-km cell, school net drift (black) and swimming (red).",
  "B: net drift of the schools over ground, of the model current at the same schools and of swimming, with bootstrap",
  "replicates; cross = daily-mean model current at the schools. C: progressive vector diagram of hourly means. D: fitted",
  "12.42-h tidal ellipses of school velocity and model current."), dy$dep[iz], grid_km_m))
add_fig("09a_clock_check.png", "Correlation between hourly mean school velocity and model current at the schools for different assumed offsets between the sonar clock and UTC.")
add_fig("09b_tide_timeseries.png", "Hourly mean school ground velocity (black) and model current at the same schools (blue), east and north components.")
add_fig("09c_vector_comparison.png", "Net drift vectors (tide and vessel terms removed) with bootstrap replicates (1-h blocks).")
add_fig("09d_current_map.png", sprintf("CMEMS mean current at %.0f m over the survey days (daily-mean product, every 3rd grid cell) and per-cell net drift of schools (black) and swimming (red).", dy$dep[iz]))
add_fig("09e_progressive_vector.png", "Progressive vector diagram: hourly mean school ground velocity, model current at the schools and their difference (swimming), integrated over the observed daytime hours. The curves still contain the tide, sampled unevenly (daytime only), so their end points differ from the fitted net drift vectors.")
add_fig("09f_estimates.png", "Net drift speed and direction with 95% bootstrap intervals: schools and current (blue), swimming variants (red) and swimming leaving out one day at a time (grey).")
est_tab <- est[, .(Estimate = model, n = n, `km/day` = f1(km_day), `95% CI km/day` = paste0(f1(speed_lo * 86.4), "–", f1(speed_hi * 86.4)),
                   `Direction (°)` = round(dir), `95% CI (°)` = paste0(round(dir_lo), "–", round(dir_hi)))]
doc <- body_add_fpar(doc, fpar(ftext("Table 1. ", fp_cap_lab), ftext("Net drift estimates (block bootstrap, 1-h blocks).", fp_cap),
                               fp_p = fp_par(keep_with_next = TRUE, padding.top = 6)), style = "Table Caption")
ft <- flextable(est_tab); ft <- theme_booktabs(ft); ft <- fontsize(ft, size = 8.5, part = "all"); ft <- bold(ft, part = "header")
ft <- width(ft, width = c(2.6, 0.5, 0.6, 1.0, 0.8, 0.9))
doc <- body_add_flextable(doc, ft, align = "center")

h1("Caveats")
lead("Model skill.", "CMEMS NWS is a 1.5-km model with data assimilation of temperature, not of currents. Tidal currents are usually well reproduced, residual currents of a few cm/s much less so; local errors of this size are comparable to the swimming estimate.")
lead("School depth.", "Profos depth is the beam-centre depth (range × sin tilt), not a measured school depth. The result was therefore repeated with the surface current and the 0–75 m depth-mean current.")
lead("Sampling.", "Model and schools are compared only where and when schools were seen (daytime, near the vessel), so the model 'net current' is the current at these samples, not the long-term mean.")
lead("Orientation.", "Swimming here is the mean displacement relative to the water. Individual schools still swim in all directions; the net value is a population-level bias.")
print(doc, target = file.path(out, sprintf("SCH302_%s_cmems_comparison.docx", if (sonar_profile == "ST90") "ST90" else paste(cfg$label, cfg$trip, sep = "_"))))
message("09: done -> ", normalizePath(out))
