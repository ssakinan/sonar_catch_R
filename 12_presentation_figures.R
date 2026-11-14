################################################################################
# 12_presentation_figures.R
#
# Simple, large-text figures and one animation for presenting the results to
# the skipper (non-scientist audience). Speeds in knots, distances in nautical
# miles. Run after 00b-10 (ST90, SU90 2025 and SU90_2026 outputs must exist):
#   source("12_presentation_figures.R")
# then build the slides with:  python 12b_presentation.py
# Output: output_presentation/
################################################################################

sonar_profile <- "ST90"
source("00_setup_functions.R")
library(sf)

out <- "output_presentation"
dir.create(out, showWarnings = FALSE)
kn  <- 1.943844                 # m/s -> knots
nmd <- 86400 / 1852             # m/s -> nautical miles per day
col <- c(ground = "#16324F", current = "#2A9DB5", swim = "#E8833A", grey = "#9AA5B1", land = "#E3E1D9")
theme_set(theme_minimal(base_size = 20) +
            theme(panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white", colour = NA),
                  axis.title = element_text(colour = "grey25"), axis.text = element_text(colour = "grey35")))
save <- function(p, f, w = 12, h = 6.75) ggsave(file.path(out, f), p, width = w, height = h, dpi = 200)

st  <- fread("output_ST90/school_movement.csv")
su  <- fread("output/school_movement.csv")             # SU90 2025
qst <- fread("output_ST90/school_qc.csv")

# ---- A: where ---------------------------------------------------------------------------------------
land <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf")
land <- suppressWarnings(st_crop(st_make_valid(land), xmin = -4.5, xmax = 2.5, ymin = 57.6, ymax = 61.2))
pA <- ggplot() +
  geom_sf(data = land, fill = col["land"], colour = "grey60", linewidth = 0.3) +
  geom_point(data = su, aes(lon, lat), colour = col["current"], size = 0.8, alpha = 0.5) +
  geom_point(data = st, aes(lon, lat), colour = col["swim"], size = 0.8, alpha = 0.4) +
  annotate("text", x = -1.35, y = 60.35, label = "Shetland", size = 6.5, fontface = "bold", colour = "grey30") +
  annotate("text", x = -3.0, y = 59.15, label = "Orkney", size = 6.5, fontface = "bold", colour = "grey30") +
  annotate("text", x = 1.35, y = 60.55, label = "Oct 2025\nmackerel", size = 6, colour = col["current"], fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 1.35, y = 59.25, label = "Jun-Jul 2026\nherring", size = 6, colour = col["swim"], fontface = "bold", lineheight = 0.9) +
  coord_sf(xlim = c(-4, 2), ylim = c(58.4, 61), expand = FALSE) +
  labs(x = NULL, y = NULL) + theme(panel.background = element_rect(fill = "#EEF5F8", colour = NA), panel.grid = element_line(colour = "white"))
save(pA, "A_where.png", w = 9, h = 7.5)

# ---- B: cleaning: real schools vs echoes travelling with the ship ---------------------------------------
qt <- qst[tracked == TRUE & is.finite(speed_lm)]
qt[, what := ifelse(excluded, "Removed", "Kept")]
pB <- ggplot(qt[speed_lm * kn < 18], aes(speed_lm * kn, fill = what)) +
  geom_histogram(binwidth = 0.5, boundary = 0, colour = "white") +
  scale_fill_manual(values = c(Kept = col[["ground"]], Removed = col[["swim"]]), name = NULL) +
  annotate("label", x = 3.5, y = 1150, label = "Real schools:\nmostly under 2 knots", size = 6.5, hjust = 0, fill = "white", linewidth = 0) +
  annotate("label", x = 8.5, y = 330, label = "'Schools' moving at the ship's\nsteaming speed (12-14 knots):\nechoes from the ship itself", size = 6.5, hjust = 0, fill = "white", linewidth = 0) +
  annotate("segment", x = 12.5, xend = 12.5, y = 230, yend = 80, colour = col[["swim"]], linewidth = 1,
           arrow = arrow(length = unit(3, "mm"))) +
  scale_x_continuous(breaks = seq(0, 18, 2)) +
  labs(x = "Speed of the tracked 'school' (knots)", y = "Number of schools") +
  theme(legend.position = c(0.88, 0.85), legend.text = element_text(size = 18))
save(pB, "B_cleaning.png")

# ---- C: one school, raw vs smoothed --------------------------------------------------------------------
trk <- readRDS("output_ST90/school_tracks.rds")
cand <- st[moving == TRUE & n_pings >= 40 & straightness > 0.6 & net_m > 150][order(-n_pings)]
sid0 <- cand$sid[1]
k <- trk[sid == sid0]
xy  <- to_xy(k$lon, k$lat); xyk <- to_xy(k$lon_kf, k$lat_kf, mean(k$lon), mean(k$lat))
k[, `:=`(x = xy[, 1], y = xy[, 2], xk = xyk[, 1], yk = xyk[, 2])]
s0 <- st[sid == sid0]
pC <- ggplot(k) +
  geom_path(aes(x, y), colour = "grey75", linewidth = 0.4) +
  geom_point(aes(x, y), colour = col[["grey"]], size = 2.3) +
  geom_path(data = k[!is.na(xk)], aes(xk, yk), colour = col[["swim"]], linewidth = 2) +
  annotate("point", x = k$xk[which(!is.na(k$xk))[1]], y = k$yk[which(!is.na(k$xk))[1]], size = 5, colour = col[["ground"]]) +
  coord_equal() +
  labs(x = "metres east", y = "metres north",
       subtitle = sprintf("Grey: raw position at each ping (%d pings, %.0f min)
Orange: smoothed track (speed %.1f knots)",
                          nrow(k), s0$duration_s / 60, s0$speed_lm * kn)) +
  theme(plot.subtitle = element_text(size = 15, colour = "grey30"))
save(pC, "C_smoothing.png", w = 10, h = 6.75)

# ---- D: why sonar depth is not school depth: beam geometry ----------------------------------------------
beam <- rbindlist(lapply(list(c("ST90", -4, col[["swim"]]), c("SU90", -9, col[["current"]])), function(b) {
  r <- seq(0, 600, 5); tl <- as.numeric(b[2])
  data.table(sonar = b[1], colr = b[3], r = r, mid = r * sin(-tl * pi / 180),        # Profos depth = range x sin(tilt)
             lo = r * sin((-tl - 4) * pi / 180), hi = r * sin((-tl + 4) * pi / 180))
}))
pD <- ggplot(beam) +
  geom_ribbon(aes(r, ymin = hi, ymax = lo, fill = sonar), alpha = 0.18) +
  geom_line(aes(r, mid, colour = sonar), linewidth = 1.2, linetype = 2) +
  annotate("point", x = 400, y = 45, size = 28, shape = 16, colour = "grey30", alpha = 0.8) +
  annotate("text", x = 400, y = 45, label = "school", colour = "white", size = 4.8, fontface = "bold") +
  annotate("label", x = 425, y = 18, label = "ST90 reports 28 m", colour = col[["swim"]], size = 6.5, hjust = 0, linewidth = 0, fill = "white") +
  annotate("label", x = 425, y = 78, label = "SU90 reports 63 m", colour = col[["current"]], size = 6.5, hjust = 0, linewidth = 0, fill = "white") +
  annotate("point", x = 0, y = 0, shape = 17, size = 6, colour = "grey20") +
  annotate("text", x = 12, y = -7, label = "vessel", hjust = 0, size = 6, colour = "grey30") +
  scale_y_reverse(limits = c(130, -12)) + scale_x_continuous(limits = c(0, 640), breaks = seq(0, 600, 100)) +
  scale_colour_manual(values = c(ST90 = col[["swim"]], SU90 = col[["current"]]), name = NULL,
                      labels = c(ST90 = "ST90 beam (tilt 4°)", SU90 = "SU90 beam (tilt 9°)")) +
  scale_fill_manual(values = c(ST90 = col[["swim"]], SU90 = col[["current"]]), guide = "none") +
  labs(x = "Distance from the vessel (m)", y = "Depth (m)") +
  theme(legend.position = "top", legend.text = element_text(size = 18))
save(pD, "D_beams.png")

# ---- E: the tide moves schools back and forth ----------------------------------------------------------
w <- 2 * pi / (12.4206 * 3600)
d <- st[sun == "day"]
d[, t := as.POSIXct(t_start, tz = "UTC") + duration_s / 2]
d[, `:=`(tt = as.numeric(t), cw = cos(w * as.numeric(t)), sw = sin(w * as.numeric(t)), hb = floor(as.numeric(t) / 3600))]
fy <- coef(lm(vy ~ cw + sw + ship_vx + ship_vy, d)); fx <- coef(lm(vx ~ cw + sw + ship_vx + ship_vy, d))
win <- as.POSIXct(c("2026-07-01 00:00", "2026-07-03 00:00"), tz = "UTC")
hh <- d[t >= win[1] & t < win[2], .(t = mean(t), n = .N, vy = mean(vy)), by = hb]
cv <- data.table(t = seq(win[1], win[2], by = 600))
cv[, vy := fy[1] + fy["cw"] * cos(w * as.numeric(t)) + fy["sw"] * sin(w * as.numeric(t))]
pE <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(data = cv, aes(t, vy * kn), colour = col[["current"]], linewidth = 1.6) +
  geom_point(data = hh, aes(t, vy * kn, size = n), colour = col[["ground"]], alpha = 0.8) +
  annotate("text", x = win[1] + 3600, y = 0.95, label = "moving north", hjust = 0, size = 6.5, colour = "grey30") +
  annotate("text", x = win[1] + 3600, y = -0.95, label = "moving south", hjust = 0, size = 6.5, colour = "grey30") +
  scale_size_area(max_size = 6, guide = "none") +
  scale_x_datetime(date_breaks = "6 hours", date_labels = "%d %b\n%H:%M") +
  coord_cartesian(ylim = c(-1.1, 1.1)) +
  labs(x = NULL, y = "North-south speed (knots)")
save(pE, "E_tide.png")

# ---- F: current or swimming? ---------------------------------------------------------------------------
ce  <- fread("output_ST90/migration/cmems_estimates.csv"); ce25 <- fread("output/migration/cmems_estimates.csv")
pick <- function(e, lab) rbind(e[grepl("^Schools: ground", model)][, what := "Movement over the seabed"],
                               e[grepl("^CMEMS current at schools", model)][, what := "Current"],
                               e[grepl("^Swimming = ", model)][, what := "Swimming through the water"])[, panel := lab]
ar <- rbind(pick(ce, "Herring, Jun-Jul 2026"), pick(ce25, "Mackerel, Oct 2025"))
ar[, `:=`(ex = vx * nmd, ey = vy * nmd, lab = sprintf("%.1f nm/day", sqrt(vx^2 + vy^2) * nmd))]
ar[, what := factor(what, c("Movement over the seabed", "Current", "Swimming through the water"))]
ring <- rbindlist(lapply(c(3, 6), function(r) data.table(r = r, a = seq(0, 2 * pi, length = 200))[, `:=`(x = r * sin(a), y = r * cos(a))]))
pF <- ggplot(ar) +
  geom_path(data = ring, aes(x, y, group = r), colour = "grey85") +
  geom_text(data = data.table(x = c(0, 7.3, 0, -7.3), y = c(7.3, 0, -7.3, 0), l = c("N", "E", "S", "W")),   # annotate() breaks with facets
            aes(x, y, label = l), size = 7, colour = "grey45", fontface = "bold") +
  geom_text(data = data.table(x = 0.25, y = c(3, 6) + 0.35, l = c("3 nm/day", "6 nm/day")),
            aes(x, y, label = l), size = 4.5, colour = "grey55", hjust = 0) +
  geom_segment(aes(0, 0, xend = ex, yend = ey, colour = what), linewidth = 2.6, arrow = arrow(length = unit(5, "mm"), type = "closed")) +
  ggrepel::geom_label_repel(aes(ex, ey, label = lab, colour = what), size = 5.5, linewidth = 0, fill = alpha("white", 0.85),
                            min.segment.length = 0, box.padding = 0.7, seed = 2, show.legend = FALSE) +
  facet_wrap(~ panel) + coord_equal(xlim = c(-8, 8), ylim = c(-8, 8)) +
  scale_colour_manual(values = c("Movement over the seabed" = col[["ground"]], "Current" = col[["current"]],
                                 "Swimming through the water" = col[["swim"]]), name = NULL) +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(), panel.grid = element_blank(), legend.position = "bottom",
        strip.text = element_text(size = 21, face = "bold"), legend.text = element_text(size = 17))
save(pF, "F_current_or_swimming.png", w = 13, h = 7.3)

# ---- G: sonar versus catch -------------------------------------------------------------------------------
ind <- fread("output_ST90/haul_sonar_indices.csv", colClasses = list(character = "haul_id"))[phase == "tow"]
pG <- ggplot(ind, aes(schools_per_h, catch)) +
  geom_point(size = 5, colour = col[["ground"]], alpha = 0.85) +
  geom_text(aes(label = haul), vjust = -1, size = 4.5, colour = "grey40") +
  labs(x = "Schools seen per hour during the tow (ST90)", y = "Catch (tonnes)") +
  scale_y_continuous(limits = c(0, max(ind$catch) * 1.1))
save(pG, "G_catch.png", w = 10, h = 6.75)

# ---- H: animation: the average daylight drift of herring, drawn on the map --------------------------------
X <- cbind(1, as.matrix(d[, .(cw, sw, ship_vx, ship_vy)]))
d[, `:=`(rx = vx - X[, -1] %*% fx[-1], ry = vy - X[, -1] %*% fy[-1])]
hr <- d[is.finite(rx) & is.finite(ry), .(t = mean(t), rx = mean(rx), ry = mean(ry)), by = hb][order(hb)]
hr[, `:=`(px = cumsum(rx) * 3600, py = cumsum(ry) * 3600)]
lon0 <- median(st$lon); lat0 <- median(st$lat)
pth <- from_xy(c(0, hr$px), c(0, hr$py), lon0, lat0)
path <- data.table(hours = 0:nrow(hr), lon = pth[, 1], lat = pth[, 2], t = c(min(hr$t) - 1800, hr$t),
                   km = sqrt(c(0, hr$px)^2 + c(0, hr$py)^2) / 1852)
land2 <- suppressWarnings(st_crop(land, xmin = -3.8, xmax = 1.6, ymin = 58.5, ymax = 60.6))
fdir <- file.path(out, "frames_drift"); unlink(fdir, recursive = TRUE); dir.create(fdir)
draw_drift <- function(i) {
  pp <- path[seq_len(min(i - 1, nrow(hr)) + 1)]          # frame 1 = start; extra frames hold the end
  p <- ggplot() +
    geom_sf(data = land2, fill = col["land"], colour = "grey60", linewidth = 0.3) +
    geom_point(data = st, aes(lon, lat), colour = "grey70", size = 0.5, alpha = 0.3) +
    geom_path(data = pp, aes(lon, lat), colour = col[["swim"]], linewidth = 1.8) +
    geom_point(data = pp[.N], aes(lon, lat), colour = col[["swim"]], size = 6) +
    annotate("point", x = lon0, y = lat0, size = 4, colour = col[["ground"]]) +
    annotate("text", x = -1.35, y = 60.35, label = "Shetland", size = 6, fontface = "bold", colour = "grey35") +
    annotate("text", x = -3.0, y = 59.05, label = "Orkney", size = 6, fontface = "bold", colour = "grey35") +
    annotate("text", x = -1.62, y = 59.47, label = "Fair Isle", size = 4.5, colour = "grey35") +
    annotate("label", x = 1.55, y = 60.55, hjust = 1, vjust = 1, size = 6.5, linewidth = 0, fill = "white",
             label = sprintf("%s\n%d daylight hours added up\n%.0f nm from the start", format(pp$t[nrow(pp)], "%d %b %H:00"), pp$hours[nrow(pp)], pp$km[nrow(pp)])) +
    coord_sf(xlim = c(-3.8, 1.6), ylim = c(58.5, 60.6), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme(panel.background = element_rect(fill = "#EEF5F8", colour = NA), panel.grid = element_line(colour = "white"))
  p
}
n_fr <- nrow(hr) + 1
render_frames(draw_drift, n_fr + 12, fdir, 1600, 1200, res = 120)   # last 12 frames hold the end
encode_mp4(fdir, file.path(out, "H_herring_drift.mp4"), fps = 6)
message("12: figures and animation in ", normalizePath(out))
