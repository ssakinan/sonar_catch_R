################################################################################
# 07_vessel_school_animation.R
#
# Real-time animation: the vessel sails its recorded track and every school
# appears at its REAL detection time and REAL start position, then moves along
# its smoothed (Kalman) track from script 03.
#
# A school is tracked for only ~20-230 s while the day lasts ~13 h, so two
# visibility aids are applied to the schools (stated on every frame):
#   school_stretch : each school's own track plays this many times slower than
#                    the vessel clock (its START time stays real)
#   exaggerate_*   : school displacement is multiplied by this factor
# The vessel always moves at true speed relative to the clock.
#
# Views (choose with `views`):
#   map    : fixed map of the whole area; vessel, its track, schools
#   follow : camera follows the vessel (local metres around it, +/- follow_half_m),
#            i.e. schools seen relative to the vessel, like a sonar display;
#            set heading_up = TRUE to rotate the picture with the vessel heading
#
# Sonar-off gaps longer than cfg$gap_s are skipped. Frames are rendered in
# parallel (ragg PNG) and joined with ffmpeg into output/animation/*.mp4.
# Requires: script 03 has been run, ffmpeg on PATH.
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root

views          <- c("map", "follow")
clock_step_s   <- 20      # real seconds of vessel clock per frame
fps            <- 30      # -> video speed = clock_step_s * fps x real time (600x)
school_stretch <- 10      # school tracks play 10x slower than the clock
fade_s         <- 300     # clock seconds a finished school takes to fade out
trail_track_s  <- 30      # school trail length, in school-track seconds
exaggerate_map    <- 10   # school displacement multiplier, map view
exaggerate_follow <- 5    # school displacement multiplier, follow view
follow_half_m  <- 1500    # half-width of the follow window (m)
heading_up     <- FALSE   # follow view: rotate so the bow points up
show_history   <- TRUE    # map view: keep earlier schools as small grey dots
sonar_range_m  <- 600     # radius of the sonar-range circle drawn around the vessel
n_cores        <- max(1, min(12, parallel::detectCores(logical = FALSE) - 2))
show_legend <- FALSE      # TRUE: show colour (speed) and size (area) legends
skip_no_school_m <- NA    # e.g. 1000: jump over periods with no school shown within this many m of the vessel
skip_pad_frames  <- 5     # with skipping: frames kept before/after each school encounter
ffmpeg   <- "ffmpeg"

# Override any setting without editing this file, e.g.
#   anim_settings <- list(views = "follow", skip_no_school_m = 1000)
#   source("07_vessel_school_animation.R")
if (exists("anim_settings")) list2env(anim_settings, envir = environment())
width_px <- if (show_legend) 1400 else 1150; height_px <- 1200; res <- 150

anim_dir <- file.path(cfg$out_dir, "animation")
dir.create(anim_dir, showWarnings = FALSE, recursive = TRUE)

# ---- data ---------------------------------------------------------------------------------
hauls  <- read_catch(cfg$catch_file, trips = cfg$trip)
ship   <- ship_track(read_profos(cfg$profos_dir, "ppe"))
tracks <- readRDS(file.path(cfg$out_dir, "school_tracks.rds"))[keep == TRUE]
mov    <- fread(file.path(cfg$out_dir, "school_movement.csv"))

# per school: real start time, anchor (first smoothed position), track on a 1-s grid
setorder(tracks, sid, t)
sch <- tracks[, {
  tg <- seq(0, max(tt), by = 1)
  .(tt = tg,
    dx = approx(tt, x_kf - x_kf[1], tg, rule = 2)$y,
    dy = approx(tt, y_kf - y_kf[1], tg, rule = 2)$y,
    t0 = as.numeric(t[1]), lon0 = lon_kf[1], lat0 = lat_kf[1], dur = max(tt))
}, by = sid]
sch <- merge(sch, mov[, .(sid, Area, speed = speed_lm)], by = "sid")
starts <- unique(sch[, .(sid, t0, dur, lon0, lat0, Area, speed)])
starts[, t_vis_end := t0 + dur * school_stretch + fade_s]

# frame clock: recorded segments only (gaps skipped)
frames <- ship[, .(T = seq(as.numeric(min(t)), as.numeric(max(t)), by = clock_step_s)), by = seg]$T
ship[, tn := as.numeric(t)]
# unwrap heading for interpolation
ship[, hd := { h <- heading; d <- c(0, diff(h)); d <- (d + 180) %% 360 - 180; h[1] + cumsum(d) }]
vessel_at <- function(T) {
  i <- findInterval(T, ship$tn, all.inside = TRUE)
  w <- (T - ship$tn[i]) / (ship$tn[i + 1] - ship$tn[i])
  list(lon = ship$lon[i] + w * (ship$lon[i + 1] - ship$lon[i]),
       lat = ship$lat[i] + w * (ship$lat[i + 1] - ship$lat[i]),
       heading = (ship$hd[i] + w * (ship$hd[i + 1] - ship$hd[i])) %% 360,
       speed = ship$speed[i])
}
ph <- label_phase(as.POSIXct(frames, origin = "1970-01-01", tz = "UTC"), hauls)

# ---- per-frame school state --------------------------------------------------------------------
school_state <- function(T) {
  act <- starts[t0 <= T & t_vis_end >= T]
  if (!nrow(act)) return(list(heads = NULL, trail = NULL))
  act[, tau := pmin((T - t0) / school_stretch, dur)]
  act[, alpha := ifelse((T - t0) / school_stretch <= dur, 1,
                        pmax(0.1, 1 - (T - t0 - dur * school_stretch) / fade_s))]
  s <- sch[act[, .(sid, tau, alpha)], on = "sid"][tt <= tau & tt >= tau - trail_track_s]
  heads <- s[, .SD[.N], by = sid]
  s[, a := alpha * (1 - (tau - tt) / trail_track_s)]
  trail <- s[, if (.N > 1) .SD, by = sid]
  list(heads = heads, trail = if (nrow(trail)) trail else NULL)
}

# helpers to place exaggerated school positions
place_lonlat <- function(d, ex) {
  d[, `:=`(lon = lon0 + ex * dx / (111320 * cos(lat0 * pi / 180)),
           lat = lat0 + ex * dy / 111320)]
}
circle_xy <- function(r, n = 120) data.table(x = r * cos(seq(0, 2 * pi, length = n)),
                                             y = r * sin(seq(0, 2 * pi, length = n)))
vessel_poly <- function(L, heading) {          # simple ship outline, bow at +y
  p <- data.table(x = c(0, 0.35, 0.3, -0.3, -0.35) * L / 2.2,
                  y = c(0.6, 0.15, -0.5, -0.5, 0.15) * L)
  a <- -heading * pi / 180
  data.table(x = p$x * cos(a) - p$y * sin(a), y = p$x * sin(a) + p$y * cos(a))
}
spd_scale  <- scale_colour_viridis_c(option = "plasma", limits = c(0, quantile(mov$speed_lm, 0.98)),
                                     oob = scales::squish, name = "School speed (m/s)")
size_scale <- scale_size_area(max_size = 5, limits = range(mov$Area), name = "School area (m2)")
phase_cols <- c(search = "orange", tow = "red", other = "grey30")

frame_label <- function(i, v) {
  st <- if (ph$phase[i] == "other") "" else sprintf("   |   %s %s",
          toupper(ph$phase[i]), sub(".*-", "haul ", ph$haul_id[i]))
  sprintf("%s UTC   |   vessel %.1f kn%s%s",
          format(as.POSIXct(frames[i], origin = "1970-01-01", tz = "UTC"), "%d %b %H:%M"),
          v$speed / 0.514444, st, skip_note[i])
}

# ---- optional: jump over periods without schools near the vessel ---------------------------------
skip_note <- rep("", length(frames)); file_suffix <- ""
if (!is.na(skip_no_school_m)) {
  near <- vapply(seq_along(frames), function(i) {
    ss <- school_state(frames[i])
    if (is.null(ss$heads)) return(FALSE)
    v <- vessel_at(frames[i])
    h <- place_lonlat(copy(ss$heads), exaggerate_follow)   # positions as drawn in the follow view
    xy <- to_xy(h$lon, h$lat, v$lon, v$lat)
    any(sqrt(xy[, 1]^2 + xy[, 2]^2) <= skip_no_school_m)
  }, logical(1))
  keep <- near
  for (k in seq_len(skip_pad_frames)) keep <- keep | shift(near, k, fill = FALSE) | shift(near, -k, fill = FALSE)
  gap <- c(0, diff(frames[keep]))
  frames <- frames[keep]; ph <- ph[keep]
  skip_note <- rep("", length(frames))
  for (j in which(gap > 1.5 * clock_step_s)) {
    skip_note[j:min(j + 2 * fps - 1, length(frames))] <-          # show for 2 s of video
      sprintf("   |   >> skipped %s", if (gap[j] >= 3600) sprintf("%.1f h", gap[j] / 3600)
                                     else sprintf("%.0f min", gap[j] / 60))
  }
  file_suffix <- sprintf("_skip%dm", skip_no_school_m)
  message(sprintf("Skipping: %d of %d frames kept (school within %d m), %d jumps",
                  sum(keep), length(keep), skip_no_school_m, sum(gap > 1.5 * clock_step_s)))
}

# ---- map view ----------------------------------------------------------------------------------
map_box <- map_coord(c(ship$lon, starts$lon0), c(ship$lat, starts$lat0))
ship_bg <- ship[seq(1, .N, 5)]
draw_map <- function(i) {
  T <- frames[i]; v <- vessel_at(T); ss <- school_state(T)
  done <- ship_bg[tn <= T]
  vp <- vessel_poly(500, v$heading); ll <- from_xy(vp$x, vp$y, v$lon, v$lat)
  vp[, `:=`(lon = ll[, 1], lat = ll[, 2])]
  rc <- circle_xy(sonar_range_m); ll <- from_xy(rc$x, rc$y, v$lon, v$lat)
  rc[, `:=`(lon = ll[, 1], lat = ll[, 2])]
  p <- ggplot() +
    geom_path(data = ship_bg, aes(lon, lat, group = seg), colour = "grey90") +
    geom_path(data = done, aes(lon, lat, group = seg), colour = "grey60", linewidth = 0.4)
  if (show_history) p <- p + geom_point(data = starts[t_vis_end < T], aes(lon0, lat0),
                                        colour = "grey70", size = 0.6)
  if (!is.null(ss$trail)) p <- p + geom_path(data = place_lonlat(ss$trail, exaggerate_map),
                                             aes(lon, lat, group = sid, colour = speed, alpha = a),
                                             linewidth = 0.6)
  if (!is.null(ss$heads)) p <- p + geom_point(data = place_lonlat(ss$heads, exaggerate_map),
                                              aes(lon, lat, colour = speed, size = Area, alpha = alpha))
  p + geom_path(data = rc, aes(lon, lat), colour = phase_cols[[ph$phase[i]]], linetype = 2, linewidth = 0.4) +
    geom_polygon(data = vp, aes(lon, lat), fill = phase_cols[[ph$phase[i]]], colour = "black", linewidth = 0.3) +
    map_box + spd_scale + size_scale + scale_alpha_identity() +
    labs(x = "Longitude", y = "Latitude", title = frame_label(i, v),
         subtitle = sprintf(paste0("Schools appear at real detection time/position; their tracks play %dx slower ",
                                   "and movement is exaggerated %dx. Dashed circle = %d m sonar range."),
                            school_stretch, exaggerate_map, sonar_range_m)) +
    theme(legend.position = if (show_legend) "right" else "none",
          plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8))
}

# ---- follow view ---------------------------------------------------------------------------------
draw_follow <- function(i) {
  T <- frames[i]; v <- vessel_at(T); ss <- school_state(T)
  rot <- if (heading_up) v$heading * pi / 180 else 0
  to_local <- function(d) {                    # exaggerated school position -> metres around vessel
    place_lonlat(d, exaggerate_follow)
    xy <- to_xy(d$lon, d$lat, v$lon, v$lat)
    d[, `:=`(x = xy[, 1] * cos(rot) - xy[, 2] * sin(rot), y = xy[, 1] * sin(rot) + xy[, 2] * cos(rot))]
  }
  win <- ship_bg_f[abs(tn - T) < 3600]
  xy <- to_xy(win$lon, win$lat, v$lon, v$lat)
  win[, `:=`(x = xy[, 1] * cos(rot) - xy[, 2] * sin(rot), y = xy[, 1] * sin(rot) + xy[, 2] * cos(rot))]
  vp <- vessel_poly(120, if (heading_up) 0 else v$heading)
  rings <- rbindlist(lapply(c(200, 400, sonar_range_m), function(r) circle_xy(r)[, r := r]))
  p <- ggplot() +
    geom_path(data = rings, aes(x, y, group = r), colour = "grey85", linetype = 3) +
    geom_path(data = win[tn <= T], aes(x, y, group = seg), colour = "grey55", linewidth = 0.5) +
    geom_path(data = win[tn > T], aes(x, y, group = seg), colour = "grey88", linewidth = 0.5)
  if (!is.null(ss$trail)) p <- p + geom_path(data = to_local(ss$trail),
                                             aes(x, y, group = sid, colour = speed, alpha = a), linewidth = 0.7)
  if (!is.null(ss$heads)) p <- p + geom_point(data = to_local(ss$heads),
                                              aes(x, y, colour = speed, size = Area, alpha = alpha))
  p + geom_polygon(data = vp, aes(x, y), fill = phase_cols[[ph$phase[i]]], colour = "black", linewidth = 0.3) +
    coord_equal(xlim = c(-1, 1) * follow_half_m, ylim = c(-1, 1) * follow_half_m, expand = FALSE) +
    spd_scale + size_scale + scale_alpha_identity() +
    labs(x = if (heading_up) "Starboard (m)" else "East of vessel (m)",
         y = if (heading_up) "Ahead (m)" else "North of vessel (m)",
         title = frame_label(i, v),
         subtitle = sprintf(paste0("Vessel-centred view%s. School tracks play %dx slower, movement exaggerated %dx. ",
                                   "Rings: 200, 400, %d m.%s"),
                            if (heading_up) " (heading up)" else " (north up)",
                            school_stretch, exaggerate_follow, sonar_range_m,
                            if (is.na(skip_no_school_m)) "" else
                              sprintf("\nPeriods without schools within %d m of the vessel are skipped.",
                                      skip_no_school_m))) +
    theme(legend.position = if (show_legend) "right" else "none",
          plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8))
}
ship_bg_f <- ship[seq(1, .N, 3)]

# ---- render (parallel) + encode ------------------------------------------------------------------
for (view in views) {
  fdir <- file.path(tempdir(), paste0("frames_vessel_", view))
  unlink(fdir, recursive = TRUE); dir.create(fdir)
  draw <- if (view == "map") draw_map else draw_follow
  render_chunk <- function(idx) {
    for (i in idx) {
      ragg::agg_png(file.path(fdir, sprintf("f_%05d.png", i)), width = width_px, height = height_px, res = res)
      print(draw(i)); dev.off()
    }
    length(idx)
  }
  message(sprintf("Rendering %d frames for view '%s' on %d cores ...", length(frames), view, n_cores))
  chunks <- split(seq_along(frames), cut(seq_along(frames), n_cores, labels = FALSE))
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(data.table); library(ggplot2) }))
    parallel::clusterExport(cl, setdiff(ls(globalenv()), "cl"), envir = globalenv())
    plot_theme <- theme_get()                       # workers start with ggplot's default theme
    parallel::clusterExport(cl, "plot_theme", envir = environment())
    parallel::clusterEvalQ(cl, theme_set(plot_theme))
    invisible(parallel::parLapply(cl, chunks, render_chunk))
    parallel::stopCluster(cl)
  } else invisible(lapply(chunks, render_chunk))

  mp4 <- normalizePath(file.path(anim_dir, sprintf("vessel_schools_%s%s.mp4", view, file_suffix)), mustWork = FALSE)
  status <- system2(ffmpeg, c("-y", "-loglevel", "error", "-framerate", fps,
                              "-i", shQuote(file.path(fdir, "f_%05d.png")),
                              "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "21",
                              "-vf", shQuote("scale=trunc(iw/2)*2:trunc(ih/2)*2"), shQuote(mp4)))
  if (status != 0) stop("ffmpeg failed for view ", view)
  message("Written: ", mp4)
}
