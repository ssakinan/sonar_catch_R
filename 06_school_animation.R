################################################################################
# 06_school_animation.R
#
# Animation of all tracked schools moving along their filtered + smoothed
# (Kalman) trajectories from script 03. Time is "time since first detection":
# every school starts at t = 0 together, regardless of its real time stamp.
#
# Two views (choose with `views`):
#   origin : every school starts at (0, 0) -> compares speeds and directions
#   map    : schools start at their real positions on the ship track;
#            displacement is exaggerated (`map_exaggerate`) because 50-100 m
#            of movement is invisible on a ~25 km map.
#
# Frames are drawn with ggplot2 (ragg PNG device) and joined with ffmpeg into
# MP4 (and optionally GIF) in output/animation/.
#
# Requires: script 03 has been run (output/school_tracks.rds), ffmpeg on PATH.
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root

views          <- c("origin", "map")
step_s         <- 1          # track seconds per frame
fps            <- 20         # frames per second -> playback speed = fps * step_s x real
trail_s        <- 30         # length of the fading trail (s); Inf = full path
moving_only    <- FALSE      # TRUE: only schools with significant movement
colour_by      <- "speed"    # "speed", "sun" or "phase"
map_exaggerate <- 25         # displacement multiplier in the map view
make_gif       <- FALSE
width_px <- 1400; height_px <- 1200; res <- 150
ffmpeg   <- "ffmpeg"         # or full path, e.g. "C:/ffmpeg/bin/ffmpeg.exe"

anim_dir <- file.path(cfg$out_dir, "animation")
dir.create(anim_dir, showWarnings = FALSE, recursive = TRUE)

# ---- data -----------------------------------------------------------------------------
tracks <- readRDS(file.path(cfg$out_dir, "school_tracks.rds"))[keep == TRUE]
mov    <- fread(file.path(cfg$out_dir, "school_movement.csv"))
hauls  <- read_catch(cfg$catch_file, trips = cfg$trip)
mov[, t_start := as.POSIXct(t_start, tz = "UTC")]
mov <- cbind(mov, label_phase(mov$t_start, hauls))
if (moving_only) mov <- mov[moving == TRUE]
tracks <- tracks[sid %in% mov$sid]

# smoothed position relative to the first smoothed point of each school
tracks[, `:=`(dx = x_kf - x_kf[1], dy = y_kf - y_kf[1]), by = sid]

# interpolate every school onto a common clock (seconds since first detection)
t_max  <- ceiling(max(tracks$tt))
frames <- seq(0, t_max + 10, by = step_s)   # +10 s so the last tracks fade out
grid <- tracks[, {
  tf <- frames[frames <= max(tt)]
  .(t = tf, dx = approx(tt, dx, tf, rule = 2)$y, dy = approx(tt, dy, tf, rule = 2)$y,
    lon0 = mean(lon_kf), lat0 = mean(lat_kf), t_end = max(tt))
}, by = sid]
grid <- merge(grid, mov[, .(sid, Area, speed = speed_lm, sun, phase, moving)], by = "sid")
grid[, `:=`(lon0 = lon0[1], lat0 = lat0[1]), by = sid]

colour_scale <- switch(colour_by,
  speed = scale_colour_viridis_c(option = "plasma", limits = c(0, quantile(mov$speed_lm, 0.98)),
                                 oob = scales::squish, name = "Speed (m/s)"),
  sun   = scale_colour_manual(values = c(day = "gold3", twilight = "orchid", night = "navy"),
                              name = "Light"),
  phase = scale_colour_manual(values = c(search = "orange", tow = "red", other = "grey40"),
                              name = "Ship activity"))
size_scale <- scale_size_area(max_size = 5, limits = range(mov$Area), name = "School area (m2)")

# ---- frame builder ---------------------------------------------------------------------------
draw_frame <- function(tnow, view) {
  trail <- grid[t <= tnow & t >= tnow - trail_s][, if (.N > 1) .SD, by = sid]
  heads <- grid[t <= tnow, .SD[.N], by = sid]          # current (or final) position
  # schools whose track has ended stay as faded points at their final position
  heads[, alpha := ifelse(tnow > t_end, pmax(0.15, 1 - (tnow - t_end) / 10), 1)]
  trail[, a := 1 - (tnow - t) / trail_s]

  if (view == "origin") {
    lim <- max(abs(c(grid$dx, grid$dy))) * 1.05
    p <- ggplot() +
      geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
      geom_path(data = trail, aes(dx, dy, group = sid, colour = .data[[colour_by]], alpha = a),
                linewidth = 0.5) +
      geom_point(data = heads, aes(dx, dy, colour = .data[[colour_by]], size = Area, alpha = alpha)) +
      coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
      labs(x = "East displacement (m)", y = "North displacement (m)",
           title = "All schools from a common start point",
           subtitle = sprintf("t = %3.0f s since first detection   |   %d schools active",
                              tnow, sum(heads$alpha == 1)))
  } else {
    ex <- map_exaggerate
    trail[, `:=`(lon = lon0 + ex * dx / (111320 * cos(lat0 * pi / 180)), lat = lat0 + ex * dy / 111320)]
    heads[, `:=`(lon = lon0 + ex * dx / (111320 * cos(lat0 * pi / 180)), lat = lat0 + ex * dy / 111320)]
    p <- ggplot() +
      geom_path(data = ship_bg, aes(lon, lat, group = seg), colour = "grey85") +
      geom_path(data = trail, aes(lon, lat, group = sid, colour = .data[[colour_by]], alpha = a),
                linewidth = 0.5) +
      geom_point(data = heads, aes(lon, lat, colour = .data[[colour_by]], size = Area, alpha = alpha)) +
      map_box +
      labs(x = "Longitude", y = "Latitude",
           title = sprintf("Schools at their observed positions (movement exaggerated %dx)", ex),
           subtitle = sprintf("t = %3.0f s since first detection   |   %d schools active",
                              tnow, sum(heads$alpha == 1)))
  }
  p + colour_scale + size_scale + scale_alpha_identity() +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
}

# map background + fixed extent (all positions incl. exaggerated movement)
if ("map" %in% views) {
  ship_bg <- ship_track(read_profos(cfg$profos_dir, "ppe"))[seq(1, .N, 5)]
  ex <- map_exaggerate
  all_lon <- c(grid$lon0 + ex * grid$dx / (111320 * cos(grid$lat0 * pi / 180)), ship_bg$lon)
  all_lat <- c(grid$lat0 + ex * grid$dy / 111320, ship_bg$lat)
  map_box <- map_coord(all_lon, all_lat)
}

# ---- render ---------------------------------------------------------------------------------
for (view in views) {
  fdir <- file.path(tempdir(), paste0("frames_", view))
  unlink(fdir, recursive = TRUE); dir.create(fdir)
  message(sprintf("Rendering %d frames for view '%s' ...", length(frames), view))
  for (i in seq_along(frames)) {
    ragg::agg_png(file.path(fdir, sprintf("f_%05d.png", i)),
                  width = width_px, height = height_px, res = res)
    print(draw_frame(frames[i], view))
    dev.off()
  }
  mp4 <- normalizePath(file.path(anim_dir, sprintf("schools_%s.mp4", view)), mustWork = FALSE)
  args <- c("-y", "-loglevel", "error", "-framerate", fps,
            "-i", shQuote(file.path(fdir, "f_%05d.png")),
            "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
            "-vf", shQuote("scale=trunc(iw/2)*2:trunc(ih/2)*2"), shQuote(mp4))
  status <- system2(ffmpeg, args)
  if (status != 0) stop("ffmpeg failed for view ", view)
  message("Written: ", mp4)

  if (make_gif) {
    gif <- sub("\\.mp4$", ".gif", mp4)
    system2(ffmpeg, c("-y", "-loglevel", "error", "-i", shQuote(mp4), "-vf",
                      shQuote("fps=10,scale=700:-1:flags=lanczos,split[a][b];[a]palettegen[p];[b][p]paletteuse"),
                      shQuote(gif)))
    message("Written: ", gif)
  }
}
