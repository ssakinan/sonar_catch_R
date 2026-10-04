################################################################################
# 02_sonar_catch_link.R
#
# How relevant is what the sonar saw for what ended up in the net?
#
#   Fig 02a  Echogram-like timeline: school detections (time x depth) with
#            tow periods and the trawl depth layer (headline .. headline +
#            vertical opening). Lower panel: schools per 10 min normalised by
#            sonar recording time, and ship speed.
#   Fig 02b  Depth coverage of the sonar: detected school depth vs range,
#            with beam-centre lines for each tilt. Detected depth is strongly
#            constrained by tilt and range -- keep this in mind before
#            comparing sonar school depth with headline depth.
#   Fig 02c  Per-haul zoom: ship track (search / tow), schools within buffer,
#            schools coloured by phase.
#   Fig 02d  Catch vs sonar indices across hauls (scales to the full trip).
#   Table    output/haul_sonar_indices.csv
#
# Sonar indices per haul and phase (search = pre_shoot_h before shoot,
# tow = shoot..haul), schools within buffer_km of the towing track:
#   n_schools          number of distinct schools
#   schools_per_10km   effort-normalised (distance sailed while recording)
#   sbi_db             10*log10( sum(sv_lin * Area) / pings recorded )
#                      = mean school backscatter per ping (relative index)
#   depth_w            backscatter-weighted school depth
#   frac_trawl_layer   fraction of school backscatter inside the trawl layer
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root

hauls <- read_catch(cfg$catch_file, trips = cfg$trip)
ppe   <- read_profos(cfg$profos_dir, "ppe")
pp    <- merge_pp(read_profos(cfg$profos_dir, "pp"))
ship  <- ship_track(ppe)
ship  <- cbind(ship, label_phase(ship$t, hauls))
pp    <- cbind(pp, label_phase(pp$t, hauls))

# Only hauls overlapping the sonar record
hz <- hauls[haul_time >= min(ship$t) & shoot_time - cfg$pre_shoot_h * 3600 <= max(ship$t)]
if (!nrow(hz)) stop("No hauls overlap the sonar data")

# ---- distance of each school ping to the tow track of its haul -------------------
pp[, dist_tow_km := NA_real_]
for (h in hz$haul_id) {
  trk <- ship[haul_id == h & phase == "tow"]
  if (!nrow(trk)) {   # no sonar during tow: fall back to logged tow line
    hh <- hz[haul_id == h]
    trk <- data.table(lon = seq(hh$shootlon, hh$haullon, length = 50),
                      lat = seq(hh$shootlat, hh$haullat, length = 50))
  }
  trk <- trk[seq(1, .N, length.out = min(.N, 300))]
  idx <- pp[haul_id == h, which = TRUE]
  if (length(idx)) {
    d <- distm(cbind(pp$lon[idx], pp$lat[idx]), cbind(trk$lon, trk$lat))
    pp[idx, dist_tow_km := apply(d, 1, min) / 1000]
  }
}
pp[, in_buffer := !is.na(dist_tow_km) & dist_tow_km <= cfg$buffer_km]

# ---- per haul x phase indices --------------------------------------------------------
eff <- ship[phase != "other",
            .(pings = .N, rec_h = sum(dt, na.rm = TRUE) / 3600,
              dist_km = sum(dist_m) / 1000, ship_speed = mean(speed)),
            by = .(haul_id, phase)]

idx <- pp[in_buffer == TRUE][hz, on = "haul_id", nomatch = 0][
  , .(n_schools = uniqueN(sid),
      sb_sum = sum(backscatter),
      depth_w = weighted.mean(depth, backscatter),
      depth_med = median(depth),
      frac_trawl_layer = sum(backscatter[depth >= headline_depth &
                                           depth <= headline_depth + vert_opening]) /
        sum(backscatter),
      area_mean = mean(Area), Sv_mean = 10 * log10(mean(sv_lin))),
  by = .(haul_id, phase)]

ind <- merge(eff, idx, by = c("haul_id", "phase"), all.x = TRUE)
ind[is.na(n_schools), `:=`(n_schools = 0L, sb_sum = 0)]
ind[, schools_per_10km := 10 * n_schools / dist_km]
ind[, sbi_db := 10 * log10(sb_sum / pings)]
ind <- merge(hz[, .(haul_id, haul, shoot_time, haul_time, headline_depth,
                    vert_opening, catch, dom_species, dom_frac)],
             ind, by = "haul_id")
fwrite(ind, file.path(cfg$out_dir, "haul_sonar_indices.csv"))
print(ind)

# ---- 02a  timeline ---------------------------------------------------------------------
hz[, layer_bot := headline_depth + vert_opening]
tl_pp <- pp[seq(1, .N, by = max(1, .N %/% 20000))]   # thin for plotting

p_top <- ggplot() +
  geom_rect(data = hz, aes(xmin = shoot_time - cfg$pre_shoot_h * 3600, xmax = shoot_time,
                           ymin = -Inf, ymax = Inf), fill = "orange", alpha = 0.12) +
  geom_rect(data = hz, aes(xmin = shoot_time, xmax = haul_time, ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 0.10) +
  geom_rect(data = hz, aes(xmin = shoot_time, xmax = haul_time,
                           ymin = headline_depth, ymax = layer_bot),
            fill = NA, colour = "red", linewidth = 0.6, linetype = 2) +
  geom_point(data = tl_pp, aes(t, depth, colour = Sv, size = Area), alpha = 0.5) +
  geom_label(data = hz, aes(x = shoot_time + (haul_time - shoot_time) / 2, y = 0,
                            label = sprintf("H%s: %.0f t %s (%.0f%%)", haul, catch,
                                            dom_species, 100 * dom_frac)),
             size = 2.8, vjust = 0) +
  scale_y_reverse(limits = c(NA, -5)) +
  scale_colour_viridis_c(option = "inferno", name = "Sv (dB)") +
  scale_size_area(max_size = 4, name = "Area (m2)") +
  labs(x = NULL, y = "School depth (m)",
       title = "School detections vs tows (red = tow, orange = search before shoot)",
       subtitle = "Dashed box = trawl layer (headline depth .. + vertical opening)")

# schools per 10 min, normalised by recording time
bin_s <- 600
ship[, bin := as.POSIXct(floor(as.numeric(t) / bin_s) * bin_s, origin = "1970-01-01", tz = "UTC")]
pp[,   bin := as.POSIXct(floor(as.numeric(t) / bin_s) * bin_s, origin = "1970-01-01", tz = "UTC")]
binned <- merge(ship[, .(rec_min = sum(dt, na.rm = TRUE) / 60, speed = mean(speed)), by = bin],
                pp[, .(n_sch = uniqueN(sid), sb = sum(backscatter), pings = uniqueN(t)), by = bin],
                by = "bin", all.x = TRUE)
binned[is.na(n_sch), `:=`(n_sch = 0L, sb = 0)]
binned[, sch_per_rec_h := ifelse(rec_min > 2, n_sch / rec_min * 60, NA)]

p_bot <- ggplot(binned, aes(bin)) +
  geom_rect(data = hz, aes(xmin = shoot_time, xmax = haul_time, ymin = -Inf, ymax = Inf),
            fill = "red", alpha = 0.10, inherit.aes = FALSE) +
  geom_col(aes(y = sch_per_rec_h), fill = "steelblue", width = bin_s * 0.9) +
  geom_line(aes(y = speed * 10), colour = "black") +
  scale_y_continuous(name = "Schools per\nrecording hour",
                     sec.axis = sec_axis(~ . / 10, name = "Ship speed (m/s)")) +
  labs(x = "Time (UTC)")

p2a <- p_top / p_bot + plot_layout(heights = c(3, 1.3), axes = "collect_x") &
  scale_x_datetime(limits = range(ship$t), date_labels = "%H:%M")
ggsave(file.path(cfg$out_dir, "02a_timeline_schools_vs_tows.png"), p2a,
       width = 12, height = 8, dpi = 200)

# ---- 02b  depth coverage of the sonar --------------------------------------------------
# NB: in the Oct-2025 export Center.dep == range * sin(|tilt|) exactly
# (R^2 = 1, intercept 0): Profos reports the BEAM-CENTRE depth at the school's
# range, not a measured depth. The real school can be anywhere within the
# vertical beam (band below). Depth indices in this script therefore mostly
# reflect range x tilt, not school depth -- use a vertical echosounder for that.
tilts <- sort(unique(pp$tilt))
fits <- pp[, {f <- lm(depth ~ range_m); .(intercept = coef(f)[1], slope = coef(f)[2],
                                          r2 = summary(f)$r.squared)}, by = tilt]
print(fits)
beam <- CJ(tilt = tilts, range_m = seq(0, max(pp$range_m), 10))
beam[, `:=`(centre = cfg$transducer_depth - range_m * sin(tilt * pi / 180),
            upper  = cfg$transducer_depth - range_m * sin((tilt + cfg$vert_beam_half) * pi / 180),
            lower  = cfg$transducer_depth - range_m * sin((tilt - cfg$vert_beam_half) * pi / 180))]
p2b <- ggplot(pp[seq(1, .N, by = max(1, .N %/% 20000))], aes(range_m, depth)) +
  geom_ribbon(data = beam, aes(range_m, ymin = upper, ymax = lower, fill = factor(tilt)),
              alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = beam, aes(range_m, centre, colour = factor(tilt)), linewidth = 0.8) +
  geom_point(aes(colour = factor(tilt)), size = 0.6, alpha = 0.4) +
  geom_hline(data = hz, aes(yintercept = headline_depth), linetype = 2) +
  scale_y_reverse() +
  labs(x = "Horizontal range from ship (m)", y = "School depth (m)",
       colour = "Tilt (deg)", fill = "Tilt (deg)",
       title = "Profos 'Center.dep' is range x sin(tilt): a beam-centre depth, not a measured depth",
       subtitle = sprintf(paste0("Fit per tilt: R2 = %s. Lines = beam centre incl. %s m transducer depth,",
                                 " bands = +/-%s deg beam; dashed = headline depth"),
                          paste(round(fits$r2, 4), collapse = " / "),
                          cfg$transducer_depth, cfg$vert_beam_half))
ggsave(file.path(cfg$out_dir, "02b_depth_vs_range_tilt.png"), p2b, width = 8, height = 6, dpi = 200)

# ---- 02c  per-haul zoom maps -------------------------------------------------------------
sch_h <- pp[phase != "other" & in_buffer == TRUE,
            .(lon = mean(lon), lat = mean(lat), depth = mean(depth), Area = mean(Area)),
            by = .(haul_id, phase, sid)]
p2c <- ggplot() +
  geom_path(data = ship[phase != "other"], aes(lon, lat, colour = phase, group = seg)) +
  geom_point(data = sch_h, aes(lon, lat, size = Area, fill = depth), shape = 21, alpha = 0.7) +
  geom_segment(data = hz, aes(shootlon, shootlat, xend = haullon, yend = haullat),
               linetype = 2, arrow = arrow(length = unit(2, "mm"))) +
  facet_wrap(~ haul_id, scales = "free") +
  scale_colour_manual(values = c(search = "orange", tow = "red")) +
  scale_fill_viridis_c(direction = -1, name = "Depth (m)") +
  scale_size_area(max_size = 5) +
  coord_quickmap() +
  labs(x = NULL, y = NULL,
       title = sprintf("Schools within %s km of the tow track (dashed = logged tow line)",
                       cfg$buffer_km))
ggsave(file.path(cfg$out_dir, "02c_haul_zoom_maps.png"), p2c, width = 11, height = 6, dpi = 200)

# ---- 02d  catch vs sonar indices ----------------------------------------------------------
long <- melt(ind, id.vars = c("haul", "phase", "catch", "dom_species"),
             measure.vars = c("schools_per_10km", "sbi_db", "depth_w", "frac_trawl_layer"))
p2d <- ggplot(long, aes(value, catch, colour = dom_species, shape = phase)) +
  geom_point(size = 3) +
  geom_text(aes(label = haul), vjust = -1, size = 2.5, show.legend = FALSE) +
  facet_wrap(~ variable, scales = "free_x",
             labeller = as_labeller(c(schools_per_10km = "Schools per 10 km sailed",
                                      sbi_db = "School backscatter index (dB)",
                                      depth_w = "Backscatter-weighted school depth (m)",
                                      frac_trawl_layer = "Fraction of backscatter in trawl layer"))) +
  scale_colour_manual(values = species_cols) +
  labs(x = NULL, y = "Catch (t)", title = "Catch vs sonar indices per haul",
       subtitle = "Only meaningful with many hauls: export Profos for the whole trip")
if (nrow(ind) >= 6) p2d <- p2d + geom_smooth(aes(group = phase), method = "lm",
                                             se = TRUE, colour = "grey30", linewidth = 0.5)
ggsave(file.path(cfg$out_dir, "02d_catch_vs_sonar_indices.png"), p2d, width = 10, height = 7, dpi = 200)

message("02: figures + haul_sonar_indices.csv written to ", cfg$out_dir)
