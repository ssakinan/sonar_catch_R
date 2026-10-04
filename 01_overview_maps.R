################################################################################
# 01_overview_maps.R
#
# Where were we fishing, what did we catch, and where/when do we have sonar?
#   Fig 01a  Trip map: hauls (shoot -> haul) coloured by dominant species,
#            sized by catch; sonar ship track overlaid.
#   Fig 01b  Trip timeline: hauls vs periods with sonar recording -> shows
#            which hauls can be linked to sonar once all days are exported.
#   Fig 01c  Sonar zoom map: ship track coloured by activity (search / tow /
#            other), school detections sized by school area, tow lines.
################################################################################

source("00_setup_functions.R")   # run from the RStudio project root
library(sf)

hauls <- read_catch(cfg$catch_file, trips = cfg$trip)
ppe   <- read_profos(cfg$profos_dir, "ppe")
pp    <- merge_pp(read_profos(cfg$profos_dir, "pp"))
ship  <- ship_track(ppe)
ship  <- cbind(ship, label_phase(ship$t, hauls))
schools <- school_summary(pp)

# ---- 01a  trip map --------------------------------------------------------------
land <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
bb   <- hauls[, .(xmin = min(shootlon, haullon) - 1, xmax = max(shootlon, haullon) + 1,
                  ymin = min(shootlat, haullat) - 0.5, ymax = max(shootlat, haullat) + 0.5)]

p1a <- ggplot() +
  geom_sf(data = land, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  geom_path(data = ship[seq(1, .N, 20)], aes(lon, lat, group = seg),
            colour = "black", linewidth = 0.4) +
  geom_segment(data = hauls, aes(shootlon, shootlat, xend = haullon, yend = haullat,
                                 colour = dom_species), linewidth = 0.8) +
  geom_point(data = hauls, aes(haullon, haullat, colour = dom_species, size = catch),
             alpha = 0.7) +
  geom_text(data = hauls, aes(haullon, haullat, label = haul), size = 2.3,
            vjust = -1.1) +
  scale_colour_manual(values = species_cols, name = "Dominant\nspecies") +
  scale_size_area(max_size = 6, name = "Catch (t)") +
  coord_sf(xlim = c(bb$xmin, bb$xmax), ylim = c(bb$ymin, bb$ymax)) +
  labs(x = NULL, y = NULL,
       title = paste("Trip", cfg$trip, "- hauls and sonar ship track (black)"))
ggsave(file.path(cfg$out_dir, "01a_trip_map.png"), p1a, width = 9, height = 7, dpi = 200)

# ---- 01b  trip timeline: hauls vs sonar coverage ---------------------------------
cov <- ship[, .(start = min(t), end = max(t)), by = seg]
p1b <- ggplot() +
  geom_rect(data = cov, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "gold", alpha = 0.6) +
  geom_segment(data = hauls, aes(x = shoot_time, xend = haul_time, y = catch,
                                 colour = dom_species), linewidth = 3) +
  geom_text(data = hauls, aes(x = shoot_time, y = catch, label = haul),
            size = 2.5, vjust = -1) +
  scale_colour_manual(values = species_cols, name = "Dominant\nspecies") +
  scale_x_datetime(date_breaks = "1 day", date_labels = "%d %b") +
  labs(x = NULL, y = "Catch (t)",
       title = "Hauls (bars = tow duration) and sonar recording periods (yellow)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(file.path(cfg$out_dir, "01b_trip_timeline.png"), p1b, width = 11, height = 4, dpi = 200)

# ---- 01c  sonar zoom map -----------------------------------------------------------
hz <- hauls[haul_time >= min(ship$t) & shoot_time <= max(ship$t)]
p1c <- ggplot() +
  geom_path(data = ship, aes(lon, lat, group = seg, colour = phase), linewidth = 0.6) +
  geom_point(data = schools, aes(lon, lat, size = Area, fill = depth),
             shape = 21, alpha = 0.7, stroke = 0.2) +
  geom_segment(data = hz, aes(shootlon, shootlat, xend = haullon, yend = haullat),
               arrow = arrow(length = unit(2, "mm")), linewidth = 1, colour = "red") +
  geom_label(data = hz, aes(haullon, haullat,
                            label = sprintf("H%s %s %.0f t", haul, dom_species, catch)),
             size = 2.8, nudge_y = 0.01) +
  scale_colour_manual(values = c(search = "orange", tow = "red", other = "grey50"),
                      name = "Ship activity") +
  scale_fill_viridis_c(direction = -1, name = "School\ndepth (m)") +
  scale_size_area(max_size = 6, name = "School\narea (m2)") +
  coord_quickmap() +
  labs(x = "Longitude", y = "Latitude",
       title = "Sonar track and detected schools (school mean positions)",
       subtitle = paste("search =", cfg$pre_shoot_h,
                        "h before shoot; tow = shoot to haul; red arrows = logged tow lines"))
ggsave(file.path(cfg$out_dir, "01c_sonar_zoom_map.png"), p1c, width = 9, height = 8, dpi = 200)

message("01: figures written to ", cfg$out_dir)
