################################################################################
# 13_interpolated_map.R
#
# Interpolated map of school backscatter for one sonar profile.
#   sonar_profile <- "ST90"; source("13_interpolated_map.R")
#
# Uses the tracked schools after QC (school_tracks.rds, so run 00b and 03 first).
# Pings are binned in grid cells; the value of a cell is the summed school
# backscatter (Sv x area, linear) per ping with schools, in relative dB, so that
# cells the vessel visited often do not stand out just for that. Cell values are
# interpolated with a Gaussian kernel weighted by pings (a plain inverse distance
# weighting gives a bullseye around every cell) and masked where data are sparse.
# The colour scale is viridis; values outside the 5-95% range are squished to the
# ends of the scale (oob = squish) to bring out the contrast in the middle.
#
# Output: <out_dir>/13_interpolated_backscatter_map.png, 13_track_bubble_map.png
################################################################################

if (!exists("sonar_profile")) sonar_profile <- "ST90"
source("00_setup_functions.R")   # run from the RStudio project root

cell_km   <- 2      # grid cell for binning
bw_km     <- 4      # Gaussian kernel bandwidth
ping_cap  <- 200    # a cell counts for at most this many pings, so no cell dominates
min_near  <- 30     # pings needed within mask_km of a grid point
mask_km   <- 9      # blank beyond this distance from the nearest cell with data
min_pings <- 5      # pings per cell to use it
squish_q  <- c(0.05, 0.95)
map_lon   <- c(-4.2, 1.6)   # map borders: north of Scotland, Orkney and Shetland
map_lat   <- c(58.2, 60.95)
grid_n    <- c(lon = 0.015, lat = 0.008)   # output grid (deg), about 0.9 km

base <- sub("/report_figs$", "", cfg$out_dir)
m <- readRDS(file.path(base, "school_tracks.rds"))[keep == TRUE]   # tracked schools after QC, one row per ping
lon0 <- mean(m$lon); lat0 <- mean(m$lat)
xy <- to_xy(m$lon, m$lat, lon0, lat0)
m[, `:=`(x = xy[, 1] / 1000, y = xy[, 2] / 1000)]
m[, `:=`(gx = floor(x / cell_km), gy = floor(y / cell_km))]

cells <- m[, .(pings = uniqueN(t), sb = sum(backscatter, na.rm = TRUE)), by = .(gx, gy)][pings >= min_pings]
cells[, `:=`(x = (gx + 0.5) * cell_km, y = (gy + 0.5) * cell_km, db = 10 * log10(sb / pings))]
cells <- cells[is.finite(db)]

# ---- ping-weighted Gaussian kernel on a lon/lat grid --------------------------------------------------
grid <- CJ(lon = seq(min(m$lon) - 0.2, max(m$lon) + 0.2, grid_n[["lon"]]),
           lat = seq(min(m$lat) - 0.1, max(m$lat) + 0.1, grid_n[["lat"]]))
gxy <- to_xy(grid$lon, grid$lat, lon0, lat0) / 1000
D <- sqrt(outer(gxy[, 1], cells$x, "-")^2 + outer(gxy[, 2], cells$y, "-")^2)
W <- exp(-D^2 / (2 * bw_km^2)) * rep(pmin(cells$pings, ping_cap), each = nrow(D))
grid[, `:=`(db = as.vector(W %*% cells$db) / rowSums(W), near = apply(D, 1, min),
            support = as.vector((D <= mask_km) %*% cells$pings))]
grid <- grid[near <= mask_km & support >= min_near & is.finite(db)]

lims <- quantile(grid$db, squish_q)
message(sprintf("%d cells, %d pings with schools; colour range %.1f to %.1f dB (values %.1f to %.1f squished)",
                nrow(cells), sum(cells$pings), lims[1], lims[2], min(grid$db), max(grid$db)))

# ---- map -------------------------------------------------------------------------------------------
land <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
theme_map <- theme_bw(base_size = 8, base_family = "Verdana") +
  theme(panel.background = element_rect(fill = "grey97"), panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "right")
save_map <- function(p, file, report_file) {      # no titles: the captions go below the figures in the report
  ggsave(file.path(base, file), p, device = ragg::agg_png, width = 15.9 / 2.54, height = 12.6 / 2.54, dpi = 300)
  rep_dir <- "output_compare/report_chapter4"
  if (sonar_profile == "ST90" && dir.exists(rep_dir)) file.copy(file.path(base, file), file.path(rep_dir, report_file), overwrite = TRUE)
}
# ship track: all ship positions stored with the school pings of this sonar, plus the second sonar of the
# same trip (no ppe export, so the track has gaps without schools; the line joins them with straight sections)
track_files <- unique(c(file.path(base, "school_tracks.rds"), if (cfg$trip == "2026E") c("output_ST90/school_tracks.rds", "output_SU90_2026/school_tracks.rds")))
ship <- rbindlist(lapply(track_files, function(f) readRDS(f)[, .(t, lon = ship_lon, lat = ship_lat)]))
ship <- unique(ship, by = "t")[order(t)][seq(1, .N, 5)]
p <- ggplot() +
  geom_tile(data = grid, aes(lon, lat, fill = db), width = grid_n[["lon"]] * 1.02, height = grid_n[["lat"]] * 1.02) +
  geom_path(data = ship, aes(lon, lat), colour = "white", linewidth = 0.15, alpha = 0.7) +
  geom_sf(data = land, fill = "grey85", colour = "grey55", linewidth = 0.2) +
  scale_fill_viridis_c(option = "viridis", limits = lims, oob = scales::squish, name = "School backscatter\nper ping (relative dB)") +
  labs(x = NULL, y = NULL) + theme_map
# ---- bubble plot: same grid resolution, along the ship's track, not interpolated ---------------------
sxy <- to_xy(m$ship_lon, m$ship_lat, lon0, lat0) / 1000
m[, `:=`(sgx = floor(sxy[, 1] / cell_km), sgy = floor(sxy[, 2] / cell_km))]
bub <- m[, .(pings = uniqueN(t), sb = sum(backscatter, na.rm = TRUE), lon = mean(ship_lon), lat = mean(ship_lat)),
         by = .(sgx, sgy)][pings >= min_pings]
bub[, db := 10 * log10(sb / pings)]
bub <- bub[is.finite(db)][order(db)]          # strongest cells drawn last
bub[, db_s := scales::squish(db, lims)]       # same limits as the map, for colour and size
# same extent for both plots
ext <- coord_sf(xlim = map_lon, ylim = map_lat, expand = FALSE)
save_map(p + ext, "13_interpolated_backscatter_map.png", "Figure_45_interpolated_map.png")
db_breaks <- pretty(lims, 4); db_breaks <- db_breaks[db_breaks >= lims[1] & db_breaks <= lims[2]]
pb <- ggplot() +
  geom_path(data = ship, aes(lon, lat), colour = "grey30", linewidth = 0.25) +
  geom_sf(data = land, fill = "grey85", colour = "grey55", linewidth = 0.2) +
  geom_point(data = bub, aes(lon, lat, size = db_s, fill = db_s), shape = 21, colour = "grey25", stroke = 0.15, alpha = 0.8) +
  scale_fill_viridis_c(option = "viridis", limits = lims, breaks = db_breaks,
                       name = "School
backscatter
per ping
(relative dB)", guide = "legend") +
  scale_radius(range = c(0.15, 3.6), limits = lims, breaks = db_breaks,
                        name = "School
backscatter
per ping
(relative dB)") +
  ext +
  labs(x = NULL, y = NULL) + theme_map
save_map(pb, "13_track_bubble_map.png", "Figure_46_track_bubbles.png")
message(sprintf("13: done; %d bubbles", nrow(bub)))
