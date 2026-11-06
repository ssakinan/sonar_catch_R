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
# Output: <out_dir>/13_interpolated_backscatter_map.png
################################################################################

if (!exists("sonar_profile")) sonar_profile <- "ST90"
source("00_setup_functions.R")   # run from the RStudio project root

cell_km   <- 2      # grid cell for binning
bw_km     <- 4      # Gaussian kernel bandwidth
ping_cap  <- 200    # a cell counts for at most this many pings, so no cell dominates
min_near  <- 100    # pings needed within mask_km of a grid point
mask_km   <- 6      # blank beyond this distance from the nearest cell with data
min_pings <- 5      # pings per cell to use it
squish_q  <- c(0.05, 0.95)
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
grid <- CJ(lon = seq(min(m$lon) - 0.1, max(m$lon) + 0.1, grid_n[["lon"]]),
           lat = seq(min(m$lat) - 0.05, max(m$lat) + 0.05, grid_n[["lat"]]))
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
ship <- unique(m[, .(t, lon = ship_lon, lat = ship_lat)], by = "t")[order(t)]
ship[, seg := cumsum(c(0, diff(as.numeric(t))) > cfg$gap_s)]
ship <- ship[seq(1, .N, 5)]
p <- ggplot() +
  geom_tile(data = grid, aes(lon, lat, fill = db), width = grid_n[["lon"]] * 1.02, height = grid_n[["lat"]] * 1.02) +
  geom_path(data = ship, aes(lon, lat, group = seg), colour = "white", linewidth = 0.15, alpha = 0.6) +
  geom_sf(data = land, fill = "grey85", colour = "grey55", linewidth = 0.2) +
  scale_fill_viridis_c(option = "viridis", limits = lims, oob = scales::squish, name = "School backscatter\nper ping (relative dB)") +
  coord_sf(xlim = range(grid$lon) + c(-0.1, 0.1), ylim = range(grid$lat) + c(-0.05, 0.05), expand = FALSE) +
  labs(x = NULL, y = NULL,
       title = sprintf("%s school backscatter, %s - %s", cfg$label, format(min(m$t), "%d %b"), format(max(m$t), "%d %b %Y")),
       subtitle = sprintf("%d-km cells smoothed with a %d-km Gaussian kernel, weighted by pings; blank = too few data within %d km; white = ship track; colours squished to the 5-95%% range",
                          cell_km, bw_km, mask_km)) +
  theme(plot.subtitle = element_text(size = 8), panel.background = element_rect(fill = "grey97"))
ggsave(file.path(base, "13_interpolated_backscatter_map.png"), p, width = 9, height = 7, dpi = 250)
message("13: done")
