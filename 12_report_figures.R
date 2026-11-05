################################################################################
# 12_report_figures.R
#
# Figures for chapter 4 of the pelAcoustics 2026 report ("Use of omnidirectional
# SONAR data"). Sized for the report's text width (15.9 cm) in its font (Verdana);
# no titles, as the report captions go below each figure.
# Run the analysis scripts for all three profiles first (00b, 03, 08, 09, 10, 10b).
#   source("12_report_figures.R")
# Output: output_compare/report_chapter4/
################################################################################

sonar_profile <- "ST90"
source("00_setup_functions.R")   # run from the RStudio project root

out <- "output_compare/report_chapter4"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
w_in <- 15.9 / 2.54
save_fig <- function(p, file, height_cm)
  ggsave(file.path(out, file), p, device = ragg::agg_png, width = w_in, height = height_cm / 2.54, dpi = 300)
theme_rep <- theme_bw(base_size = 8, base_family = "Verdana") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92"),
        strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 8),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
cols_d <- c("Schools over ground" = "grey15", "Model current" = "#0072B2", "Swimming" = "#D55E00")
cols_s <- c(ST90 = "#0072B2", SU90 = "#D55E00")
sets <- c(`(a) SU90, October 2025` = "output", `(b) ST90, June-July 2026` = "output_ST90", `(c) SU90, June-July 2026` = "output_SU90_2026")

# ---- Figure 41: where the data come from ---------------------------------------------------------
land <- sf::st_as_sf(rnaturalearth::ne_countries(scale = "large", returnclass = "sf"))
pos <- rbindlist(lapply(names(sets), function(s) fread(file.path(sets[[s]], "school_movement.csv"))[, .(lon, lat)][, set := s]))
pos[, set := factor(set, names(sets))]
hl <- rbindlist(lapply(c("2025K", "2026E"), function(tr) read_catch(cfg$catch_file, trips = tr)[, .(trip, shoot_time, lon = shootlon, lat = shootlat)]))
ind <- rbindlist(lapply(names(sets), function(s) fread(file.path(sets[[s]], "haul_sonar_indices.csv"), colClasses = list(character = "haul_id"))[, .(haul_id, haul)][, set := s]))
all_h <- rbindlist(lapply(c("2025K", "2026E"), function(tr) read_catch(cfg$catch_file, trips = tr)[, .(haul_id, lon = shootlon, lat = shootlat)]))
hs <- unique(ind[, .(haul_id, set)])[all_h, on = "haul_id", nomatch = NULL]
hs[, set := factor(set, names(sets))]
p41 <- ggplot() +
  geom_sf(data = land, fill = "grey85", colour = "grey60", linewidth = 0.2) +
  geom_point(data = pos, aes(lon, lat, colour = "Tracked school"), size = 0.25, alpha = 0.5) +
  geom_point(data = hs, aes(lon, lat, shape = "Start of tow"), size = 1.4, colour = "black", fill = "white", stroke = 0.4) +
  facet_wrap(~ set) +
  coord_sf(xlim = c(-2.6, 2.1), ylim = c(58.6, 61), expand = FALSE) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1)) +
  scale_colour_manual(values = c("Tracked school" = "#0072B2"), name = NULL, guide = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  scale_shape_manual(values = c("Start of tow" = 21), name = NULL) +
  labs(x = NULL, y = NULL) + theme_rep
save_fig(p41, "Figure_41_data_map.png", 7.2)

# ---- Figure 42: why tracks are smoothed ----------------------------------------------------------
tr <- readRDS("output_ST90/school_tracks.rds")
mv <- fread("output_ST90/school_movement.csv")
ex <- mv[moving == TRUE & duration_s > 180 & net_m > 100 & net_m < 170 & range_m > 250][order(-n_pings)][1:3, sid]
tx <- tr[sid %in% ex]
tx[, lab := paste0("School ", match(sid, ex))]
sm <- tx[keep == TRUE & !is.na(x_kf)]
span <- tx[, .(cx = (max(x) + min(x)) / 2, cy = (max(y) + min(y)) / 2, h = max(diff(range(x)), diff(range(y))) / 2 * 1.08), by = lab]
span[, h := max(h)]                                  # same scale in every panel
ext <- rbind(span[, .(lab, x = cx - h, y = cy - h)], span[, .(lab, x = cx + h, y = cy + h)])
p42 <- ggplot() +
  geom_blank(data = ext, aes(x, y)) +
  geom_point(data = tx, aes(x, y, colour = "Raw positions from each ping"), size = 0.5) +
  geom_path(data = sm, aes(x_kf, y_kf, colour = "Smoothed track"), linewidth = 0.8) +
  facet_wrap(~ lab, scales = "free") +
  scale_colour_manual(values = c("Raw positions from each ping" = "grey60", "Smoothed track" = "#D55E00"), name = NULL) +
  labs(x = "East (m)", y = "North (m)") + theme_rep + theme(aspect.ratio = 1)
save_fig(p42, "Figure_42_smoothing.png", 6.8)

# ---- Figure 43: drift, current and swimming --------------------------------------------------------
est <- rbindlist(lapply(names(sets), function(s) {
  x <- fread(file.path(sets[[s]], "migration/cmems_estimates.csv"))
  rbind(x[grepl("^Schools", model)][, what := "Schools over ground"], x[grepl("^CMEMS current at schools", model)][, what := "Model current"],
        x[grepl("^Swimming = ", model)][, what := "Swimming"])[, set := s]
}))
est[, `:=`(set = factor(set, names(sets)), what = factor(what, names(cols_d)), e = vx * 86.4, n = vy * 86.4)]
lim <- 13.5
rings <- CJ(r = c(5, 10), a = seq(0, 2 * pi, length = 120))[, `:=`(x = r * sin(a), y = r * cos(a))]
p43 <- ggplot(est) +
  geom_path(data = rings, aes(x, y, group = r), colour = "grey85", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_text(data = data.table(y = c(5.7, 10.7), l = c("5", "10 km/day")), aes(0.4, y, label = l), hjust = 0, size = 2.2,
            family = "Verdana", colour = "grey40") +
  geom_segment(aes(0, 0, xend = e, yend = n, colour = what), arrow = arrow(length = unit(1.8, "mm"), type = "closed"), linewidth = 0.7) +
  facet_wrap(~ set) + coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  scale_colour_manual(values = cols_d, name = NULL) +
  labs(x = "East (km/day)", y = "North (km/day)") + theme_rep
save_fig(p43, "Figure_43_drift_current_swimming.png", 7.6)

# ---- Figure 44: same school, two beams -------------------------------------------------------------
pr <- fread("output_compare/compare_matched_schools.csv")
hb <- cfg$vert_beam_half; tilts <- c(ST90 = -4, SU90 = -9)
rr <- seq(0, 600, 10)
beam <- rbindlist(lapply(names(tilts), function(s) data.table(sonar = s, r = rr, centre = rr * sin(-tilts[[s]] * pi / 180),
  lo = pmax(rr * sin(-(tilts[[s]] + hb) * pi / 180), 0), hi = rr * sin(-(tilts[[s]] - hb) * pi / 180))))
pg <- rbind(pr[, .(id = .I, r = range_m_st, depth = depth_st, sonar = "ST90")], pr[, .(id = .I, r = range_m_su, depth = depth_su, sonar = "SU90")])
p44 <- ggplot() +
  geom_ribbon(data = beam, aes(r, ymin = lo, ymax = hi, fill = sonar), alpha = 0.15) +
  geom_line(data = pg, aes(r, depth, group = id), colour = "grey75", linewidth = 0.12) +
  geom_point(data = pg, aes(r, depth, colour = sonar, shape = sonar), size = 0.7) +
  geom_line(data = beam, aes(r, centre, colour = sonar), linetype = 2, linewidth = 0.4) +
  scale_y_reverse() + scale_colour_manual(values = cols_s, name = NULL) + scale_fill_manual(values = cols_s, name = NULL) +
  scale_shape_manual(values = c(ST90 = 16, SU90 = 17), name = NULL) +
  labs(x = "Range from the vessel (m)", y = "Depth (m)") + theme_rep
save_fig(p44, "Figure_44_two_beams.png", 7)
message("12: figures written to ", normalizePath(out))
