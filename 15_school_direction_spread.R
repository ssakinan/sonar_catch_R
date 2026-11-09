################################################################################
# 15_school_direction_spread.R
#
# Spread of movement directions of individual daytime schools (report comments
# by B. Berges): displacement over ground (SOG) and through water (STW = SOG
# minus model current), one rose per dataset, with the mean direction and the
# mean resultant length R (0 = no common direction, 1 = all schools the same).
# Uses migration/cmems_school_velocities.csv written by 09_cmems_comparison.R.
#   source("15_school_direction_spread.R")
# Output: output_compare/report_chapter4/Figure_school_direction_spread.png and
#         school_direction_stats.csv
################################################################################

sonar_profile <- "ST90"
source("00_setup_functions.R")   # run from the RStudio project root

moving_only <- TRUE   # only schools whose speed is clearly above its uncertainty
out <- "output_compare/report_chapter4"
sets <- c(`(a) SU90, October 2025` = "output", `(b) ST90, June-July 2026` = "output_ST90",
          `(c) SU90, June-July 2026` = "output_SU90_2026")
kinds <- c(SOG = "Displacement over ground (SOG)", STW = "Displacement through water (STW)")

d <- rbindlist(lapply(names(sets), function(s)
  fread(file.path(sets[[s]], "migration", "cmems_school_velocities.csv"))[, set := s]))
if (moving_only) d <- d[moving == TRUE]
long <- rbind(d[, .(set, kind = kinds[["SOG"]], vx, vy)], d[, .(set, kind = kinds[["STW"]], vx = sx, vy = sy)])
long[, `:=`(dir = bearing_of(vx, vy), set = factor(set, names(sets)), kind = factor(kind, kinds))]

stats <- long[, {
  cs <- circ_summary(dir)
  .(n = .N, mean_dir = cs$mean_dir, Rbar = cs$Rbar, p = cs$rayleigh_p,
    net_dir = bearing_of(mean(vx), mean(vy)), net_speed = sqrt(mean(vx)^2 + mean(vy)^2))
}, by = .(set, kind)]
fwrite(stats, file.path(out, "school_direction_stats.csv"))
print(as.data.frame(stats[, .(set, kind, n, mean_dir = round(mean_dir), Rbar = round(Rbar, 2), p = signif(p, 2),
                             net_dir = round(net_dir), net_kmday = round(net_speed * 86.4, 1))]))

theme_rep <- theme_bw(base_size = 8, base_family = "Verdana") +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 6.8, lineheight = 0.95), axis.text.y = element_text(size = 5.5, colour = "grey50"))
# one panel per dataset and kind; title carries the statistics
stats[, panel := sprintf("%s %s\nn = %d, mean %.0f°, R = %.2f", sub(",.*", "", set), sub(".*\\((.*)\\)", "\\1", kind),
                         n, mean_dir, Rbar)]
stats[, panel := paste0(panel, ifelse(grepl("2025", set), " (2025)", " (2026)"))]
lv <- stats[order(kind, set), panel]
long <- stats[, .(set, kind, panel)][long, on = .(set, kind)]
long[, panel := factor(panel, lv)]
bins <- long[, .(n = .N), by = .(panel, b = floor(dir / 15))][, pct := 100 * n / sum(n), by = panel][, mid := b * 15 + 7.5]
ylim <- max(bins$pct) * 1.05
arr <- stats[, .(panel = factor(panel, lv), mean_dir, len = Rbar * ylim)]   # arrow length = R x outer radius
p <- ggplot(bins, aes(mid, pct)) +
  geom_col(width = 15, fill = "#0072B2", colour = "white", linewidth = 0.2) +
  geom_segment(data = arr, aes(x = mean_dir, xend = mean_dir, y = 0, yend = len), colour = "#D55E00",
               linewidth = 0.8, arrow = arrow(length = unit(1.5, "mm"))) +
  coord_polar(start = 0) +
  scale_x_continuous(limits = c(0, 360), breaks = seq(0, 270, 90), labels = c("0°", "90°", "180°", "270°")) +
  scale_y_continuous(limits = c(0, ylim), breaks = pretty(c(0, ylim), 3), labels = function(b) paste0(b, "%")) +
  facet_wrap(~ panel, ncol = 3) +
  labs(x = NULL, y = NULL) + theme_rep
ggsave(file.path(out, "Figure_school_direction_spread.png"), p, device = ragg::agg_png,
       width = 15.9 / 2.54, height = 12 / 2.54, dpi = 300)
message("15: done")
