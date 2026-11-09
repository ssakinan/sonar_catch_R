################################################################################
# 14_catch_vs_sonar.R
#
# Catch rate against sonar indices per tow, for the report (comment by B. Berges):
#   (a) number of schools in three 30-min windows: before shooting, the first
#       30 min of the tow, and after hauling
#   (b) backscatter of the individual schools seen during each tow, as one
#       boxplot per tow placed at its catch rate
# for the three datasets (SU90 2025, ST90 2026, SU90 2026).
# Uses the tracked schools after QC (school_tracks.rds; run 00b and 03 first for
# each profile). Haul times are taken from the catch log, corrected to the sonar
# clock as in the analysis.
#   source("14_catch_vs_sonar.R")
# Output: output_compare/report_chapter4/Figure_catch_vs_nschools.png,
#         Figure_catch_vs_backscatter.png and catch_vs_sonar_windows.csv
################################################################################

sonar_profile <- "ST90"
source("00_setup_functions.R")   # run from the RStudio project root

win_min   <- 30     # window length (min)
guard_min <- 60     # a window counts only if the sonar recorded schools within this time of it
out <- "output_compare/report_chapter4"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
sets <- list(list(lab = "(a) SU90, October 2025", dir = "output", trip = "2025K", offset = profiles$SU90$catch_offset_h %||% 0),
             list(lab = "(b) ST90, June-July 2026", dir = "output_ST90", trip = "2026E", offset = profiles$ST90$catch_offset_h),
             list(lab = "(c) SU90, June-July 2026", dir = "output_SU90_2026", trip = "2026E", offset = profiles$SU90_2026$catch_offset_h))
wins <- c("30 min before shooting", "First 30 min of tow", "30 min after hauling")

res <- list(); sch <- list()
for (s in sets) {
  k <- readRDS(file.path(s$dir, "school_tracks.rds"))[keep == TRUE, .(sid, t = as.numeric(t), backscatter)]
  cfg$catch_offset_h <- s$offset
  h <- read_catch(cfg$catch_file, trips = s$trip)
  h <- h[shoot_time >= as.POSIXct(min(k$t), origin = "1970-01-01", tz = "UTC") - 3600 &
         shoot_time <= as.POSIXct(max(k$t), origin = "1970-01-01", tz = "UTC") + 3600]
  h[, `:=`(t0 = as.numeric(shoot_time), t1 = as.numeric(haul_time))]
  h[, rate := catch / ((t1 - t0) / 3600)]
  pt <- sort(unique(k$t))                              # times with schools: the only sign the sonar was on
  on_near <- function(a, b) any(pt >= a - guard_min * 60 & pt <= b + guard_min * 60)
  for (i in seq_len(nrow(h))) {
    w <- rbind(c(h$t0[i] - win_min * 60, h$t0[i]), c(h$t0[i], min(h$t0[i] + win_min * 60, h$t1[i])),
               c(h$t1[i], h$t1[i] + win_min * 60))
    for (j in 1:3) {
      ok <- on_near(w[j, 1], w[j, 2])
      res[[length(res) + 1]] <- data.table(set = s$lab, haul = h$haul[i], catch = h$catch[i], rate = h$rate[i],
        window = wins[j], valid = ok, n_schools = if (ok) uniqueN(k[t >= w[j, 1] & t < w[j, 2], sid]) else NA_integer_)
    }
    tow <- k[t >= h$t0[i] & t <= h$t1[i], .(bs = mean(backscatter)), by = sid]
    if (nrow(tow)) sch[[length(sch) + 1]] <- tow[, .(set = s$lab, haul = h$haul[i], rate = h$rate[i], db = 10 * log10(bs))]
  }
}
res <- rbindlist(res); sch <- rbindlist(sch)
res[, `:=`(set = factor(set, sapply(sets, `[[`, "lab")), window = factor(window, wins))]
sch[, set := factor(set, levels(res$set))]
fwrite(res, file.path(out, "catch_vs_sonar_windows.csv"))

rho_a <- res[valid == TRUE, {
  ct <- suppressWarnings(cor.test(rate, n_schools, method = "spearman", exact = FALSE))
  .(n = .N, rho = ct$estimate, p = ct$p.value)
}, by = .(set, window)]
rho_b <- sch[, .(m = median(db)), by = .(set, haul, rate)][, {
  ct <- suppressWarnings(cor.test(rate, m, method = "spearman", exact = FALSE))
  .(n = .N, rho = ct$estimate, p = ct$p.value)
}, by = set]
print(rho_a[, .(set, window, n, rho = round(rho, 2), p = signif(p, 2))]); print(rho_b)

theme_rep <- theme_bw(base_size = 8, base_family = "Verdana") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92"),
        strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 7.5))
lab_rho <- function(d) d[, lab := sprintf("rho = %.2f\nn = %d", rho, n)]

pa <- ggplot(res[valid == TRUE], aes(rate, n_schools)) +
  geom_point(size = 1.3, alpha = 0.7, colour = "#0072B2") +
  geom_text(data = lab_rho(copy(rho_a)), aes(x = Inf, y = Inf, label = lab), hjust = 1.1, vjust = 1.2,
            size = 2.3, family = "Verdana", colour = "grey30") +
  facet_grid(set ~ window, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = "Catch rate (t per hour of towing)", y = "Number of tracked schools") + theme_rep
ggsave(file.path(out, "Figure_catch_vs_nschools.png"), pa, device = ragg::agg_png, width = 15.9 / 2.54, height = 13 / 2.54, dpi = 300)

# tows side by side, ordered by catch rate (labels = catch rate), so that tows with similar rates do not overlap
sch[, tow := sprintf("%d (H%d)", round(rate), haul)]
lv <- unique(sch[order(set, rate, haul), .(set, tow)])[, paste(set, tow)]
sch[, x := factor(paste(set, tow), levels = lv)]
pb <- ggplot(sch, aes(x, db)) +
  geom_boxplot(aes(group = x), outlier.shape = NA, colour = "grey40", fill = "grey90", linewidth = 0.3, width = 0.7) +
  geom_point(position = position_jitter(width = 0.18, height = 0, seed = 1), size = 0.45, alpha = 0.35, colour = "#0072B2") +
  geom_text(data = lab_rho(copy(rho_b)), aes(x = Inf, y = Inf, label = lab), hjust = 1.1, vjust = 1.2,
            size = 2.3, family = "Verdana", colour = "grey30") +
  facet_wrap(~ set, ncol = 1, scales = "free") +
  scale_x_discrete(labels = function(b) sub("^\\([abc]\\) [^,]+, [^ ]+ [0-9]{4} ", "", b)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "Catch rate of the tow, t per hour of towing (haul number); tows ordered by catch rate",
       y = "Mean backscatter of each school\nduring the tow (relative dB)") +
  theme_rep + theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(out, "Figure_catch_vs_backscatter.png"), pb, device = ragg::agg_png, width = 15.9 / 2.54, height = 14 / 2.54, dpi = 300)
message("14: done")
