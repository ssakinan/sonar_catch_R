################################################################################
# 11_project_summary.R
#
# Two-page summary of all sonar work on SCH302 (SU90 2025, ST90 and SU90 2026).
# Numbers are read from the outputs of scripts 00b-10, so run those first.
#   source("11_project_summary.R")
# Output: SCH302_sonar_summary.docx in the project root.
################################################################################

library(data.table); library(officer); library(flextable)

rd  <- function(f, ...) fread(f, ...)
m25 <- rd("output/school_movement.csv")
h25 <- rd("output/haul_sonar_indices.csv", colClasses = list(character = "haul_id"))
c25 <- rd("output/movement_circular_stats.csv")
q25 <- rd("output/school_qc.csv")
g25 <- rd("output/migration/migration_estimates.csv")[grepl("\\(main\\)", model)]
w25 <- rd("output/migration/cmems_estimates.csv")
w25_sw <- w25[grepl("^Swimming = ", model)]; w25_c <- w25[grepl("^CMEMS current at schools", model)]
cor25 <- h25[, rbindlist(lapply(c("schools_per_h", "sbi_h_db"), function(v) rbindlist(lapply(c("catch", "rate"), function(y) {
  x <- .SD[is.finite(get(v))]; if (y == "rate") x[, rate := catch / as.numeric(difftime(haul_time, shoot_time, units = "hours"))]
  data.table(rho = suppressWarnings(cor(x[[y]], x[[v]], method = "spearman")))
})))), by = phase]
n25 <- uniqueN(h25$haul_id)
cmp <- rd("output_compare/compare_summary.csv")
cdr <- rd("output_compare/compare_drift.csv")
ccr <- rd("output_compare/compare_catch_correlations.csv")
prs <- rd("output_compare/compare_matched_schools.csv")
mst <- rd("output_compare/compare_matched_stats.csv")
M <- function(i) mst[item == i]$value
mig <- rd("output_ST90/migration/migration_estimates.csv")
loo <- mig[grepl("^Without", model)]

f0 <- function(x) formatC(round(x), format = "d", big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
ang <- function(a, b) abs(((a - b + 180) %% 360) - 180)
unwrap <- function(x, ref) ref + ((x - ref + 180) %% 360) - 180
S  <- function(s, v) cmp[sonar == s][[v]]
D  <- function(s, pat) cdr[sonar == s & grepl(pat, what)]
st_g <- D("ST90", "over ground"); su_g <- D("SU90", "over ground")
st_w <- D("ST90", "Swimming");   su_w <- D("SU90", "Swimming")
st_c <- D("ST90", "CMEMS current"); su_c <- D("SU90", "CMEMS current")
loo_rng <- range(unwrap(loo$dir, st_g$dir)) %% 360
hz25 <- c25[grepl("^haul", group)]
day25 <- c25[group == "sun: day"]
cv <- c(x = cor(prs$vx_st, prs$vx_su), y = cor(prs$vy_st, prs$vy_su))

fp_b <- fp_text(bold = TRUE, font.size = 10); fp_n <- fp_text(font.size = 10)
doc <- read_docx()
doc <- body_set_default_section(doc, prop_section(
  page_size = page_size(orient = "portrait", width = 21 / 2.54, height = 29.7 / 2.54),
  page_margins = page_mar(top = 0.7, bottom = 0.7, left = 0.8, right = 0.8)))
par_fmt <- fp_par(padding.bottom = 4, line_spacing = 1.05)
h  <- function(x) doc <<- body_add_fpar(doc, fpar(ftext(x, fp_text(bold = TRUE, font.size = 12, color = "#1b4f72")),
                                                  fp_p = fp_par(padding.top = 6, padding.bottom = 2, keep_with_next = TRUE)))
p  <- function(...) doc <<- body_add_fpar(doc, fpar(ftext(paste0(...), fp_n), fp_p = par_fmt))
bl <- function(head, ...) doc <<- body_add_fpar(doc, fpar(ftext("•  ", fp_n), ftext(paste0(head, " "), fp_b), ftext(paste0(...), fp_n),
                                                        fp_p = fp_par(padding.bottom = 3, padding.left = 12, line_spacing = 1.05)))

doc <- body_add_fpar(doc, fpar(ftext("Omnidirectional sonar on SCH302: summary of the analyses", fp_text(bold = TRUE, font.size = 16))))
doc <- body_add_fpar(doc, fpar(ftext(sprintf("SU90 23–29 October 2025 (mackerel) and ST90 + SU90 June–July 2026 (herring). Code: sonar_catch_R (github.com/ssakinan/sonar_catch_R). %s.",
                                             format(Sys.Date(), "%d %B %Y")), fp_text(font.size = 9, color = "grey40")),
                               fp_p = fp_par(padding.bottom = 6)))

h("Aim and data")
p("The aim was to see what the vessel's own omnidirectional sonar adds to the interpretation of commercial catches, and what it ",
  "shows about school movement. Data are LSSS Profos school exports (per ping and per school) and the trip catch log (",
  "SCH302_2025K_2026E.csv). Three datasets were analysed with the same R workflow and quality control: SU90, 23–29 October 2025, ",
  sprintf("%d tows east of Shetland (mainly mackerel, %s t); ST90, 28 June – 5 July 2026, %d herring tows east of Shetland and Orkney; and SU90 for the same 2026 period.",
          n25, f0(sum(unique(h25[, .(haul_id, catch)])$catch)), ccr[sonar == "ST90" & phase == "tow"]$n[1]))

h("Method")
bl("Quality control.", sprintf("Pings with Sv < −70 dB or near-horizontal tilt were dropped. Schools were excluded that were faster than 2.5 m/s, moved with the ship, showed tracker jumps, were near-field (< 150 m), weak (< −60 dB) or extreme in size. On the ST90 the excluded set formed a second speed mode at the ship's steaming speed (echoes tied to the vessel). Excluded: %s of %s schools (SU90 2025), %s of %s (ST90) and %s of %s (SU90 2026).",
   f0(sum(q25$excluded)), f0(nrow(q25)), f0(S("ST90", "excluded")), f0(S("ST90", "detected")), f0(S("SU90", "excluded")), f0(S("SU90", "detected"))))
bl("Track smoothing.", "Raw ping-to-ping positions jitter by 3–5 m per axis, which makes raw speeds meaningless. Tracks were cleaned (merging duplicates, outlier filter), smoothed with a spline and a Kalman smoother, and a school counted as moving only when its speed exceeded twice its standard error.")
bl("Movement.", "Directions were summarised with circular statistics. Daytime net drift was estimated as the intercept of a model with a 12.42-h tidal harmonic and ship-velocity terms, with a bootstrap over 1-h blocks. CMEMS North-West Shelf model currents (1.5 km, hourly) were interpolated to each school's position, time and depth, and swimming = ground velocity − current.")

h("Results")
bl("SU90, October 2025 (mackerel).", sprintf(paste(
   "%s schools were tracked after quality control (median smoothed speed %s m/s, %.0f%% moving significantly). All %d tows fished",
   "within the depths the sonar could see, but school numbers and backscatter still did not predict catch (|ρ| ≤ %s).",
   "Overall there was little common direction (R̄ = %s), but daytime schools went south (%.0f°, R̄ = %s). Their net drift,",
   "%s km/day towards %.0f°, closely matches the CMEMS current at the schools (%s km/day towards %.0f°): the mackerel were",
   "largely carried by the current, and the swimming that remains (%s km/day) has no clear direction."),
   f0(nrow(m25)), f2(median(m25$speed_lm)), 100 * mean(m25$moving), n25, f2(max(abs(cor25$rho), na.rm = TRUE)),
   f2(c25[group == "all"]$Rbar), day25$mean_dir, f2(day25$Rbar), f1(g25$km_day), g25$dir, f1(w25_c$km_day), w25_c$dir, f1(w25_sw$km_day)))
bl("ST90, 2026: sonar versus catch.", sprintf("After quality control %s schools were tracked (median speed %s m/s). School numbers and backscatter before or during a tow did not predict catch or catch rate (|ρ| ≤ %s over %d tows). Most tows fished deeper than the reported (beam-centre) school depths.",
   f0(S("ST90", "schools")), f2(S("ST90", "speed_med")), f2(max(abs(ccr[sonar == "ST90"]$rho))), ccr[sonar == "ST90" & phase == "tow"]$n[1]))
bl("ST90, 2026: daytime migration.", sprintf("Daytime schools show a net drift of %s km/day towards %.0f° (WSW; 95%% CI %.0f–%.0f°), stable when any day is left out (%.0f–%.0f°). The semidiurnal tide moves schools back and forth at 0.13–0.19 m/s, more than the drift, and individual schools swim in all directions; the drift is a population-level bias.",
   f1(st_g$km_day), st_g$dir, st_g$dir_lo, st_g$dir_hi, loo_rng[1], loo_rng[2]))
bl("Current or swimming?", sprintf("The CMEMS current at the schools is weak and points south-east (%s km/day towards %.0f°), so it does not explain the drift. Swimming relative to the water is %s km/day towards %.0f° (95%% CI %.0f–%.0f°), robust to the depth of the current used and to the clock offset. Schools follow the model tide only partly (60–80%%).",
   f1(st_c$km_day), st_c$dir, f1(st_w$km_day), st_w$dir, st_w$dir_lo, st_w$dir_hi))
bl("ST90 versus SU90, 2026.", sprintf("The SU90 clock runs 1 h 03 min 38 s behind the ST90 (found from the ship positions in both exports) and was aligned before comparison. Hourly school counts are correlated (ρ = 0.56). The SU90 independently gives the same picture: drift %s km/day towards %.0f° (%.0f° from the ST90), swimming %s km/day towards %.0f°, current %s km/day towards %.0f°. Neither sonar predicts catch (|ρ| ≤ %s).",
   f1(su_g$km_day), su_g$dir, ang(su_g$dir, st_g$dir), f1(su_w$km_day), su_w$dir, f1(su_c$km_day), su_c$dir, f2(max(abs(ccr$rho)))))
bl("Same schools, two sonars.", sprintf(paste(
   "%s schools were tracked by both sonars at the same time and place (%.0f%% of ST90 and %.0f%% of SU90 schools). The overlap",
   "is low mainly because the beams were tilted differently (ST90 %.0f°, SU90 %.0f°) and look at different depths; large, dense",
   "schools are the ones seen by both. A school in both beams must reach the band where they overlap, about %.0f–%.0f m, while",
   "Profos reports %.0f m and %.0f m for the same school. The SU90 reads %.1f dB higher Sv (+%.1f dB per 100 m range) and %.0f%%",
   "of the ST90 school area, so the two are not interchangeable without inter-calibration. Velocities agree without bias;",
   "the paired differences put the velocity error of one sonar at about %.2f m/s per component (%.2f m/s for tracks over",
   "2 min), so single-school directions are uncertain while averages over many schools are reliable."),
   f0(M("pairs")), 100 * M("ST90 share matched"), 100 * M("SU90 share matched"), M("ST90 modal tilt"), M("SU90 modal tilt"),
   M("median depth bound low (m)"), M("median depth bound high (m)"), median(prs$depth_st), median(prs$depth_su),
   M("median Sv SU90 - ST90 (dB)"), M("Sv difference per 100 m range (dB)"), 100 * 10^M("median log10 area ratio SU90/ST90"),
   M("velocity error per sonar, robust (m/s)"), prs[pmin(duration_s_st, duration_s_su) > 120, mad(c(vx_su - vx_st, vy_su - vy_st)) / sqrt(2)]))

tab <- data.table(
  ` ` = c("Period", "Tracked schools (after QC)", "Median smoothed speed (m/s)", "Daytime drift over ground", "Swimming relative to CMEMS", "Max |ρ| school count or backscatter vs catch"),
  `SU90 2025` = c("23–29 Oct 2025", f0(nrow(m25)), f2(median(m25$speed_lm)), sprintf("%s km/d, %.0f°", f1(g25$km_day), g25$dir),
                  sprintf("%s km/d, unclear", f1(w25_sw$km_day)), f2(max(abs(cor25$rho), na.rm = TRUE))),
  `ST90 2026` = c("28 Jun – 5 Jul 2026", f0(S("ST90", "schools")), f2(S("ST90", "speed_med")),
                  sprintf("%s km/d, %.0f°", f1(st_g$km_day), st_g$dir), sprintf("%s km/d, %.0f°", f1(st_w$km_day), st_w$dir),
                  f2(max(abs(ccr[sonar == "ST90"]$rho)))),
  `SU90 2026` = c("same", f0(S("SU90", "schools")), f2(S("SU90", "speed_med")),
                  sprintf("%s km/d, %.0f°", f1(su_g$km_day), su_g$dir), sprintf("%s km/d, %.0f°", f1(su_w$km_day), su_w$dir),
                  f2(max(abs(ccr[sonar == "SU90"]$rho)))))
ft <- flextable(tab); ft <- theme_booktabs(ft); ft <- fontsize(ft, size = 9, part = "all"); ft <- bold(ft, part = "header")
ft <- bold(ft, j = 1); ft <- width(ft, width = c(2.3, 1.4, 1.6, 1.6)); ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all")
doc <- body_add_fpar(doc, fpar(ftext("Key numbers", fp_b), fp_p = fp_par(padding.top = 6, keep_with_next = TRUE)))
doc <- body_add_flextable(doc, ft, align = "left")

h("Conclusions")
bl("Catch.", "In none of the three datasets did the sonar indices explain catch, even in 2025 when all tows fished within the depths the sonar could see. The forward ±90° sector and ~600 m range limit what is seen, catch per haul is capped by the skipper and the vessel, and counts within 5 km describe the area more than the fish that entered the net.")
bl("Movement.", "Smoothing is essential: smoothed tracks give realistic speeds (medians 0.6–0.75 m/s). In 2026 two independent sonars show the daytime herring moving west to south-west at about 10 km/day relative to the water against a weak south-easterly current, consistent with a movement towards the Orkney–Shetland spawning grounds; this is suggestive of migration, not proof. In 2025 the southward drift of the mackerel is explained by the current, which shows the value of the current-model check.")
bl("Data quality.", "Profos depth is the beam-centre depth (range × sin tilt): the same school is reported at ~29 m by the ST90 and ~63 m by the SU90. Clocks differ: the ST90 is about 1 h ahead of the catch log and the SU90.")

h("Next steps")
bl("Export settings.", "Re-export with the full sector and longer range, and the per-ping effort file (ppe), which the current exports lack. Set the two sonars' tilts on purpose (same tilt to count the same schools, different tilts to cover two layers) and inter-calibrate them.")
bl("Clocks and depth.", "Fix all clocks to UTC at source. Use the vertical fan or echosounder for true school depth to compare with headline depth.")
bl("More data.", "Repeat the migration analysis on other trips and years to see whether the westward movement recurs at this time of year, and validate the swimming estimate against current measurements (e.g. ADCP).")
bl("Outputs.", "Reports: output/SCH302_sonar_catch_report.docx, output_ST90/SCH302_ST90_sonar_catch_report.docx, output_ST90/migration/ (migration and CMEMS documents), output_compare/SCH302_ST90_vs_SU90.docx; animations in output*/animation.")

print(doc, target = "SCH302_sonar_summary.docx")
message("11: done -> ", normalizePath("SCH302_sonar_summary.docx"))
