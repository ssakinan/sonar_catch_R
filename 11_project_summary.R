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
cmp <- rd("output_compare/compare_summary.csv")
cdr <- rd("output_compare/compare_drift.csv")
ccr <- rd("output_compare/compare_catch_correlations.csv")
prs <- rd("output_compare/compare_matched_schools.csv")
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
doc <- body_add_fpar(doc, fpar(ftext(sprintf("SU90 October 2025 (mackerel) and ST90 + SU90 June–July 2026 (herring). Code: sonar_catch_R (github.com/ssakinan/sonar_catch_R). %s.",
                                             format(Sys.Date(), "%d %B %Y")), fp_text(font.size = 9, color = "grey40")),
                               fp_p = fp_par(padding.bottom = 6)))

h("Aim and data")
p("The aim was to see what the vessel's own omnidirectional sonar adds to the interpretation of commercial catches, and what it ",
  "shows about school movement. Data are LSSS Profos school exports (per ping and per school) and the trip catch log (",
  "SCH302_2025K_2026E.csv). Three datasets were analysed with the same R workflow: SU90, 23 October 2025, around two mackerel tows ",
  sprintf("(%s and %s t); ST90, 28 June – 5 July 2026, %d herring tows east of Shetland and Orkney; and SU90 for the same 2026 period.",
          f0(h25[haul_id == "2025K-7"]$catch[1]), f0(h25[haul_id == "2025K-8"]$catch[1]), ccr[sonar == "ST90" & phase == "tow"]$n[1]))

h("Method")
bl("Quality control.", sprintf("Pings with Sv < −70 dB or near-horizontal tilt were dropped. Schools were excluded that were faster than 2.5 m/s, moved with the ship, showed tracker jumps, were near-field (< 150 m), weak (< −60 dB) or extreme in size. On the ST90 the excluded set formed a second speed mode at the ship's steaming speed (echoes tied to the vessel). Excluded: %s of %s ST90 schools and %s of %s SU90 schools.",
   f0(S("ST90", "excluded")), f0(S("ST90", "detected")), f0(S("SU90", "excluded")), f0(S("SU90", "detected"))))
bl("Track smoothing.", "Raw ping-to-ping positions jitter by 3–5 m per axis, which makes raw speeds meaningless. Tracks were cleaned (merging duplicates, outlier filter), smoothed with a spline and a Kalman smoother, and a school counted as moving only when its speed exceeded twice its standard error.")
bl("Movement.", "Directions were summarised with circular statistics. Daytime net drift was estimated as the intercept of a model with a 12.42-h tidal harmonic and ship-velocity terms, with a bootstrap over 1-h blocks. CMEMS North-West Shelf model currents (1.5 km, hourly) were interpolated to each school's position, time and depth, and swimming = ground velocity − current.")

h("Results")
bl("SU90, October 2025.", sprintf("%s of 154 schools were tracked; positioning noise was about 4 m and the median smoothed speed %s m/s, with %.0f%% of schools moving significantly. Schools around each tow moved coherently (haul 7: %.0f°, R̄ = %s; haul 8: %.0f°, R̄ = %s), but there was no common direction overall (R̄ = %s). Both tows fished with the headline at %d–%d m, mostly below the (beam-centre) depths of the schools the sonar reported.",
   f0(nrow(m25)), f2(median(m25$speed_lm)), 100 * mean(m25$moving), hz25[group == "haul: 2025K-7"]$mean_dir, f2(hz25[group == "haul: 2025K-7"]$Rbar),
   hz25[group == "haul: 2025K-8"]$mean_dir, f2(hz25[group == "haul: 2025K-8"]$Rbar), f2(c25[group == "all"]$Rbar),
   min(h25$headline_depth), max(h25$headline_depth)))
bl("ST90, 2026: sonar versus catch.", sprintf("After quality control %s schools were tracked (median speed %s m/s). School numbers and backscatter before or during a tow did not predict catch or catch rate (|ρ| ≤ %s over %d tows). Most tows fished deeper than the reported (beam-centre) school depths.",
   f0(S("ST90", "schools")), f2(S("ST90", "speed_med")), f2(max(abs(ccr[sonar == "ST90"]$rho))), ccr[sonar == "ST90" & phase == "tow"]$n[1]))
bl("ST90, 2026: daytime migration.", sprintf("Daytime schools show a net drift of %s km/day towards %.0f° (WSW; 95%% CI %.0f–%.0f°), stable when any day is left out (%.0f–%.0f°). The semidiurnal tide moves schools back and forth at 0.13–0.19 m/s, more than the drift, and individual schools swim in all directions; the drift is a population-level bias.",
   f1(st_g$km_day), st_g$dir, st_g$dir_lo, st_g$dir_hi, loo_rng[1], loo_rng[2]))
bl("Current or swimming?", sprintf("The CMEMS current at the schools is weak and points south-east (%s km/day towards %.0f°), so it does not explain the drift. Swimming relative to the water is %s km/day towards %.0f° (95%% CI %.0f–%.0f°), robust to the depth of the current used and to the clock offset. Schools follow the model tide only partly (60–80%%).",
   f1(st_c$km_day), st_c$dir, f1(st_w$km_day), st_w$dir, st_w$dir_lo, st_w$dir_hi))
bl("ST90 versus SU90, 2026.", sprintf("The SU90 clock runs 1 h 03 min 38 s behind the ST90 (found from the ship positions in both exports) and was aligned before comparison. Hourly school counts are correlated (ρ = 0.56) and %s schools were tracked by both sonars; their velocities agree moderately (r = %s east, %s north). The SU90 independently gives the same picture: drift %s km/day towards %.0f° (%.0f° from the ST90), swimming %s km/day towards %.0f°, current %s km/day towards %.0f°. Neither sonar predicts catch (|ρ| ≤ %s).",
   f0(nrow(prs)), f2(cv["x"]), f2(cv["y"]), f1(su_g$km_day), su_g$dir, ang(su_g$dir, st_g$dir), f1(su_w$km_day), su_w$dir,
   f1(su_c$km_day), su_c$dir, f2(max(abs(ccr$rho)))))

tab <- data.table(
  ` ` = c("Period", "Tracked schools (after QC)", "Median smoothed speed (m/s)", "Daytime drift over ground", "Swimming relative to CMEMS", "Max |ρ| sonar index vs catch"),
  `SU90 2025` = c("23 Oct 2025", f0(nrow(m25)), f2(median(m25$speed_lm)), "–", "–", "2 tows only"),
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
bl("Catch.", "With the current exports the sonar indices do not explain catch: tows fish deeper than the reported schools, the forward ±90° sector and ~550–600 m range limit what is seen, and catch also depends on where and how the skipper chose to fish.")
bl("Movement.", "Smoothing is essential: smoothed tracks give realistic speeds (median ~0.6 m/s) and coherent local movement. Two independent sonars show the daytime herring moving west to south-west at about 10 km/day relative to the water, while the model current is weak and south-easterly, consistent with a movement towards the Orkney–Shetland spawning grounds; this is suggestive of migration, not proof.")
bl("Data quality.", "Profos depth is the beam-centre depth (range × sin tilt): the same school is reported at ~29 m by the ST90 and ~63 m by the SU90. Clocks differ: the ST90 is about 1 h ahead of the catch log and the SU90.")

h("Next steps")
bl("Export settings.", "Re-export with the full sector and longer range, and the per-ping effort file (ppe), which the 2026 exports lack.")
bl("Clocks and depth.", "Fix all clocks to UTC at source. Use the vertical fan or echosounder for true school depth to compare with headline depth.")
bl("More data.", "Repeat the migration analysis on other trips and years to see whether the westward movement recurs at this time of year, and validate the swimming estimate against current measurements (e.g. ADCP).")
bl("Outputs.", "Reports: output/SCH302_sonar_catch_report.docx, output_ST90/SCH302_ST90_sonar_catch_report.docx, output_ST90/migration/ (migration and CMEMS documents), output_compare/SCH302_ST90_vs_SU90.docx; animations in output*/animation.")

print(doc, target = "SCH302_sonar_summary.docx")
message("11: done -> ", normalizePath("SCH302_sonar_summary.docx"))
