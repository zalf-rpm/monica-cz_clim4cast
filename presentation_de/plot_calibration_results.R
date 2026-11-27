#### Figures of the German calibration and cross-validation (presentation) ####
#
# Run from the repository root (monica-cz_clim4cast), e.g.
#   Rscript presentation_de/plot_calibration_results.R
#
# Input (repository root):
#   national_cv_WW.csv          cross-validated national means (evaluate-cv_de.py, national folds)
#   regional_cv_WW.csv          cross-validated regional yields (evaluate-cv_de.py, regional grid)
#   regional_cv_summary_WW.csv  CV RMSE per region
#   regional_best_WW.csv        calibrated FieldConditionModifier per region (collect-calibration-results_de.py)
#   calibration_WW_de.csv, Weights_WW_de.csv
# Simulated values are the leave-one-year-out cross-validated predictions throughout.

require("sf")
require("ggplot2")
require("RColorBrewer")


#### Settings ####

myCrop       <- "WW"
myCropName   <- "Winter wheat"
myDryYears   <- c(2003, 2018)
firstYear    <- 1995
minObsYears  <- 2        # regions with fewer observed years use the national FieldConditionModifier
outFolder    <- "presentation_de/figures"

dir.create(outFolder, recursive = TRUE, showWarnings = FALSE)


#### Load data ####

national <- read.csv(paste0("national_cv_", myCrop, ".csv"))
regional <- read.csv(paste0("regional_cv_", myCrop, ".csv"))
cvSummary <- read.csv(paste0("regional_cv_summary_", myCrop, ".csv"))
best <- read.csv(paste0("regional_best_", myCrop, ".csv"))
weights <- read.csv(paste0("Weights_", myCrop, "_de.csv"))

# a region-year predicted twice (e.g. old and re-run grid evaluated together): keep the last one
dup <- duplicated(regional[, c("region_id", "year")], fromLast = TRUE)
if (any(dup)) {
  warning(sum(dup), " duplicated region-years in regional_cv_", myCrop, ".csv; the last one is used")
  regional <- regional[!dup, ]
}

# yields in t/ha
national$obs <- national$obs_kg_ha / 1000
national$sim <- national$sim_kg_ha / 1000
regional$obs <- regional$obs_kg_ha / 1000
regional$sim <- regional$sim_kg_ha / 1000

# number of observed years per region in the calibration period
obsFile <- read.csv(paste0("calibration_", myCrop, "_de.csv"), check.names = FALSE)
yearCols <- names(obsFile)[suppressWarnings(as.numeric(names(obsFile))) >= firstYear]
yearCols <- yearCols[!is.na(yearCols)]
noOfObs <- data.frame(region_id = obsFile$region_id,
                      no_of_obs = rowSums(obsFile[, yearCols] >= 0))

# district polygons (AGS) linked to the region ids via the weights file
counties <- st_read("data/de/DE_counties_25832.shp", quiet = TRUE)
weights$AGS <- sprintf("%05d", as.integer(weights$NUTS_NAME))
counties <- merge(counties, weights[, c("AGS", "region_id", "NUTS_ID")], by = "AGS", all.x = TRUE)

states  <- st_geometry(st_read("R_scripts_plots/NUTS250_N1.shp", quiet = TRUE))
Germany <- st_geometry(st_read("R_scripts_plots/DE_genau.shp", quiet = TRUE))
Germany <- st_transform(Germany, crs = st_crs(states))


#### Helper functions ####

stats_label <- function(obs, sim, r2 = TRUE) {
  d <- sim - obs
  lines <- c(sprintf("RMSE = %.2f t/ha", sqrt(mean(d^2))),
             sprintf("Bias = %+.2f t/ha", mean(d)))
  if (r2) lines <- c(sprintf("R² = %.2f", summary(lm(sim ~ obs))$r.squared), lines)
  paste(lines, collapse = "\n")
}

presentation_theme <- theme_bw() +
  theme(
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 14),
    legend.position   = c(0.15, 0.87),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title  = element_blank(),
    legend.text   = element_text(size = 14),
    legend.key.size = unit(1.5, "lines")
  )

save_plot <- function(p, name, width = 8, height = 6) {
  ggsave(file.path(outFolder, paste0(name, ".jpeg")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(outFolder, paste0(name, ".pdf")), p, width = width, height = height)
}

scatter_plot <- function(df, title, pointColor, pointAlpha, pointSize) {
  lims <- range(c(df$obs, df$sim), na.rm = TRUE)
  label <- stats_label(df$obs, df$sim, r2 = FALSE)
  ggplot(df, aes(x = obs, y = sim)) +
    geom_point(color = pointColor, alpha = pointAlpha, size = pointSize) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                color = "black", linewidth = 0.6, linetype = "solid") +
    geom_abline(slope = 1, intercept = 0, linewidth = 1, linetype = "dashed") +
    annotate("text", x = Inf, y = -Inf, label = label,
             hjust = 1.05, vjust = -0.3, size = 5, lineheight = 1.1) +
    scale_x_continuous(limits = lims) +
    scale_y_continuous(limits = lims) +
    coord_fixed(ratio = 1) +
    labs(title = title, subtitle = "Leave-one-year-out cross-validation, 1995-2024",
         x = "Observed yield (t/ha)", y = "Simulated yield (t/ha)") +
    presentation_theme
}

# colour of every value: 101 colours between zlim, values outside clipped, NA grey
value_colors <- function(values, zlim, palette, naColor = "grey85") {
  cols <- colorRampPalette(palette)(101)
  v <- pmin(pmax(values, min(zlim)), max(zlim))
  idx <- round((v - min(zlim)) / diff(zlim) * 100) + 1
  out <- cols[idx]
  out[is.na(values)] <- naColor
  out
}

# horizontal colour bar below a map panel; barWidth: share of the panel width, barHeight: share of the
# legend row height (cf. legend.shrink / legend.width in plot_absolute_maps_Claas.R)
color_bar <- function(zlim, palette, breaks, title, barWidth = 0.5, barHeight = 0.22) {
  cols <- colorRampPalette(palette)(101)
  par(plt = c(0.5 - barWidth / 2, 0.5 + barWidth / 2, 0.55 - barHeight, 0.55))
  image(x = seq(min(zlim), max(zlim), length.out = 101), y = 1, z = matrix(seq_len(101), ncol = 1),
        col = cols, axes = FALSE, xlab = "", ylab = "")
  axis(1, at = breaks, labels = breaks, cex.axis = 1.3)
  box(col = "grey30")
  mtext(title, side = 3, line = 0.4, cex = 1.1)
}

map_panel <- function(values, title, zlim, palette) {
  par(mar = c(0.1, 0.1, 3, 0.1))
  plot(st_geometry(counties), col = value_colors(values, zlim, palette),
       border = "grey70", lwd = 0.2, main = "")
  plot(states, col = NA, border = "grey50", lwd = 0.5, add = TRUE)
  plot(Germany, col = NA, border = "grey30", lwd = 1, add = TRUE)
  title(main = title, cex.main = 2)
}

# panels: list of list(values, title, zlim, palette, breaks, legend)
plot_maps <- function(panels, file, mainTitle, width = 17, height = 9) {
  png(filename = file.path(outFolder, file), width = width, height = height, units = "in", res = 300)
  n <- length(panels)
  layout(rbind(seq_len(n), n + seq_len(n)), heights = c(5, 1))
  par(oma = c(0.5, 1, 4, 1), xaxs = "i", yaxs = "i")
  for (p in panels) map_panel(p$values, p$title, p$zlim, p$palette)
  for (p in panels) color_bar(p$zlim, p$palette, p$breaks, p$legend)
  mtext(mainTitle, side = 3, line = 1, outer = TRUE, cex = 2, font = 2, adj = 0)
  dev.off()
}

yieldPalette <- brewer.pal(11, "RdYlGn")
diffPalette  <- brewer.pal(11, "RdBu")
yieldLimits  <- c(3, 9)
yieldBreaks  <- seq(3, 9, 1)
diffLimits   <- c(-2, 2)
diffBreaks   <- seq(-2, 2, 1)

# observed, simulated and difference per region for the given years (same region-years for both)
yield_panels <- function(years, label) {
  sel <- regional[regional$year %in% years, ]
  agg <- aggregate(cbind(obs, sim) ~ region_id, data = sel, FUN = mean)
  m <- match(counties$region_id, agg$region_id)
  list(
    list(values = agg$obs[m], title = paste("Observed", label), zlim = yieldLimits,
         palette = yieldPalette, breaks = yieldBreaks, legend = "Yield (t/ha)"),
    list(values = agg$sim[m], title = paste("Simulated", label), zlim = yieldLimits,
         palette = yieldPalette, breaks = yieldBreaks, legend = "Yield (t/ha)"),
    list(values = agg$sim[m] - agg$obs[m], title = "Simulated - observed", zlim = diffLimits,
         palette = diffPalette, breaks = diffBreaks, legend = "Difference (t/ha)")
  )
}


#### 1 Scatter plot (districts) ####

p <- scatter_plot(regional, paste(myCropName, "- district yields"),
                  pointColor = "#366C62", pointAlpha = 0.15, pointSize = 1.5)
save_plot(p, paste0(myCrop, "_scatter_regional"))


#### 2 National time series ####

lineData <- rbind(data.frame(year = national$year, yield = national$obs, type = "Observed"),
                  data.frame(year = national$year, yield = national$sim, type = "Simulated"))
p <- ggplot(lineData, aes(x = year, y = yield, color = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  annotate("text", x = Inf, y = Inf, label = stats_label(national$obs, national$sim, r2 = FALSE),
           hjust = 1.05, vjust = 1.3, size = 5, lineheight = 1.1) +
  scale_color_manual(values = c("Observed" = "grey30", "Simulated" = "#366C62")) +
  scale_x_continuous(breaks = seq(1995, 2025, by = 5)) +
  labs(title = paste(myCropName, "- national mean yields"),
       subtitle = "Simulated: leave-one-year-out cross-validation",
       x = "", y = "Yield (t/ha)") +
  presentation_theme +
  theme(legend.position = "top", legend.justification = "left",
        legend.background = element_blank())
save_plot(p, paste0(myCrop, "_timeseries_national"), width = 10, height = 5.5)


#### 3 Maps: mean over all years ####

plot_maps(yield_panels(sort(unique(regional$year)), "(mean)"),
          paste0(myCrop, "_map_mean.png"),
          paste(myCropName, "- mean yield", firstYear, "-", max(regional$year)))


#### 4 Maps: dry years ####

for (y in myDryYears) {
  plot_maps(yield_panels(y, y), paste0(myCrop, "_map_", y, ".png"), paste(myCropName, "-", y))
}


#### 5 Maps: FieldConditionModifier and CV RMSE ####

best <- merge(best, noOfObs, by = "region_id", all.x = TRUE)
regionalFcm <- ifelse(!is.na(best$rmse_kg_ha) & best$no_of_obs >= minObsYears, best$FieldConditionModifier, NA)
fcm <- regionalFcm[match(counties$region_id, best$region_id)]
rmse <- cvSummary$rmse_kg_ha[match(counties$region_id, cvSummary$region_id)] / 1000

plot_maps(list(
  list(values = fcm, title = "FieldConditionModifier", zlim = c(0.74, 1.12),
       palette = brewer.pal(9, "YlGnBu"), breaks = seq(0.8, 1.1, 0.1), legend = "FieldConditionModifier (-)"),
  list(values = rmse, title = "Cross-validated RMSE", zlim = c(0, 1.5),
       palette = brewer.pal(9, "YlOrRd"), breaks = seq(0, 1.5, 0.5), legend = "RMSE (t/ha)")
), paste0(myCrop, "_map_fcm_rmse.png"), paste(myCropName, "- regional calibration"), width = 12)
