# =============================================================================
# FD_sacriv_fig2_3_4.R
# =============================================================================
#
# Title:
#   Reproducible workflow for Sacramento River Basin Sr isoscape manuscript
#   Figures 2, 3, and 4
#
# Purpose:
#   This script reproduces the figure workflows used to characterize
#   87Sr/86Sr diversity and uncertainty-informed isotopic suites across the
#   Sacramento River Basin isoscapes described in the manuscript.
#
#   Specifically, this script reproduces:
#     - Figure 2: histogram comparison of modeled 87Sr/86Sr values across
#                 present-day anadromous and historical (pre-dam) stream networks
#     - Figure 3: uncertainty-informed isotopic suites for the present-day
#                 anadromous network
#     - Figure 4: uncertainty-informed isotopic suites for the historical
#                 full basin network
#
# Data inputs expected in ./data:
#   - srbiso_full.RDS
#       Full historical / pre-dam isoscape as an sf object
#   - intersected_lines.RDS
#       Present-day network isoscape as an sf object
#
# Notes:
#   - The present-day anadromous network is created by removing stream segments
#     upstream of impassable dams from intersected_lines.RDS.
#   - The root mean squared prediction error (RMSE) of the SSNM prediction
#     standard errors is used both to:
#       (i) threshold the histogram bin width in Figure 2
#       (ii) estimate the number of isotopically resolvable suites (K)
#            for Figures 3 and 4
#
# Output:
#   By default, figures are generated in the R session.
#   An optional save section is included at the end.
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load packages
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(scales)
  library(tibble)
  library(tmap)
})

# Optional interactive map mode for inspecting sf objects
tmap_mode("view")


# -----------------------------------------------------------------------------
# 2. Define paths
# -----------------------------------------------------------------------------
# Assumes script is run from the repo root.
data_dir <- "./data"
out_dir  <- "./out"

full_isoscape_file <- file.path(data_dir, "srbiso_full.RDS")
anadromous_file    <- file.path(data_dir, "intersected_lines.RDS")

stopifnot(file.exists(full_isoscape_file))
stopifnot(file.exists(anadromous_file))


# -----------------------------------------------------------------------------
# 3. Read isoscape data
# -----------------------------------------------------------------------------
# srbiso_full: full historical / pre-dam network isoscape
# intersected_lines: present-day network, later filtered to anadromous reaches
srbiso_full <- readRDS(full_isoscape_file)
intersected_lines <- readRDS(anadromous_file)

if (!inherits(srbiso_full, "sf")) {
  stop("srbiso_full.RDS must be an sf object.")
}
if (!inherits(intersected_lines, "sf")) {
  stop("intersected_lines.RDS must be an sf object.")
}

required_cols <- c("rid", "Sr8786", "geometry")
if (!all(required_cols %in% names(srbiso_full))) {
  stop("srbiso_full.RDS must contain: rid, Sr8786, geometry")
}
if (!all(required_cols %in% names(intersected_lines))) {
  stop("intersected_lines.RDS must contain: rid, Sr8786, geometry")
}

# The historical object may have predSE instead of Sr8786.predSE
if ("predSE" %in% names(srbiso_full) && !"Sr8786.predSE" %in% names(srbiso_full)) {
  srbiso_full <- srbiso_full %>% rename(Sr8786.predSE = predSE)
}

if (!"Sr8786.predSE" %in% names(srbiso_full)) {
  stop("srbiso_full.RDS must contain Sr8786.predSE (or predSE).")
}
if (!"Sr8786.predSE" %in% names(intersected_lines)) {
  stop("intersected_lines.RDS must contain Sr8786.predSE.")
}


# -----------------------------------------------------------------------------
# 4. Build the present-day anadromous isoscape
# -----------------------------------------------------------------------------
# Following the manuscript workflow, remove stream segments above dams or
# otherwise inaccessible to present-day anadromous salmon.
exclude_rids <- c(
  536, 511, 487, 477, 659,   # above Indian Valley Dam
  655, 617, 595, 616,        # above Keswick
  668, 753, 589, 613, 652, 583, 551, 540, 550, 584, 560, 569, 537, 539
)

anadiso <- intersected_lines %>%
  filter(!rid %in% exclude_rids) %>%
  dplyr::select(rid, Sr8786, Sr8786.predSE, geometry)

# Optional quick inspection
# tm_shape(anadiso) + tm_lines(col = "Sr8786")
# tm_shape(srbiso_full) + tm_lines(col = "Sr8786")


# -----------------------------------------------------------------------------
# 5. Compute RMSE for each isoscape
# -----------------------------------------------------------------------------
# The manuscript uses the RMSE derived from the prediction standard errors
# to define a practical isotopic resolution threshold.
rmse_from_predse <- function(x) {
  sqrt(mean(x^2, na.rm = TRUE))
}

rmse_anad <- rmse_from_predse(anadiso$Sr8786.predSE)
rmse_pre  <- rmse_from_predse(srbiso_full$Sr8786.predSE)

message("Present-day anadromous RMSE: ", signif(rmse_anad, 6))
message("Historical pre-dam RMSE: ", signif(rmse_pre, 6))


# -----------------------------------------------------------------------------
# 6. Prepare histogram data for Figure 2
# -----------------------------------------------------------------------------
# To match the manuscript logic, stream segments are thresholded by each
# isoscape's RMSE before plotting distributions.
anad_df <- anadiso %>%
  st_drop_geometry() %>%
  filter(Sr8786.predSE <= rmse_anad) %>%
  mutate(Isoscape = "Anadromous")

pre_df <- srbiso_full %>%
  st_drop_geometry() %>%
  filter(Sr8786.predSE <= rmse_pre) %>%
  mutate(Isoscape = "Pre-Dam")


# -----------------------------------------------------------------------------
# 7. Define the previously used discrete Sr groups for Figure 2
# -----------------------------------------------------------------------------
# These are the legacy discrete ranges shown for comparison against the
# continuous SSNM-derived isoscape distributions.
sr_ranges <- tibble::tribble(
  ~Region, ~xmin,     ~xmax,     ~fill_color,
  "AME",   0.71021,   0.71031,   "darkgoldenrod1",
  "BAT",   0.70382,   0.70431,   "lightblue",
  "BUT",   0.70466,   0.70489,   "green3",
  "CNH",   0.70481,   0.70615,   "gold",
  "DEE",   0.70407,   0.70413,   "orchid",
  "FEA",   0.705702,  0.706526,  "skyblue2",
  "FRH",   0.7068,    0.70739,   "orange",
  "MIL",   0.70392,   0.7042,    "darkseagreen3",
  "NIH",   0.70968,   0.70989,   "tan",
  "THE",   0.70569,   0.70599,   "slateblue3",
  "YUB",   0.70762,   0.70852,   "plum2"
)


# -----------------------------------------------------------------------------
# 8. Figure 2: histogram comparison of isoscape diversity
# -----------------------------------------------------------------------------
fig2 <- ggplot() +
  geom_rect(
    data = sr_ranges,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, color = Region),
    fill = NA,
    linewidth = 1.1,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  geom_histogram(
    data = pre_df,
    aes(x = Sr8786, fill = "Pre-Dam"),
    binwidth = rmse_pre,
    position = "identity",
    color = "black"
  ) +
  geom_histogram(
    data = anad_df,
    aes(x = Sr8786, fill = "Anadromous"),
    binwidth = rmse_anad,
    position = "identity",
    color = "black"
  ) +
  scale_fill_manual(
    name = "Sr River Isoscapes",
    values = c("Anadromous" = "indianred2", "Pre-Dam" = "royalblue")
  ) +
  scale_color_manual(
    name = expression("Discrete " ^87 * "Sr/" ^86 * "Sr Groups"),
    values = setNames(sr_ranges$fill_color, sr_ranges$Region)
  ) +
  labs(
    title = "Thresholded Sr Isotope Diversity Across Stream Segments",
    subtitle = paste0(
      "RMSE (Historical) = ", signif(rmse_pre, 6),
      " | RMSE (Present-day Anadromous) = ", signif(rmse_anad, 6)
    ),
    x = expression("" ^87 * "Sr/" ^86 * "Sr Ratio"),
    y = "Count of Stream Segments"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(fig2)


# -----------------------------------------------------------------------------
# 9. Helper functions for Figures 3 and 4
# -----------------------------------------------------------------------------

# Build a long palette by concatenating several Brewer palettes.
get_brewer_colors <- function() {
  c(
    brewer.pal(12, "Set3"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Pastel1"),
    brewer.pal(8, "Set2"),
    brewer.pal(12, "Paired")
  )
}

# Given an sf isoscape and RMSE, estimate number of resolvable clusters.
estimate_cluster_count <- function(sf_obj, rmse_value) {
  sr_range <- range(sf_obj$Sr8786, na.rm = TRUE)
  k <- round(diff(sr_range) / rmse_value)
  max(k, 1)
}

# Apply k-means clustering to the Sr8786 field and create ordered cluster labels.
cluster_isoscape <- function(sf_obj, rmse_value, seed = 42) {
  
  n_clusters <- estimate_cluster_count(sf_obj, rmse_value)
  
  set.seed(seed)
  sf_obj$cluster_kmeans <- kmeans(sf_obj$Sr8786, centers = n_clusters)$cluster
  
  cluster_labels <- sf_obj %>%
    st_drop_geometry() %>%
    group_by(cluster_kmeans) %>%
    summarise(
      min_sr = round(min(Sr8786, na.rm = TRUE), 5),
      max_sr = round(max(Sr8786, na.rm = TRUE), 5),
      .groups = "drop"
    ) %>%
    arrange(min_sr) %>%
    mutate(
      cluster_ordered = row_number(),
      cluster_label = paste0(cluster_ordered, " (", min_sr, "–", max_sr, ")")
    )
  
  sf_obj <- sf_obj %>%
    left_join(
      cluster_labels %>% dplyr::select(cluster_kmeans, cluster_label, cluster_ordered),
      by = "cluster_kmeans"
    )
  
  color_vec <- get_brewer_colors()[1:n_clusters]
  color_map <- setNames(color_vec, sort(unique(sf_obj$cluster_label)))
  
  list(
    sf = sf_obj,
    labels = cluster_labels,
    n_clusters = n_clusters,
    color_map = color_map
  )
}

# Build the combined map + dot plot figure used for Figures 3 and 4.
plot_cluster_figure <- function(clustered_obj, rmse_value, map_title, xbreaks = NULL) {
  
  sf_obj <- clustered_obj$sf
  color_map <- clustered_obj$color_map
  n_clusters <- clustered_obj$n_clusters
  
  p_map <- ggplot(sf_obj) +
    geom_sf(aes(color = factor(cluster_label)), linewidth = 0.6) +
    scale_color_manual(values = color_map, guide = "none") +
    labs(
      title = map_title,
      subtitle = paste0("RMSE = ", signif(rmse_value, 6), ", K = ", n_clusters),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11)
  
  p_dot <- ggplot(
    sf_obj,
    aes(x = Sr8786, y = factor(cluster_label), color = cluster_label)
  ) +
    geom_jitter(height = 0.2, alpha = 0.8, size = 1.7) +
    scale_color_manual(values = color_map, guide = "none") +
    labs(
      x = expression("" ^87 * "Sr/" ^86 * "Sr Ratio"),
      y = "Isoscape-Derived Suite"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9)
    )
  
  if (!is.null(xbreaks)) {
    p_dot <- p_dot +
      scale_x_continuous(
        breaks = xbreaks,
        labels = scales::number_format(accuracy = 0.001)
      )
  }
  
  p_map + p_dot + patchwork::plot_layout(ncol = 2)
}


# -----------------------------------------------------------------------------
# 10. Figure 3: anadromous uncertainty-informed isotopic suites
# -----------------------------------------------------------------------------
clustered_anad <- cluster_isoscape(anadiso, rmse_anad)

fig3 <- plot_cluster_figure(
  clustered_obj = clustered_anad,
  rmse_value = rmse_anad,
  map_title = "Present-day Anadromous Stream Segments Clustered by Sr Isotope Ratio",
  xbreaks = seq(0.701, 0.713, by = 0.003)
)

print(fig3)


# -----------------------------------------------------------------------------
# 11. Figure 4: historical / pre-dam uncertainty-informed isotopic suites
# -----------------------------------------------------------------------------
clustered_pre <- cluster_isoscape(srbiso_full, rmse_pre)

fig4 <- plot_cluster_figure(
  clustered_obj = clustered_pre,
  rmse_value = rmse_pre,
  map_title = "Historical Stream Segments Clustered by Sr Isotope Ratio",
  xbreaks = seq(0.701, 0.713, by = 0.003)
)

print(fig4)


# -----------------------------------------------------------------------------
# 12. Optional quick inspection maps
# -----------------------------------------------------------------------------
# Uncomment to inspect the clustered isoscapes interactively.
#
# tm_shape(clustered_anad$sf) + tm_lines(col = "cluster_label")
# tm_shape(clustered_pre$sf) + tm_lines(col = "cluster_label")


# -----------------------------------------------------------------------------
# 13. Optional save section
# -----------------------------------------------------------------------------
# By default, nothing is saved. Set save_figures <- TRUE to write PDFs.
save_figures <- FALSE

if (save_figures) {
  
  fig_dir <- file.path(out_dir, "figures_2_3_4")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  ggsave(
    filename = file.path(fig_dir, "figure2_histogram_isoscape_diversity.pdf"),
    plot = fig2,
    width = 8,
    height = 6
  )
  
  ggsave(
    filename = file.path(fig_dir, "figure3_anadromous_isotope_suites.pdf"),
    plot = fig3,
    width = 7,
    height = 5
  )
  
  ggsave(
    filename = file.path(fig_dir, "figure4_historical_isotope_suites.pdf"),
    plot = fig4,
    width = 7,
    height = 5
  )
  
  message("Saved figures to: ", fig_dir)
}


# -----------------------------------------------------------------------------
# 14. End of workflow
# -----------------------------------------------------------------------------
message("Figures 2, 3, and 4 workflow completed successfully.")

