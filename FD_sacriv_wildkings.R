# ──────────────────────────────────────────────────────────────────────────────
# Manuscript title: A Strontium River Isoscape for Fisheries Applications 
#                   in the Sacramento Basin, California.
# Background: Case Study — Wild Chinook Natal Assignments (3 adults) via Sr Isoscape
# Authors: Kyle G. Brennan; with data from Willmes et al. (2021) and collaborators
# Submitted to: North American Journal for Fisheries Management 
#
# Brief:
#   Generates probabilistic natal-origin maps for three wild Adult Chinook
#   using a Sacramento River 87Sr/86Sr isoscape (SSN-derived) and a Gaussian,
#   isotope-only likelihood. Outputs interactive tmap objects (and optional HTML)
#   for quick visualization of likely natal tributaries.
#
# Inputs (examples):
#   • Wild-only otolith CSV (e.g., mwilmes2021_otodata.csv)
#   • SSN isoscape lines with predicted Sr and SE (intersected_lines.RDS)
#   • (Optional) hatchery reference CSV to estimate within-site SD
#
# Output:
#   • Per-fish probability-of-origin maps; optional 1×3 panel figure
#
# Notes:
#   • Error model combines within-site + analytical (+ optional isoscape SE).
#   • No spatial priors applied here (isotope-only demonstration).
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# Repro setup (Zenodo record: 10.5281/zenodo.17595030)
# 1) Download and unzip the Zenodo archive to a folder on your machine.
# 2) Set base_dir below to that folder path (the unzipped folder).
# 3) Run the rest of the script; all inputs are read from base_dir.
# ──────────────────────────────────────────────────────────────────────────────

# Example: base_dir <- "/Users/you/Downloads/zenodo_17595030"
base_dir <- "/PATH/TO/YOUR/zenodo_17595030"

paths <- list(
  hatch_csv   = file.path(base_dir, "otolith_ref_data_sac_hatcheries(all).csv"),
  otolith_csv = file.path(base_dir, "mwilmes2021_otodata.csv"),
  lines_rds   = file.path(base_dir, "intersected_lines.RDS"),
  out_dir     = file.path(base_dir, "outputs")
)

# Fail fast if inputs are missing; create outputs dir
req  <- c(paths$hatch_csv, paths$otolith_csv, paths$lines_rds)
miss <- req[!file.exists(req)]
if (length(miss)) {
  stop("Missing required input(s):\n- ",
       paste(basename(miss), collapse = "\n- "),
       "\nCheck that base_dir points to the unzipped Zenodo folder.")
}
dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)

# ──────────────────────────────────────────────────────────────────────────────
# 0) Packages
# ──────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(tmap)
  library(ggplot2)
  library(grid)          # for unit()
  library(RColorBrewer)  # ensures the palette is available
})

# ──────────────────────────────────────────────────────────────────────────────
# 1) Inputs
# ──────────────────────────────────────────────────────────────────────────────

# Hatchery reference data (used for SD mean -> within_site)
hatch_srsr <- readr::read_csv(paths$hatch_csv)

grouped_df <- hatch_srsr %>%
  dplyr::group_by(NAT_LOC) %>%
  dplyr::summarise(
    Count      = n(),
    Avg_Sr8786 = mean(Sr8786, na.rm = TRUE),
    SD_Sr8786  = sd(Sr8786,  na.rm = TRUE),
    Avg_Lat    = mean(Lat_dd,  na.rm = TRUE),
    Avg_Long   = mean(Long_dd, na.rm = TRUE),
    .groups    = "drop"
  )

# Willmes 2021 otolith data (Wild only, Adults)
wilmsotos <- readr::read_csv(paths$otolith_csv) %>%
  dplyr::filter(`Simplified classification` == "Wild",
                `Life stage` == "Adult")

# Keep needed columns and standardize names
gdt <- wilmsotos %>%
  dplyr::select(
    `Life stage`,
    `Grouped classification`,
    `Natal mean 87Sr/86Sr`,
    `2sd`,
    `Simplified classification`,
    dplyr::any_of("Sample ID")
  ) %>%
  dplyr::rename(
    life_stage = `Life stage`,
    group      = `Grouped classification`,
    Sr8786     = `Natal mean 87Sr/86Sr`,
    SE         = `2sd`,
    sample_id0 = `Sample ID`
  ) %>%
  dplyr::mutate(
    sample_id = if (!is.null(sample_id0)) sample_id0
    else paste0(gsub("\\s+", "", life_stage), "_", dplyr::row_number())
  ) %>%
  dplyr::select(-sample_id0)

# Isoscape lines (already filtered upstream in your workflow)
intersected_lines <- readRDS(paths$lines_rds)

# Remove specified RIDs (residual segments above dams)
intersected_lines2 <- intersected_lines %>%
  dplyr::filter(!rid %in% c(
    536, 511, 487, 477, 659,                # above Indian Valley Dam
    655, 617, 595, 616,                     # above Keswick
    668, 753, 589, 613, 652, 583, 551, 540, 550, 584, 560, 569, 537, 539
  ))

# Keep only fields needed (must remain sf)
isoscape <- intersected_lines2 %>%
  dplyr::select(Sr8786, Sr8786.predSE, rid, geometry)

# ──────────────────────────────────────────────────────────────────────────────
# 2) Assignment math setup
# ──────────────────────────────────────────────────────────────────────────────
pid_iso   <- isoscape$Sr8786
pid_isose <- isoscape$Sr8786.predSE

# If you later add a spatial/stream-order prior, plug it into this vector:
pid_isose_mod <- rep(0, length(pid_iso))  # placeholder

# Error terms (same as original)
otointer    <- mean(grouped_df$SD_Sr8786, na.rm = TRUE)  # hatchery pop SD mean (optional context)
# pooled within-site SD
pooled_var <- with(grouped_df,
                   sum((Count - 1) * SD_Sr8786^2, na.rm = TRUE) / sum(Count - 1, na.rm = TRUE))
sigma_within_pop <- sqrt(pooled_var)
within_site <- sigma_within_pop
analyt      <- 0.00011 / 2                               # ICPMS machine analytical (2SD -> SD)
within_pop  <- within_site - analyt                      # retained for completeness

# Total error per segment: (placeholder) + within-site + analytical
error <- sqrt((pid_isose_mod)^2 + (within_site)^2 + (analyt)^2)

# Streams geometry for joining results back to lines (keep sf)
streams <- isoscape %>% dplyr::select(geometry, rid)

# ──────────────────────────────────────────────────────────────────────────────
# 3) Core functions (preserve sf with mutate)
# ──────────────────────────────────────────────────────────────────────────────

# Per-sample isotope-only likelihood -> returns sf with Origin_Prob + sample_id
assign_one_sample <- function(sample_sr, sample_id) {
  ll   <- (1 / sqrt(2 * pi * error^2)) * exp(-1 * (sample_sr - pid_iso)^2 / (2 * error^2))
  prob <- ll / sum(ll)
  
  # Preserve sf class by mutating on the sf object
  streams %>%
    dplyr::mutate(
      Origin_Prob = prob,
      sample_id   = sample_id
    )
}

# Map for one sample (expects sf)
map_one_sample <- function(df_for_sample) {
  tmap::tm_shape(df_for_sample) +
    tmap::tm_lines(
      col = "Origin_Prob",
      palette = "-RdYlBu",
      style = "cont",
      lwd = 4,
      title.col = "Probability of Natal Origin"
    ) +
    tmap::tm_layout(
      legend.outside = TRUE,
      title = paste0(unique(df_for_sample$sample_id), " — Probability of Natal Origin"),
      title.size = 1.0
    )
}

# Build maps for a data.frame of samples (returns named list of tmap objects)
build_maps_for_df <- function(df) {
  stopifnot(all(c("Sr8786", "sample_id") %in% names(df)))
  assignments <- purrr::map2(df$Sr8786, df$sample_id, assign_one_sample)  # list of sf
  maps        <- purrr::map(assignments, map_one_sample)                  # list of tmap
  names(maps) <- df$sample_id
  maps
}

# ──────────────────────────────────────────────────────────────────────────────
# 4) Split by life_stage and build nested outputs
# ──────────────────────────────────────────────────────────────────────────────
gdt_sub <- gdt %>%
  dplyr::filter(life_stage %in% c("Adult", "Juvenile"))

df_adult    <- gdt_sub %>% dplyr::filter(life_stage == "Adult")
df_juvenile <- gdt_sub %>% dplyr::filter(life_stage == "Juvenile")

otoassmaps <- list(
  Adult    = build_maps_for_df(df_adult),
  Juvenile = build_maps_for_df(df_juvenile)
)

# ──────────────────────────────────────────────────────────────────────────────
# 5) Preview interactively
# ──────────────────────────────────────────────────────────────────────────────
tmap::tmap_mode("view")

# Adults — print each map with simple titles
for (i in seq_along(otoassmaps$Adult)) {
  print(otoassmaps$Adult[[i]] + tmap::tm_layout(title = paste("Adult", i)))
}

# ──────────────────────────────────────────────────────────────────────────────
# 6) (Optional) Save each map to HTML per group
# ──────────────────────────────────────────────────────────────────────────────
save_group_maps <- function(maps_list, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  purrr::walk(names(maps_list), function(id) {
    f <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_\\-]+", "_", id), "_natal_origin.html"))
    tmap::tmap_save(maps_list[[id]], filename = f)
  })
}

# Examples (uncomment to use):
# save_group_maps(otoassmaps$Adult,    file.path(paths$out_dir, "Adult"))
# save_group_maps(otoassmaps$Juvenile, file.path(paths$out_dir, "Juvenile"))

# ──────────────────────────────────────────────────────────────────────────────
# 7) Case study: pick 3 Adults and make a 1×3 panel figure
# ──────────────────────────────────────────────────────────────────────────────

# 7.1) Pick the 3 Adults you want to assign
samples_adult3 <- gdt %>%
  dplyr::filter(life_stage == "Adult") %>%
  dplyr::slice(1:3) %>%                  # change which 3 you want here
  dplyr::select(sample_id, Sr8786)

# 7.2) Assignment for each of the 3 Adults (0–1 rescaled probabilities)
assfun <- function(x, id){
  assignIsotopeOnly      <- (1/sqrt((2*pi*error^2))) * exp(-1*(x - pid_iso)^2 / (2*error^2))
  assignIsotopeOnly_norm <- assignIsotopeOnly / sum(assignIsotopeOnly)
  assignIsotopeOnly_max  <- assignIsotopeOnly_norm / max(assignIsotopeOnly_norm)
  
  streams %>%
    dplyr::mutate(
      # Origin_Prob = assignIsotopeOnly_norm,   # <-- use to have probs sum to 1
      Origin_Prob = assignIsotopeOnly_max,      # <-- 0–1 rescaled (max = 1)
      sample_id   = id
    )
}

results_adult3 <- mapply(
  assfun,
  samples_adult3$Sr8786,
  samples_adult3$sample_id,
  SIMPLIFY = FALSE
)

# 7.3) Make tmap objects for the 3 Adults (optional, not used in panel)
mapfun <- function(x){
  tmap::tm_shape(x) +
    tmap::tm_lines(
      col     = "Origin_Prob",
      palette = "-RdYlBu",
      style   = "cont",
      lwd     = 4,
      title.col = "Probability of Natal Origin"
    )
}
otoassmaps_adult3 <- lapply(results_adult3, mapfun)

# 7.4) Combine the three Adult results into one sf with a panel label (ggplot)
panels <- dplyr::bind_rows(
  dplyr::mutate(results_adult3[[1]], Panel = "Adult 1"),
  dplyr::mutate(results_adult3[[2]], Panel = "Adult 2"),
  dplyr::mutate(results_adult3[[3]], Panel = "Adult 3")
)

# Ensure lon/lat for degree axes
if (!sf::st_is_longlat(panels)) panels <- sf::st_transform(panels, 4326)

# Tight bounding box with a small buffer
buffer_km <- 10
crs_ll   <- sf::st_crs(panels)
crs_proj <- 3310  # California Albers (meters)

net_union <- sf::st_union(sf::st_geometry(panels))
bbox_poly <- sf::st_as_sfc(sf::st_bbox(net_union), crs = crs_ll)
bbox_buf  <- bbox_poly |>
  sf::st_transform(crs_proj) |>
  sf::st_buffer(buffer_km * 1000) |>
  sf::st_transform(crs_ll)

bb   <- sf::st_bbox(bbox_buf)
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])

# Plot with a single shared legend and Brewer "-RdYlBu" (direction = -1)
p <- ggplot(panels) +
  geom_sf(aes(color = Origin_Prob), linewidth = 0.5) +
  facet_grid(. ~ Panel) +
  scale_color_distiller(
    palette   = "RdYlBu",
    direction = -1,                 # reverses to "-RdYlBu"
    limits    = c(0, 1),
    breaks    = seq(0, 1, by = 0.1),
    name      = "Scaled probability (0–1)",
    guide     = guide_colorbar(barheight = grid::unit(80, "pt"))
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_line(linewidth = 0.2),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.spacing    = grid::unit(0.6, "lines"),
    axis.title       = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    strip.text       = element_text(face = "bold")
  )

print(p)

# Save panel figure to outputs dir
ggsave(file.path(paths$out_dir, "adult_assignment_1x3_RdYlBu_rev.pdf"),
       plot = p, width = 12, height = 6, units = "in", device = "pdf")

# ──────────────────────────────────────────────────────────────────────────────
# Notes / toggles:
# • To include per-segment isoscape SE directly in the error term, replace:
#     pid_isose_mod <- rep(0, length(pid_iso))
#   with:
#     pid_isose_mod <- pid_isose
#   then recompute:
#     error <- sqrt(pid_isose_mod^2 + within_site^2 + analyt^2)
# • If you prefer map titles to use `group` instead of `sample_id`, swap in the
#   map/build helpers accordingly.
# ──────────────────────────────────────────────────────────────────────────────

