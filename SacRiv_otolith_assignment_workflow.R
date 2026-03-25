# =============================================================================
# Author: Dr. Kyle Gerard Brennan (2026) kyle.brennan@utah.edu
# Manuscript- A Dendritic Strontium River Isoscape for Fisheries 
#            Applications in the Sacramento River Basin, California, USA. 

# PURPOSE
# -------
# Interactive workflow for assigning otolith 87Sr/86Sr transect intervals to a
# Sacramento River Basin dendritic isoscape.
#
# This script is designed to:
#   1) load otolith transect data and SSNM-derived isoscape data from ./data
#   2) display interactive otolith transects in Plotly
#   3) let users define natal / rearing intervals in R
#   4) generate continuous assignment maps in interactive tmap
#   5) reproduce manuscript example fish: 2007_11 and 2008_69
#
# DEFAULT BEHAVIOR
# ----------------
# Plots are displayed in R and NOT saved by default.
#
# OPTIONAL
# --------
# A final section at the end shows how to save the manuscript example figures.
#
# REPO STRUCTURE ASSUMED
# ----------------------
# repo/
# ├─ data/
# │  ├─ WR_allSr_forR.csv
# │  ├─ wr_sub.csv
# │  ├─ intersected_lines.RDS
# │  ├─ otolith_ref_data_sac_hatcheri...csv
# │  └─ sac_sites_clean.RDS
# └─ out/
#    └─ otolith_assignment_workflow.R
#
# REQUIRED COLUMNS IN OTOLITH TRANSECT DATA
# -----------------------------------------
#   ID, Year, Fish_no, radius_in_microns, otoSr, SrV
#
# REQUIRED FIELDS IN intersected_lines.RDS (sf object)
# ----------------------------------------------------
#   rid, Sr8786, Sr8786.predSE, geometry
#
# =============================================================================


# -----------------------------------------------------------------------------
# 0) Packages
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(plotly)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(tmap)
  library(rlang)
  library(scales)
  library(grid)
})

# Interactive tmap mode for RStudio Viewer / browser
tmap_mode("view")


# -----------------------------------------------------------------------------
# 1) Paths
# -----------------------------------------------------------------------------
# Assumes you run this script from the repo root.
data_dir <- "./data"
out_dir  <- "./out"

wr_file       <- file.path(data_dir, "WR_allSr_forR.csv")
wr_sub_file   <- file.path(data_dir, "wr_sub.csv")
isoscape_file <- file.path(data_dir, "intersected_lines.RDS")

# Attempt to locate the hatchery reference file automatically
hatchery_candidates <- list.files(
  data_dir,
  pattern = "^otolith_ref_data_sac_hatcheri.*\\.(csv|CSV)$",
  full.names = TRUE
)

if (length(hatchery_candidates) == 0) {
  stop(
    paste0(
      "Could not find hatchery reference csv in ./data.\n",
      "Expected something like: otolith_ref_data_sac_hatcheries_all.csv"
    ),
    call. = FALSE
  )
}
hatchery_file <- hatchery_candidates[1]


# -----------------------------------------------------------------------------
# 2) Helper checks and utilities
# -----------------------------------------------------------------------------
check_required_cols <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        object_name, " is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

safe_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x)
}

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

rescale_within_group <- function(x) {
  xmax <- max(x, na.rm = TRUE)
  if (!is.finite(xmax) || xmax == 0) return(rep(0, length(x)))
  x / xmax
}


# -----------------------------------------------------------------------------
# 3) Load data
# -----------------------------------------------------------------------------
wr <- read_csv(wr_file, show_col_types = FALSE)
check_required_cols(
  wr,
  c("ID", "Year", "Fish_no", "radius_in_microns", "otoSr", "SrV"),
  "WR_allSr_forR.csv"
)

wr_sub <- NULL
if (file.exists(wr_sub_file)) {
  wr_sub <- read_csv(wr_sub_file, show_col_types = FALSE)
  check_required_cols(
    wr_sub,
    c("ID", "Year", "Fish_no", "radius_in_microns", "otoSr", "SrV"),
    "wr_sub.csv"
  )
}

hatch_srsr <- read_csv(hatchery_file, show_col_types = FALSE)
check_required_cols(
  hatch_srsr,
  c("NAT_LOC", "Sr8786"),
  "hatchery reference csv"
)

intersected_lines <- readRDS(isoscape_file)
if (!inherits(intersected_lines, "sf")) {
  stop("intersected_lines.RDS must be an sf object.", call. = FALSE)
}
check_required_cols(
  intersected_lines,
  c("rid", "Sr8786", "Sr8786.predSE", "geometry"),
  "intersected_lines.RDS"
)


# -----------------------------------------------------------------------------
# 4) Prepare the isoscape used for present-day anadromous assignments
# -----------------------------------------------------------------------------
# Remove residual segments above dams, consistent with manuscript workflow
exclude_rids <- c(
  536, 511, 487, 477, 659,                # above Indian Valley Dam
  655, 617, 595, 616,                     # above Keswick
  668, 753, 589, 613, 652, 583, 551, 540, 550, 584, 560, 569, 537, 539
)

isoscape <- intersected_lines %>%
  filter(!rid %in% exclude_rids) %>%
  dplyr::select(rid, Sr8786, Sr8786.predSE, geometry)

# Ensure valid CRS
if (is.na(st_crs(isoscape))) {
  warning("isoscape has no CRS; assigning EPSG:5070 as fallback.")
  st_crs(isoscape) <- 5070
}


# -----------------------------------------------------------------------------
# 5) Estimate uncertainty terms
# -----------------------------------------------------------------------------
# Pooled within-population SD from hatchery reference groups
grouped_df <- hatch_srsr %>%
  group_by(NAT_LOC) %>%
  summarise(
    Count      = n(),
    Avg_Sr8786 = mean(Sr8786, na.rm = TRUE),
    SD_Sr8786  = safe_sd(Sr8786),
    .groups    = "drop"
  )

pooled_var <- with(
  grouped_df,
  sum((Count - 1) * SD_Sr8786^2, na.rm = TRUE) / sum(Count - 1, na.rm = TRUE)
)
sigma_within_pop <- sqrt(pooled_var)

# Analytical SD from your existing workflow
analytical_sd <- 0.00011 / 2

message("Pooled within-population SD: ", signif(sigma_within_pop, 4))
message("Analytical SD: ", signif(analytical_sd, 4))


# -----------------------------------------------------------------------------
# 6) Interactive Plotly otolith transect viewer
# -----------------------------------------------------------------------------
# This displays otoSr and SrV against radius for a given fish ID.
# Use this to visually choose natal / rearing intervals.

plot_oto_transect <- function(data, fish_id) {
  
  fish_dat <- data %>%
    filter(ID == fish_id) %>%
    arrange(radius_in_microns)
  
  if (nrow(fish_dat) == 0) {
    stop(paste("No rows found for fish ID:", fish_id), call. = FALSE)
  }
  
  yr <- unique(fish_dat$Year)
  fish_no <- unique(fish_dat$Fish_no)
  
  plot_ly(fish_dat) %>%
    add_trace(
      x = ~radius_in_microns,
      y = ~otoSr,
      type = "scatter",
      mode = "lines+markers",
      name = "otoSr (87Sr/86Sr)",
      line = list(width = 1.4),
      marker = list(size = 5),
      text = ~ID,
      hovertemplate = paste(
        "Fish ID: %{text}<br>",
        "Radius: %{x} µm<br>",
        "otoSr: %{y:.6f}<extra></extra>"
      )
    ) %>%
    add_trace(
      x = ~radius_in_microns,
      y = ~SrV,
      type = "scatter",
      mode = "lines",
      name = "SrV",
      yaxis = "y2",
      line = list(dash = "dot", width = 1.1),
      text = ~ID,
      hovertemplate = paste(
        "Fish ID: %{text}<br>",
        "Radius: %{x} µm<br>",
        "SrV: %{y:.3f}<extra></extra>"
      )
    ) %>%
    layout(
      title = paste0(
        "Interactive otolith transect: ", fish_id,
        " (Year ", yr[1], ", Fish ", fish_no[1], ")"
      ),
      xaxis = list(title = "Radius from core (µm)"),
      yaxis = list(title = "otoSr (87Sr/86Sr)"),
      yaxis2 = list(
        title = "SrV",
        overlaying = "y",
        side = "right"
      ),
      hovermode = "closest"
    )
}


# -----------------------------------------------------------------------------
# 7) Static ggplot transect helper (optional, useful for final figures)
# -----------------------------------------------------------------------------
plot_transect_static <- function(data, fish_ids) {
  
  df <- data %>%
    filter(ID %in% fish_ids) %>%
    arrange(ID, radius_in_microns)
  
  if (nrow(df) == 0) {
    stop("No rows found for requested fish IDs.", call. = FALSE)
  }
  
  range_oto <- range(df$otoSr, na.rm = TRUE)
  range_SrV <- range(df$SrV,   na.rm = TRUE)
  
  scale_factor <- diff(range_oto) / diff(range_SrV)
  offset <- range_oto[1] - range_SrV[1] * scale_factor
  
  df <- df %>%
    mutate(
      fish_label = str_replace_all(ID, "_", "-"),
      SrV_scaled = SrV * scale_factor + offset
    )
  
  ggplot(df, aes(x = radius_in_microns, color = fish_label)) +
    geom_line(aes(y = otoSr), linewidth = 0.8) +
    geom_point(aes(y = otoSr), size = 1.5) +
    geom_line(aes(y = SrV_scaled), linetype = "dashed", linewidth = 0.7) +
    scale_y_continuous(
      name = expression(otolith~{}^{87}*Sr/{}^{86}*Sr),
      sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "SrV")
    ) +
    labs(
      x = "Radius from core (µm)",
      color = "Fish ID"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
}


# -----------------------------------------------------------------------------
# 8) Interval summarization
# -----------------------------------------------------------------------------
# Converts a selected radius interval into one observed otolith Sr value
# suitable for continuous assignment.

summarise_interval <- function(data, fish_id, rmin, rmax, interval_label = "interval") {
  
  sub <- data %>%
    filter(
      ID == fish_id,
      radius_in_microns >= rmin,
      radius_in_microns <= rmax
    ) %>%
    filter(!is.na(otoSr))
  
  if (nrow(sub) == 0) {
    warning(
      paste0("No rows for fish ", fish_id, " in interval ", rmin, "-", rmax, " µm")
    )
    return(NULL)
  }
  
  tibble(
    ID = fish_id,
    interval_label = interval_label,
    rmin = rmin,
    rmax = rmax,
    n_spots = nrow(sub),
    Sr8786_obs = mean(sub$otoSr, na.rm = TRUE),
    Sr8786_obs_sd = safe_sd(sub$otoSr)
  )
}


# -----------------------------------------------------------------------------
# 9) Continuous assignment function
# -----------------------------------------------------------------------------
# Computes posterior probabilities over the isoscape.
#
# Total uncertainty combines:
#   - isoscape prediction SE
#   - pooled within-population SD
#   - analytical SD
#   - SD within the selected otolith interval (if >1 point)

assign_interval_sf <- function(
    isoscape_sf,
    Sr8786_obs,
    sample_sd = NA_real_,
    sigma_within_pop,
    analytical_sd,
    fish_id = "sample",
    interval_label = "interval",
    rmin = NA_real_,
    rmax = NA_real_
) {
  
  sample_sd_use <- ifelse(is.na(sample_sd), 0, sample_sd)
  
  total_sd <- sqrt(
    isoscape_sf$Sr8786.predSE^2 +
      sigma_within_pop^2 +
      analytical_sd^2 +
      sample_sd_use^2
  )
  
  likelihood <- dnorm(
    x = Sr8786_obs,
    mean = isoscape_sf$Sr8786,
    sd = total_sd
  )
  
  posterior <- likelihood / sum(likelihood, na.rm = TRUE)
  
  isoscape_sf %>%
    mutate(
      fish_id = fish_id,
      interval_label = interval_label,
      rmin = rmin,
      rmax = rmax,
      Sr8786_obs = Sr8786_obs,
      Sr8786_obs_sd = sample_sd,
      total_sd = total_sd,
      Origin_Prob = posterior
    )
}


# -----------------------------------------------------------------------------
# 10) Wrapper to run multiple intervals
# -----------------------------------------------------------------------------
# intervals_tbl must contain:
#   ID, interval_label, rmin, rmax

run_assignments_for_intervals <- function(
    data,
    isoscape_sf,
    intervals_tbl,
    sigma_within_pop,
    analytical_sd
) {
  
  check_required_cols(
    intervals_tbl,
    c("ID", "interval_label", "rmin", "rmax"),
    "intervals_tbl"
  )
  
  interval_summaries <- purrr::pmap_dfr(
    intervals_tbl,
    function(ID, interval_label, rmin, rmax, ...) {
      summarise_interval(
        data = data,
        fish_id = ID,
        rmin = rmin,
        rmax = rmax,
        interval_label = interval_label
      )
    }
  )
  
  assignment_list <- purrr::pmap(
    interval_summaries,
    function(ID, interval_label, rmin, rmax, n_spots, Sr8786_obs, Sr8786_obs_sd) {
      assign_interval_sf(
        isoscape_sf = isoscape_sf,
        Sr8786_obs = Sr8786_obs,
        sample_sd = Sr8786_obs_sd,
        sigma_within_pop = sigma_within_pop,
        analytical_sd = analytical_sd,
        fish_id = ID,
        interval_label = interval_label,
        rmin = rmin,
        rmax = rmax
      )
    }
  )
  
  assignments_sf <- bind_rows(assignment_list) %>%
    mutate(
      panel_label = paste0(
        fish_id, ": ", interval_label,
        " (", rmin, "-", rmax, " µm)"
      )
    )
  
  list(
    interval_summaries = interval_summaries,
    assignments_sf = assignments_sf
  )
}


# -----------------------------------------------------------------------------
# 11) Interactive tmap assignment viewer
# -----------------------------------------------------------------------------
# By default, probabilities are rescaled 0-1 within each panel to make visual
# comparison easier.

plot_assignment_map <- function(assignments_sf, scale_within_panel = TRUE) {
  
  map_sf <- assignments_sf
  
  if (scale_within_panel) {
    map_sf <- map_sf %>%
      group_by(panel_label) %>%
      mutate(Origin_Prob_plot = rescale_within_group(Origin_Prob)) %>%
      ungroup()
    title_name <- "Scaled probability"
  } else {
    map_sf <- map_sf %>%
      mutate(Origin_Prob_plot = Origin_Prob)
    title_name <- "Posterior probability"
  }
  
  tm_shape(map_sf) +
    tm_lines(
      col = "Origin_Prob_plot",
      palette = "-RdYlBu",
      style = "cont",
      lwd = 4,
      title.col = title_name,
      popup.vars = c(
        "Fish ID"      = "fish_id",
        "Interval"     = "interval_label",
        "rmin"         = "rmin",
        "rmax"         = "rmax",
        "Obs 87Sr/86Sr"= "Sr8786_obs",
        "Prob"         = "Origin_Prob"
      )
    ) +
    tm_facets(by = "panel_label", ncol = 2) +
    tm_view(set.view = c(-121.5, 39.5, 7)) +
    tm_layout(
      legend.outside = TRUE,
      frame = FALSE
    )
}


# -----------------------------------------------------------------------------
# 12) Manuscript example intervals
# -----------------------------------------------------------------------------
# Edit here if your final chosen intervals differ.
# These are currently set to match your latest workflow logic.

example_intervals <- tibble::tribble(
  ~ID,       ~interval_label, ~rmin, ~rmax,
  "2007_11", "natal",          230,   300,
  "2007_11", "rearing",        420,   540,
  "2008_69", "natal",          230,   300,
  "2008_69", "rearing",        431,   556
)


# -----------------------------------------------------------------------------
# 13) USER INSTRUCTIONS
# -----------------------------------------------------------------------------
#
# To use this workflow on your own fish:
#
# Step 1:
#   Ensure your otolith transect data are in ./data/WR_allSr_forR.csv
#   with required columns:
#     ID, Year, Fish_no, radius_in_microns, otoSr, SrV
#
# Step 2:
#   Display the interactive transect for your fish:
#     plot_oto_transect(wr, "YOUR_FISH_ID")
#
# Step 3:
#   Choose one or more radius intervals based on the transect.
#
# Step 4:
#   Define an interval table like this:
#
#   my_intervals <- tibble::tribble(
#     ~ID,            ~interval_label, ~rmin, ~rmax,
#     "YOUR_FISH_ID", "natal",          230,   300,
#     "YOUR_FISH_ID", "rearing",        420,   540
#   )
#
# Step 5:
#   Run assignments:
#
#   my_results <- run_assignments_for_intervals(
#     data = wr,
#     isoscape_sf = isoscape,
#     intervals_tbl = my_intervals,
#     sigma_within_pop = sigma_within_pop,
#     analytical_sd = analytical_sd
#   )
#
# Step 6:
#   View interval summaries:
#     my_results$interval_summaries
#
# Step 7:
#   View interactive assignment maps:
#     plot_assignment_map(my_results$assignments_sf)
#
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# 14) Interactive manuscript examples in R
# -----------------------------------------------------------------------------
# These objects are created but not saved by default.

# Interactive Plotly transects
p_transect_2007_11 <- plot_oto_transect(wr, "2007_11")
p_transect_2008_69 <- plot_oto_transect(wr, "2008_69")

# Print these in the R console / RStudio to view interactively:
 p_transect_2007_11
 p_transect_2008_69

# Run manuscript interval assignments
example_results <- run_assignments_for_intervals(
  data = wr,
  isoscape_sf = isoscape,
  intervals_tbl = example_intervals,
  sigma_within_pop = sigma_within_pop,
  analytical_sd = analytical_sd
)

example_interval_summaries <- example_results$interval_summaries
example_assignments_sf <- example_results$assignments_sf

# Inspect interval summaries in R
 example_interval_summaries

# Interactive tmap map
example_assignment_map <- plot_assignment_map(example_assignments_sf)

# Print in R to view interactively:
example_assignment_map


# -----------------------------------------------------------------------------
# 15) Fish-specific manuscript example maps
# -----------------------------------------------------------------------------
# Useful if you want to inspect one fish at a time interactively.

assign_2007_11 <- example_assignments_sf %>%
  filter(fish_id == "2007_11")

assign_2008_69 <- example_assignments_sf %>%
  filter(fish_id == "2008_69")

map_2007_11 <- plot_assignment_map(assign_2007_11)
map_2008_69 <- plot_assignment_map(assign_2008_69)

# Print in R to view:
 map_2007_11
 map_2008_69


# -----------------------------------------------------------------------------
# 16) Optional: helper function to build a one-off fish workflow
# -----------------------------------------------------------------------------
# This is convenient for users working on one fish at a time.

run_single_fish_workflow <- function(data, fish_id, intervals_tbl) {
  
  message("Displaying interactive transect for ", fish_id)
  print(plot_oto_transect(data, fish_id))
  
  fish_results <- run_assignments_for_intervals(
    data = data,
    isoscape_sf = isoscape,
    intervals_tbl = intervals_tbl,
    sigma_within_pop = sigma_within_pop,
    analytical_sd = analytical_sd
  )
  
  message("Interval summaries:")
  print(fish_results$interval_summaries)
  
  message("Displaying interactive assignment map...")
  fish_map <- plot_assignment_map(fish_results$assignments_sf)
  print(fish_map)
  
  invisible(fish_results)
}

# Example usage:
# my_intervals <- tibble::tribble(
#   ~ID,       ~interval_label, ~rmin, ~rmax,
#   "2007_11", "natal",          230,   300,
#   "2007_11", "rearing",        420,   540
# )
# run_single_fish_workflow(wr, "2007_11", my_intervals)


# -----------------------------------------------------------------------------
# 17) Optional save section for manuscript examples
# -----------------------------------------------------------------------------
# Nothing below runs unless you set save_example_figures <- TRUE

save_example_figures <- FALSE

if (save_example_figures) {
  
  # Create output subfolder if desired
  fig_dir <- file.path(out_dir, "manuscript_example_figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Static transect figure for manuscript examples
  p_static_examples <- plot_transect_static(wr, c("2007_11", "2008_69"))
  
  ggsave(
    filename = file.path(fig_dir, "transects_2007_11_2008_69.pdf"),
    plot = p_static_examples,
    width = 7,
    height = 4,
    units = "in"
  )
  
  # Save interactive tmap widgets as html if desired
  # Requires htmlwidgets; tmap view objects can be saved with tmap_save
  tmap_save(
    example_assignment_map,
    filename = file.path(fig_dir, "example_assignment_map.html")
  )
  
  tmap_save(
    map_2007_11,
    filename = file.path(fig_dir, "assignment_map_2007_11.html")
  )
  
  tmap_save(
    map_2008_69,
    filename = file.path(fig_dir, "assignment_map_2008_69.html")
  )
  
  message("Saved manuscript example figures to: ", fig_dir)
}


# -----------------------------------------------------------------------------
# 18) End of script
# -----------------------------------------------------------------------------
message("Otolith assignment workflow loaded successfully.")
message("To view the manuscript example transects, print:")
message("  p_transect_2007_11")
message("  p_transect_2008_69")
message("To view the combined manuscript example assignment map, print:")
message("  example_assignment_map")