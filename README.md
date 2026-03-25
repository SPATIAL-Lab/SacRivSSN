# SacRivSSN

Code and data products associated with the manuscript:  
**A Dendritic Strontium River Isoscape for Fisheries Applications in the Sacramento River Basin, California, USA**

---

## Repository structure

### `out/`
R scripts used to reproduce the analyses and figures in the manuscript.

- **FD_sacriv_SSNM.R**  
  Builds the Sacramento River Basin spatial stream network model (SSNM) for river-water 87Sr/86Sr.

- **FD_sacriv_wildkings.R**  
  Generates probabilistic natal-origin and early-rearing habitat maps for example wild adult Chinook salmon using the Sacramento River 87Sr/86Sr isoscape.

- **SacRiv_otolith_assignment_workflow.R**  
  Interactive workflow for assigning otolith 87Sr/86Sr transect intervals to the Sacramento River Basin dendritic isoscape.

- **FD_sacriv_fig2_3_4.R**  
  Reproduces the workflows for manuscript Figures 2, 3, and 4, including 87Sr/86Sr diversity histograms and uncertainty-informed isotopic suite maps.

---

### `data/`
Input data used by the scripts in this repository.

#### Data archived on Zenodo

Most full data products supporting the manuscript are archived at:  
**DOI: https://doi.org/10.5281/zenodo.17595030**

Additional archived data products are available at:  
**DOI: https://doi.org/10.5281/zenodo.17545916**

These include:

- **sacriv32224.ssn.zip**  
  Sacramento River spatial stream network object.

---

#### Data included directly in this GitHub repository

- **intersected_lines.RDS**  
  Sacramento River 87Sr/86Sr isoscape stream-line object for present-day assignments.

- **otolith_ref_data_sac_hatcheries(all).csv**  
  Hatchery reference otolith dataset used to estimate within-population / within-site error for Bayesian probabilistic assignments.

- **sac_sites_clean.RDS**  
  Cleaned river-water sampling sites used for SSNM development. Duplicate site observations were summarized, yielding 106 sites used in `FD_sacriv_SSNM.R`.

- **wr_sub.csv**  
  Subset of endangered winter-run Chinook salmon otolith transects used in example assignments from Phillis et al. (2018).

- **WR_allSr_forR.csv**  
  Full otolith transect dataset used in the assignment workflow examples.

- **srbiso_full.RDS**  
  Full Sacramento River Basin strontium river isoscape, including above-dam river environments.

---

## How to use this repository

A typical workflow is:

1. Run **FD_sacriv_SSNM.R**  
   to build or inspect the Sacramento River 87Sr/86Sr spatial stream network model.

2. Run **FD_sacriv_fig2_3_4.R**  
   to reproduce manuscript Figures 2, 3, and 4.

3. Run **SacRiv_otolith_assignment_workflow.R**  
   to interactively inspect otolith transects and generate natal and rearing assignment maps.

4. Run **FD_sacriv_wildkings.R**  
   for the manuscript example fish assignment analyses.

---

## Required R packages

The scripts in this repository use the following R packages: sf
dplyr
readr
ggplot2
tmap
plotly
patchwork
RColorBrewer
purrr
tibble
tidyr
---

## Notes

- All scripts use **relative paths** assuming the repository structure shown above.
- Larger data products are archived on **Zenodo** and not stored directly in GitHub.
- The `out/` directory contains analysis scripts (not just output figures).

---

## Citation

If you use this repository or associated data products, please cite the manuscript and Zenodo archive:

- Manuscript: *A Dendritic Strontium River Isoscape for Fisheries Applications in the Sacramento River Basin, California, USA*
- Publisher: North American Journal of Fisheries Management (2026)
- Data archive: https://doi.org/10.5281/zenodo.17595030
scales

