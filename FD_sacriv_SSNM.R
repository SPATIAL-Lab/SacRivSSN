# ──────────────────────────────────────────────────────────────────────────────
# Sacramento River SSNM — Clean Workflow (Part 1/2)
# Author: Kyle Gerard Brennan  |  kyle.brennan@utah.edu  |  2023–2025
# Repro DOI: 10.5281/zenodo.17595030
#
# Manuscript title: A Strontium River Isoscape for Fisheries Applications 
#                   in the Sacramento Basin, California.
# Submitted to: North American Journal for Fisheries Management 
# ──────────────────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────────────────
# Packages (load with checks)
# ──────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  ensure_pkg <- function(pkgs) {
    to_get <- setdiff(pkgs, rownames(installed.packages()))
    if (length(to_get)) install.packages(to_get)
    invisible(lapply(pkgs, require, character.only = TRUE))
  }
  
  # CRAN / commonly available
  cran_pkgs <- c(
    "dplyr","moments","sf","terra","writexl","plotly","combinat","car","tmap"
  )
  ensure_pkg(cran_pkgs)
  
  # SSN2 is on CRAN; SSN (legacy) may require the GitHub release tarball
  if (!requireNamespace("SSN2", quietly = TRUE)) install.packages("SSN2")
  suppressMessages(library(SSN2))
  
  if (!requireNamespace("SSN", quietly = TRUE)) {
    message("\nSSN package not found. If install from CRAN fails, use the release tarball:")
    message("install.packages(\"https://github.com/jayverhoef/SSN/releases/download/v1.1.18/SSN_1.1.18.tar.gz\",",
            " repos = NULL, type = \"source\")\n")
    try(install.packages(
      "https://github.com/jayverhoef/SSN/releases/download/v1.1.18/SSN_1.1.18.tar.gz",
      repos = NULL, type = "source"
    ), silent = TRUE)
  }
  suppressMessages(library(SSN))
  
  # set tmap to interactive view mode
  tmap::tmap_mode("view")
})

# Sanity check: confirm the SSN object folder exists
if (!dir.exists(paths$ssn_dir)) {
  stop("SSN directory not found at: ", paths$ssn_dir,
       "\nCheck that base_dir points to your unzipped Zenodo folder.\n",
       "Tip: use find_in_base(\"^sacriv32224\\.ssn$\") to auto-discover.")
}
# ──────────────────────────────────────────────────────────────────────────────
# Repro setup (Zenodo: 10.5281/zenodo.17595030)
# 1) Download and unzip the Zenodo archive to a folder on your machine.
# 2) Set base_dir to that unzipped folder.
# 3) Run the script; file paths below point into base_dir.
# ──────────────────────────────────────────────────────────────────────────────

# Example:
# base_dir <- "/Users/you/Downloads/zenodo_17595030"
base_dir <- "/PATH/TO/YOUR/zenodo_17595030"

# (Simple) explicit path(s) within the Zenodo folder:
paths <- list(
  ssn_dir = file.path(base_dir, "sacriv32224.ssn"),  # SSN object folder
  out_dir = file.path(base_dir, "outputs")
)

# Optional: auto-find (useful if folder structure changes)
find_in_base <- function(pattern) {
  hits <- list.files(base_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(hits)) stop("File not found under base_dir: ", pattern)
  hits[[1]]
}
# If you prefer auto-discovery, uncomment this line:
# paths$ssn_dir <- find_in_base("^sacriv32224\\.ssn$")

# Create outputs dir if you write anything later
dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# PART 1 — Quick run of selected models, outlier handling, predictions
#         (builds the SSN objects used again in Part 2)
# =============================================================================

## Step 1. Import SSN object and prediction points ----------------------------
sacssn <- importSSN(paths$ssn_dir, o.write = TRUE)

# Attach mid-network prediction points
sacssn <- importPredpts(sacssn, "midpred", "ssn")

## Step 2. Clean raw fields and log-transform [Sr] ----------------------------
ObsDF <- getSSNdata.frame(sacssn)
#remove incorrectly snapped data point
ObsDF <- ObsDF %>%
  mutate(
    Sr_mgL = ifelse(pid == 213, NA, sr_mgL),
    Sr8786 = ifelse(pid == 213, NA, Sr8786)
  )
# REQUIRED: normality for [Sr]
ObsDF <- ObsDF %>% mutate(ln_sr = log(sr_mgL))
sacssn <- putSSNdata.frame(ObsDF, sacssn)

## Step 3. Add required additive function (flow weight precursor) -------------
# REQUIRED for [Sr] SSNM: precip accumulation as avf
sacssn <- additive.function(sacssn, "ppt_c", "afv_fasac")

## Step 4. Build distance matrices (REQUIRED for SSN models) ------------------
createDistMat(sacssn, predpts = "midpred", o.write = T)

## Step 5. Fit initial spatial [Sr] model (best AIC LM already known) ---------
# Fixed to tail-up (solute-like downstream transport)
firstmod <- glmssn(ln_sr ~ ign_ap+ign_bp+scp+agemin+agemax+x5p+x7p,
                   EstMeth = "ML",
                   ssn.object = sacssn,
                   CorModels = c("Exponential.tailup"),
                   addfunccol = "afv_fasac")
firstmod

## Step 6. Identify residual outliers using Tukey whiskers --------------------
# Derive residuals
resm1  <- residuals(firstmod)
ObsDFr <- getSSNdata.frame(resm1)

# Interactive boxplot to visualize limits (optional)
plot_ly(data = ObsDFr, y = ~`_resid_`, type = "box", quartilemethod = "exclusive")

# Compute lower/upper whisker bounds (Q1 ± 1.5*IQR; use whisker extrema)
x <- ObsDFr$`_resid_`; x <- x[is.finite(x)]
qtype <- c(linear=7, exclusive=6, inclusive=8)[["exclusive"]]
q     <- quantile(x, c(.25,.75), type=qtype); iqr <- diff(q)
lower <- min(x[x >= q[1] - 1.5*iqr], na.rm = TRUE)
upper <- max(x[x <= q[2] + 1.5*iqr], na.rm = TRUE)

# Flag and set outside-whisker residuals to NA in ln_sr
ObsDF <- getSSNdata.frame(sacssn)
indOutlier <- ObsDFr["_resid_"] > upper | ObsDFr["_resid_"] < lower
indOutlier[is.na(indOutlier)] <- FALSE
ObsDF[indOutlier, "ln_sr"] <- NA
ObsDF$ln_sr  # quick check
sacssn1 <- putSSNdata.frame(ObsDF, sacssn)

## Step 7. Refit [Sr] spatial model after outlier removal ---------------------
mod2 <- glmssn(ln_sr ~ ign_ap+ign_bp+scp+agemin+agemu+x5p+x7p,
               EstMeth = "ML",
               ssn.object = sacssn1,
               CorModels = c("LinearSill.tailup" , "LinearSill.taildown" , "Gaussian.Euclid"),
               addfunccol = "afv_fasac")
mod2  # “clean” ML ln[Sr] model

## Step 8. LOOCV diagnostics for the cleaned [Sr] model -----------------------
CrossValidationStatsSSN(mod2)
cv.out <- CrossValidationSSN(mod2)

par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(mod2$sampinfo$z, cv.out[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction")
abline(0, 1)

plot(na.omit(getSSNdata.frame(sacssn1)[, "ln_sr"]),
     cv.out[, "cv.se"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction SE")

# (Optional) Re-inspect residuals visually
resm2  <- residuals(mod2)
hist(resm2)
ObsDFr <- getSSNdata.frame(resm2)
ObsDF  <- getSSNdata.frame(sacssn1)
plot_ly(data = ObsDFr, y = ~`_resid_`, type = "box", quartilemethod = "exclusive")

# Re-save post-visualization (no additional removal here by design)
sacssn2 <- putSSNdata.frame(ObsDF, sacssn1)

## Step 9. Assign final [Sr] model (ML + chosen correlation set) --------------
srfmod <- glmssn(ln_sr ~ ign_ap+ign_bp+scp+agemin+agemu+x5p+x7p,
                 EstMeth = "ML",
                 ssn.object = sacssn1,
                 CorModels = c("LinearSill.tailup" , "LinearSill.taildown" , "Gaussian.Euclid"),
                 addfunccol = "afv_fasac")
summary(srfmod)
varcomp(srfmod)

# LOOCV vs observed for [Sr]
srfmodcv.out <- CrossValidationSSN(srfmod)
srfmodcv.lm  <- lm(srfmodcv.out[, "cv.pred"] ~ srfmod$sampinfo$z)
summary(srfmodcv.lm)

## Step 10. Predict [Sr] at mid-network points and attach to SSN --------------
srfmod.pred <- predict(srfmod, "midpred")

# Back-transform predictions and SE
predDF <- getSSNdata.frame(srfmod.pred, Name = 'midpred')
predDF <- predDF %>% mutate(predsr_c    = exp(ln_sr))
predDF <- predDF %>% mutate(predsrSE_c  = exp(ln_sr.predSE))

# Save prediction attributes onto SSN (pred layer + observed layer join)
sacssn_a <- putSSNdata.frame(predDF, sacssn1, "midpred")
pd.id    <- predDF %>% dplyr::select(predsr_c)
sacssn_a@data <- cbind(sacssn_a@data, pd.id)
pd.id    <- predDF %>% dplyr::select(predsrSE_c)
sacssn_a@data <- cbind(sacssn_a@data, pd.id)

obsDFa   <- getSSNdata.frame(sacssn_a)
ob.id    <- obsDFa
pd.id    <- predDF %>% dplyr::select(rid, predsr_c)
obsDF_a  <- left_join(ob.id, pd.id, by="rid")
sacssn_a <- putSSNdata.frame(obsDF_a, sacssn_a)

## Step 11. Build hydrologic weight for isotope model -------------------------
# REQUIRED: pptSr = predicted [Sr] * precipitation accumulation
sacssn_a@data$pptSr <- sacssn_a@data$predsr_c * sacssn_a@data$ppt_c
sacssn_a <- additive.function(sacssn_a, 'pptSr', 'afv_pptSr')

## Step 12. Initial isotope model, residual check, and outlier handling -------
srsrezmod <- glmssn(Sr8786 ~ ign_ap+sedsp+sup+agemin+accomp,
                    EstMeth = "ML",
                    ssn.object = sacssn_a,
                    CorModels = c("Exponential.tailup", "Exponential.Euclid"),
                    addfunccol = "afv_pptSr")
srsrezmod

# Residuals and boxplot
srsr.res      <- residuals(srsrezmod)
hist(srsr.res)
ObsDFres.srsr <- getSSNdata.frame(srsr.res)
box_plot      <- plot_ly(data = ObsDFres.srsr, y = ~`_resid_`, type = "box"); box_plot

# Tukey whiskers (exclusive quartiles)
x <- ObsDFres.srsr$`_resid_`; x <- x[is.finite(x)]
qtype <- c(linear=7, exclusive=6, inclusive=8)[["exclusive"]]
q     <- quantile(x, c(.25,.75), type=qtype); iqr <- diff(q)
lower <- min(x[x >= q[1] - 1.5*iqr], na.rm = TRUE)
upper <- max(x[x <= q[2] + 1.5*iqr], na.rm = TRUE)

# Remove residual outliers from Sr8786 and save
ObsDF <- getSSNdata.frame(sacssn_a)
indOutlier <- ObsDFres.srsr["_resid_"] > upper | ObsDFres.srsr["_resid_"] < lower
indOutlier[is.na(indOutlier)] <- FALSE
ObsDF[indOutlier, "Sr8786"] <- NA
ObsDF$Sr8786
sacssn_f <- putSSNdata.frame(ObsDF, sacssn_a)

## Step 13. Refit isotope model after outlier removal + quick checks ----------
srsrRRmod <- glmssn(Sr8786 ~ ign_ap+sedsp+sup+agemin+accomp,
                    EstMeth = "ML",
                    ssn.object = sacssn_f,
                    CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
                    addfunccol = "afv_pptSr")
srsrRRmod

# Visual residual check (optional)
srsr.res <- residuals(srsrRRmod)
hist(sacssn_f, "Sr8786")
ObsDFr <- getSSNdata.frame(srsr.res)
box_plot <- plot_ly(data = ObsDFr, y = ~`_resid_`, type = "box"); box_plot

## Step 14. Optional “final” ML isotope model + LOOCV visual summary ----------
srsrfmod <- glmssn(Sr8786 ~ ign_ap+sedsp+sup+agemin+accomp,
                   EstMeth = "ML",
                   ssn.object = sacssn_f,
                   CorModels = c("Epanech.tailup", "Epanech.taildown", "Spherical.Euclid" ),
                   addfunccol = "afv_pptSr")

summary(srsrfmod)
varcomp(srsrfmod)

srsrcv.out <- CrossValidationSSN(srsrfmod)
srsr.lm    <- lm(srsrcv.out[, "cv.pred"] ~ srsrfmod$sampinfo$z)
summary(srsr.lm)

# LOOCV scatter with SE-scaled points (visual only)
loocv_data <- data.frame(
  Observed  = srsrfmod$sampinfo$z,
  Predicted = srsrcv.out[, "cv.pred"],
  SE        = srsrcv.out[, "cv.se"]
)

isomodplot <- ggplot(loocv_data, aes(x = Observed, y = Predicted)) +
  geom_point(aes(size = SE)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  scale_size_continuous(range = c(1, 6)) +
  labs(
    x = expression("Observed Water " ^87 * "Sr/" ^86 * "Sr"),
    y = expression("Predicted Water " ^87 * "Sr/" ^86 * "Sr"),
    title = "LOOCV Results: Predicted vs Observed"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
print(isomodplot)




# =============================================================================
# PART 2 — Exhaustive model searches (ML predictors → REML structures)
#         Uses the SSN objects created in Part 1 (sacssn, sacssn1, sacssn2, sacssn_a, sacssn_f)
# =============================================================================

## Step 15. Non-spatial screening for ln[Sr] (nugget-only ML) -----------------
all_predictors <- c("ign_ap", "ign_bp", "ign_ip", "mtp", "pyp", "scp", "sedsp", "sup",
                    "agemax", "agemin", "agemu", "accomp", "arcp", "riftp", "x5p", "x7p")

# Generate all predictor combinations
srML_generate_predictor_combinations <- function(predictors) {
  predictor_combinations <- list()
  for (i in 1:length(predictors)) {
    predictor_combinations[[i]] <- combn(predictors, i, simplify = FALSE)
  }
  return(do.call(c, predictor_combinations))
}
srML_predictor_combinations <- srML_generate_predictor_combinations(all_predictors)

# Fit nugget-only ML models with significance screen
srML_fit_models <- function(predictor_combinations, ssn_object, significance_threshold = 0.1) {
  models <- list()
  pb <- txtProgressBar(min = 0, max = length(predictor_combinations), style = 3)
  for (i in 1:length(predictor_combinations)) {
    model_formula <- as.formula(paste("ln_sr", "~", paste(predictor_combinations[[i]], collapse = "+")))
    model_attempt <- tryCatch({
      glmssn(model_formula,
             EstMeth = "ML",
             ssn.object = ssn_object,
             CorModels = NULL,
             use.nugget = T)
    }, error = function(e) { message("Error encountered with predictors: ", paste(predictor_combinations[[i]], collapse = "_"), ". Skipping..."); NULL })
    if (!is.null(model_attempt)) {
      model_summary <- summary(model_attempt)
      p_values <- model_summary$fixed.effects.estimates[-1, "prob.t"]
      if (all(p_values < significance_threshold)) {
        models[[paste(predictor_combinations[[i]], collapse = "_")]] <- model_attempt
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(models)
}
srML_models <- srML_fit_models(srML_predictor_combinations, sacssn)

# Compare and extract summaries
srML_best_model <- find_best_model(srML_models)
srML_best_model

srML.model_list  <- as.list(srML_models)
srML.comp.models <- InfoCritCompare(srML.model_list)
r_squared        <- sapply(srML.model_list, function(srML.models) GR2(srML.models))
srML.comp.models$R_squared <- r_squared
srML.comp.models <- srML.comp.models[, c("formula","EstMethod","Variance_Components","neg2LogL",
                                         "bias","std.bias","RAV","RMSPE","std.MSPE","AIC",
                                         "cov.80","cov.90","cov.95","R_squared")]
View(srML.comp.models)

## Step 16. Spatial ML sweep for ln[Sr] (predictors search under fixed TU/TD/E)
# (Re-define srML_fit_models for spatial ML with fixed correlation structures)
srML_fit_models <- function(predictor_combinations, ssn_object, significance_threshold = 0.1) {
  models <- list()
  pb <- txtProgressBar(min = 0, max = length(predictor_combinations), style = 3)
  for (i in 1:length(predictor_combinations)) {
    model_formula <- as.formula(paste("ln_sr", "~", paste(predictor_combinations[[i]], collapse = "+")))
    model_attempt <- tryCatch({
      glmssn(model_formula,
             EstMeth = "ML",
             ssn.object = ssn_object,
             CorModels = c("Exponential.tailup", "Exponential.taildown", "Exponential.Euclid"),
             addfunccol = "afv_fasac")
    }, error = function(e) { message("Error encountered with predictors: ", paste(predictor_combinations[[i]], collapse = "_"), ". Skipping..."); NULL })
    if (!is.null(model_attempt)) {
      model_summary <- summary(model_attempt)
      p_values <- model_summary$fixed.effects.estimates[-1, "prob.t"]
      if (all(p_values < significance_threshold)) {
        models[[paste(predictor_combinations[[i]], collapse = "_")]] <- model_attempt
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(models)
}

# Run spatial ML on the cleaned object (sacssn1)
srML_predictor_combinations <- srML_generate_predictor_combinations(all_predictors)
str(srML_predictor_combinations)
srML_models <- srML_fit_models(srML_predictor_combinations, sacssn1)

srML.model_list  <- as.list(srML_models)
srML.comp.models <- InfoCritCompare(srML.model_list)
r_squared        <- sapply(srML.model_list, function(srML.models) GR2(srML.models))
srML.comp.models$R_squared <- r_squared
srML.comp.models <- srML.comp.models[, c("formula","EstMethod","Variance_Components","neg2LogL",
                                         "bias","std.bias","RAV","RMSPE","std.MSPE","AIC",
                                         "cov.80","cov.90","cov.95","R_squared")]
View(srML.comp.models)

## Step 17. Spatial REML sweep for ln[Sr] (structures search with predictors fixed)
lnsr_REML_generate_models <- function(ssn_object, response_var, predictor_vars) {
  tail_up_models  <- c("Exponential.tailup", "LinearSill.tailup", "Spherical.tailup", "Mariah.tailup", "Epanech.tailup")
  tail_down_models<- c("Exponential.taildown","LinearSill.taildown","Spherical.taildown","Mariah.taildown","Epanech.taildown")
  euclidean_models<- c("Spherical.Euclid", "Gaussian.Euclid", "Exponential.Euclid", "Cauchy.Euclid")
  model_list <- list()
  total_iters <- length(tail_up_models) * length(tail_down_models) * length(euclidean_models)
  pb <- txtProgressBar(min = 0, max = total_iters, style = 3); current_iter <- 0
  for (tail_up in tail_up_models) {
    for (tail_down in tail_down_models) {
      for (euclidean in euclidean_models) {
        model_attempt <- try({
          glmssn(as.formula(paste(response_var, "~", paste(predictor_vars, collapse = "+"))),
                 EstMeth = "REML",
                 ssn.object = ssn_object,
                 CorModels = c(tail_up, tail_down, euclidean),
                 addfunccol = "afv_fasac")
        }, silent = TRUE)
        if (!inherits(model_attempt, "try-error")) {
          model_list[[paste(tail_up, tail_down, euclidean, sep = "_")]] <- model_attempt
        } else {
          cat("Failed to generate model for:", tail_up, tail_down, euclidean, "\n")
        }
        current_iter <- current_iter + 1
        setTxtProgressBar(pb, current_iter)
      }
    }
  }
  close(pb)
  return(model_list)
}

# Run REML grid with the predictors you fixed (from ML step)
sr_REML.models <- lnsr_REML_generate_models(
  ssn_object   = sacssn2,
  response_var = "ln_sr",
  predictor_vars = c("ign_ap","ign_bp","scp","agemin","agemu","x5p","x7p")
)

sr_REML.model_list  <- as.list(sr_REML.models)
sr_REML.comp.models <- InfoCritCompare(sr_REML.model_list)
r_squared           <- sapply(sr_REML.model_list, function(sr_REML.models) GR2(sr_REML.models))
sr_REML.comp.models$R_squared <- r_squared
sr_REML.comp.models <- sr_REML.comp.models[, c("formula","EstMethod","Variance_Components","neg2LogL",
                                               "bias","std.bias","RAV","RMSPE","std.MSPE","AIC",
                                               "cov.80","cov.90","cov.95","R_squared")]
View(sr_REML.comp.models)

## Step 18. Predictor pre-screen for isotope model (VIF + LM significance) ----
# Candidate pool (same as earlier)
all_predictors <- c("ign_ap","ign_bp","ign_ip","mtp","pyp","scp","sedsp","sup",
                    "accomp","arcp","riftp","agemin","agemax","agemu","x5p","x7p")
subset_obsDF_a <- ObsDF[, c(all_predictors, "ln_sr")]

# VIF-based reduction helper
reduce_collinearity_with_response <- function(data, response, threshold = 10) {
  clean_data <- data[!is.na(data[[response]]), ]
  predictors <- setdiff(names(clean_data), response)
  repeat {
    fit_formula <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))
    fit <- lm(fit_formula, data = clean_data)
    vif_values <- car::vif(fit)
    max_vif <- max(vif_values, na.rm = TRUE)
    if (max_vif < threshold) break
    drop_predictor <- names(which.max(vif_values))
    predictors <- setdiff(predictors, drop_predictor)
  }
  return(predictors)
}

cleaned_data   <- subset_obsDF_a[!is.na(subset_obsDF_a$ln_sr), ] %>% st_drop_geometry()
vif_predictors <- reduce_collinearity_with_response(data = cleaned_data, "ln_sr"); vif_predictors

# (Your chosen screen)
vif_predictors <- c("ign_ap","ign_ip","mtp","scp","sedsp","sup","agemin","accomp","riftp","x7p")

# Significance-only screen on LM
assess_significance <- function(data_frame, response_variable, predictors) {
  data <- getSSNdata.frame(data_frame)
  significant_formulas <- list()
  n <- length(predictors)
  total_models <- sum(2^(1:n) - 1)
  counter <- 0
  for (i in 1:n) {
    combi <- combn(predictors, i, simplify = FALSE)
    for (j in 1:length(combi)) {
      counter <- counter + 1
      print(paste("Processing model", counter, "of", total_models))
      formula_str <- as.formula(paste(response_variable, "~", paste(combi[[j]], collapse = "+")))
      model <- lm(formula_str, data = data)
      p_values <- summary(model)$coefficients[,4]
      if (all(p_values[-1] < 0.05)) {
        significant_formulas <- append(significant_formulas, list(formula_str))
      }
    }
  }
  return(significant_formulas)
}

srsr_assesment <- assess_significance(sacssn_a, "Sr8786", vif_predictors)
srsr_aic_values <- sapply(srsr_assesment, function(formula) {
  model <- lm(formula, data = getSSNdata.frame(sacssn_a))
  return(AIC(model))
})
srsr_best_formula_index <- which.min(srsr_aic_values)
srsr_best_formula <- srsr_assesment[[srsr_best_formula_index]]
srsr_lmbest_model <- lm(srsr_best_formula, data = getSSNdata.frame(sacssn_a))
summary(srsr_lmbest_model)

# Your final screened predictor pool
vif_signif_predictors <- c("ign_ap","sedsp","sup","agemin","accomp","riftp")

## Step 19. Spatial ML sweep for isotope model (predictors under fixed TU) -----
srsrML_generate_predictor_combinations <- function(predictors) {
  predictor_combinations <- list()
  for (i in 1:length(predictors)) {
    predictor_combinations[[i]] <- combn(predictors, i, simplify = FALSE)
  }
  return(do.call(c, predictor_combinations))
}
srsr_predictor_combinations <- srsrML_generate_predictor_combinations(vif_signif_predictors)

srsrML_fit_models <- function(predictor_combinations, ssn_object, significance_threshold = 0.1) {
  models <- list()
  pb <- txtProgressBar(min = 0, max = length(predictor_combinations), style = 3)
  for (i in 1:length(predictor_combinations)) {
    model_formula <- as.formula(paste("Sr8786", "~", paste(predictor_combinations[[i]], collapse = "+")))
    model_attempt <- try({
      model <- glmssn(model_formula,
                      EstMeth = "ML",
                      ssn.object = ssn_object,
                      CorModels = "Exponential.tailup",
                      addfunccol = "afv_pptSr")
      model_summary <- summary(model)
      p_values <- model_summary$fixed.effects.estimates[-1, "prob.t"]
      if (all(p_values < significance_threshold)) {
        models[[paste(predictor_combinations[[i]], collapse = "_")]] <- model
      }
    }, silent = TRUE)
    if (inherits(model_attempt, "try-error")) {
      cat("Error encountered with predictors:", predictor_combinations[[i]], "\n")
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(models)
}
srsrML.models <- srsrML_fit_models(srsr_predictor_combinations, sacssn_a)

srsrML.model_list  <- as.list(srsrML.models)
srsrML.comp.models <- InfoCritCompare(srsrML.model_list)
r_squared          <- sapply(srsrML.model_list, function(m) GR2(m))
srsrML.comp.models$R_squared <- r_squared
srsrML.comp.models <- srsrML.comp.models[, c("formula","EstMethod","Variance_Components","neg2LogL",
                                             "bias","std.bias","RAV","RMSPE","std.MSPE","AIC",
                                             "cov.80","cov.90","cov.95","R_squared")]
View(srsrML.comp.models)

## Step 20. Spatial REML sweep for isotope model (structures with predictors fixed)
srsrREML_generate_models <- function(ssn_object, response_var, predictor_vars) {
  tail_up_models   <- c("Exponential.tailup","LinearSill.tailup","Spherical.tailup","Mariah.tailup","Epanech.tailup")
  tail_down_models <- c("Exponential.taildown","LinearSill.taildown","Spherical.taildown","Mariah.taildown","Epanech.taildown")
  euclidean_models <- c("Spherical.Euclid","Gaussian.Euclid","Exponential.Euclid","Cauchy.Euclid")
  model_list <- list()
  total_iterations <- length(tail_up_models)*length(tail_down_models)*length(euclidean_models)
  pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
  iteration_count <- 0
  for (tail_up in tail_up_models) {
    for (tail_down in tail_down_models) {
      for (euclidean in euclidean_models) {
        iteration_count <- iteration_count + 1
        setTxtProgressBar(pb, iteration_count)
        model_attempt <- tryCatch({
          glmssn(as.formula(paste(response_var, "~", paste(predictor_vars, collapse = "+"))),
                 EstMeth = "REML",
                 ssn.object = ssn_object,
                 CorModels = c(tail_up, tail_down, euclidean),
                 addfunccol = "afv_pptSr")
        }, error = function(e) { return(paste("Error for CorModels:", tail_up, tail_down, euclidean, "\nError message:", e$message)) })
        if (!inherits(model_attempt, "character")) {
          model_list[[paste(tail_up, tail_down, euclidean, sep = "_")]] <- model_attempt
        } else {
          cat(model_attempt, "\n")
        }
      }
    }
  }
  close(pb)
  return(model_list)
}

srsrREML.models <- srsrREML_generate_models(
  ssn_object   = sacssn_f,
  response_var = "Sr8786",
  predictor_vars = c("ign_ap","sedsp","sup","agemin","accomp")  # from LM screen best set
)

srsrREML.model_list  <- as.list(srsrREML.models)
srsrREML.comp_model  <- InfoCritCompare(srsrREML.model_list)
r_squared            <- sapply(srsrREML.model_list, function(m) GR2(m))
srsrREML.comp_model$R_squared <- r_squared
srsrREML.comp_model <- srsrREML.comp_model[, c("formula","EstMethod","Variance_Components","neg2LogL",
                                               "bias","std.bias","RAV","RMSPE","std.MSPE","AIC",
                                               "cov.80","cov.90","cov.95","R_squared")]
View(srsrREML.comp_model)
  