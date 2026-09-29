#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Export Nature Climate Change "Source Data" files.
#
# Writes the numerical data underlying every plotted panel of the main and
# Extended Data figures as CSV files (one folder + one zip per figure) under
# source_data/.
#
# This script only *reads* the data, fitted models and precomputed posterior
# draws used by the figure code; it does not refit models or change any
# modelling choice. Data preprocessing is copied verbatim from
# book/code_figures.qmd, book/SI_individual_mismatch_trajectories.qmd,
# book/estimated-fleds.qmd, book/sensitivity_drop_enso.qmd and
# book/1-bloom-onset.qmd.
#
# Model-prediction panels are summarised exactly as ggdist::stat_lineribbon
# draws them: per x value, posterior median and 95% quantile interval
# (ggdist::median_qi, .width = 0.95), computed on the same draw data frame
# passed to the plot. A ggplot_build() comparison against the actual
# stat_lineribbon layer (with the figure's axis limits) is run for every
# model panel and written to source_data/_export_checks.csv.
#
# Usage (from the repository root, or from book/):
#   Rscript Rscripts/export_source_data.R [repo_root] [output_dir] [limit_mode]
#   limit_mode = "scale" (default; matches current figure code) or "coord"
#   (use after the figures switch to coord_cartesian zooming).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidyr)
  library(tidybayes)
  library(ggdist)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1) args[[1]] else if (dir.exists("Rdata")) "." else ".."
repo_root <- normalizePath(repo_root, winslash = "/")
out_root  <- if (length(args) >= 2) args[[2]] else file.path(repo_root, "source_data")
# How the figure code applies axis limits when rebuilding the rendered layers:
#   "scale" = scale_[xy]_continuous(limits = ...) as in the current figure code
#             (ggplot2 drops out-of-range draws BEFORE stat_lineribbon summarises);
#   "coord" = coord_cartesian(xlim, ylim) zoom (no draws dropped).
# Switch to "coord" (3rd argument or env var SOURCE_DATA_LIMIT_MODE) once the
# figures are changed to zoom rather than censor; the full-posterior columns
# are unaffected by this setting.
limit_mode <- if (length(args) >= 3) args[[3]] else Sys.getenv("SOURCE_DATA_LIMIT_MODE", "scale")
stopifnot(limit_mode %in% c("scale", "coord"))
message("Axis-limit mode for rendered-layer rebuild: ", limit_mode)

data_dir  <- file.path(repo_root, "data")
model_dir <- file.path(repo_root, "Rdata")
draw_dir  <- file.path(model_dir, "precomputed_draws")
stopifnot(dir.exists(data_dir), dir.exists(model_dir), dir.exists(draw_dir))

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
read_precomputed_draw <- function(file) readRDS(file.path(draw_dir, file))

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
file_log  <- list()
check_log <- list()

write_panel <- function(df, fig, file) {
  dir.create(file.path(out_root, fig), showWarnings = FALSE, recursive = TRUE)
  path <- file.path(out_root, fig, file)
  utils::write.csv(as.data.frame(df), path, row.names = FALSE, na = "")
  file_log[[length(file_log) + 1]] <<- data.frame(
    figure = fig, file = file, rows = nrow(df), cols = ncol(df)
  )
  message("  wrote ", fig, "/", file, " (", nrow(df), " rows)")
  invisible(path)
}

write_readme <- function(fig, lines) {
  dir.create(file.path(out_root, fig), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, file.path(out_root, fig, "README.txt"), useBytes = TRUE)
}

within_lim <- function(v, lim) {
  if (is.null(lim)) return(rep(TRUE, length(v)))
  !is.na(v) & v >= lim[1] & v <= lim[2]
}

# Posterior median + 95% quantile interval per x (and optional group),
# i.e. what stat_lineribbon(.width = 0.95) draws.
summarise_ribbon <- function(draws, x, y, group = NULL, x_name = x, xlim = NULL) {
  d <- draws %>%
    ungroup() %>%
    transmute(.x = .data[[x]], .value = .data[[y]],
              .g = if (is.null(group)) NA_character_ else as.character(.data[[group]]))
  s <- d %>%
    group_by(.g, .x) %>%
    median_qi(.value, .width = 0.95) %>%
    ungroup() %>%
    arrange(.g, .x)
  out <- tibble(!!x_name := s$.x)
  if (!is.null(group)) out <- tibble(!!group := s$.g, out)
  out %>%
    mutate(
      posterior_median = s$.value,
      lower_95 = s$.lower,
      upper_95 = s$.upper,
      within_x_axis_limits = within_lim(.data[[x_name]], xlim)
    )
}

# Rebuild the stat_lineribbon layer with the figure's axis limits and compare
# the drawn line/ribbon to the exported summary.
check_ribbon <- function(label, draws, x, y, export, x_name = x, group = NULL,
                         xlim = NULL, ylim = NULL) {
  d <- ungroup(draws)
  map <- if (is.null(group)) aes(x = .data[[x]], y = .data[[y]]) else
    aes(x = .data[[x]], y = .data[[y]], fill = .data[[group]], group = .data[[group]])
  p <- ggplot() + stat_lineribbon(data = d, mapping = map, .width = 0.95)
  if (limit_mode == "scale") {
    if (!is.null(xlim)) p <- p + scale_x_continuous(limits = xlim)
    if (!is.null(ylim)) p <- p + scale_y_continuous(limits = ylim)
  } else {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  n_oob <- if (limit_mode == "coord") 0L else sum(!within_lim(d[[y]], ylim) & !is.na(d[[y]])) +
    sum(!within_lim(d[[x]], xlim) & !is.na(d[[x]]))
  built <- suppressWarnings(ggplot_build(p))$data[[1]]
  built <- built[built$.width == 0.95 | is.na(built$.width), ]
  ex <- export %>% filter(within_x_axis_limits)
  if (!is.null(group)) {
    lev <- sort(unique(as.character(d[[group]])))
    built$grp_label <- lev[built$group]
    ex$grp_label <- as.character(ex[[group]])
  } else {
    built$grp_label <- "all"; ex$grp_label <- "all"
  }
  j <- inner_join(
    ex %>% transmute(grp_label, xk = round(.data[[x_name]], 9),
                     posterior_median, lower_95, upper_95),
    built %>% transmute(grp_label, xk = round(x, 9), y, ymin, ymax),
    by = c("grp_label", "xk")
  )
  # Draws per x that survive the axis-limit censoring (what the stat actually used)
  kept <- d %>%
    mutate(in_lim = if (limit_mode == "coord") !is.na(.data[[y]]) else within_lim(.data[[y]], ylim)) %>%
    group_by(grp_label = if (is.null(group)) "all" else as.character(.data[[group]]),
             xk = round(.data[[x]], 9)) %>%
    summarise(n_draws_within_axis_limits = sum(in_lim), .groups = "drop")
  res <- data.frame(
    panel = label,
    exported_rows_in_x_limits = nrow(ex),
    drawn_rows = nrow(built),
    matched_rows = nrow(j),
    draws_outside_axis_limits = n_oob,
    max_abs_diff_full_vs_drawn_median = max(abs(j$posterior_median - j$y), na.rm = TRUE),
    max_abs_diff_full_vs_drawn_lower = max(abs(j$lower_95 - j$ymin), na.rm = TRUE),
    max_abs_diff_full_vs_drawn_upper = max(abs(j$upper_95 - j$ymax), na.rm = TRUE),
    n_drawn_NA = sum(is.na(j$y) | is.na(j$ymin) | is.na(j$ymax))
  )
  check_log[[length(check_log) + 1]] <<- res
  message(sprintf("  check %-8s matched %d/%d, max|diff| median=%.3g lower=%.3g upper=%.3g, OOB draws=%d",
                  label, res$matched_rows, res$exported_rows_in_x_limits,
                  res$max_abs_diff_full_vs_drawn_median, res$max_abs_diff_full_vs_drawn_lower,
                  res$max_abs_diff_full_vs_drawn_upper, n_oob))
  # Return the export augmented with the values actually rendered by
  # stat_lineribbon under the figure's axis limits.
  export %>%
    mutate(grp_label = if (is.null(group)) "all" else as.character(.data[[group]]),
           xk = round(.data[[x_name]], 9)) %>%
    left_join(built %>% transmute(grp_label, xk = round(x, 9),
                                  drawn_median = y, drawn_lower_95 = ymin, drawn_upper_95 = ymax),
              by = c("grp_label", "xk")) %>%
    left_join(kept, by = c("grp_label", "xk")) %>%
    select(-grp_label, -xk)
}


flag_points <- function(df, x, y, xlim = NULL, ylim = NULL) {
  df %>% mutate(within_axis_limits = within_lim(.data[[x]], xlim) & within_lim(.data[[y]], ylim))
}

# -----------------------------------------------------------------------------
# Data (preprocessing copied from book/code_figures.qmd)
# -----------------------------------------------------------------------------
message("Loading and preprocessing data")
bloom_dat_raw  <- read.csv(file.path(data_dir, "bloom_dates.csv"))
female_dat_raw <- read.csv(file.path(data_dir, "Lay_date_mismatches_250702.csv"))
male_dat_raw   <- female_dat_raw
stopifnot("YEAR" %in% names(female_dat_raw))

bloom_dat <- bloom_dat_raw %>%
  mutate(
    bloom_start_thr5_as_date = as.Date(bloom_start_thr5, format = "%d/%m/%Y"),
    winter_start_date = as.Date(paste0(season_year - 1, "-12-21")),
    delta5 = as.numeric(difftime(bloom_start_thr5_as_date, winter_start_date, units = "days")),
    season_year = as.numeric(season_year),
    TIME = as.numeric(scale(season_year, center = TRUE, scale = FALSE))
  )

prep_mismatch <- function(raw, sex) raw %>%
  filter(SEX == sex) %>%
  mutate(
    ori_age = as.numeric(AGE),
    ori_time = as.numeric(YEAR),
    TIME = scale(YEAR, center = TRUE, scale = FALSE),
    AGE = scale(as.numeric(AGE), center = TRUE, scale = FALSE),
    AFR = scale(as.numeric(AFR), center = TRUE, scale = FALSE),
    LONGEVITY = scale(as.numeric(LONGEVITY), center = TRUE, scale = FALSE),
    YEAR = as.factor(YEAR),
    RING = as.factor(RING)
  )

prep_fledging <- function(raw, sex) raw %>%
  filter(SEX == sex) %>%
  mutate(
    ori_age = as.numeric(AGE),
    ori_time = as.numeric(YEAR),
    TIME = scale(YEAR, center = TRUE, scale = FALSE),
    ZX5MISMATCH = scale(as.numeric(X5MISMATCH), center = TRUE, scale = FALSE),
    AGE = scale(as.numeric(AGE), center = TRUE, scale = FALSE),
    AFR = scale(as.numeric(AFR), center = TRUE, scale = FALSE),
    LONGEVITY = scale(as.numeric(LONGEVITY), center = TRUE, scale = FALSE),
    FLEDS = factor(FLEDS, levels = c(0, 1, 2, 3), ordered = TRUE)
  )

dat_female_mismatch <- prep_mismatch(female_dat_raw, 1)
dat_male_mismatch   <- prep_mismatch(male_dat_raw, 0)
dat_female_fledging <- prep_fledging(female_dat_raw, 1)
dat_male_fledging   <- prep_fledging(male_dat_raw, 0)

# -----------------------------------------------------------------------------
# Absolute model residuals (the sigma-panel points), exactly as in the figure
# code: residuals(model, summary = TRUE)[, "Estimate"]. Models are loaded one
# at a time to limit memory use.
# -----------------------------------------------------------------------------
message("Computing absolute residuals for sigma panels")
bloom_model <- readRDS(file.path(model_dir, "model_delta5.rds"))
res_bloom <- residuals(bloom_model, summary = TRUE)[, "Estimate"]
stopifnot(length(res_bloom) == nrow(bloom_dat))
dat_with_residuals_bloom <- bloom_dat %>% mutate(raw_sigma = abs(res_bloom))
rm(bloom_model); invisible(gc())

abs_resid <- function(model_file, dat) {
  m <- readRDS(file.path(model_dir, model_file))
  r <- residuals(m, summary = TRUE)[, "Estimate"]
  stopifnot(length(r) == nrow(dat))
  rm(m); invisible(gc())
  dat %>% mutate(raw_sigma = abs(r))
}
dat_with_residuals_female <- abs_resid("female_fit_m11.rds", dat_female_mismatch)
dat_with_residuals_male   <- abs_resid("male_fit_m11.rds", dat_male_mismatch)

# -----------------------------------------------------------------------------
# Observed-point summaries (as in book/code_figures.qmd)
# -----------------------------------------------------------------------------
count_mismatch <- function(dat, xvar) dat %>%
  group_by(.data[[xvar]], X5MISMATCH) %>% summarize(n = n(), .groups = "drop")
count_sigma <- function(dat, xvar) dat %>%
  group_by(.data[[xvar]], raw_sigma) %>% summarize(n = n(), .groups = "drop")

success_props <- function(dat, xvar) dat %>%
  mutate(SUCCESS = ifelse(FLEDS == "0", "Failure (0)", "Success (1+)")) %>%
  count(.data[[xvar]], SUCCESS) %>%
  rename(!!xvar := 1) %>%
  tidyr::complete(!!rlang::sym(xvar), SUCCESS = c("Failure (0)", "Success (1+)"),
                  fill = list(n = 0)) %>%
  group_by(.data[[xvar]]) %>%
  mutate(total_n_year = sum(n),
         proportion = ifelse(total_n_year > 0, n / total_n_year, 0)) %>%
  filter(SUCCESS == "Success (1+)") %>%
  ungroup()

grouped_props <- function(dat) dat %>%
  filter(!is.na(X5MISMATCH)) %>%
  mutate(fled_group = case_when(
    FLEDS == "0" ~ "0",
    FLEDS %in% c("1", "2") ~ "1-2",
    FLEDS == "3" ~ "3",
    TRUE ~ "Other")) %>%
  filter(fled_group %in% c("0", "1-2", "3")) %>%
  count(X5MISMATCH, fled_group) %>%
  tidyr::complete(X5MISMATCH, fled_group = c("0", "1-2", "3"), fill = list(n = 0)) %>%
  group_by(X5MISMATCH) %>%
  mutate(total_n_year = sum(n),
         proportion = ifelse(total_n_year > 0, n / total_n_year, 0)) %>%
  ungroup()

# Axis limits used in the figure code (points/lines outside are not drawn).
LIM_YEAR <- NULL
LIM_AGE  <- c(1, 23)
LIM_MM   <- c(-90, 210)
LIM_SIG  <- c(0, 150)
LIM_P    <- c(0, 1)
LIM_BLOOM_MEAN  <- c(0, 60)   # Fig. 2A y-axis
LIM_BLOOM_SIGMA <- c(0, 100)  # Fig. 2B y-axis
LIM_ZETA        <- c(0, 0.5)  # Fig. 2F / ED Fig. 4D y-axis

# =============================================================================
# Temporal panels (Fig. 2 females; Extended Data Fig. 4 males)
# =============================================================================
export_temporal <- function(sex, fig, labs) {
  sfx  <- if (sex == "female") "f" else "m"
  dmis <- if (sex == "female") dat_female_mismatch else dat_male_mismatch
  dres <- if (sex == "female") dat_with_residuals_female else dat_with_residuals_male
  dfl  <- if (sex == "female") dat_female_fledging else dat_male_fledging

  mt <- read_precomputed_draw(paste0("mismatch_time_draws", sfx, ".rds"))
  s <- summarise_ribbon(mt, "actual_time", ".epred", x_name = "year")
  s <- check_ribbon(paste(fig, labs["mm"]), mt, "actual_time", ".epred", s, "year", ylim = LIM_MM)
  write_panel(count_mismatch(dmis, "ori_time") %>%
                rename(year = ori_time, mismatch_days = X5MISMATCH, n_records = n) %>%
                flag_points("year", "mismatch_days", ylim = LIM_MM),
              fig, paste0(labs["mm"], "_points.csv"))
  write_panel(s, fig, paste0(labs["mm"], "_mismatch_mean.csv"))

  s <- summarise_ribbon(mt, "actual_time", "sigma", x_name = "year")
  s <- check_ribbon(paste(fig, labs["sd"]), mt, "actual_time", "sigma", s, "year", ylim = LIM_SIG)
  write_panel(count_sigma(dres, "ori_time") %>%
                rename(year = ori_time, abs_residual_days = raw_sigma, n_records = n) %>%
                flag_points("year", "abs_residual_days", ylim = LIM_SIG),
              fig, paste0(labs["sd"], "_points.csv"))
  write_panel(s, fig, paste0(labs["sd"], "_mismatch_sigma.csv"))
  rm(mt); invisible(gc())

  fs <- read_precomputed_draw(paste0("fledging_success_draws", sfx, ".rds"))
  s <- summarise_ribbon(fs, "actual_time", "prob_success", x_name = "year")
  s <- check_ribbon(paste(fig, labs["ns"]), fs, "actual_time", "prob_success", s, "year", ylim = LIM_P)
  write_panel(success_props(dfl, "ori_time") %>%
                transmute(year = ori_time, n_successful = n, n_records_year = total_n_year,
                          proportion_successful = proportion) %>%
                flag_points("year", "proportion_successful", ylim = LIM_P),
              fig, paste0(labs["ns"], "_points.csv"))
  write_panel(s, fig, paste0(labs["ns"], "_nest_success.csv"))
  rm(fs); invisible(gc())

  # Zeta panel: the figure passes the category-level draw frame (4 identical
  # disc values per draw); summarised on the same frame so values match.
  ft <- read_precomputed_draw(paste0("fledging_time_draws", sfx, ".rds"))
  s <- summarise_ribbon(ft, "actual_time", "disc", x_name = "year")
  s <- check_ribbon(paste(fig, labs["zeta"]), ft, "actual_time", "disc", s, "year", ylim = LIM_ZETA)
  write_panel(s, fig, paste0(labs["zeta"], "_zeta.csv"))
  rm(ft); invisible(gc())
}

# ---- Fig. 2 -----------------------------------------------------------------
message("Fig. 2")
fig <- "Fig2"
bt <- read_precomputed_draw("bloom_time_draws.rds")
s <- summarise_ribbon(bt, "actual_year", ".epred", x_name = "year")
s <- check_ribbon("Fig2 Fig2A", bt, "actual_year", ".epred", s, "year", ylim = LIM_BLOOM_MEAN)
write_panel(s, fig, "Fig2A_bloom_onset_mean.csv")
write_panel(bloom_dat %>% transmute(year = season_year, bloom_onset_days = delta5) %>%
              flag_points("year", "bloom_onset_days", ylim = LIM_BLOOM_MEAN),
            fig, "Fig2A_points.csv")
s <- summarise_ribbon(bt, "actual_year", "sigma", x_name = "year")
s <- check_ribbon("Fig2 Fig2B", bt, "actual_year", "sigma", s, "year", ylim = LIM_BLOOM_SIGMA)
write_panel(s, fig, "Fig2B_bloom_timing_sigma.csv")
write_panel(dat_with_residuals_bloom %>% transmute(year = season_year, abs_residual_days = raw_sigma) %>%
              flag_points("year", "abs_residual_days", ylim = LIM_BLOOM_SIGMA),
            fig, "Fig2B_points.csv")
rm(bt); invisible(gc())

export_temporal("female", "Fig2", c(mm = "Fig2C", sd = "Fig2D", ns = "Fig2E", zeta = "Fig2F"))

# ---- Extended Data Fig. 4 ---------------------------------------------------
message("Extended Data Fig. 4")
export_temporal("male", "EDFig4", c(mm = "EDFig4A", sd = "EDFig4B", ns = "EDFig4C", zeta = "EDFig4D"))

# =============================================================================
# Age / mismatch panels (Fig. 3 females; Extended Data Fig. 5 males)
# =============================================================================
export_age <- function(sex, fig, labs) {
  sfx_a <- if (sex == "female") "female" else "male"
  sfx   <- if (sex == "female") "f" else "m"
  dmis  <- if (sex == "female") dat_female_mismatch else dat_male_mismatch
  dres  <- if (sex == "female") dat_with_residuals_female else dat_with_residuals_male
  dfl   <- if (sex == "female") dat_female_fledging else dat_male_fledging

  ad <- read_precomputed_draw(paste0("age_draws_", sfx_a, ".rds"))
  s <- summarise_ribbon(ad, "actual_age", ".epred", x_name = "age_years", xlim = LIM_AGE)
  s <- check_ribbon(paste(fig, labs["mm"]), ad, "actual_age", ".epred", s, "age_years",
               xlim = LIM_AGE, ylim = LIM_MM)
  write_panel(count_mismatch(dmis, "ori_age") %>%
                rename(age_years = ori_age, mismatch_days = X5MISMATCH, n_records = n) %>%
                flag_points("age_years", "mismatch_days", LIM_AGE, LIM_MM),
              fig, paste0(labs["mm"], "_points.csv"))
  write_panel(s, fig, paste0(labs["mm"], "_mismatch_mean_age.csv"))

  s <- summarise_ribbon(ad, "actual_age", "sigma", x_name = "age_years", xlim = LIM_AGE)
  s <- check_ribbon(paste(fig, labs["sd"]), ad, "actual_age", "sigma", s, "age_years",
               xlim = LIM_AGE, ylim = LIM_SIG)
  write_panel(count_sigma(dres, "ori_age") %>%
                rename(age_years = ori_age, abs_residual_days = raw_sigma, n_records = n) %>%
                flag_points("age_years", "abs_residual_days", LIM_AGE, LIM_SIG),
              fig, paste0(labs["sd"], "_points.csv"))
  write_panel(s, fig, paste0(labs["sd"], "_mismatch_sigma_age.csv"))
  rm(ad); invisible(gc())

  fa <- read_precomputed_draw(paste0("fledging_success_draws_age", sfx, ".rds"))
  s <- summarise_ribbon(fa, "actual_age", "prob_success", x_name = "age_years")
  s <- check_ribbon(paste(fig, labs["ns"]), fa, "actual_age", "prob_success", s, "age_years", ylim = LIM_P)
  write_panel(success_props(dfl, "ori_age") %>%
                transmute(age_years = ori_age, n_successful = n, n_records_age = total_n_year,
                          proportion_successful = proportion) %>%
                flag_points("age_years", "proportion_successful", ylim = LIM_P),
              fig, paste0(labs["ns"], "_points.csv"))
  write_panel(s, fig, paste0(labs["ns"], "_nest_success_age.csv"))
  rm(fa); invisible(gc())

  fg <- read_precomputed_draw(paste0("fledging_grouped_draws", sfx, ".rds"))
  s <- summarise_ribbon(fg, "actual_mismatch", "prob_val", group = "fled_group",
                        x_name = "mismatch_days", xlim = LIM_MM) %>%
    rename(fledgling_category = fled_group)
  s <- check_ribbon(paste(fig, labs["cat"]), fg, "actual_mismatch", "prob_val",
               s %>% rename(fled_group = fledgling_category), "mismatch_days",
               group = "fled_group", xlim = LIM_MM, ylim = LIM_P) %>%
    rename(fledgling_category = fled_group)
  write_panel(grouped_props(dfl) %>%
                transmute(mismatch_days = X5MISMATCH, fledgling_category = fled_group,
                          n_records_category = n, n_records_mismatch_value = total_n_year,
                          proportion = proportion) %>%
                flag_points("mismatch_days", "proportion", LIM_MM, LIM_P),
              fig, paste0(labs["cat"], "_points.csv"))
  write_panel(s, fig, paste0(labs["cat"], "_fledging_category_probs.csv"))
  rm(fg); invisible(gc())
}

message("Fig. 3")
export_age("female", "Fig3", c(mm = "Fig3A", sd = "Fig3B", ns = "Fig3C", cat = "Fig3D"))
message("Extended Data Fig. 5")
export_age("male", "EDFig5", c(mm = "EDFig5A", sd = "EDFig5B", ns = "EDFig5C", cat = "EDFig5D"))

# =============================================================================
# Extended Data Fig. 2: individual mismatch trajectories
# (book/SI_individual_mismatch_trajectories.qmd)
# =============================================================================
message("Extended Data Fig. 2")
fig <- "EDFig2"
for (sx in c("female", "male")) {
  pn <- if (sx == "female") "EDFig2A" else "EDFig2B"
  ad <- read_precomputed_draw(paste0("age_draws_", sx, ".rds"))
  s <- summarise_ribbon(ad, "actual_age", ".epred", x_name = "age_years", xlim = LIM_AGE)
  s <- check_ribbon(paste(fig, pn), ad, "actual_age", ".epred", s, "age_years",
               xlim = LIM_AGE, ylim = LIM_MM)
  rm(ad); invisible(gc())
  write_panel(s, fig, paste0(pn, "_population_", sx, ".csv"))
  ind <- read_precomputed_draw(paste0("individual_predictions_", sx, "_si.rds")) %>%
    ungroup() %>%
    transmute(ring = as.character(RING), age_years = actual_age,
              posterior_median_mismatch_days = .value) %>%
    arrange(ring, age_years) %>%
    flag_points("age_years", "posterior_median_mismatch_days", LIM_AGE, LIM_MM)
  stopifnot(n_distinct(ind$ring) == 50)
  write_panel(ind, fig, paste0(pn, "_individuals_", sx, ".csv"))
}

# =============================================================================
# Extended Data Fig. 3: expected fledglings, 2003 vs 2023 (book/estimated-fleds.qmd)
# =============================================================================
message("Extended Data Fig. 3")
fig <- "EDFig3"
dat_female_fl3 <- female_dat_raw %>% filter(SEX == 1) %>%
  mutate(ori_time = as.numeric(YEAR), TIME = scale(YEAR, center = TRUE, scale = FALSE)[, 1])
dat_male_fl3 <- male_dat_raw %>% filter(SEX == 0) %>%
  mutate(ori_time = as.numeric(YEAR), TIME = scale(YEAR, center = TRUE, scale = FALSE)[, 1])

compute_expected_fleds_endpoints <- function(fit, dat) {
  mean_year <- mean(dat$ori_time)
  draws <- fit %>%
    epred_draws(
      newdata = expand.grid(TIME = c(min(dat$TIME), max(dat$TIME)),
                            ZX5MISMATCH = 0, AGE = 0, AFR = 0, LONGEVITY = 0),
      re_formula = NA
    ) %>%
    mutate(actual_time = round(TIME + mean_year),
           cat_num = as.numeric(as.character(.category)))
  draws %>%
    group_by(.draw, actual_time) %>%
    summarise(expected_fleds = sum(cat_num * .epred), .groups = "drop")
}

fit <- readRDS(file.path(model_dir, "female_fit_Q17.rds"))
endpoints_female <- compute_expected_fleds_endpoints(fit, dat_female_fl3) %>% mutate(sex = "Females")
rm(fit); invisible(gc())
fit <- readRDS(file.path(model_dir, "male_fit_Q17.rds"))
endpoints_male <- compute_expected_fleds_endpoints(fit, dat_male_fl3) %>% mutate(sex = "Males")
rm(fit); invisible(gc())

endpoints_all <- bind_rows(endpoints_female, endpoints_male)
ed3_draws <- endpoints_all %>%
  transmute(sex, year = actual_time, draw = .draw, expected_fledglings = expected_fleds) %>%
  arrange(sex, year, draw)
write_panel(ed3_draws, fig, "EDFig3_expected_fledglings_draws.csv")

ed3_summary <- ed3_draws %>%
  group_by(sex, year) %>%
  summarise(n_draws = n(),
            posterior_mean = mean(expected_fledglings),
            posterior_median = median(expected_fledglings),
            lower_2.5 = quantile(expected_fledglings, 0.025, names = FALSE),
            upper_97.5 = quantile(expected_fledglings, 0.975, names = FALSE),
            .groups = "drop")
write_panel(ed3_summary, fig, "EDFig3_summary.csv")

ed3_pct <- ed3_draws %>%
  group_by(sex, draw) %>%
  mutate(period = if_else(year == min(year), "first", "last")) %>%
  ungroup() %>%
  select(sex, draw, period, expected_fledglings) %>%
  pivot_wider(names_from = period, values_from = expected_fledglings) %>%
  mutate(pct_change = 100 * (last - first) / first) %>%
  group_by(sex) %>%
  summarise(median_pct_change = median(pct_change),
            lower_2.5 = quantile(pct_change, 0.025, names = FALSE),
            upper_97.5 = quantile(pct_change, 0.975, names = FALSE),
            .groups = "drop")
write_panel(ed3_pct, fig, "EDFig3_pct_change_2003_2023.csv")

# Check the violin overlay (stat_summary mean + 95% quantiles) against ggplot_build.
p_ed3 <- ggplot(endpoints_all %>% mutate(year = factor(actual_time)),
                aes(x = year, y = expected_fleds)) +
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = quantile(x, 0.025),
                                                 ymax = quantile(x, 0.975)),
               geom = "errorbar") +
  facet_wrap(~sex)
b <- ggplot_build(p_ed3)$data[[1]]
b_cmp <- b %>% arrange(PANEL, x)
e_cmp <- ed3_summary %>% arrange(sex, year)
check_log[[length(check_log) + 1]] <- data.frame(
  panel = "EDFig3 errorbars", exported_rows_in_x_limits = nrow(e_cmp), drawn_rows = nrow(b_cmp),
  matched_rows = nrow(b_cmp), draws_outside_axis_limits = 0,
  max_abs_diff_full_vs_drawn_median = max(abs(b_cmp$y - e_cmp$posterior_mean)),
  max_abs_diff_full_vs_drawn_lower = max(abs(b_cmp$ymin - e_cmp$lower_2.5)),
  max_abs_diff_full_vs_drawn_upper = max(abs(b_cmp$ymax - e_cmp$upper_97.5)),
  n_drawn_NA = 0
)

# =============================================================================
# Extended Data Fig. 6: breeding records by fledgling number (book/sensitivity_drop_enso.qmd)
# =============================================================================
message("Extended Data Fig. 6")
fig <- "EDFig6"
ed6 <- female_dat_raw %>%
  filter(SEX %in% c(0, 1), !is.na(YEAR), !is.na(FLEDS)) %>%
  mutate(Year = as.integer(YEAR),
         Sex_Label = factor(if_else(SEX == 0, "Males", "Females"), levels = c("Females", "Males")),
         Fledglings = factor(FLEDS, levels = c(0, 1, 2, 3),
                             labels = c("0 fledglings", "1 fledgling", "2 fledglings", "3 fledglings"))) %>%
  count(Sex_Label, Year, Fledglings, name = "records") %>%
  transmute(sex = as.character(Sex_Label), year = Year,
            fledglings = as.integer(sub(" .*", "", as.character(Fledglings))),
            breeding_records = records)
write_panel(ed6, fig, "EDFig6_breeding_records_by_fledglings.csv")
write_panel(data.frame(shaded_period = "ENSO 2015-2016", xmin_year = 2014.5, xmax_year = 2016.5),
            fig, "EDFig6_shaded_period.csv")

# =============================================================================
# Extended Data Fig. 1: study-site map geometry (book/1-bloom-onset.qmd)
# =============================================================================
message("Extended Data Fig. 1")
fig <- "EDFig1"
marker_lat_deg <- 21 + 50/60 + 59.0/3600
marker_lon_deg <- -(105 + 52/60 + 54.0/3600)
write_panel(data.frame(feature = "Isla Isabel (colony marker)",
                       latitude_deg = marker_lat_deg, longitude_deg = marker_lon_deg,
                       latitude_dms = "21 50 59.0 N", longitude_dms = "105 52 54.0 W"),
            fig, "EDFig1_isla_isabel_location.csv")
write_panel(data.frame(vertex = 1:5,
                       latitude_deg  = c(22.1203, 22.1203, 21.5797, 21.5797, 22.1203),
                       longitude_deg = c(-106.1763, -105.5937, -105.5937, -106.1763, -106.1763)),
            fig, "EDFig1_chlorophyll_extraction_polygon.csv")
write_panel(data.frame(
  map = c("main", "inset"),
  lon_min = c(-119, marker_lon_deg - 0.58), lon_max = c(-88, marker_lon_deg + 0.58),
  lat_min = c(13, marker_lat_deg - 0.48),  lat_max = c(33.5, marker_lat_deg + 0.48)),
  fig, "EDFig1_map_extents.csv")

# =============================================================================
# READMEs
# =============================================================================
ribbon_cols <- c(
  "  posterior_median      full-posterior median of the plotted quantity at that x (all 16,000 draws)",
  "  lower_95, upper_95    full-posterior 2.5% and 97.5% quantiles (ggdist::median_qi, .width = 0.95, as used by stat_lineribbon)",
  "  within_x_axis_limits  TRUE if the x value lies inside the figure's x-axis limits (rows with FALSE are not drawn)",
  "  drawn_median, drawn_lower_95, drawn_upper_95  line and ribbon values as actually rendered in the figure",
  "                        (ggplot_build of the stat_lineribbon layer with the figure's axis limits). ggplot2 removes",
  "                        posterior draws outside the y-axis limits before the stat is computed, so these can differ",
  "                        from the full-posterior columns where n_draws_within_axis_limits < total draws.",
  "  n_draws_within_axis_limits  number of posterior draw rows at that x that fall inside the y-axis limits"
)
censor_note_for <- function(fig) {
  ck <- bind_rows(check_log)
  ck <- ck[startsWith(ck$panel, paste0(fig, " ")) &
           pmax(ck$max_abs_diff_full_vs_drawn_median, ck$max_abs_diff_full_vs_drawn_lower,
                ck$max_abs_diff_full_vs_drawn_upper) > 1e-8, ]
  if (nrow(ck) == 0) return(character(0))
  censor_note(paste(sub("^[^ ]+ ", "", ck$panel), collapse = ", "))
}
censor_note <- function(panels) c(
  "",
  paste0("NOTE: in ", panels, " some posterior draws fall outside the plotted y-axis range; the rendered line/ribbon"),
  "(drawn_* columns) is therefore computed from the retained draws only and differs from the full-posterior",
  "summary (posterior_median, lower_95, upper_95). Both are provided. In all other panels the two are identical",
  "to machine precision."
)
point_note <- "  within_axis_limits    TRUE if the point lies inside the figure's axis limits (FALSE points are outside the plotted range and not drawn)"
general <- c(
  "Posterior summaries are computed from the precomputed posterior draws in Rdata/precomputed_draws",
  "(16,000 posterior draws; 100-point prediction grid; random effects excluded, re_formula = NA;",
  "other covariates held at their sample means). Generated by Rscripts/export_source_data.R.",
  ""
)

temporal_readme <- function(fig, sexlab, labs, bloom = FALSE) {
  c(paste0(fig, " source data (", sexlab, ", temporal trends)"), "", general,
    if (bloom) c(
      "Fig2A_bloom_onset_mean.csv : predicted mean bloom onset vs year (bloom model, model_delta5).",
      "  year = calendar year; values in days since 21 December of the preceding year.",
      ribbon_cols,
      "Fig2A_points.csv : observed annual bloom onset (5% threshold). year; bloom_onset_days (days since 21 Dec); within_axis_limits.",
      "Fig2B_bloom_timing_sigma.csv : predicted residual SD (sigma, unpredictability, days) of bloom onset vs year.",
      ribbon_cols,
      "Fig2B_points.csv : absolute residuals of observed bloom onset from the bloom model (posterior-mean residual). year; abs_residual_days.",
      ""),
    paste0(labs["mm"], "_mismatch_mean.csv : predicted mean phenological mismatch (days; laying date minus bloom onset, 5% threshold) vs year."),
    ribbon_cols,
    paste0(labs["mm"], "_points.csv : observed mismatch values. year; mismatch_days; n_records = number of breeding records with that",
           " year x mismatch value (mapped to point size). Points were horizontally jittered (width 0.2 yr) in the figure",
           if (sexlab == "males") " code but drawn with geom_point (no jitter applied)" else "", "; unjittered years are given."),
    point_note,
    paste0(labs["sd"], "_mismatch_sigma.csv : predicted residual SD of mismatch (sigma, inconsistency, days) vs year."),
    ribbon_cols,
    paste0(labs["sd"], "_points.csv : absolute residuals (days) of observed mismatch from the mismatch model (m11, posterior-mean residual),",
           " one row per unique year x residual value; n_records = number of records with that value (not mapped; fixed point size).",
           " Points were horizontally jittered (width 0.2 yr); unjittered years given."),
    point_note,
    paste0(labs["ns"], "_nest_success.csv : predicted probability of nest success (>=1 fledgling) vs year (ordinal model Q17)."),
    ribbon_cols,
    paste0(labs["ns"], "_points.csv : observed annual proportion of successful breeding records. year; n_successful;",
           " n_records_year (mapped to point size); proportion_successful."),
    point_note,
    paste0(labs["zeta"], "_zeta.csv : predicted discrimination parameter (zeta; 'consistency', unitless) of the ordinal model vs year.",
           " The plotted frame holds one row per draw x outcome category (4 identical zeta values per draw), so",
           " n_draws_within_axis_limits counts 4 x 16,000 rows; medians and quantiles are computed on that same frame as in the figure."),
    ribbon_cols
  )
}
write_readme("Fig2", c(temporal_readme("Fig. 2", "females", c(mm = "Fig2C", sd = "Fig2D", ns = "Fig2E", zeta = "Fig2F"), bloom = TRUE), censor_note_for("Fig2")))
write_readme("EDFig4", c(temporal_readme("Extended Data Fig. 4", "males", c(mm = "EDFig4A", sd = "EDFig4B", ns = "EDFig4C", zeta = "EDFig4D")), censor_note_for("EDFig4")))

age_readme <- function(fig, sexlab, labs) {
  c(paste0(fig, " source data (", sexlab, ", age and mismatch effects)"), "", general,
    paste0(labs["mm"], "_mismatch_mean_age.csv : predicted mean phenological mismatch (days) vs age (years). x-axis limits 1-23 yr."),
    ribbon_cols,
    paste0(labs["mm"], "_points.csv : observed mismatch values by age. age_years; mismatch_days; n_records (mapped to point size)."),
    point_note,
    paste0(labs["sd"], "_mismatch_sigma_age.csv : predicted residual SD of mismatch (sigma, days) vs age."),
    ribbon_cols,
    paste0(labs["sd"], "_points.csv : absolute residuals (days) from the mismatch model by age, one row per unique age x residual value;",
           " n_records not mapped. Points were horizontally jittered (width 0.2 yr); unjittered ages given."),
    point_note,
    paste0(labs["ns"], "_nest_success_age.csv : predicted probability of nest success (>=1 fledgling) vs age."),
    ribbon_cols,
    paste0(labs["ns"], "_points.csv : observed proportion of successful records by age. age_years; n_successful; n_records_age (point size); proportion_successful."),
    point_note,
    paste0(labs["cat"], "_fledging_category_probs.csv : predicted probability of each fledgling category (0, 1-2, 3 fledglings) vs mismatch (days)."),
    "  fledgling_category    0, 1-2 or 3 fledglings",
    ribbon_cols,
    paste0(labs["cat"], "_points.csv : observed proportion of records in each fledgling category at each observed mismatch value."),
    "  n_records_category; n_records_mismatch_value (mapped to point size); proportion.",
    point_note
  )
}
write_readme("Fig3", c(age_readme("Fig. 3", "females", c(mm = "Fig3A", sd = "Fig3B", ns = "Fig3C", cat = "Fig3D")), censor_note_for("Fig3")))
write_readme("EDFig5", c(age_readme("Extended Data Fig. 5", "males", c(mm = "EDFig5A", sd = "EDFig5B", ns = "EDFig5C", cat = "EDFig5D")), censor_note_for("EDFig5")))

write_readme("EDFig2", c(
  "Extended Data Fig. 2 source data (individual mismatch trajectories with age)", "", general,
  "EDFig2A_population_female.csv / EDFig2B_population_male.csv : population-level (fixed-effect) predicted mismatch (days) vs age (years).",
  ribbon_cols,
  "EDFig2A_individuals_female.csv / EDFig2B_individuals_male.csv : model-implied trajectories of 50 randomly sampled individuals",
  "  (set.seed(42)), including each individual's random intercept and age slope (re_formula = ~(1 + AGE | RING)).",
  "  ring = individual band ID; age_years; posterior_median_mismatch_days = posterior median of the linear predictor (days).",
  point_note,
  censor_note_for("EDFig2")
))
write_readme("EDFig3", c(
  "Extended Data Fig. 3 source data (posterior expected fledglings, first vs last study year)", "",
  "Computed from the ordinal nest-success model Q17 (female_fit_Q17.rds / male_fit_Q17.rds) at the first (2003)",
  "and last (2023) study year, other predictors at their centred means (0), random effects excluded.",
  "Expected fledglings = sum over categories (0,1,2,3) of category x posterior category probability, per draw.", "",
  "EDFig3_expected_fledglings_draws.csv : all posterior draws plotted as violins. sex; year; draw (posterior draw index); expected_fledglings.",
  "EDFig3_summary.csv : per sex x year: n_draws; posterior_mean (plotted point); posterior_median; lower_2.5, upper_97.5 (plotted error bar).",
  "EDFig3_pct_change_2003_2023.csv : per-draw percent change 100*(2023 - 2003)/2003, summarised as median and 2.5/97.5% quantiles."
))
write_readme("EDFig6", c(
  "Extended Data Fig. 6 source data (annual breeding records by fledgling number)", "",
  "EDFig6_breeding_records_by_fledglings.csv : sex; year; fledglings (0-3 fledglings per record); breeding_records (count, bar height).",
  "EDFig6_shaded_period.csv : x-range of the shaded 2015-2016 ENSO period (year units)."
))
write_readme("EDFig1", c(
  "Extended Data Fig. 1 source data (study-site map)", "",
  "The map was drawn in book/1-bloom-onset.qmd with Natural Earth basemaps (rnaturalearth, scale = 'medium'; public data, not reproduced here).",
  "EDFig1_isla_isabel_location.csv : Isla Isabel marker, decimal degrees (WGS84) and DMS.",
  "EDFig1_chlorophyll_extraction_polygon.csv : closed polygon vertices (WGS84) of the ~60 x 60 km chlorophyll-a extraction box.",
  "EDFig1_map_extents.csv : longitude/latitude limits of the main map and the inset."
))

# =============================================================================
# Logs and zips
# =============================================================================
checks <- bind_rows(check_log)
utils::write.csv(checks, file.path(out_root, "_export_checks.csv"), row.names = FALSE)
files <- bind_rows(file_log)
files$bytes <- file.size(file.path(out_root, files$figure, files$file))
utils::write.csv(files, file.path(out_root, "_file_manifest.csv"), row.names = FALSE)

old_wd <- setwd(out_root)
for (f in sort(unique(files$figure))) {
  zf <- paste0(f, "_source_data.zip")
  if (file.exists(zf)) file.remove(zf)
  zip::zipr(zf, files = f)
  message("  zipped ", zf, " (", file.size(zf), " bytes)")
}
setwd(old_wd)

print(checks)
message("Done.")
