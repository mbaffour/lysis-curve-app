# Lysis Curve OD Visualization App  v3.0
# Michael Baffour Awuah / Ramsey Lab
#
# v3.0 -- Complete phenotyping platform:
#   BUG FIXES          -- Fixed case_when crash in stats (dplyr compatibility).
#                         Fixed if() vs ifelse() in infection metrics (vectorization).
#   SIGNIFICANCE       -- Bar plots now show significance brackets when stats
#                         have been run for the selected metric.
#   ERROR BAR OPTIONS  -- Choose SD, SEM, 95% CI, or None for bar plot error bars.
#   CUSTOM BAR COLORS  -- Auto-match from Plot tab colors, OR override with
#                         custom HEX per bar via color pickers.
#   HEATMAP SELECTION  -- Checkboxes to pick which metrics appear in the
#                         phenotype heatmap.
#   OD TIME HEATMAP    -- New: samples × time heatmap colored by OD value
#                         (inferno palette). Shows full dynamics in one image.
#   EXPORT CONTROLS    -- Every Analysis plot now has its own format selector
#                         (PDF/PNG/SVG/TIFF/JPEG) and export dimensions.
#   SMOOTHING HELP     -- Expandable explanation of what smoothing does.
#   PLASMID NOTE       -- Infection strength section notes applicability to
#                         plasmid expression experiments (use empty-vector as ref).
#
# v2.9 -- Full Phenotyping Dashboard:
#   15 METRICS         -- Initial OD, max OD, time to max OD, final OD, AUC,
#                         μmax, doubling time, lag phase, stationary phase
#                         duration, lysis time, lysis rate, OD drop magnitude,
#                         residual OD, recovery slope. Plus 5 infection metrics
#                         when a reference is selected.
#   DERIVATIVE PLOT    -- dOD/dt vs time shows growth/lysis rates visually.
#                         Positive = growth, negative = lysis. Zero-crossing
#                         marks lysis onset.
#   INFECTION METRICS  -- Select a reference (control) sample to compute:
#                         infection strength, relative growth, relative μmax,
#                         relative max OD, lysis onset Δ from lag.
#   ANNOTATED CURVES   -- Growth curves with toggleable phase markers:
#                         lag end (green ▲), max OD (blue ◆), lysis onset
#                         (red ▼).
#   HEATMAP            -- Z-score normalised phenotype heatmap (samples ×
#                         metrics) with actual values printed in each cell.
#   BAR PLOTS          -- Now supports all 19 metrics including infection
#                         metrics. Colors inherited from Plot tab.
#   STATS              -- Same t-test/ANOVA/pairwise engine, now covering all
#                         metrics. "Download All Stats" applies BH correction.
#
# v2.8 -- Analysis Dashboard:
#   METRICS ENGINE     -- New Analysis tab with automatic growth curve metric
#                         extraction: max OD, μmax, doubling time, AUC, lag
#                         phase, lysis onset, lysis rate per sample.
#   BAR PLOTS          -- Interactive bar plot comparison of any metric across
#                         samples, with jittered data points and SEM error bars.
#                         Inherits sample colors from the Plot tab.
#   STATISTICS         -- Pairwise statistical comparisons (t-test for 2 groups,
#                         ANOVA + pairwise t with BH correction for ≥3 groups).
#   EXPORT             -- Download metrics table, stats results, and bar plots.
#                         "Download All Stats" runs all metrics at once with
#                         BH multiple-testing correction.
#   INTERACTIVE TABLE  -- DT-powered sortable/searchable metrics table.
#
# v2.7 -- Bug fixes:
#   SHADOW FIX         -- Shadow/ribbon error display no longer crashes when
#                         points are also enabled. The fill scale conflict
#                         between geom_ribbon (scale_fill_manual) and
#                         geom_point (scale_fill_identity) is resolved by
#                         using a unified scale_fill_manual when shadow mode
#                         is active.
#   NA ERROR FIX       -- Single-replicate samples (where sd() returns NA)
#                         no longer crash error bars or ribbons. NA sd values
#                         are replaced with 0 before computing error bounds.
#   OPERATOR ORDER FIX -- The %||% (null-coalescing) operator is now defined
#                         before its first use in renderUI.
#   ERROR GUARD        -- Error column values are explicitly guarded against
#                         NA before computing err_lo / err_hi bounds.
#
# v2.6 -- Stability & performance improvements for large replicate counts:
#   OBSERVER LEAK FIX  -- Color observeEvent handlers are registered at most
#                         once per sample ID (tracked in registered_color_obs).
#   POINT LOOP -> O(1) -- build_plot() no longer adds one geom_point() layer
#                         per sample; shape/fill are vectorised into the data
#                         frame and a single geom_point() call covers all.
#   PIVOT GUARD        -- prepare_plot_data() pre-selects only the requested
#                         sample columns before pivot_longer().
#   VECTORISED FILTER  -- apply_time_filters() uses outer() matrix comparison
#                         instead of per-row vapply for time-point exclusions.
#   DEBOUNCE           -- Main plot render waits 400 ms after the last input
#                         change before rebuilding, preventing render pile-up.
#   CRASH GUARD        -- build_plot() wrapped in tryCatch(); errors display
#                         as an in-plot message rather than crashing the session.
#   PAGINATION         -- Variable Styling panel shows at most 20 samples at a
#                         time with Prev/Next navigation; keeps the DOM small.
# Michael Baffour Awuah / Ramsey Lab
#
# v2.5 changes:
#   SETTINGS  – Save/load now stores BOTH per-sample aesthetics (color, shape,
#               linetype keyed by sample name) AND all global visual settings.
#               On import: matched samples get their saved aesthetics applied
#               automatically; a status panel shows exactly which names matched
#               and which didn't.  "Clear Imported" fully revokes everything and
#               restores defaults — per-sample overrides included.
#
#   GIF       – Rebuilt without gganimate. Each frame is a real ggplot rendered
#               to PNG via ragg/png, then stitched with gifski. Full fidelity:
#               same theme, axes, error bars, shapes, colors as the main plot.
#               Each frame adds one more line (cumulative reveal).
#
#   PPTX      – Each slide is a pixel-perfect render of build_plot() with the
#               cumulative subset + grey-out of non-active samples. Final slide
#               shows all samples in full color. Slide dimensions driven by the
#               export width/height inputs.

# ── Packages ──────────────────────────────────────────────────────────────────
library(shiny)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggrepel)
library(ggprism)
library(svglite)
library(jsonlite)
library(zoo)    # rolling calculations for growth metrics
library(DT)     # interactive tables in Analysis tab

has_officer   <- requireNamespace("officer",   quietly = TRUE)
has_rvg       <- requireNamespace("rvg",       quietly = TRUE)
has_gifski    <- requireNamespace("gifski",    quietly = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────

safe_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

normalize_hex_color <- function(hex) {
  if (is.null(hex) || is.na(hex) || nchar(trimws(hex)) == 0) return("#000000")
  hex <- trimws(hex)
  if (!grepl("^#", hex)) hex <- paste0("#", hex)
  if (grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", hex)) return(toupper(hex))
  "#000000"
}

rgb_to_hex <- function(r, g, b, a = 1) {
  rgb(min(max(as.integer(r), 0), 255) / 255,
      min(max(as.integer(g), 0), 255) / 255,
      min(max(as.integer(b), 0), 255) / 255,
      min(max(a, 0), 1))
}

parse_excluded_timepoints <- function(text, all_timepoints) {
  if (is.null(text) || nchar(trimws(text)) == 0) return(numeric(0))
  parts <- strsplit(trimws(text), "[,; ]+")[[1]]
  vals  <- suppressWarnings(as.numeric(parts))
  vals  <- vals[!is.na(vals)]
  if (length(vals) == 0 || is.null(all_timepoints)) return(numeric(0))
  matched <- sapply(vals, function(v) {
    diffs <- abs(all_timepoints - v)
    tol   <- max(abs(all_timepoints)) * 0.001 + 0.001
    if (min(diffs) <= tol) all_timepoints[which.min(diffs)] else NA_real_
  })
  unique(matched[!is.na(matched)])
}

detect_replicate_column <- function(data, time_col, group_col, value_col) {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  ignore_cols <- c(time_col, group_col, value_col)
  names_lower <- tolower(colnames(data))
  rep_hits <- names_lower %in% c("replicate", "rep", "rep_id", "replicate_id",
                                 "well", "well_id", "wellid",
                                 "technical_replicate", "biological_replicate", "repeat")
  rep_cols <- colnames(data)[rep_hits]
  rep_cols <- setdiff(rep_cols, ignore_cols)
  if (length(rep_cols) > 0) return(rep_cols[1])
  NULL
}

# Detect replicates in wide format by finding repeated time blocks
infer_wide_replicates <- function(data, time_col) {
  if (is.null(data) || nrow(data) == 0) return(list(reps = NULL, block_id = NULL))
  t_vec <- suppressWarnings(as.numeric(as.character(data[[time_col]])))
  if (all(is.na(t_vec))) return(list(reps = NULL, block_id = NULL))
  resets <- which(diff(t_vec) < 0)
  if (length(resets) == 0) return(list(reps = 1L, block_id = rep(1L, nrow(data))))
  block_starts <- c(1L, resets + 1L)
  block_ends   <- c(resets, nrow(data))
  block_id <- integer(nrow(data))
  for (i in seq_along(block_starts)) {
    block_id[block_starts[i]:block_ends[i]] <- i
  }
  list(reps = length(block_starts), block_id = block_id)
}

# ── Growth Metrics Engine ─────────────────────────────────────────────────────
# Calculates per-sample growth curve metrics from the summarised plot data.
# Works with the app's existing prepare_plot_data() output (time, variable,
# mean_value columns).  Returns one row per sample.

calculate_growth_metrics <- function(plot_data, smooth_window = 5) {
  if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)
  
  has_rep <- "replicate" %in% names(plot_data)
  grp_cols <- if (has_rep) c("variable", "replicate") else "variable"
  
  results <- plot_data %>%
    group_by(across(all_of(grp_cols))) %>%
    arrange(time) %>%
    do({
      d <- .
      t_vec  <- d$time
      od_vec <- d$mean_value
      n_pts  <- length(t_vec)
      
      # ── Basic metrics ──────────────────────────────────────────────────────
      initial_od  <- od_vec[1]
      max_od      <- max(od_vec, na.rm = TRUE)
      time_max_od <- t_vec[which.max(od_vec)]
      final_od    <- tail(od_vec, 1)
      
      # ── AUC (trapezoidal rule) ──────────────────────────────────────────────
      auc <- if (n_pts >= 2) {
        sum(diff(t_vec) * (head(od_vec, -1) + tail(od_vec, -1)) / 2)
      } else NA_real_
      
      # ── Growth rate (μmax) via rolling slope on log(OD) ─────────────────────
      mu_max     <- NA_real_
      lag_phase  <- NA_real_
      win        <- min(smooth_window, n_pts)
      
      if (n_pts >= win && win >= 3) {
        log_od <- log(pmax(od_vec, 1e-9))

        roll_slopes <- tryCatch({
          rollapply(seq_len(n_pts), width = win, FUN = function(idx) {
            tt <- t_vec[idx]
            ll <- log_od[idx]
            if (diff(range(tt)) == 0) return(0)
            coef(lm(ll ~ tt))[2]
          }, align = "center", fill = NA)
        }, error = function(e) rep(NA_real_, n_pts))

        
        valid_slopes <- which(is.finite(roll_slopes) & roll_slopes > 0)
        if (length(valid_slopes) > 0) {
          mu_max <- max(roll_slopes[valid_slopes], na.rm = TRUE)
          
          threshold  <- 0.10 * mu_max
          lag_idx    <- which(roll_slopes >= threshold & is.finite(roll_slopes))
          lag_phase  <- if (length(lag_idx) > 0) t_vec[lag_idx[1]] else NA_real_
        }
      }
      
      # ── Doubling time ───────────────────────────────────────────────────────
      doubling_time <- if (is.finite(mu_max) && mu_max > 0) log(2) / mu_max else NA_real_
      
      # ── Stationary phase duration ───────────────────────────────────────────
      # Time spent within 95% of max_od before any decline
      stat_phase_dur <- NA_real_
      max_idx <- which.max(od_vec)
      if (max_idx > 1) {
        near_max <- which(od_vec >= 0.95 * max_od)
        if (length(near_max) >= 2) {
          stat_phase_dur <- t_vec[max(near_max)] - t_vec[min(near_max)]
        }
      }
      
      # ── Lysis detection ─────────────────────────────────────────────────────
      lysis_time    <- NA_real_
      lysis_rate    <- NA_real_
      od_drop       <- NA_real_   # absolute drop from peak
      residual_od   <- NA_real_   # minimum OD after lysis
      recovery_slope <- NA_real_  # regrowth after lysis trough
      
      if (max_idx < n_pts) {
        post_peak   <- (max_idx + 1):n_pts
        od_decline  <- max_od - od_vec[post_peak]
        lysis_idx   <- which(od_decline > 0.05 * max_od)
        if (length(lysis_idx) > 0) {
          abs_idx    <- post_peak[lysis_idx[1]]
          lysis_time <- t_vec[abs_idx]
          
          # Lysis rate: slope from peak to deepest point after peak
          trough_idx <- max_idx + which.min(od_vec[post_peak])
          if (trough_idx > max_idx) {
            lysis_window <- max_idx:trough_idx
            if (length(lysis_window) >= 2) {
              lysis_fit <- tryCatch(
                coef(lm(od_vec[lysis_window] ~ t_vec[lysis_window]))[2],
                error = function(e) NA_real_
              )
              lysis_rate <- unname(lysis_fit)
            }
          }
          
          # OD drop magnitude
          post_min   <- min(od_vec[post_peak], na.rm = TRUE)
          od_drop    <- max_od - post_min
          residual_od <- post_min
          
          # Recovery slope: positive slope after the trough
          if (trough_idx < n_pts) {
            recovery_window <- trough_idx:n_pts
            if (length(recovery_window) >= 3) {
              rec_fit <- tryCatch(
                coef(lm(od_vec[recovery_window] ~ t_vec[recovery_window]))[2],
                error = function(e) NA_real_
              )
              rec_val <- unname(rec_fit)
              recovery_slope <- if (is.finite(rec_val) && rec_val > 0) rec_val else NA_real_
            }
          }
        }
      }
      
      out <- tibble(
        initial_od     = initial_od,
        max_od         = max_od,
        time_max_od    = time_max_od,
        final_od       = final_od,
        auc            = auc,
        mu_max         = mu_max,
        doubling_time  = doubling_time,
        lag_phase      = lag_phase,
        stat_phase_dur = stat_phase_dur,
        lysis_time     = lysis_time,
        lysis_rate     = lysis_rate,
        od_drop        = od_drop,
        residual_od    = residual_od,
        recovery_slope = recovery_slope
      )
      
      if (has_rep) {
        out$replicate <- d$replicate[1]
      }
      out
    }) %>%
    ungroup() %>%
    rename(sample = variable)
  
  results
}

# ── Derivative data (dOD/dt) for plotting ─────────────────────────────────────
# Returns a tibble with time, variable, dODdt columns for the derivative plot.
calculate_derivative <- function(plot_data, smooth_window = 5) {
  if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)
  
  plot_data %>%
    group_by(variable) %>%
    arrange(time) %>%
    mutate(
      dODdt = {
        od  <- mean_value
        tt  <- time
        n   <- length(tt)
        win <- min(smooth_window, n)
        if (n >= win && win >= 3) {
          tryCatch({
            rollapply(seq_len(n), width = win, FUN = function(idx) {
              if (diff(range(tt[idx])) == 0) return(0)
              coef(lm(od[idx] ~ tt[idx]))[2]
            }, align = "center", fill = NA)
          }, error = function(e) rep(NA_real_, n))
        } else rep(NA_real_, n)
      }
    ) %>%
    ungroup() %>%
    filter(is.finite(dODdt))
}

# ── Infection strength (requires reference condition) ─────────────────────────
calculate_infection_metrics <- function(metrics_df, ref_sample) {
  if (is.null(metrics_df) || nrow(metrics_df) == 0 || is.null(ref_sample)) return(NULL)
  if (!ref_sample %in% metrics_df$sample) return(NULL)
  
  ref_row <- metrics_df %>% filter(sample == ref_sample)
  ref_auc    <- mean(ref_row$auc,    na.rm = TRUE)
  ref_mu_max <- mean(ref_row$mu_max, na.rm = TRUE)
  ref_max_od <- mean(ref_row$max_od, na.rm = TRUE)
  
  metrics_df %>%
    mutate(
      relative_growth    = ifelse(is.finite(ref_auc) & ref_auc > 0, auc / ref_auc, NA_real_),
      infection_strength = ifelse(is.finite(ref_auc) & ref_auc > 0, 1 - (auc / ref_auc), NA_real_),
      relative_mu_max    = ifelse(is.finite(ref_mu_max) & ref_mu_max > 0, mu_max / ref_mu_max, NA_real_),
      relative_max_od    = ifelse(is.finite(ref_max_od) & ref_max_od > 0, max_od / ref_max_od, NA_real_),
      lysis_onset_delta  = lysis_time - lag_phase
    )
}

# ── Replicate Summary ─────────────────────────────────────────────────────────

summarise_metrics <- function(metrics_df) {
  if (is.null(metrics_df) || nrow(metrics_df) == 0) return(NULL)
  
  metric_cols <- c("initial_od", "max_od", "time_max_od", "final_od", "auc",
                   "mu_max", "doubling_time", "lag_phase", "stat_phase_dur",
                   "lysis_time", "lysis_rate", "od_drop", "residual_od",
                   "recovery_slope",
                   "relative_growth", "infection_strength", "relative_mu_max",
                   "relative_max_od", "lysis_onset_delta")
  present_cols <- intersect(metric_cols, names(metrics_df))
  
  long_df <- metrics_df %>%
    pivot_longer(cols = all_of(present_cols),
                 names_to = "metric", values_to = "value") %>%
    filter(is.finite(value))
  
  has_rep <- "replicate" %in% names(metrics_df)
  has_dup <- any(duplicated(long_df[, c("sample", "metric")]))
  
  if (has_rep || has_dup) return(long_df)
  
  long_df %>%
    group_by(sample, metric) %>%
    summarise(value = mean(value), .groups = "drop")
}

# ── Statistical Comparisons ───────────────────────────────────────────────────
# Compares a chosen metric across all samples using t-test (2 groups) or
# ANOVA (>2 groups) + pairwise post-hoc.

compare_metrics <- function(metrics_long, metric_name) {
  if (is.null(metrics_long) || nrow(metrics_long) == 0) return(NULL)
  
  d <- metrics_long %>% filter(metric == metric_name, is.finite(value))
  groups <- unique(d$sample)
  if (length(groups) < 2) return(NULL)
  
  counts <- d %>% group_by(sample) %>% summarise(n = n(), .groups = "drop")
  if (any(counts$n < 2)) return(NULL)
  
  results <- list()
  
  if (length(groups) == 2) {
    # Pairwise t-test
    g1 <- d$value[d$sample == groups[1]]
    g2 <- d$value[d$sample == groups[2]]
    if (length(g1) >= 2 && length(g2) >= 2) {
      tt <- tryCatch(t.test(g1, g2), error = function(e) NULL)
      wt <- tryCatch(wilcox.test(g1, g2, exact = FALSE), error = function(e) NULL)
      if (!is.null(tt))
        results[[length(results) + 1]] <- tibble(
          metric = metric_name,
          comparison = paste(groups[1], "vs", groups[2]),
          test = "t.test", p_value = tt$p.value)
      if (!is.null(wt))
        results[[length(results) + 1]] <- tibble(
          metric = metric_name,
          comparison = paste(groups[1], "vs", groups[2]),
          test = "wilcox.test", p_value = wt$p.value)
    }
  } else {
    # ANOVA + pairwise t-tests
    aov_res <- tryCatch({
      fit <- aov(value ~ sample, data = d)
      s   <- summary(fit)
      tibble(metric = metric_name, comparison = "Overall (ANOVA)",
             test = "anova", p_value = s[[1]]$`Pr(>F)`[1])
    }, error = function(e) NULL)
    if (!is.null(aov_res)) results[[length(results) + 1]] <- aov_res
    
    # Pairwise comparisons
    pw <- tryCatch({
      pt <- pairwise.t.test(d$value, d$sample, p.adjust.method = "BH")
      pmat <- pt$p.value
      pairs <- which(!is.na(pmat), arr.ind = TRUE)
      if (nrow(pairs) > 0) {
        tibble(
          metric     = metric_name,
          comparison = paste(rownames(pmat)[pairs[,1]], "vs", colnames(pmat)[pairs[,2]]),
          test       = "pairwise.t (BH)",
          p_value    = pmat[pairs]
        )
      } else NULL
    }, error = function(e) NULL)
    if (!is.null(pw)) results[[length(results) + 1]] <- pw
  }
  
  if (length(results) == 0) return(NULL)
  out <- bind_rows(results)
  pv  <- out$p_value
  out$significance <- ifelse(is.na(pv), "ns",
                             ifelse(pv < 0.001, "***",
                                    ifelse(pv < 0.01,  "**",
                                           ifelse(pv < 0.05,  "*",
                                                  ifelse(pv < 0.1,   ".", "ns")))))
  out
}

# Friendly labels for metrics in UI dropdowns and plot titles
metric_labels <- c(
  initial_od         = "Initial OD",
  max_od             = "Max OD",
  time_max_od        = "Time to Max OD",
  final_od           = "Final OD",
  auc                = "Area Under Curve (AUC)",
  mu_max             = "Max Growth Rate (μmax, h⁻¹)",
  doubling_time      = "Doubling Time",
  lag_phase          = "Lag Phase",
  stat_phase_dur     = "Stationary Phase Duration",
  lysis_time         = "Lysis Onset Time",
  lysis_rate         = "Lysis Rate (h⁻¹)",
  od_drop            = "OD Drop Magnitude",
  residual_od        = "Residual OD (post-lysis)",
  recovery_slope     = "Recovery Slope (post-lysis)",
  relative_growth    = "Relative Growth (vs reference)",
  infection_strength = "Infection Strength",
  relative_mu_max    = "Relative μmax (vs reference)",
  relative_max_od    = "Relative Max OD (vs reference)",
  lysis_onset_delta  = "Lysis Onset Δ from Lag"
)

core_metric_choices <- c(
  "Initial OD"             = "initial_od",
  "Max OD"                 = "max_od",
  "Time to Max OD"         = "time_max_od",
  "Final OD"               = "final_od",
  "AUC"                    = "auc",
  "Max Growth Rate"        = "mu_max",
  "Doubling Time"          = "doubling_time",
  "Lag Phase"              = "lag_phase",
  "Stationary Phase Dur."  = "stat_phase_dur"
)

lysis_metric_choices <- c(
  "Lysis Onset Time"       = "lysis_time",
  "Lysis Rate"             = "lysis_rate",
  "OD Drop Magnitude"      = "od_drop",
  "Residual OD"            = "residual_od",
  "Recovery Slope"         = "recovery_slope"
)

infection_metric_choices <- c(
  "Infection Strength"     = "infection_strength",
  "Relative Growth"        = "relative_growth",
  "Relative μmax"          = "relative_mu_max",
  "Relative Max OD"        = "relative_max_od",
  "Lysis Onset Δ Lag"       = "lysis_onset_delta"
)

all_metric_choices <- c(
  "── Growth ──" = "",
  core_metric_choices,
  "── Lysis ──"  = "",
  lysis_metric_choices,
  "── Infection (needs reference) ──" = "",
  infection_metric_choices
)

# ── Optional UI widgets — must be pre-computed before fluidPage() ─────────────

pptx_ui <- if (has_officer && has_rvg) {
  downloadButton("downloadPPTX", "Download PowerPoint (.pptx)",
                 style = "background:#C0392B;color:white;border:none;width:100%;")
} else {
  p(style = "color:#999;font-size:.85em;",
    "Install 'officer' + 'rvg' to enable PPTX export.")
}

gif_ui <- if (has_gifski) {
  downloadButton("downloadGIF", "Download Animated GIF",
                 style = "background:#27AE60;color:white;border:none;width:100%;")
} else {
  p(style = "color:#999;font-size:.85em;",
    "Install 'gifski' to enable GIF export.")
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel(div(
    h3("OD Growth Curve Analyzer", style="margin:0;color:#1b2838;font-weight:700;"),
    p("Bacteria vs. Phage \u2014 who wins? Let the data decide.",
      style="margin:0;font-size:.8em;color:#555;font-style:italic;")
  )),
  
  tags$head(
    tags$script(HTML(
      "Shiny.addCustomMessageHandler('updateColorPreview', function(msg) {
         var el = document.getElementById(msg.id);
         if (el) el.style.backgroundColor = msg.color;
       });"
    )),
    uiOutput("night_mode_css"),
    tags$style(HTML("
      /* ── Base layout ─────────────────────────────────────────── */
      body { font-family: 'Segoe UI', Arial, sans-serif; }
      .panel-section {
        background:#f8f9fa; padding:12px; border-radius:5px;
        margin:10px 0; border:1px solid #e9ecef;
      }
      .panel-title { margin-top:0; margin-bottom:10px; font-weight:bold; color:#495057; }

      /* ── Collapsible details elements ─────────────────────────── */
      details > summary {
        cursor:pointer; font-weight:bold; padding:6px 10px;
        background:#e9ecef; border-radius:3px; list-style:none; user-select:none;
      }
      details > summary::-webkit-details-marker { display:none; }
      details[open] > summary::before { content:'\\25BC  '; font-size:.75em; }
      details > summary::before       { content:'\\25B6  '; font-size:.75em; }
      details { margin-bottom:6px; }

      /* ── Dark sidebar overrides ───────────────────────────────── */
      .well {
        background:#1b2838 !important;
        border-color:#2c3e50 !important;
        color:#e8f4f8 !important;
      }
      .well label, .well .control-label,
      .well p, .well small, .well span { color:#c8d8e8 !important; }
      .well h4.panel-title { color:#4ecdc4 !important; }
      .well .panel-section {
        background:#243447 !important;
        border-color:#3a4f6e !important;
        color:#e8f4f8 !important;
      }
      .well details > summary {
        background:#2c3e50 !important;
        color:#e8f4f8 !important;
      }
      .well .form-control {
        background:#1e3044 !important;
        color:#e8f4f8 !important;
        border-color:#3a4f6e !important;
      }
      .well .shiny-input-container { color:#c8d8e8 !important; }
      .well .checkbox label, .well .radio label { color:#c8d8e8 !important; }

      /* ── Color preview swatch ─────────────────────────────────── */
      .color-preview {
        display:inline-block; width:28px; height:28px;
        border:1px solid #ccc; border-radius:4px; flex-shrink:0;
      }

      /* ── Settings status banners ──────────────────────────────── */
      .settings-status { padding:8px 12px; border-radius:4px; font-size:.85em; margin-top:6px; }
      .settings-status.ok   { background:#d4edda; color:#155724; border:1px solid #c3e6cb; }
      .settings-status.warn { background:#fff3cd; color:#856404; border:1px solid #ffeeba; }
      .match-list { margin:4px 0 0 0; padding-left:16px; }
      .excl-preview {
        background:#fff3cd; color:#856404; padding:6px 10px;
        border-radius:4px; font-size:.82em; margin-top:4px;
      }

      /* ── Analysis sub-tabs ────────────────────────────────────── */
      #analysis_subtabs > li > a {
        font-weight:600; color:#2c3e50;
      }
      #analysis_subtabs > li.active > a {
        color:#4ecdc4 !important; border-bottom:2px solid #4ecdc4;
      }

      /* ── Witty label style ────────────────────────────────────── */
      .wit-label {
        font-size:.75em; color:#888; font-style:italic; margin-top:-4px; margin-bottom:6px;
      }
      .well .wit-label { color:#7a9bb5 !important; }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 
                 fileInput("file", "Choose CSV File", accept = c("text/csv", ".csv")),
                 div(style = "text-align:right; margin-top:-8px; margin-bottom:4px;",
                     checkboxInput("night_mode", "\U1F319 Night Mode", value = FALSE)),

                 # ── Settings ─────────────────────────────────────────────────────────────
                 div(class = "panel-section",
                     h4("Save / Load Settings", class = "panel-title"),
                     p(style = "font-size:.82em;color:#555;margin-bottom:8px;",
                       "Saves all visual settings: axis scales, labels, error bars, fonts, line/point options, ",
                       "AND per-sample aesthetics (color, shape, linetype) keyed by sample name. ",
                       "Loading restores all those settings — but never touches your data, ",
                       "sample selection, or time filter."),
                     fluidRow(
                       column(6, downloadButton("saveSettings",  "Save Settings",  style = "width:100%")),
                       column(6, actionButton("clearSettings", "Clear & Reset",
                                              icon  = icon("times"),
                                              style = "width:100%;background:#dc3545;color:white;border:none;"))
                     ),
                     br(),
                     fileInput("loadSettings", "Load Settings File",
                               accept = c("application/json", ".json")),
                     uiOutput("settings_status_ui")
                 ),
                 
                 # ── Axis Settings ─────────────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Axis Settings"),
                   div(class = "panel-section",
                       h4("Axis Scales", class = "panel-title"),
                       selectInput("x_scale_type", "X-axis Scale:",
                                   choices  = c("Linear" = "linear", "Logarithmic" = "log",
                                                "Square Root" = "sqrt", "Reverse" = "reverse"),
                                   selected = "linear"),
                       selectInput("y_scale_type", "Y-axis Scale:",
                                   choices  = c("Linear" = "linear", "Logarithmic" = "log",
                                                "Square Root" = "sqrt", "Reverse" = "reverse"),
                                   selected = "log"),
                       checkboxInput("use_advanced_ticks", "Advanced Tick Customization", value = TRUE),
                       conditionalPanel(
                         condition = "input.use_advanced_ticks == true && input.y_scale_type == 'log'",
                         h5("Y-axis Log Tick Range"),
                         fluidRow(
                           column(6, numericInput("y_log_min_exponent", "Min Exponent:", value = -2)),
                           column(6, numericInput("y_log_max_exponent", "Max Exponent:", value =  0))
                         )
                       ),
                       conditionalPanel(
                         condition = "input.use_advanced_ticks == true && input.x_scale_type == 'linear'",
                         h5("X-axis Tick Control", style = "margin:10px 0 4px;"),
                         fluidRow(
                           column(6,
                                  numericInput("x_tick_interval", "Major Interval:",
                                               value = NA, min = 0.001, step = 1)
                           ),
                           column(6,
                                  p(style = "font-size:.8em;color:#777;margin-top:28px;",
                                    "Leave blank = auto")
                           )
                         ),
                         textInput("x_extra_ticks", "Extra tick positions:",
                                   value = "", placeholder = "e.g. 270, 450"),
                         p(style = "font-size:.8em;color:#777;margin-top:-6px;",
                           "Values always shown regardless of interval. Comma-separated.")
                       ),
                       checkboxInput("custom_x_limits", "Custom X-axis limits", value = FALSE),
                       conditionalPanel(condition = "input.custom_x_limits == true",
                                        fluidRow(
                                          column(6, numericInput("x_min", "X-min:", value = 0)),
                                          column(6, numericInput("x_max", "X-max:", value = 100))
                                        )
                       ),
                       checkboxInput("custom_y_limits", "Custom Y-axis limits", value = FALSE),
                       conditionalPanel(condition = "input.custom_y_limits == true",
                                        fluidRow(
                                          column(6, numericInput("y_min", "Y-min:", value = 0.01)),
                                          column(6, numericInput("y_max", "Y-max:", value = 2))
                                        )
                       ),
                       sliderInput("x_expand_left",   "X Left Expand:",   0, 0.1, 0,    0.01),
                       sliderInput("x_expand_right",  "X Right Expand:",  0, 0.1, 0.05, 0.01),
                       sliderInput("y_expand_bottom", "Y Bottom Expand:", 0, 0.1, 0,    0.01),
                       sliderInput("y_expand_top",    "Y Top Expand:",    0, 0.1, 0.05, 0.01)
                   ),
                   div(class = "panel-section",
                       h4("Labels & Formatting", class = "panel-title"),
                       textInput("x_axis_label",  "X-axis Label:",  "Time (minutes)"),
                       textInput("y_axis_label",  "Y-axis Label:",  "A550"),
                       textInput("plot_title",    "Plot Title:",    ""),
                       textInput("plot_subtitle", "Plot Subtitle:", ""),
                       h5("Gridlines", style = "margin-top:12px;"),
                       checkboxInput("show_major_gridlines", "Show Major Gridlines", FALSE),
                       checkboxInput("show_minor_gridlines", "Show Minor Gridlines", FALSE),
                       conditionalPanel(condition = "input.show_major_gridlines == true",
                                        selectInput("major_gridline_color", "Major Gridline Color:",
                                                    choices  = c("Light Gray" = "#E6E6E6", "Medium Gray" = "#CCCCCC",
                                                                 "Dark Gray"  = "#999999", "Light Blue"  = "#DEEBF7",
                                                                 "Custom"     = "custom"),
                                                    selected = "#E6E6E6"),
                                        conditionalPanel(condition = "input.major_gridline_color == 'custom'",
                                                         textInput("major_gridline_color_custom", "HEX:", "#E6E6E6")),
                                        sliderInput("major_gridline_size", "Width:", 0.1, 1, 0.5, 0.1)
                       ),
                       conditionalPanel(condition = "input.show_minor_gridlines == true",
                                        selectInput("minor_gridline_color", "Minor Gridline Color:",
                                                    choices  = c("Light Gray" = "#F2F2F2", "Medium Gray" = "#E6E6E6",
                                                                 "Dark Gray"  = "#CCCCCC", "Light Blue"  = "#EFF6FB",
                                                                 "Custom"     = "custom"),
                                                    selected = "#F2F2F2"),
                                        conditionalPanel(condition = "input.minor_gridline_color == 'custom'",
                                                         textInput("minor_gridline_color_custom", "HEX:", "#F2F2F2")),
                                        sliderInput("minor_gridline_size", "Width:", 0.1, 1, 0.3, 0.1)
                       ),
                       h5("Font Settings", style = "margin-top:12px;"),
                       selectInput("font_family", "Font Family:",
                                   choices  = c("Sans Serif (Arial/Helvetica)" = "sans",
                                                "Serif (Times New Roman)"      = "serif",
                                                "Monospace (Courier New)"      = "mono",
                                                "Helvetica Neue"               = "Helvetica Neue",
                                                "Garamond"                     = "Garamond",
                                                "Palatino"                     = "Palatino"),
                                   selected = "sans"),
                       numericInput("title_font_size",      "Title Font Size:",  20, 8, 48),
                       numericInput("axis_label_font_size", "Axis Label Size:",  20, 8, 36),
                       numericInput("axis_text_font_size",  "Axis Text Size:",   16, 6, 32),
                       checkboxInput("bold_title",         "Bold Title",         TRUE),
                       checkboxInput("italic_axis_labels", "Italic Axis Labels", FALSE),
                       selectInput("axis_text_angle", "X-axis Text Angle:",
                                   choices  = c("Horizontal (0°)" = "0",
                                                "Angled (45°)"    = "45",
                                                "Vertical (90°)"  = "90"),
                                   selected = "0")
                   )
                 ),
                 
                 # ── Time Point Filtering ──────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Time Point Filtering"),
                   div(class = "panel-section",
                       h4("Restrict Displayed Time Range", class = "panel-title"),
                       p(style = "font-size:.85em;color:#555;",
                         "Filters data before plotting. All statistics, error bars, and exports reflect the filtered data."),
                       checkboxInput("enable_time_filter", "Enable Time Filtering", value = FALSE),
                       conditionalPanel(condition = "input.enable_time_filter == true",
                                        h5("Time Range", style = "margin:8px 0 4px;"),
                                        uiOutput("time_range_slider_ui"),
                                        hr(style = "margin:10px 0;"),
                                        h5("Exclude Specific Time Points", style = "margin:8px 0 4px;"),
                                        p(style = "font-size:.82em;color:#666;", "Comma-separated values to remove (e.g. 0, 5, 180)."),
                                        textInput("exclude_timepoints", label = NULL, value = "", placeholder = "e.g. 0, 5, 180"),
                                        uiOutput("excluded_timepoints_preview")
                       )
                   )
                 ),
                 
                 # ── Region Highlighting ───────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Region Highlighting"),
                   div(class = "panel-section",
                       checkboxInput("enable_highlighting", "Enable Region Highlighting", FALSE),
                       conditionalPanel(condition = "input.enable_highlighting == true",
                                        numericInput("region_count", "Number of Regions:", 1, 1, 5),
                                        uiOutput("region_settings")
                       )
                   )
                 ),
                 
                 # ── Time Point Markers ────────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Time Point Markers"),
                   div(class = "panel-section",
                       checkboxInput("enable_time_markers", "Enable Time Markers", FALSE),
                       conditionalPanel(condition = "input.enable_time_markers == true",
                                        numericInput("marker_count", "Number of Markers:", 1, 1, 10),
                                        uiOutput("time_marker_settings")
                       )
                   )
                 ),
                 
                 # ── Color Palettes ────────────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Color Palettes"),
                   div(class = "panel-section",
                       selectInput("color_palette", "Color Palette:",
                                   choices  = c("Custom"              = "custom",
                                                "Viridis"             = "viridis",
                                                "Plasma"              = "plasma",
                                                "Colorblind-friendly" = "colorblind",
                                                "Publication"         = "publication",
                                                "Rainbow"             = "rainbow",
                                                "Grayscale"           = "gray"),
                                   selected = "custom"),
                       htmlOutput("palette_preview")
                   )
                 ),
                 
                 # ── Line & Point Settings ─────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Line & Point Settings"),
                   div(class = "panel-section",
                       h4("Line Options", class = "panel-title"),
                       sliderInput("line_thickness", "Line Thickness:", 0.1, 3, 1, 0.1)
                   ),
                   div(class = "panel-section",
                       h4("Point Options", class = "panel-title"),
                       checkboxInput("show_points", "Show Points", TRUE),
                       conditionalPanel(condition = "input.show_points == true",
                                        sliderInput("shape_size",   "Point Size:",   0.5, 8, 3,   0.1),
                                        sliderInput("point_stroke", "Stroke Width:", 0,   2, 0.5, 0.1)
                       ),
                   )
                 ),

                 # ── Legend ────────────────────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Legend"),
                   div(class = "panel-section",
                       fluidRow(
                         column(6, selectInput("legend_position", "Position:",
                                               choices = c("Right" = "right", "Left" = "left",
                                                           "Top" = "top", "Bottom" = "bottom",
                                                           "Inside (custom)" = "inside",
                                                           "None" = "none"),
                                               selected = "right")),
                         column(6, selectInput("legend_text_face", "Text Style:",
                                               choices = c("Plain" = "plain", "Bold" = "bold",
                                                           "Italic" = "italic", "Bold Italic" = "bold.italic"),
                                               selected = "plain"))
                       ),
                       # Inside-plot positioning controls
                       conditionalPanel(
                         condition = "input.legend_position == 'inside'",
                         div(style = "margin-top:6px;",
                           p(style = "font-size:.82em;color:#555;margin-bottom:4px;font-weight:600;",
                             "Corner Presets:"),
                           fluidRow(
                             column(3, actionButton("leg_tl", "Top-L",  style = "font-size:11px;padding:3px 6px;width:100%;")),
                             column(3, actionButton("leg_tr", "Top-R",  style = "font-size:11px;padding:3px 6px;width:100%;")),
                             column(3, actionButton("leg_bl", "Bot-L",  style = "font-size:11px;padding:3px 6px;width:100%;")),
                             column(3, actionButton("leg_br", "Bot-R",  style = "font-size:11px;padding:3px 6px;width:100%;"))
                           ),
                           fluidRow(style = "margin-top:6px;",
                             column(6, numericInput("legend_x", "X position (0-1):", 0.85, 0, 1, 0.05)),
                             column(6, numericInput("legend_y", "Y position (0-1):", 0.95, 0, 1, 0.05))
                           ),
                           checkboxInput("legend_click_mode", "Click plot to place legend", FALSE),
                           p(style = "font-size:.80em;color:#888;margin-top:2px;",
                             "When checked: click anywhere on the plot to move the legend there.")
                         )
                       ),
                       fluidRow(style = "margin-top:4px;",
                         column(12, checkboxInput("legend_no_box", "No legend background/border", FALSE))
                       ),
                       # Reserved space control — shown when legend won't occupy external space
                       conditionalPanel(
                         condition = "input.legend_position == 'none' || input.legend_position == 'inside' || input.show_end_labels == true",
                         fluidRow(style = "margin-top:6px;",
                           column(12, numericInput("legend_reserve_space",
                                                   "Reserved right space (pt):", 130, 0, 400, 10)),
                           column(12, p(style = "font-size:.80em;color:#888;margin-top:0;",
                                        "Keeps graph panel width constant. Increase if legend was wide, decrease if it was narrow."))
                         )
                       ),
                       fluidRow(style = "margin-top:6px;",
                         column(12, sliderInput("legend_wrap_width",
                                                "Wrap labels at (chars):", 8, 60, 20, 1)),
                         column(12, p(style = "font-size:.80em;color:#888;margin-top:0;",
                                      "Long names are wrapped at this width so the graph panel area stays constant. Set to 60 to disable wrapping."))
                       ),
                       p(style = "font-size:.82em;color:#888;margin-top:4px;margin-bottom:0;",
                         "Font size follows Axis Text size in Axis Settings.")
                   )
                 ),

                 # ── Label Options ─────────────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Label Options"),
                   div(class = "panel-section",
                       checkboxInput("show_end_labels", "Show End-of-Line Labels", FALSE),
                       conditionalPanel(condition = "input.show_end_labels == true",
                                        numericInput("label_font_size", "Label Font Size (pt):", 12, 3, 36),
                                        checkboxInput("label_bold",  "Bold Labels",              TRUE),
                                        numericInput("label_offset", "Label Offset (% x-axis):", 3.5, 0, 20)
                       )
                   )
                 ),
                 
                 # ── Error Bars / Shadow ───────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Variability / Error Display"),
                   div(class = "panel-section",
                       selectInput("error_display_mode", "Display Mode:",
                                   choices = c(
                                     "None"                       = "none",
                                     "Error Bars"                 = "bars",
                                     "Shadow / Ribbon"            = "shadow",
                                     "Spaghetti (Replicates)"     = "spaghetti",
                                     "Quantile Bands (IQR + 95%)" = "quantile_bands",
                                     "Jitter Points"              = "jitter",
                                     "Combo (Traces + CI Band)"   = "combo"
                                   ),
                                   selected = "bars"),
                       # ── Stat-based controls (bars / shadow / combo) ─────────────
                       conditionalPanel(
                         condition = "['bars','shadow','combo'].indexOf(input.error_display_mode) !== -1",
                         selectInput("error_type", "Error Statistic:",
                                     choices  = c("SD" = "sd", "SEM" = "sem", "95% CI" = "ci95"),
                                     selected = "sem"),
                         numericInput("error_multiplier", "Error Multiplier:", 1, 0.1, 5, 0.1),
                         checkboxInput("asymmetric_error", "Floor errors at zero", FALSE)
                       ),
                       # ── Error bars specific ──────────────────────────────────────
                       conditionalPanel(
                         condition = "input.error_display_mode == 'bars'",
                         selectInput("error_bar_style", "Bar Style:",
                                     choices  = c("T-shaped"     = "T",
                                                  "Solid lines"  = "solid",
                                                  "Dashed lines" = "dashed"),
                                     selected = "T"),
                         sliderInput("error_bar_width",     "Bar Width:",     0.01, 2, 1,   0.01),
                         sliderInput("error_bar_thickness", "Bar Thickness:", 0.1,  2, 0.8, 0.1),
                         selectInput("error_bar_position", "Position:",
                                     choices  = c("Middle" = "middle", "Dodge" = "dodge"),
                                     selected = "middle")
                       ),
                       # ── Ribbon alpha (shadow + combo) ────────────────────────────
                       conditionalPanel(
                         condition = "['shadow','combo'].indexOf(input.error_display_mode) !== -1",
                         sliderInput("shadow_alpha", "Ribbon Alpha:", 0.05, 0.6, 0.2, 0.05)
                       ),
                       # ── Spaghetti / combo replicate trace alpha ──────────────────
                       conditionalPanel(
                         condition = "['spaghetti','combo'].indexOf(input.error_display_mode) !== -1",
                         sliderInput("spaghetti_alpha", "Replicate Trace Alpha:", 0.05, 0.9, 0.35, 0.05)
                       ),
                       # ── Quantile bands ───────────────────────────────────────────
                       conditionalPanel(
                         condition = "input.error_display_mode == 'quantile_bands'",
                         sliderInput("qband_inner_alpha", "Inner Band Alpha (IQR):",  0.05, 0.8, 0.40, 0.05),
                         sliderInput("qband_outer_alpha", "Outer Band Alpha (95%):",  0.05, 0.5, 0.15, 0.05),
                         p(style = "font-size:.80em;color:#888;margin-top:0;",
                           "Inner band = 25th–75th percentile; outer band = 2.5th–97.5th percentile of replicates.")
                       ),
                       # ── Jitter ───────────────────────────────────────────────────
                       conditionalPanel(
                         condition = "input.error_display_mode == 'jitter'",
                         sliderInput("jitter_alpha", "Point Alpha:",  0.1, 1.0, 0.55, 0.05),
                         sliderInput("jitter_size",  "Point Size:",   0.3, 5.0, 1.5,  0.1),
                         sliderInput("jitter_width", "Jitter Width:", 0.0, 2.0, 0.0,  0.05),
                         p(style = "font-size:.80em;color:#888;margin-top:0;",
                           "Plots each raw replicate observation. Width > 0 adds horizontal jitter.")
                       )
                   )
                 ),
                 
                 # ── Plot Dimensions & Export ──────────────────────────────────────────────
                 tags$details(
                   tags$summary("Plot Dimensions & Export"),
                   div(class = "panel-section",
                       h4("Canvas Size", class = "panel-title"),
                       fluidRow(
                         column(6, numericInput("plot_width",  "Width (px):",  820, 400, 2000)),
                         column(6, numericInput("plot_height", "Height (px):", 600, 300, 1500))
                       ),
                       checkboxInput("custom_aspect_ratio", "Custom Aspect Ratio", FALSE),
                       conditionalPanel(condition = "input.custom_aspect_ratio == true",
                                        numericInput("aspect_ratio", "Width/Height:", 1, 0.3, 4, 0.1)
                       ),
                       h4("Image Export", class = "panel-title"),
                       selectInput("export_format", "Format:",
                                   choices  = c("PDF" = "pdf", "PNG" = "png", "JPEG" = "jpeg",
                                                "SVG" = "svg", "TIFF" = "tiff"),
                                   selected = "pdf"),
                       fluidRow(
                         column(6, numericInput("export_width",  "Width (in):",  10, 2, 30)),
                         column(6, numericInput("export_height", "Height (in):", 8,  2, 20))
                       ),
                       numericInput("export_dpi", "DPI:", 300, 72, 1200),
                       downloadButton("downloadPlot", "Download Image"),
                       hr(),
                       h4("PowerPoint Export", class = "panel-title"),
                       p(style = "font-size:.85em;color:#666;",
                         "Cumulative build: slide 1 shows sample 1, slide 2 adds sample 2, etc. Each slide title shows the newly-added sample name. Final slide = all samples, identical to the main plot."),
                       pptx_ui,
                       hr(),
                       h4("Animated GIF Export", class = "panel-title"),
                       p(style = "font-size:.85em;color:#666;",
                         "Cumulative build: each frame adds one more line. No grey masks — lines already drawn stay fully visible. Axes, theme, and error bars match the main plot exactly."),
                       numericInput("gif_fps",    "Frames per second:", 1,   0.2, 10,   0.2),
                       numericInput("gif_width",  "GIF Width (px):",    800, 300, 2000, 50),
                       numericInput("gif_height", "GIF Height (px):",   600, 200, 1500, 50),
                       gif_ui
                   )
                 ),
                 
                 # ── Variable Styling ──────────────────────────────────────────────────────
                 tags$details(
                   tags$summary("Variable Styling"),
                   div(class = "panel-section",
                       h4("Sample Selection", class = "panel-title"),
                       fluidRow(
                         column(8, selectInput("selected_samples", "Select Samples (click order = legend order):",
                                              choices = character(0), multiple = TRUE)),
                         column(2, actionButton("samples_select_all", "All",
                                                style = "margin-top:25px;width:100%;font-size:.8em;")),
                         column(2, actionButton("samples_deselect_all", "None",
                                                style = "margin-top:25px;width:100%;font-size:.8em;background:#dc3545;color:white;border:none;"))
                       ),
                       uiOutput("var_settings")
                   )
                 )
                 
    ), # end sidebarPanel
    
    mainPanel(width = 9,
              tabsetPanel(id = "main_tabs", type = "tabs",
                          
                          # ── Tab 1: Plot (existing functionality) ──────────────────────────
                          tabPanel("Plot", value = "plot_tab",
                                   uiOutput("plot_container")
                          ),
                          
                          # ── Tab 2: Analysis ───────────────────────────────────────────────
                          tabPanel("Analysis", value = "analysis_tab",
                                   br(),
                                   tabsetPanel(id = "analysis_subtabs", type = "tabs",

                                     # ── Sub-tab 1: Metrics & Stats ─────────────────────────
                                     tabPanel("Metrics & Stats",
                                       br(),
                                       p(class = "wit-label", style="padding-left:4px;",
                                         "Step 1: Calculate metrics. Step 2: Run stats. Step 3: Publish."),

                                       # Growth Curve Metrics
                                       div(class = "panel-section",
                                           h4("Growth Curve Metrics", class = "panel-title"),
                                           p(style = "font-size:.85em;color:#555;",
                                             "Extracts quantitative metrics per sample: growth physiology, ",
                                             "biomass production, lysis dynamics, and recovery."),
                                           tags$details(
                                             tags$summary("What is smoothing?"),
                                             p(style = "font-size:.82em;color:#555;padding:8px;",
                                               "The smoothing window controls how many consecutive time points are used ",
                                               "when calculating slopes (growth rate, lysis rate). A window of 5 means each ",
                                               "slope is computed from 5 adjacent points using linear regression. ",
                                               "Larger windows (7\u201310) reduce noise but can blur rapid transitions like lysis onset. ",
                                               "Smaller windows (3\u20134) capture fast changes but amplify noise. ",
                                               "Default of 5 is good for most plate-reader data. If your data has > 100 time points, try 7\u201310. ",
                                               "If you have < 20 time points, try 3.")
                                           ),
                                           fluidRow(
                                             column(3,
                                                    numericInput("metrics_smooth_window", "Smoothing Window:",
                                                                 value = 5, min = 3, max = 20, step = 1)),
                                             column(3,
                                                    actionButton("calc_metrics", "Calculate Metrics",
                                                                 icon = icon("calculator"),
                                                                 style = "margin-top:25px;background:#2C3E50;color:white;border:none;")),
                                             column(3,
                                                    downloadButton("download_metrics_csv", "Download Metrics CSV",
                                                                   style = "margin-top:25px;width:100%;")),
                                             column(3,
                                                    downloadButton("download_all_stats", "Download All Stats CSV",
                                                                   style = "margin-top:25px;width:100%;background:#E67E22;color:white;border:none;"))
                                           ),
                                           br(),
                                           DT::DTOutput("metrics_table")
                                       ),

                                       # Infection Strength
                                       div(class = "panel-section",
                                           h4("Infection Strength Metrics", class = "panel-title"),
                                           p(style = "font-size:.85em;color:#555;",
                                             "Compare each sample against a reference (uninfected control). ",
                                             "Also works for plasmid expression experiments: use the empty-vector ",
                                             "control as reference to quantify the effect of your expressed gene. ",
                                             "Infection Strength = 1 \u2212 (AUC\u209b\u2090\u2098\u209a\u2097\u2091 / AUC\u1d63\u2091f). ",
                                             "Values > 0 mean growth suppression; < 0 means enhanced growth."),
                                           fluidRow(
                                             column(4,
                                                    selectInput("ref_sample", "Reference (Control) Sample:",
                                                                choices = character(0))),
                                             column(4,
                                                    actionButton("calc_infection", "Compute Infection Metrics",
                                                                 icon = icon("virus"),
                                                                 style = "margin-top:25px;background:#C0392B;color:white;border:none;")),
                                             column(4,
                                                    downloadButton("download_infection_csv", "Download CSV",
                                                                   style = "margin-top:25px;width:100%;"))
                                           ),
                                           br(),
                                           DT::DTOutput("infection_table")
                                       ),

                                       # Statistical Comparisons
                                       div(class = "panel-section",
                                           h4("Statistical Comparisons", class = "panel-title"),
                                           p(style = "font-size:.85em;color:#555;",
                                             "t-test (2 groups) or ANOVA + pairwise t with BH correction (\u22653 groups). ",
                                             "Results feed into significance brackets on the bar plot."),
                                           fluidRow(
                                             column(8, uiOutput("stats_sample_selector")),
                                             column(2, actionButton("stats_select_all", "All",
                                                                    style = "margin-top:25px;width:100%;font-size:.8em;")),
                                             column(2, actionButton("stats_deselect_all", "None",
                                                                    style = "margin-top:25px;width:100%;font-size:.8em;background:#dc3545;color:white;border:none;"))
                                           ),
                                           fluidRow(
                                             column(4,
                                                    selectInput("stats_metric", "Metric to Test:",
                                                                choices = c(core_metric_choices,
                                                                            lysis_metric_choices))),
                                             column(4,
                                                    actionButton("run_stats", "Run Tests",
                                                                 icon = icon("chart-bar"),
                                                                 style = "margin-top:25px;background:#8E44AD;color:white;border:none;")),
                                             column(4,
                                                    downloadButton("download_stats_csv", "Download Stats CSV",
                                                                   style = "margin-top:25px;width:100%;"))
                                           ),
                                           br(),
                                           DT::DTOutput("stats_table")
                                       )
                                     ),

                                     # ── Sub-tab 2: Plots ────────────────────────────────────
                                     tabPanel("Plots",
                                       br(),

                                       # Bar Plot
                                       div(class = "panel-section",
                                           h4("Metric Comparison (Bar Plot)", class = "panel-title"),
                                           p(style = "font-size:.85em;color:#555;",
                                             "Bar plot with individual data points and error bars. ",
                                             "Significance brackets shown when statistical tests have been run."),

                                           # Sample selector
                                           fluidRow(
                                             column(8, uiOutput("barplot_sample_selector")),
                                             column(2, actionButton("barplot_select_all", "All",
                                                                    style="margin-top:25px;width:100%;font-size:.8em;")),
                                             column(2, actionButton("barplot_deselect_all", "None",
                                                                    style="margin-top:25px;width:100%;font-size:.8em;background:#dc3545;color:white;border:none;"))
                                           ),

                                           fluidRow(
                                             column(3,
                                                    selectInput("barplot_type", "Plot Type:",
                                                                choices = c("Bar + Dots" = "bar",
                                                                            "Dot Plot" = "dot",
                                                                            "Box Plot" = "box",
                                                                            "Violin" = "violin"),
                                                                selected = "bar")),
                                             column(3,
                                                    selectInput("barplot_metric", "Select Metric:",
                                                                choices = c(core_metric_choices,
                                                                            lysis_metric_choices))),
                                             column(3,
                                                    selectInput("barplot_error_type", "Error Bars:",
                                                                choices = c("SEM" = "sem", "SD" = "sd",
                                                                            "95% CI" = "ci95", "None" = "none"),
                                                                selected = "sem")),
                                             column(3,
                                                    checkboxInput("barplot_show_sig", "Show significance brackets", TRUE))
                                           ),
                                           fluidRow(
                                             column(3,
                                                    selectInput("barplot_sort_order", "Bar Order:",
                                                                choices = c("Data order" = "default",
                                                                            "Alphabetical" = "alpha",
                                                                            "By value (asc)" = "asc",
                                                                            "By value (desc)" = "desc",
                                                                            "Custom order" = "custom"),
                                                                selected = "default")),
                                             column(3,
                                                    selectInput("barplot_sizing_mode", "Bar Sizing:",
                                                                choices = c("Shrink plot" = "shrink",
                                                                            "Fixed (gaps)" = "fixed"),
                                                                selected = "shrink")),
                                             column(3,
                                                    checkboxInput("barplot_horizontal", "Horizontal bars", FALSE))
                                           ),
                                           conditionalPanel(
                                             condition = "input.barplot_sort_order === 'custom'",
                                             div(style = "padding:4px 15px 8px;",
                                               p(style = "font-size:.82em;color:#555;margin-bottom:4px;",
                                                 "Click samples below in desired left-to-right order (deselect all first, then re-click in order):"),
                                               uiOutput("barplot_custom_order_ui")
                                             )
                                           ),
                                           fluidRow(
                                             column(3,
                                                    checkboxInput("barplot_use_sample_colors",
                                                                  "Auto-match Plot tab colors", TRUE)),
                                             column(3,
                                                    checkboxInput("barplot_custom_colors",
                                                                  "Customize bar colors", FALSE)),
                                             column(3,
                                                    numericInput("barplot_width",  "Width (px):",  700, 400, 1600)),
                                             column(3,
                                                    numericInput("barplot_height", "Height (px):", 500, 300, 1200))
                                           ),
                                           conditionalPanel(
                                             condition = "input.barplot_custom_colors == true",
                                             uiOutput("barplot_color_pickers")
                                           ),

                                           # Style controls
                                           tags$details(
                                             tags$summary("Bar Plot Style"),
                                             div(style="padding:8px 0;",
                                               fluidRow(
                                                 column(4, textInput("barplot_title", "Custom Title:", value = "")),
                                                 column(4, sliderInput("barplot_bar_width_ctrl", "Bar Width:", 0.3, 0.95, 0.7, 0.05)),
                                                 column(4, sliderInput("barplot_bar_alpha_ctrl", "Bar Opacity:", 0.3, 1.0, 0.85, 0.05))
                                               ),
                                               fluidRow(
                                                 column(3, numericInput("barplot_title_size", "Title Size:", 16, 8, 32)),
                                                 column(3, numericInput("barplot_axis_text_size", "Axis Text Size:", 12, 6, 24)),
                                                 column(3, numericInput("barplot_axis_label_size", "Axis Label Size:", 14, 6, 32)),
                                                 column(3, numericInput("barplot_sig_text_size", "Significance Text Size:", 4.5, 2, 10, 0.5))
                                               ),
                                               fluidRow(
                                                 column(3, selectInput("barplot_x_angle", "X Label Angle:",
                                                                       choices = c("0\u00b0" = "0", "30\u00b0" = "30",
                                                                                   "45\u00b0" = "45", "90\u00b0" = "90"),
                                                                       selected = "45")),
                                                 column(3, selectInput("barplot_legend_pos", "Legend:",
                                                                       choices = c("None" = "none", "Right" = "right",
                                                                                   "Top" = "top", "Bottom" = "bottom"),
                                                                       selected = "none")),
                                                 column(3, numericInput("barplot_legend_size", "Legend Font Size:", 12, 6, 24)),
                                                 column(3, selectInput("barplot_legend_face", "Legend Style:",
                                                                       choices = c("Plain" = "plain", "Bold" = "bold",
                                                                                   "Italic" = "italic", "Bold Italic" = "bold.italic"),
                                                                       selected = "plain"))
                                               )
                                             )
                                           ),

                                           fluidRow(
                                             column(3,
                                                    selectInput("barplot_export_fmt", "Export Format:",
                                                                choices = c("PDF"="pdf","PNG"="png","SVG"="svg",
                                                                            "TIFF"="tiff","JPEG"="jpeg"),
                                                                selected = "pdf")),
                                             column(3,
                                                    numericInput("barplot_export_w", "Export Width (in):", 8, 2, 20)),
                                             column(3,
                                                    numericInput("barplot_export_h", "Export Height (in):", 6, 2, 16)),
                                             column(3,
                                                    downloadButton("download_barplot", "Download Bar Plot",
                                                                   style = "margin-top:25px;width:100%;"))
                                           ),
                                           uiOutput("barplot_container")
                                       ),

                                       # Derivative Plot
                                       div(class = "panel-section",
                                           h4("Derivative Plot (dOD/dt)", class = "panel-title"),
                                           p(style = "font-size:.85em;color:#555;",
                                             "Rate of OD change over time. Positive = growth, negative = lysis/decline. ",
                                             "The zero-crossing after the peak marks lysis onset. ",
                                             "For plasmid expression: negative values indicate your gene product is killing cells."),
                                           fluidRow(
                                             column(3, checkboxInput("deriv_show_zero", "Show zero line", TRUE)),
                                             column(3, numericInput("deriv_width",  "Width (px):",  700, 400, 1600)),
                                             column(3, numericInput("deriv_height", "Height (px):", 400, 200, 1000)),
                                             column(3,
                                                    selectInput("deriv_export_fmt", "Export Format:",
                                                                choices = c("PDF"="pdf","PNG"="png","SVG"="svg",
                                                                            "TIFF"="tiff","JPEG"="jpeg"),
                                                                selected = "pdf"))
                                           ),
                                           fluidRow(
                                             column(3, numericInput("deriv_title_size", "Title Size:", 16, 8, 32)),
                                             column(3, numericInput("deriv_axis_text_size", "Axis Text Size:", 12, 6, 24)),
                                             column(3, numericInput("deriv_axis_label_size", "Axis Label Size:", 14, 6, 32)),
                                             column(3, numericInput("deriv_legend_size", "Legend Text Size:", 12, 6, 24))
                                           ),
                                           fluidRow(
                                             column(3, selectInput("deriv_title_face", "Title Style:",
                                                                   choices = c("Bold" = "bold", "Plain" = "plain",
                                                                               "Italic" = "italic", "Bold Italic" = "bold.italic"),
                                                                   selected = "bold")),
                                             column(3, selectInput("deriv_axis_face", "Axis Label Style:",
                                                                   choices = c("Plain" = "plain", "Bold" = "bold",
                                                                               "Italic" = "italic", "Bold Italic" = "bold.italic"),
                                                                   selected = "plain")),
                                             column(3, selectInput("deriv_x_angle", "X-axis Angle:",
                                                                   choices = c("0\u00b0" = "0", "30\u00b0" = "30",
                                                                               "45\u00b0" = "45", "90\u00b0" = "90"),
                                                                   selected = "0")),
                                             column(3, selectInput("deriv_legend_pos", "Legend:",
                                                                   choices = c("Right" = "right", "Left" = "left",
                                                                               "Top" = "top", "Bottom" = "bottom",
                                                                               "None" = "none"),
                                                                   selected = "right"))
                                           ),
                                           fluidRow(
                                             column(3, numericInput("deriv_export_w", "Export Width (in):", 8, 2, 20)),
                                             column(3, numericInput("deriv_export_h", "Export Height (in):", 5, 2, 16)),
                                             column(3,
                                                    downloadButton("download_deriv_plot", "Download Derivative Plot",
                                                                   style = "margin-top:25px;width:100%;"))
                                           ),
                                           uiOutput("deriv_plot_container")
                                       ),

                                       # Annotated Growth Curves
                                       div(class = "panel-section",
                                           h4("Annotated Growth Curves", class = "panel-title"),
                                           p(style = "font-size:.85em;color:#555;",
                                             "Growth curves with phase boundaries: lag end (\u25b2 green), ",
                                             "max OD (\u25c6 blue), lysis/decline onset (\u25bc red)."),
                                           fluidRow(
                                             column(2, checkboxInput("annot_lag",   "Lag end",     TRUE)),
                                             column(2, checkboxInput("annot_mumax", "Max OD",      TRUE)),
                                             column(2, checkboxInput("annot_lysis", "Lysis onset", TRUE)),
                                             column(3, numericInput("annot_width",  "Width (px):",  700, 400, 1600)),
                                             column(3, numericInput("annot_height", "Height (px):", 500, 200, 1200))
                                           ),
                                           fluidRow(
                                             column(3, numericInput("annot_title_size", "Title Size:", 16, 8, 32)),
                                             column(3, numericInput("annot_axis_text_size", "Axis Text Size:", 12, 6, 24)),
                                             column(3, numericInput("annot_axis_label_size", "Axis Label Size:", 14, 6, 32)),
                                             column(3, numericInput("annot_legend_size", "Legend Text Size:", 12, 6, 24))
                                           ),
                                           fluidRow(
                                             column(3, selectInput("annot_title_face", "Title Style:",
                                                                   choices = c("Bold" = "bold", "Plain" = "plain",
                                                                               "Italic" = "italic", "Bold Italic" = "bold.italic"),
                                                                   selected = "bold")),
                                             column(3, selectInput("annot_axis_face", "Axis Label Style:",
                                                                   choices = c("Plain" = "plain", "Bold" = "bold",
                                                                               "Italic" = "italic", "Bold Italic" = "bold.italic"),
                                                                   selected = "plain")),
                                             column(3, selectInput("annot_x_angle", "X-axis Angle:",
                                                                   choices = c("0\u00b0" = "0", "30\u00b0" = "30",
                                                                               "45\u00b0" = "45", "90\u00b0" = "90"),
                                                                   selected = "0")),
                                             column(3, selectInput("annot_legend_pos", "Legend:",
                                                                   choices = c("Right" = "right", "Left" = "left",
                                                                               "Top" = "top", "Bottom" = "bottom",
                                                                               "None" = "none"),
                                                                   selected = "right"))
                                           ),
                                           fluidRow(
                                             column(3,
                                                    selectInput("annot_export_fmt", "Export Format:",
                                                                choices = c("PDF"="pdf","PNG"="png","SVG"="svg",
                                                                            "TIFF"="tiff","JPEG"="jpeg"),
                                                                selected = "pdf")),
                                             column(3, numericInput("annot_export_w", "Export Width (in):", 8, 2, 20)),
                                             column(3, numericInput("annot_export_h", "Export Height (in):", 6, 2, 16)),
                                             column(3,
                                                    downloadButton("download_annot_plot", "Download",
                                                                   style = "margin-top:25px;width:100%;"))
                                           ),
                                           uiOutput("annot_plot_container")
                                       )
                                     ),

                                     # ── Sub-tab 3: Heatmaps ─────────────────────────────────
                                     tabPanel("Heatmaps",
                                       br(),

                                       # Phenotype Heatmap
                                       div(class = "panel-section",
                                           h4("Phenotype Heatmap", class = "panel-title"),
                                           p(style = "font-size:.85em;color:#555;",
                                             "Z-score normalised heatmap: rows = samples, columns = metrics. ",
                                             "Select which samples and metrics to include."),
                                           fluidRow(
                                             column(8, uiOutput("heatmap_sample_selector")),
                                             column(2, actionButton("heatmap_select_all", "All",
                                                                    style="margin-top:25px;width:100%;font-size:.8em;")),
                                             column(2, actionButton("heatmap_deselect_all", "None",
                                                                    style="margin-top:25px;width:100%;font-size:.8em;background:#dc3545;color:white;border:none;"))
                                           ),
                                           uiOutput("heatmap_metric_checkboxes"),
                                           fluidRow(
                                             column(3, numericInput("heatmap_width",  "Width (px):",  700, 400, 1400)),
                                             column(3, numericInput("heatmap_height", "Height (px):", 400, 200, 1000)),
                                             column(3,
                                                    selectInput("heatmap_export_fmt", "Export Format:",
                                                                choices = c("PDF"="pdf","PNG"="png","SVG"="svg",
                                                                            "TIFF"="tiff","JPEG"="jpeg"),
                                                                selected = "pdf")),
                                             column(3,
                                                    downloadButton("download_heatmap", "Download Heatmap",
                                                                   style = "margin-top:25px;width:100%;"))
                                           ),
                                           fluidRow(
                                             column(3, numericInput("heatmap_export_w", "Export Width (in):", 10, 2, 20)),
                                             column(3, numericInput("heatmap_export_h", "Export Height (in):", 6, 2, 16))
                                           ),
                                           uiOutput("heatmap_container")
                                       ),

                                       # OD Over Time Heatmap
                                       div(class = "panel-section",
                                           h4("OD Over Time Heatmap", class = "panel-title"),
                                           p(style = "font-size:.85em;color:#555;",
                                             "Samples \u00d7 time heatmap colored by OD value. Shows the full ",
                                             "growth/lysis dynamics in a single image. Bright = high OD, dark = low OD."),
                                           fluidRow(
                                             column(5, uiOutput("od_heatmap_sample_selector")),
                                             column(2, actionButton("od_heatmap_select_all", "All",
                                                                    style="margin-top:25px;width:100%;font-size:.8em;")),
                                             column(2, actionButton("od_heatmap_deselect_all", "None",
                                                                    style="margin-top:25px;width:100%;font-size:.8em;background:#dc3545;color:white;border:none;")),
                                             column(3, selectInput("od_heatmap_palette", "Palette:",
                                                                   choices = c("Inferno"="inferno","Viridis"="viridis",
                                                                               "Plasma"="plasma","Magma"="magma"),
                                                                   selected = "inferno"))
                                           ),
                                           fluidRow(
                                             column(3, numericInput("od_heatmap_width",  "Width (px):",  800, 400, 1600)),
                                             column(3, numericInput("od_heatmap_height", "Height (px):", 400, 200, 1000)),
                                             column(3,
                                                    selectInput("od_heatmap_export_fmt", "Export Format:",
                                                                choices = c("PDF"="pdf","PNG"="png","SVG"="svg",
                                                                            "TIFF"="tiff","JPEG"="jpeg"),
                                                                selected = "pdf")),
                                             column(3,
                                                    downloadButton("download_od_heatmap", "Download OD Heatmap",
                                                                   style = "margin-top:25px;width:100%;"))
                                           ),
                                           fluidRow(
                                             column(3, numericInput("od_heatmap_export_w", "Export Width (in):", 12, 2, 24)),
                                             column(3, numericInput("od_heatmap_export_h", "Export Height (in):", 5, 2, 16))
                                           ),
                                           uiOutput("od_heatmap_container")
                                       )
                                     )

                                   ) # end tabsetPanel
                          ), # end Analysis tabPanel

                          # ── Tab 3: Experiment Notes ───────────────────────────────────────
                          tabPanel("Experiment Notes", value = "notes_tab",
                            br(),
                            fluidRow(

                              # ── Left column: input fields ─────────────────────────────────
                              column(7,

                                # Experiment Identity
                                div(class = "panel-section",
                                  tags$h5(style = "margin-top:0;font-weight:700;color:#444;border-bottom:1px solid #ddd;padding-bottom:4px;",
                                          "Experiment Identity"),
                                  fluidRow(
                                    column(4, textInput("notes_exp_id",       "Experiment ID / Name", "")),
                                    column(4, dateInput("notes_date",         "Date", value = Sys.Date())),
                                    column(4, textInput("notes_experimenter", "Experimenter",         ""))
                                  ),
                                  fluidRow(
                                    column(6, textInput("notes_project",    "Project / Campaign", "")),
                                    column(6, textInput("notes_institution","Lab / Institution",  ""))
                                  )
                                ),

                                # Biology
                                div(class = "panel-section", style = "margin-top:10px;",
                                  tags$h5(style = "margin-top:0;font-weight:700;color:#444;border-bottom:1px solid #ddd;padding-bottom:4px;",
                                          "Biological Parameters"),
                                  fluidRow(
                                    column(4, selectInput("notes_exp_type", "Experiment Type",
                                                          choices = c("Phage Infection", "Plasmid Induction",
                                                                      "Growth Curve", "Lysis Curve",
                                                                      "Dose-Response", "Other"),
                                                          selected = "Phage Infection")),
                                    column(4, textInput("notes_host_strain",  "Host Strain(s)", "")),
                                    column(4, textInput("notes_phage_plasmid","Phage / Plasmid", ""))
                                  ),
                                  fluidRow(
                                    column(4, textInput("notes_moi",       "MOI",            "")),
                                    column(4, textInput("notes_replicate", "Replicate / Batch", "")),
                                    column(4, textInput("notes_passage",   "Passage #",      ""))
                                  )
                                ),

                                # Conditions
                                div(class = "panel-section", style = "margin-top:10px;",
                                  tags$h5(style = "margin-top:0;font-weight:700;color:#444;border-bottom:1px solid #ddd;padding-bottom:4px;",
                                          "Conditions"),
                                  fluidRow(
                                    column(4, textInput("notes_media",       "Growth Medium",         "")),
                                    column(4, textInput("notes_temperature", "Temperature (°C)",      "")),
                                    column(4, textInput("notes_time_inf",    "Time of Infection/Induction", ""))
                                  ),
                                  fluidRow(
                                    column(4, textInput("notes_inducer",      "Inducer",               "")),
                                    column(4, textInput("notes_inducer_conc", "Inducer Concentration", "")),
                                    column(4, textInput("notes_extra_cond",   "Other Condition",       ""))
                                  )
                                ),

                                # Free notes
                                div(class = "panel-section", style = "margin-top:10px;",
                                  tags$h5(style = "margin-top:0;font-weight:700;color:#444;border-bottom:1px solid #ddd;padding-bottom:4px;",
                                          "Notes"),
                                  textAreaInput("notes_observations", "Observations / Results",
                                                "", rows = 3, width = "100%"),
                                  textAreaInput("notes_issues",       "Issues / Anomalies",
                                                "", rows = 2, width = "100%"),
                                  textAreaInput("notes_next_steps",   "Next Steps",
                                                "", rows = 2, width = "100%"),
                                  textInput("notes_tags", "Tags (comma-separated)",
                                            placeholder = "e.g. pilot, MG1655, T4, low-MOI")
                                )
                              ), # end left column

                              # ── Right column: integration & export ────────────────────────
                              column(5,

                                # Graph caption integration
                                div(class = "panel-section",
                                  tags$h5(style = "margin-top:0;font-weight:700;color:#444;border-bottom:1px solid #ddd;padding-bottom:4px;",
                                          "Show on Graph"),
                                  checkboxInput("notes_show_caption", "Add experiment info as graph caption", FALSE),
                                  conditionalPanel(condition = "input.notes_show_caption == true",
                                    checkboxGroupInput("notes_caption_fields",
                                      "Fields to include in caption:",
                                      choices = c(
                                        "Experiment ID"         = "exp_id",
                                        "Date"                  = "date",
                                        "Experimenter"          = "experimenter",
                                        "Experiment Type"       = "exp_type",
                                        "Host Strain(s)"        = "host_strain",
                                        "Phage / Plasmid"       = "phage_plasmid",
                                        "MOI"                   = "moi",
                                        "Replicate / Batch"     = "replicate",
                                        "Inducer"               = "inducer",
                                        "Inducer Concentration" = "inducer_conc",
                                        "Growth Medium"         = "media",
                                        "Temperature"           = "temperature"
                                      ),
                                      selected = c("exp_id", "exp_type", "host_strain",
                                                   "phage_plasmid", "moi", "replicate")
                                    ),
                                    numericInput("notes_caption_size", "Caption font size (pt):", 9, 6, 16),
                                    div(style = "background:#f8f8f8;border:1px solid #ddd;border-radius:4px;padding:8px;margin-top:6px;",
                                      tags$strong(style = "font-size:.82em;color:#555;", "Caption preview:"),
                                      textOutput("notes_caption_preview")
                                    )
                                  )
                                ),

                                # PPTX integration
                                div(class = "panel-section", style = "margin-top:10px;",
                                  tags$h5(style = "margin-top:0;font-weight:700;color:#444;border-bottom:1px solid #ddd;padding-bottom:4px;",
                                          "PowerPoint Integration"),
                                  checkboxInput("notes_pptx_slide",
                                                "Add experiment notes slide to PPTX download", TRUE)
                                ),

                                # Downloads
                                div(class = "panel-section", style = "margin-top:10px;",
                                  tags$h5(style = "margin-top:0;font-weight:700;color:#444;border-bottom:1px solid #ddd;padding-bottom:4px;",
                                          "Export Notes"),
                                  p(style = "font-size:.82em;color:#666;margin-bottom:8px;",
                                    "Download your experiment record in any format."),
                                  fluidRow(
                                    column(6,
                                      downloadButton("download_notes_txt", "TXT",
                                                     style = "width:100%;margin-bottom:6px;")),
                                    column(6,
                                      downloadButton("download_notes_csv", "CSV (log row)",
                                                     style = "width:100%;margin-bottom:6px;"))
                                  ),
                                  fluidRow(
                                    column(6,
                                      downloadButton("download_notes_pdf", "PDF",
                                                     style = "width:100%;margin-bottom:6px;")),
                                    column(6,
                                      downloadButton("download_notes_html", "HTML",
                                                     style = "width:100%;margin-bottom:6px;"))
                                  ),
                                  fluidRow(
                                    column(6,
                                      downloadButton("download_notes_json", "JSON",
                                                     style = "width:100%;margin-bottom:6px;")),
                                    column(6,
                                      downloadButton("download_notes_png", "Image (PNG)",
                                                     style = "width:100%;margin-bottom:6px;"))
                                  ),
                                  p(style = "font-size:.80em;color:#888;margin-top:4px;",
                                    tags$b("TXT:"), " readable summary. ",
                                    tags$b("CSV:"), " one-row log for a master spreadsheet. ",
                                    tags$b("PDF/PNG:"), " shareable document/photo. ",
                                    tags$b("HTML:"), " open in browser & print. ",
                                    tags$b("JSON:"), " machine-readable record.")
                                ),

                                # Custom key-value pairs for extensible metadata
                                div(class = "panel-section", style = "margin-top:10px;",
                                  tags$h5(style = "margin-top:0;font-weight:700;color:#444;border-bottom:1px solid #ddd;padding-bottom:4px;",
                                          "Custom Fields"),
                                  p(style = "font-size:.82em;color:#666;margin-bottom:6px;",
                                    "Add any extra key-value metadata (e.g. 'Antibiotic: Amp 100 µg/mL')."),
                                  fluidRow(
                                    column(5, textInput("notes_custom_key1",   "Key",   "", placeholder = "Field name")),
                                    column(7, textInput("notes_custom_val1",   "Value", "", placeholder = "Value"))
                                  ),
                                  fluidRow(
                                    column(5, textInput("notes_custom_key2",   "Key",   "", placeholder = "Field name")),
                                    column(7, textInput("notes_custom_val2",   "Value", "", placeholder = "Value"))
                                  ),
                                  fluidRow(
                                    column(5, textInput("notes_custom_key3",   "Key",   "", placeholder = "Field name")),
                                    column(7, textInput("notes_custom_val3",   "Value", "", placeholder = "Value"))
                                  )
                                )
                              ) # end right column
                            ) # end fluidRow
                          ) # end Experiment Notes tabPanel
              )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  rv <- reactiveValues(
    data              = NULL,
    time_col          = NULL,
    od_vars           = NULL,
    od_vars_raw       = NULL,
    group_col         = NULL,
    value_col         = NULL,
    rep_col           = NULL,
    wide_block_id     = NULL,
    wide_rep_count    = NULL,
    is_long_format    = FALSE,
    palette_colors    = NULL,
    all_timepoints    = NULL,
    imported_settings = NULL,
    settings_status   = NULL,
    style_page        = 1L
  )
  
  # ── Palette helpers ──────────────────────────────────────────────────────────
  get_palette_colors <- function(pal, n) {
    switch(pal,
           viridis     = scales::viridis_pal(option = "D")(n),
           plasma      = scales::viridis_pal(option = "C")(n),
           colorblind  = rep_len(c("#000000","#E69F00","#56B4E9","#009E73",
                                   "#F0E442","#0072B2","#D55E00","#CC79A7"), n),
           publication = rep_len(c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                                   "#FF7F00","#FFFF33","#A65628","#F781BF"), n),
           rainbow     = rainbow(n),
           gray        = gray.colors(n, start = 0.15, end = 0.85),
           NULL
    )
  }
  
  # ── Day / Night Mode CSS ─────────────────────────────────────────────────────
  output$night_mode_css <- renderUI({
    if (!isTRUE(input$night_mode)) return(NULL)
    tags$style(HTML("
      /* ── Night mode: main panel & content area ─────────────────── */
      body { background-color: #121212 !important; color: #e0e0e0 !important; }
      .main-panel, .col-sm-9, .shiny-bound-output,
      .tab-content, .shiny-tab-content,
      .tabbable > .tab-content { background-color: #121212 !important; color: #e0e0e0 !important; }
      .panel-section {
        background: #1e2d3d !important;
        border-color: #2c3e50 !important;
        color: #e0e0e0 !important;
      }
      .panel-title { color: #4ecdc4 !important; }
      h4, h3, h2, h1, p, label, .control-label, span { color: #e0e0e0 !important; }
      .nav-tabs { border-color: #2c3e50 !important; }
      .nav-tabs > li > a {
        background-color: #1e2d3d !important;
        color: #c8d8e8 !important;
        border-color: #2c3e50 !important;
      }
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:hover {
        background-color: #121212 !important;
        color: #4ecdc4 !important;
        border-bottom-color: #121212 !important;
      }
      details > summary {
        background: #1e2d3d !important;
        color: #c8d8e8 !important;
      }
      .form-control {
        background-color: #1e2d3d !important;
        color: #e0e0e0 !important;
        border-color: #2c3e50 !important;
      }
      .dataTables_wrapper, table.dataTable { background: #1e2d3d !important; color: #e0e0e0 !important; }
      table.dataTable thead th { background: #243447 !important; color: #c8d8e8 !important; }
      table.dataTable.stripe tbody tr.odd,
      table.dataTable.stripe tbody tr.even { background-color: #1e2d3d !important; }
      .btn-default { background:#243447 !important; color:#c8d8e8 !important; border-color:#3a4f6e !important; }
      .shiny-notification { background:#1e2d3d !important; color:#e0e0e0 !important; }
      .wit-label { color: #7a9bb5 !important; }
    "))
  })

  observeEvent(input$color_palette, {
    rv$palette_colors <- if (input$color_palette != "custom" && !is.null(rv$od_vars))
      get_palette_colors(input$color_palette, length(rv$od_vars)) else NULL
  })
  
  output$palette_preview <- renderUI({
    if (input$color_palette == "custom" || is.null(rv$od_vars)) return(NULL)
    cols <- get_palette_colors(input$color_palette, min(8, length(rv$od_vars)))
    div(p("Preview:"),
        div(style = "margin-top:6px;",
            lapply(cols, function(col)
              div(style = paste0("display:inline-block;width:28px;height:20px;",
                                 "background:", col, ";margin-right:4px;",
                                 "border:1px solid #ccc;")))))
  })
  
  # ── Time filter UI ───────────────────────────────────────────────────────────
  output$time_range_slider_ui <- renderUI({
    tp <- rv$all_timepoints
    if (is.null(tp) || length(tp) < 2)
      return(p(style = "color:#999;font-size:.85em;", "Load a data file to enable the range slider."))
    t_min  <- min(tp); t_max <- max(tp)
    t_step <- max(min(diff(tp)), 0.001)
    sliderInput("time_filter_range", label = NULL,
                min = t_min, max = t_max, value = c(t_min, t_max),
                step = t_step, width = "100%")
  })
  
  output$excluded_timepoints_preview <- renderUI({
    if (is.null(rv$all_timepoints)) return(NULL)
    excl <- parse_excluded_timepoints(input$exclude_timepoints, rv$all_timepoints)
    if (length(excl) == 0) return(NULL)
    div(class = "excl-preview",
        paste0("Will exclude ", length(excl), " point(s): ", paste(sort(excl), collapse = ", ")))
  })
  
  # ── Data loading ─────────────────────────────────────────────────────────────
  is_long_format_detect <- function(data) {
    if (nrow(data) == 0) return(FALSE)
    potential <- c("variable","Variable","condition","Condition","treatment","Treatment",
                   "sample","Sample","group","Group")
    if (any(potential %in% colnames(data))) return(TRUE)
    if (ncol(data) <= 4) {
      vd <- sapply(data, function(x) length(unique(x)) / nrow(data))
      return(any(vd < 0.3 & vd > 0))
    }
    FALSE
  }
  
  detect_time_column <- function(data) {
    hits <- grep("^(time|Time|TIME|t|T)$", colnames(data), value = TRUE)
    if (length(hits)) return(hits[1])
    hits2 <- grep("time|Time", colnames(data), value = TRUE)
    if (length(hits2)) return(hits2[1])
    for (nm in colnames(data)) {
      v <- suppressWarnings(as.numeric(data[[nm]]))
      if (!any(is.na(v)) && length(v) > 1 && all(diff(v) >= 0)) return(nm)
    }
    colnames(data)[1]
  }
  
  detect_od_columns <- function(data, exclude_cols) {
    cols <- setdiff(colnames(data), exclude_cols)
    hits <- grep("od|OD|absorbance|Absorbance|value|Value", cols,
                 ignore.case = TRUE, value = TRUE)
    if (length(hits)) return(hits)
    cols[sapply(data[cols], is.numeric)]
  }
  
  observeEvent(input$file, {
    req(input$file)
    data <- tryCatch(
      read.csv(input$file$datapath, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e)
        read.csv(input$file$datapath, stringsAsFactors = FALSE, check.names = TRUE)
    )

    # ── Sanitise raw CSV ──────────────────────────────────────────────────────
    # 1. Drop phantom columns produced by trailing commas on every row
    #    (e.g. "Time,A,B," creates an unnamed "" column — discard it)
    blank_cols <- nchar(trimws(colnames(data))) == 0
    if (any(blank_cols)) data <- data[, !blank_cols, drop = FALSE]

    # 2. Strip rows whose time value is not a finite number.
    #    Catches: trailing empty rows (from Excel/CSV padding), stray "*" rows,
    #    and any other non-numeric sentinel.  MUST run before infer_wide_replicates()
    #    because that function counts rows to assign block/replicate IDs — if junk
    #    rows are present it inflates the last block, and the ID vector then
    #    mismatches after apply_time_filters() strips those rows later → crash.
    tc_guess <- grep("^(time|Time|TIME|t|T)$", colnames(data), value = TRUE)
    if (length(tc_guess) == 0) tc_guess <- colnames(data)[1]
    t_num_tmp <- suppressWarnings(as.numeric(as.character(data[[tc_guess[1]]])))
    data <- data[is.finite(t_num_tmp), , drop = FALSE]
    # ─────────────────────────────────────────────────────────────────────────

    rv$is_long_format <- is_long_format_detect(data)
    rv$time_col       <- detect_time_column(data)
    
    if (rv$is_long_format) {
      potential <- c("variable","Variable","condition","Condition","treatment","Treatment",
                     "sample","Sample","group","Group")
      found     <- potential[potential %in% colnames(data)]
      group_col <- if (length(found) > 0) found[1] else {
        dv    <- sapply(data, function(x) length(unique(x)) / nrow(data))
        cands <- setdiff(names(which(dv < 0.3 & dv > 0)), rv$time_col)
        if (length(cands)) cands[1] else setdiff(colnames(data), rv$time_col)[1]
      }
      vc <- detect_od_columns(data, c(rv$time_col, group_col))
      if (length(vc) > 1) vc <- vc[1]
      rv$group_col   <- group_col
      rv$value_col   <- vc
      rv$rep_col     <- detect_replicate_column(data, rv$time_col, group_col, vc)
      rv$od_vars     <- as.character(unique(data[[group_col]]))
      rv$od_vars_raw <- rv$od_vars
      rv$wide_block_id  <- NULL
      rv$wide_rep_count <- NULL
    } else {
      oc <- detect_od_columns(data, rv$time_col)
      rv$od_vars     <- oc
      rv$od_vars_raw <- oc
      rv$rep_col     <- NULL
      wide_rep <- infer_wide_replicates(data, rv$time_col)
      rv$wide_block_id  <- wide_rep$block_id
      rv$wide_rep_count <- wide_rep$reps
    }
    rv$data <- data
    
    tr <- range(suppressWarnings(as.numeric(as.character(data[[rv$time_col]]))), na.rm = TRUE)
    updateNumericInput(session, "x_min", value = tr[1])
    updateNumericInput(session, "x_max", value = tr[2])
    updateSelectInput(session, "selected_samples", choices = rv$od_vars, selected = rv$od_vars)

    all_t <- sort(unique(suppressWarnings(as.numeric(as.character(data[[rv$time_col]])))))
    rv$all_timepoints <- all_t[is.finite(all_t)]
    updateTextInput(session, "exclude_timepoints", value = "")
    
    if (!is.null(input$color_palette) && input$color_palette != "custom")
      rv$palette_colors <- get_palette_colors(input$color_palette, length(rv$od_vars))
    
    apply_imported_sample_aesthetics()
  })

  # ── Select All / Deselect All for Plot tab sample selector ─────────────────
  observeEvent(input$samples_select_all, {
    updateSelectInput(session, "selected_samples", selected = rv$od_vars)
  })
  observeEvent(input$samples_deselect_all, {
    updateSelectInput(session, "selected_samples", selected = character(0))
  })

  # ── Settings: Save ───────────────────────────────────────────────────────────
  output$saveSettings <- downloadHandler(
    filename = function()
      paste0("plot_settings_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"),
    content = function(file) {
      all_inputs <- reactiveValuesToList(input)
      data_inputs     <- c("file", "loadSettings", "selected_samples",
                           "time_filter_range", "exclude_timepoints")
      sample_prefixes <- c("line_type_","shape_","shape_filled_","color_selector_",
                           "color_","use_rgb_","red_","green_","blue_","alpha_",
                           "use_hex_","hex_color_","legend_label_",
                           "apply_rgb_","apply_hex_")
      is_sample_inp <- function(nm) any(vapply(sample_prefixes, startsWith, logical(1), x = nm))
      
      global_settings <- all_inputs[vapply(names(all_inputs), function(nm)
        !(nm %in% data_inputs) && !is_sample_inp(nm), logical(1))]
      
      sample_aesthetics <- list()
      samps <- isolate(input$selected_samples)
      if (!is.null(samps)) {
        for (vn in samps) {
          vid     <- safe_id(vn)
          leg     <- input[[paste0("legend_label_", vid)]]
          col_sel <- input[[paste0("color_selector_", vid)]]
          color   <- if (!is.null(col_sel) && col_sel == "custom")
            normalize_hex_color(input[[paste0("color_", vid)]])
          else if (!is.null(col_sel)) col_sel
          else "#000000"
          sample_aesthetics[[vn]] <- list(
            color        = color,
            shape        = input[[paste0("shape_",        vid)]],
            shape_filled = input[[paste0("shape_filled_", vid)]],
            line_type    = input[[paste0("line_type_",    vid)]],
            legend_label = if (!is.null(leg) && nchar(leg) > 0) leg else vn
          )
        }
      }
      
      out <- list(
        version            = "2.5",
        saved_at           = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        global_settings    = global_settings,
        sample_aesthetics  = sample_aesthetics,
        experiment_notes   = get_notes()
      )
      write_json(out, file, pretty = TRUE, auto_unbox = TRUE)
    }
  )
  
  # ── Settings: Load ───────────────────────────────────────────────────────────
  observeEvent(input$loadSettings, {
    req(input$loadSettings)
    s <- tryCatch(
      fromJSON(input$loadSettings$datapath, simplifyVector = FALSE),
      error = function(e) {
        rv$settings_status <- list(type = "warn",
                                   msg = paste("Could not parse settings file:", e$message), details = NULL)
        NULL
      }
    )
    if (is.null(s)) return()
    
    global     <- if (!is.null(s$global_settings))   s$global_settings   else list()
    sample_aes <- if (!is.null(s$sample_aesthetics)) s$sample_aesthetics
    else if (!is.null(s$samples))       s$samples
    else list()
    
    rv$imported_settings <- list(global = global, samples = sample_aes)
    
    data_inputs     <- c("file","loadSettings","selected_samples",
                         "time_filter_range","exclude_timepoints")
    sample_prefixes <- c("line_type_","shape_","shape_filled_","color_selector_",
                         "color_","use_rgb_","red_","green_","blue_","alpha_",
                         "use_hex_","hex_color_","legend_label_","apply_rgb_","apply_hex_")
    is_sample_inp <- function(nm) any(vapply(sample_prefixes, startsWith, logical(1), x = nm))
    
    n_global <- 0L
    for (nm in names(global)) {
      if (nm %in% data_inputs || is_sample_inp(nm)) next
      if (!(nm %in% names(input)))                  next
      tryCatch({
        val <- global[[nm]]
        if (is.logical(val) && length(val) == 1) {
          updateCheckboxInput(session, nm, value = val)
        } else if (is.numeric(val) && length(val) == 1) {
          # Call both — sliderInput needs updateSliderInput, numericInput needs
          # updateNumericInput; each silently ignores the wrong message type.
          updateNumericInput(session, nm, value = val)
          updateSliderInput( session, nm, value = val)
        } else if (is.character(val) && length(val) == 1) {
          updateSelectInput(session, nm, selected = val)
          updateTextInput(session,   nm, value    = val)
        } else if (is.character(val) && length(val) > 1) {
          # checkboxGroupInput stores a character vector
          updateCheckboxGroupInput(session, nm, selected = val)
        }
        n_global <- n_global + 1L
      }, error = function(e) invisible(NULL))
    }
    
    apply_imported_sample_aesthetics()

    # Restore experiment notes if present in the settings file
    if (!is.null(s$experiment_notes)) {
      en <- s$experiment_notes
      restore_text <- function(id, val) {
        if (!is.null(val) && nchar(trimws(val)) > 0)
          updateTextInput(session, id, value = trimws(val))
      }
      restore_text("notes_exp_id",       en$exp_id)
      restore_text("notes_experimenter", en$experimenter)
      restore_text("notes_project",      en$project)
      restore_text("notes_institution",  en$institution)
      restore_text("notes_host_strain",  en$host_strain)
      restore_text("notes_phage_plasmid",en$phage_plasmid)
      restore_text("notes_moi",          en$moi)
      restore_text("notes_replicate",    en$replicate)
      restore_text("notes_passage",      en$passage)
      restore_text("notes_media",        en$media)
      restore_text("notes_temperature",  en$temperature)
      restore_text("notes_time_inf",     en$time_inf)
      restore_text("notes_inducer",      en$inducer)
      restore_text("notes_inducer_conc", en$inducer_conc)
      restore_text("notes_extra_cond",   en$extra_cond)
      restore_text("notes_tags",         en$tags)
      restore_text("notes_custom_key1",  en$custom_key1)
      restore_text("notes_custom_val1",  en$custom_val1)
      restore_text("notes_custom_key2",  en$custom_key2)
      restore_text("notes_custom_val2",  en$custom_val2)
      restore_text("notes_custom_key3",  en$custom_key3)
      restore_text("notes_custom_val3",  en$custom_val3)
      if (!is.null(en$observations)) updateTextAreaInput(session, "notes_observations", value = en$observations)
      if (!is.null(en$issues))       updateTextAreaInput(session, "notes_issues",       value = en$issues)
      if (!is.null(en$next_steps))   updateTextAreaInput(session, "notes_next_steps",   value = en$next_steps)
      if (!is.null(en$exp_type))     updateSelectInput(session,   "notes_exp_type",     selected = en$exp_type)
      if (!is.null(en$date)) {
        d <- tryCatch(as.Date(en$date), error = function(e) NULL)
        if (!is.null(d)) updateDateInput(session, "notes_date", value = d)
      }
    }

    cur_samps   <- isolate(input$selected_samples)
    saved_names <- names(sample_aes)
    matched     <- intersect(cur_samps,  saved_names)
    unmatched   <- setdiff(saved_names,  cur_samps)
    
    sample_msg <- if (length(saved_names) == 0) "No sample styles in file."
    else paste0(length(matched), "/", length(saved_names), " sample style(s) applied.")
    
    details_html <- if (length(saved_names) > 0) {
      match_lines   <- if (length(matched)   > 0)
        paste0("<li style='color:#155724'>&#10003; ", matched, "</li>", collapse = "")
      else ""
      unmatch_lines <- if (length(unmatched) > 0)
        paste0("<li style='color:#856404'>&#9675; ", unmatched,
               " <em>(not in current data)</em></li>", collapse = "")
      else ""
      paste0("<ul class='match-list'>", match_lines, unmatch_lines, "</ul>")
    } else NULL
    
    rv$settings_status <- list(
      type    = "ok",
      msg     = paste0(n_global, " global setting(s) applied. ", sample_msg),
      details = details_html
    )
  })
  
  # Applies imported per-sample aesthetics to any matching currently-loaded samples.
  apply_imported_sample_aesthetics <- function() {
    imp <- rv$imported_settings
    if (is.null(imp) || is.null(imp$samples) || length(imp$samples) == 0) return()
    cur_samps <- isolate(input$selected_samples)
    if (is.null(cur_samps) || length(cur_samps) == 0) return()
    for (vn in cur_samps) {
      if (!vn %in% names(imp$samples)) next
      aes   <- imp$samples[[vn]]
      vid   <- safe_id(vn)
      if (!is.null(aes$color)) {
        updateSelectInput(session, paste0("color_selector_", vid), selected = "custom")
        updateTextInput(session,   paste0("color_",          vid), value   = aes$color)
        updateTextInput(session,   paste0("hex_color_",      vid), value   = aes$color)
      }
      if (!is.null(aes$shape))
        updateSelectInput(session, paste0("shape_", vid), selected = as.character(aes$shape))
      if (!is.null(aes$shape_filled))
        updateCheckboxInput(session, paste0("shape_filled_", vid), value = aes$shape_filled)
      if (!is.null(aes$line_type))
        updateSelectInput(session, paste0("line_type_", vid), selected = aes$line_type)
      if (!is.null(aes$legend_label))
        updateTextInput(session, paste0("legend_label_", vid), value = aes$legend_label)
    }
  }
  
  # ── Settings: Clear ──────────────────────────────────────────────────────────
  observeEvent(input$clearSettings, {
    if (is.null(rv$imported_settings)) {
      rv$settings_status <- list(type = "warn", msg = "No imported settings to clear.", details = NULL)
      return()
    }
    
    defs <- list(
      x_scale_type          = "linear",  y_scale_type           = "log",
      use_advanced_ticks    = TRUE,       y_log_min_exponent     = -2,
      y_log_max_exponent    = 0,          custom_x_limits        = FALSE,
      custom_y_limits       = FALSE,      x_expand_left          = 0,
      x_expand_right        = 0.05,       y_expand_bottom        = 0,
      y_expand_top          = 0.05,       x_axis_label           = "Time (minutes)",
      y_axis_label          = "A550",     plot_title             = "",
      plot_subtitle         = "",         show_major_gridlines   = FALSE,
      show_minor_gridlines  = FALSE,      font_family            = "sans",
      title_font_size       = 20,         axis_label_font_size   = 20,
      axis_text_font_size   = 16,         bold_title             = TRUE,
      italic_axis_labels    = FALSE,      axis_text_angle        = "0",
      enable_highlighting   = FALSE,      enable_time_markers    = FALSE,
      color_palette         = "custom",   line_thickness         = 1,
      show_points           = TRUE,       shape_size             = 3,
      point_stroke          = 0.5,        show_end_labels        = FALSE,
      label_font_size       = 12,         label_bold             = TRUE,
      label_offset          = 3.5,        legend_x               = 0.85,
      legend_y              = 0.95,       legend_no_box          = FALSE,
      legend_click_mode     = FALSE,      legend_reserve_space   = 130,
      error_type             = "sem",
      error_display_mode    = "bars",     error_multiplier       = 1,
      asymmetric_error      = FALSE,      error_bar_style        = "T",
      error_bar_width       = 1,          error_bar_thickness    = 0.8,
      error_bar_position    = "middle",   shadow_alpha           = 0.2,
      spaghetti_alpha       = 0.35,       qband_inner_alpha      = 0.40,
      qband_outer_alpha     = 0.15,       jitter_alpha           = 0.55,
      jitter_size           = 1.5,        jitter_width           = 0.0,
      legend_wrap_width     = 20L,
      custom_aspect_ratio   = FALSE,      aspect_ratio           = 1,
      export_format         = "pdf",      export_width           = 10,
      export_height         = 8,          export_dpi             = 300,
      plot_width            = 820,        plot_height            = 600,
      gif_fps               = 1,          gif_width              = 800,
      gif_height            = 600,        x_extra_ticks          = ""
    )
    for (nm in names(defs)) {
      tryCatch({
        val <- defs[[nm]]
        if (is.logical(val)) {
          updateCheckboxInput(session, nm, value = val)
        } else if (is.numeric(val)) {
          updateNumericInput(session, nm, value = val)
          updateSliderInput( session, nm, value = val)
        } else {
          updateSelectInput(session, nm, selected = val)
          updateTextInput(session,   nm, value    = val)
        }
      }, error = function(e) invisible(NULL))
    }
    
    cur_samps  <- isolate(input$selected_samples)
    def_colors <- c("#000000","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                    "#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB")
    shape_opts <- c(16, 15, 17, 18, 25, 4, 8, 23, 3, 10)
    if (!is.null(cur_samps)) {
      for (i in seq_along(cur_samps)) {
        vn        <- cur_samps[i]; vid <- safe_id(vn)
        def_col   <- def_colors[(i - 1) %% length(def_colors) + 1]
        def_shape <- as.character(shape_opts[(i - 1) %% length(shape_opts) + 1])
        tryCatch({
          updateSelectInput(session,   paste0("color_selector_", vid), selected = "custom")
          updateTextInput(session,     paste0("color_",          vid), value    = def_col)
          updateTextInput(session,     paste0("hex_color_",      vid), value    = def_col)
          updateSelectInput(session,   paste0("shape_",          vid), selected = def_shape)
          updateCheckboxInput(session, paste0("shape_filled_",   vid), value    = TRUE)
          updateSelectInput(session,   paste0("line_type_",      vid), selected = "solid")
          updateTextInput(session,     paste0("legend_label_",   vid), value    = vn)
        }, error = function(e) invisible(NULL))
      }
    }
    
    rv$imported_settings <- NULL
    rv$settings_status   <- list(
      type    = "ok",
      msg     = "All imported settings cleared. Global and per-sample settings restored to defaults.",
      details = NULL
    )
  })
  
  output$settings_status_ui <- renderUI({
    if (is.null(rv$settings_status)) return(NULL)
    cls <- if (rv$settings_status$type == "ok") "settings-status ok" else "settings-status warn"
    tagList(
      div(class = cls, rv$settings_status$msg),
      if (!is.null(rv$settings_status$details))
        HTML(paste0('<div style="font-size:.82em;margin-top:4px;">',
                    rv$settings_status$details, '</div>'))
    )
  })
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  PAGE_SIZE <- 20L
  
  observeEvent(input$selected_samples, { rv$style_page <- 1L }, ignoreNULL = FALSE)
  observeEvent(input$style_prev_page,  { rv$style_page <- max(1L, rv$style_page - 1L) })
  observeEvent(input$style_next_page,  {
    n_pages <- ceiling(length(input$selected_samples) / PAGE_SIZE)
    rv$style_page <- min(n_pages, rv$style_page + 1L)
  })
  
  output$var_settings <- renderUI({
    req(rv$od_vars, input$selected_samples)
    all_samps <- input$selected_samples
    n_total   <- length(all_samps)
    n_pages   <- ceiling(n_total / PAGE_SIZE)
    page      <- rv$style_page
    idx_start <- (page - 1L) * PAGE_SIZE + 1L
    idx_end   <- min(page * PAGE_SIZE, n_total)
    page_samps <- all_samps[idx_start:idx_end]
    
    shape_opts <- c("Circle" = 16, "Square" = 15, "Triangle Up" = 17,
                    "Diamond" = 18, "Triangle Down" = 25, "Cross" = 4,
                    "X" = 8, "Open Circle" = 21, "Open Square" = 22,
                    "Open Triangle" = 24, "Plus" = 3, "Star" = 10)
    default_colors <- c("#000000","#E41A1C","#377EB8","#4DAF4A",
                        "#984EA3","#FF7F00","#A65628","#F781BF",
                        "#999999","#66C2A5","#FC8D62","#8DA0CB")
    named_colors   <- c("#000000","#E41A1C","#377EB8","#4DAF4A",
                        "#984EA3","#FF7F00","#FFFF33","#A65628",
                        "#F781BF","#999999","#009E73")
    
    nav_bar <- if (n_pages > 1L) {
      div(style = "display:flex;align-items:center;justify-content:space-between;margin-bottom:10px;",
          actionButton("style_prev_page", "◄ Prev",
                       style = "font-size:.8em;padding:3px 8px;",
                       disabled = if (page <= 1L) "disabled" else NULL),
          span(style = "font-size:.85em;color:#555;",
               paste0("Page ", page, " of ", n_pages,
                      " (samples ", idx_start, "–", idx_end, " of ", n_total, ")")),
          actionButton("style_next_page", "Next ►",
                       style = "font-size:.8em;padding:3px 8px;",
                       disabled = if (page >= n_pages) "disabled" else NULL)
      )
    } else NULL
    
    sample_widgets <- lapply(seq_along(page_samps), function(j) {
      i   <- idx_start + j - 1L
      vn  <- page_samps[j]
      vid <- safe_id(vn)
      shape_def <- as.character(shape_opts[1 + (i - 1L) %% length(shape_opts)])
      
      imp_aes <- rv$imported_settings$samples[[vn]]
      
      if (!is.null(imp_aes)) {
        color_def <- imp_aes$color %||% default_colors[1]
        sel_def   <- "custom"
      } else if (!is.null(rv$palette_colors) && input$color_palette != "custom") {
        idx       <- match(vn, rv$od_vars)
        color_def <- if (!is.na(idx)) rv$palette_colors[idx] else default_colors[1]
        sel_def   <- "custom"
      } else {
        color_def <- default_colors[1 + (i - 1L) %% length(default_colors)]
        sel_def   <- if (color_def %in% named_colors) color_def else "custom"
      }
      if (is.na(color_def) || is.null(color_def)) color_def <- "#000000"
      
      default_shape  <- if (!is.null(imp_aes$shape))       as.character(imp_aes$shape)    else shape_def
      default_lt     <- if (!is.null(imp_aes$line_type))    imp_aes$line_type              else "solid"
      default_filled <- if (!is.null(imp_aes$shape_filled)) imp_aes$shape_filled           else TRUE
      default_leg    <- if (!is.null(imp_aes$legend_label)) imp_aes$legend_label           else vn
      
      div(
        style = "margin-bottom:12px;padding:10px;border:1px solid #ddd;border-radius:6px;background:#fafafa;",
        h4(vn, style = "margin:0 0 8px 0;font-size:1em;color:#333;"),
        selectInput(paste0("line_type_", vid), "Line Type:",
                    choices  = c("Solid" = "solid", "Dashed" = "dashed",
                                 "Dotted" = "dotted", "DotDash" = "dotdash",
                                 "LongDash" = "longdash", "TwoDash" = "twodash"),
                    selected = default_lt),
        h5("Shape", style = "margin:8px 0 4px;"),
        fluidRow(
          column(6, selectInput(paste0("shape_", vid), "Shape:",
                                choices = shape_opts, selected = default_shape)),
          column(6, checkboxInput(paste0("shape_filled_", vid), "Filled", default_filled))
        ),
        div(class = "panel-section",
            h5("Color", class = "panel-title"),
            selectInput(paste0("color_selector_", vid), "Pick Color:",
                        choices  = c("Black" = "#000000", "Red"    = "#E41A1C",
                                     "Blue"  = "#377EB8", "Green"  = "#4DAF4A",
                                     "Purple"= "#984EA3", "Orange" = "#FF7F00",
                                     "Yellow"= "#FFFF33", "Brown"  = "#A65628",
                                     "Pink"  = "#F781BF", "Gray"   = "#999999",
                                     "Teal"  = "#009E73", "Custom" = "custom"),
                        selected = sel_def),
            conditionalPanel(
              condition = paste0("input['color_selector_", vid, "'] == 'custom'"),
              div(style = "display:flex;align-items:center;gap:8px;",
                  textInput(paste0("color_", vid), "HEX:", value = color_def),
                  div(id    = paste0("color_preview_", vid),
                      class = "color-preview",
                      style = paste0("background-color:", color_def, ";"))
              )
            ),
            checkboxInput(paste0("use_rgb_", vid), "RGB Sliders", FALSE),
            conditionalPanel(
              condition = paste0("input['use_rgb_", vid, "'] == true"),
              sliderInput(paste0("red_",   vid), "R:", 0, 255, 0, 1),
              sliderInput(paste0("green_", vid), "G:", 0, 255, 0, 1),
              sliderInput(paste0("blue_",  vid), "B:", 0, 255, 0, 1),
              sliderInput(paste0("alpha_", vid), "A:", 0,   1, 1, 0.01),
              actionButton(paste0("apply_rgb_", vid), "Apply RGB",
                           style = "margin-top:4px;font-size:.85em;")
            ),
            checkboxInput(paste0("use_hex_", vid), "Direct HEX Entry", FALSE),
            conditionalPanel(
              condition = paste0("input['use_hex_", vid, "'] == true"),
              div(style = "display:flex;align-items:center;gap:8px;",
                  textInput(paste0("hex_color_", vid), "HEX Code:", value = color_def),
                  actionButton(paste0("apply_hex_", vid), "Apply",
                               style = "font-size:.85em;")
              )
            )
        ),
        textInput(paste0("legend_label_", vid), "Legend Label:", value = default_leg)
      )
    })
    
    tagList(nav_bar, do.call(tagList, sample_widgets))
  })
  
  registered_color_obs <- character(0)
  
  observe({
    req(input$selected_samples, rv$data)
    new_samps <- setdiff(safe_id(input$selected_samples), registered_color_obs)
    if (length(new_samps) == 0) return()
    
    for (vn in input$selected_samples[safe_id(input$selected_samples) %in% new_samps]) {
      local({
        vi <- safe_id(vn)
        registered_color_obs <<- c(registered_color_obs, vi)
        
        observeEvent(input[[paste0("color_selector_", vi)]], {
          sel <- input[[paste0("color_selector_", vi)]]
          if (!is.null(sel) && sel != "custom") {
            rc <- tryCatch(col2rgb(sel), error = function(e) NULL)
            if (!is.null(rc)) {
              updateSliderInput(session, paste0("red_",   vi), value = rc[1, 1])
              updateSliderInput(session, paste0("green_", vi), value = rc[2, 1])
              updateSliderInput(session, paste0("blue_",  vi), value = rc[3, 1])
              updateTextInput(session,   paste0("color_",     vi), value = sel)
              updateTextInput(session,   paste0("hex_color_", vi), value = sel)
            }
          }
        }, ignoreInit = TRUE)
        
        observeEvent(input[[paste0("color_", vi)]], {
          col <- normalize_hex_color(input[[paste0("color_", vi)]])
          session$sendCustomMessage("updateColorPreview",
                                    list(id = paste0("color_preview_", vi), color = col))
          rc <- tryCatch(col2rgb(col), error = function(e) NULL)
          if (!is.null(rc)) {
            updateSliderInput(session, paste0("red_",   vi), value = rc[1, 1])
            updateSliderInput(session, paste0("green_", vi), value = rc[2, 1])
            updateSliderInput(session, paste0("blue_",  vi), value = rc[3, 1])
          }
        }, ignoreInit = TRUE)
        
        observeEvent(input[[paste0("apply_rgb_", vi)]], {
          hex <- rgb_to_hex(input[[paste0("red_",   vi)]],
                            input[[paste0("green_", vi)]],
                            input[[paste0("blue_",  vi)]],
                            input[[paste0("alpha_", vi)]])
          updateSelectInput(session, paste0("color_selector_", vi), selected = "custom")
          updateTextInput(session,   paste0("color_",     vi), value = hex)
          updateTextInput(session,   paste0("hex_color_", vi), value = hex)
        })
        
        observeEvent(input[[paste0("apply_hex_", vi)]], {
          hex <- normalize_hex_color(input[[paste0("hex_color_", vi)]])
          updateSelectInput(session, paste0("color_selector_", vi), selected = "custom")
          updateTextInput(session, paste0("color_", vi), value = hex)
          rc <- tryCatch(col2rgb(hex), error = function(e) NULL)
          if (!is.null(rc)) {
            updateSliderInput(session, paste0("red_",   vi), value = rc[1, 1])
            updateSliderInput(session, paste0("green_", vi), value = rc[2, 1])
            updateSliderInput(session, paste0("blue_",  vi), value = rc[3, 1])
          }
        })
      })
    }
  })
  
  # ── Dynamic UI: region / marker settings ────────────────────────────────────
  output$region_settings <- renderUI({
    req(input$region_count, input$enable_highlighting)
    do.call(tagList, lapply(seq_len(input$region_count), function(i) {
      div(style = "border:1px solid #ddd;padding:10px;margin-bottom:8px;border-radius:5px;",
          h5(paste("Region", i)),
          fluidRow(
            column(6, numericInput(paste0("region_x_min_", i), "X min:", 0)),
            column(6, numericInput(paste0("region_x_max_", i), "X max:", 60))
          ),
          fluidRow(
            column(6, numericInput(paste0("region_y_min_", i), "Y min:", 0.01)),
            column(6, numericInput(paste0("region_y_max_", i), "Y max:", 1))
          ),
          selectInput(paste0("region_color_", i), "Fill Color:",
                      choices  = c("Light Gray"   = "#ebebeb", "Light Blue"  = "#e6f3ff",
                                   "Light Red"    = "#ffebeb", "Light Green" = "#ebffeb",
                                   "Light Yellow" = "#ffffeb", "Custom"      = "custom"),
                      selected = "#ebebeb"),
          conditionalPanel(
            condition = paste0("input['region_color_", i, "'] == 'custom'"),
            textInput(paste0("region_color_custom_", i), "HEX:", "#ebebeb")
          ),
          sliderInput(paste0("region_alpha_", i), "Transparency:", 0, 1, 0.3, 0.05)
      )
    }))
  })
  
  output$time_marker_settings <- renderUI({
    req(input$marker_count, input$enable_time_markers)
    do.call(tagList, lapply(seq_len(input$marker_count), function(i) {
      div(style = "border:1px solid #ddd;padding:10px;margin-bottom:8px;border-radius:5px;",
          h5(paste("Marker", i)),
          numericInput(paste0("marker_time_", i), "Time Point:", 30),
          selectInput(paste0("marker_line_type_", i), "Line Type:",
                      choices  = c("Dashed" = "dashed", "Solid" = "solid",
                                   "Dotted" = "dotted", "DotDash" = "dotdash",
                                   "LongDash" = "longdash"),
                      selected = "dashed"),
          selectInput(paste0("marker_color_", i), "Color:",
                      choices  = c("Black" = "#000000", "Red"   = "#E41A1C",
                                   "Blue"  = "#377EB8", "Green" = "#4DAF4A",
                                   "Gray"  = "#999999", "Custom" = "custom"),
                      selected = "#000000"),
          conditionalPanel(
            condition = paste0("input['marker_color_", i, "'] == 'custom'"),
            textInput(paste0("marker_color_custom_", i), "HEX:", "#000000")
          ),
          sliderInput(paste0("marker_size_", i), "Width:", 0.1, 3, 1, 0.1),
          checkboxInput(paste0("marker_label_", i), "Add Label", FALSE),
          conditionalPanel(
            condition = paste0("input['marker_label_", i, "'] == 'true'"),
            textInput(paste0("marker_text_", i),           "Label:",    paste0("t=", 30)),
            numericInput(paste0("marker_label_size_", i),  "Size:",     4, 2, 12),
            selectInput(paste0("marker_label_position_", i), "Position:",
                        choices = c("Top" = "top", "Middle" = "middle", "Bottom" = "bottom"),
                        selected = "top"),
            numericInput(paste0("marker_label_hjust_", i), "H-Offset:", 0, -5, 5, 0.1)
          )
      )
    }))
  })
  
  # ── Time filtering ───────────────────────────────────────────────────────────
  apply_time_filters <- function(d, time_col) {
    d[[time_col]] <- suppressWarnings(as.numeric(as.character(d[[time_col]])))
    if (!isTRUE(input$enable_time_filter)) return(d)
    rng <- input$time_filter_range
    if (!is.null(rng) && length(rng) == 2) {
      t_vec <- d[[time_col]]
      d <- d[!is.na(t_vec) & t_vec >= rng[1] & t_vec <= rng[2], , drop = FALSE]
    }
    excl <- parse_excluded_timepoints(input$exclude_timepoints, rv$all_timepoints)
    if (length(excl) > 0) {
      t_vec <- d[[time_col]]
      tol   <- max(abs(rv$all_timepoints), na.rm = TRUE) * 0.001 + 0.001
      hit   <- rowSums(abs(outer(t_vec, excl, "-")) <= tol) == 0L
      d     <- d[hit, , drop = FALSE]
    }
    d
  }
  
  # ── Data prep ────────────────────────────────────────────────────────────────
  prepare_plot_data <- function(samples = NULL) {
    req(rv$data, rv$time_col, input$selected_samples)
    samps <- if (!is.null(samples)) samples else input$selected_samples
    if (rv$is_long_format) {
      req(rv$group_col, rv$value_col)
      d  <- apply_time_filters(rv$data, rv$time_col)
      d  <- d[d[[rv$group_col]] %in% samps, , drop = FALSE]
      pd <- d %>%
        group_by(across(all_of(c(rv$time_col, rv$group_col)))) %>%
        summarise(mean_value = mean(.data[[rv$value_col]], na.rm = TRUE),
                  sd_value   = sd(  .data[[rv$value_col]], na.rm = TRUE),
                  n          = n(), .groups = "drop") %>%
        mutate(sd_value   = ifelse(is.na(sd_value), 0, sd_value),
               sem_value  = sd_value / sqrt(n),
               ci95_value = qt(0.975, df = pmax(n - 1, 1)) * sem_value)
      names(pd)[names(pd) == rv$group_col] <- "variable"
      names(pd)[names(pd) == rv$time_col]  <- "time"
    } else {
      d    <- apply_time_filters(rv$data, rv$time_col)
      meas <- intersect(rv$od_vars_raw, samps)
      d_sub <- d[, c(rv$time_col, meas), drop = FALSE]
      pd   <- d_sub %>%
        pivot_longer(cols = all_of(meas), names_to = "variable", values_to = "value") %>%
        group_by(time = .data[[rv$time_col]], variable) %>%
        summarise(mean_value = mean(value, na.rm = TRUE),
                  sd_value   = sd(  value, na.rm = TRUE),
                  n          = n(), .groups = "drop") %>%
        mutate(sd_value   = ifelse(is.na(sd_value), 0, sd_value),
               sem_value  = sd_value / sqrt(n),
               ci95_value = qt(0.975, df = pmax(n - 1, 1)) * sem_value)
    }
    pd %>% filter(is.finite(time) & is.finite(mean_value))
  }
  
  prepare_metrics_data <- function(samples = NULL) {
    req(rv$data, rv$time_col, input$selected_samples)
    samps <- if (!is.null(samples)) samples else input$selected_samples
    if (rv$is_long_format) {
      req(rv$group_col, rv$value_col)
      d  <- apply_time_filters(rv$data, rv$time_col)
      d  <- d[d[[rv$group_col]] %in% samps, , drop = FALSE]
      grp_cols <- c(rv$time_col, rv$group_col)
      if (!is.null(rv$rep_col) && rv$rep_col %in% names(d)) grp_cols <- c(grp_cols, rv$rep_col)
      pd <- d %>%
        group_by(across(all_of(grp_cols))) %>%
        summarise(mean_value = mean(.data[[rv$value_col]], na.rm = TRUE),
                  .groups = "drop")
      names(pd)[names(pd) == rv$group_col] <- "variable"
      names(pd)[names(pd) == rv$time_col]  <- "time"
      if (!is.null(rv$rep_col) && rv$rep_col %in% names(pd))
        names(pd)[names(pd) == rv$rep_col] <- "replicate"
    } else {
      d  <- apply_time_filters(rv$data, rv$time_col)
      meas <- intersect(rv$od_vars_raw, samps)
      if (is.null(rv$wide_block_id)) {
        wide_rep <- infer_wide_replicates(d, rv$time_col)
        rv$wide_block_id  <- wide_rep$block_id
        rv$wide_rep_count <- wide_rep$reps
      }
      if (!is.null(rv$wide_block_id)) d$replicate <- rv$wide_block_id
      grp_cols_w <- c(rv$time_col, "variable")
      if ("replicate" %in% names(d)) grp_cols_w <- c(grp_cols_w, "replicate")
      pd <- d %>%
        pivot_longer(cols = all_of(meas), names_to = "variable", values_to = "value") %>%
        group_by(across(all_of(grp_cols_w))) %>%
        summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
        rename(time = !!sym(rv$time_col))
    }
    pd %>% filter(is.finite(time) & is.finite(mean_value))
  }

  # ── Replicate-level data (for spaghetti / jitter / quantile-band displays) ──
  prepare_replicate_data <- function(samples = NULL) {
    req(rv$data, rv$time_col)
    samps <- if (!is.null(samples)) samples else input$selected_samples
    if (rv$is_long_format) {
      req(rv$group_col, rv$value_col)
      d <- apply_time_filters(rv$data, rv$time_col)
      d <- d[d[[rv$group_col]] %in% samps, , drop = FALSE]
      grp_cols <- c(rv$time_col, rv$group_col)
      if (!is.null(rv$rep_col) && rv$rep_col %in% names(d))
        grp_cols <- c(grp_cols, rv$rep_col)
      pd <- d %>%
        group_by(across(all_of(grp_cols))) %>%
        summarise(rep_value = mean(.data[[rv$value_col]], na.rm = TRUE), .groups = "drop")
      names(pd)[names(pd) == rv$group_col] <- "variable"
      names(pd)[names(pd) == rv$time_col]  <- "time"
      if (!is.null(rv$rep_col) && rv$rep_col %in% names(pd))
        names(pd)[names(pd) == rv$rep_col] <- "replicate"
      else
        pd$replicate <- 1L
    } else {
      d    <- apply_time_filters(rv$data, rv$time_col)
      meas <- intersect(rv$od_vars_raw, samps)
      if (is.null(rv$wide_block_id)) {
        wide_rep <- infer_wide_replicates(d, rv$time_col)
        rv$wide_block_id  <- wide_rep$block_id
        rv$wide_rep_count <- wide_rep$reps
      }
      d$replicate <- if (!is.null(rv$wide_block_id)) rv$wide_block_id else 1L
      pd <- d %>%
        pivot_longer(cols = all_of(meas), names_to = "variable", values_to = "rep_value") %>%
        rename(time = !!sym(rv$time_col)) %>%
        select(time, variable, replicate, rep_value)
    }
    pd %>% filter(is.finite(time) & is.finite(rep_value))
  }

  # ── Aesthetics resolver ──────────────────────────────────────────────────────
  resolve_aesthetics <- function(samples) {
    n          <- length(samples)
    shapes     <- setNames(numeric(n),   samples)
    colors     <- setNames(character(n), samples)
    line_types <- setNames(character(n), samples)
    leg_labels <- setNames(character(n), samples)
    filled_map <- setNames(logical(n),   samples)
    def_colors <- c("#000000","#E41A1C","#377EB8","#4DAF4A",
                    "#984EA3","#FF7F00","#A65628","#F781BF",
                    "#999999","#66C2A5","#FC8D62","#8DA0CB")
    def_shapes <- c(16, 15, 17, 18, 25, 4, 8, 23, 3, 10)
    
    for (i in seq_len(n)) {
      vn  <- samples[i]; vid <- safe_id(vn)
      si  <- input[[paste0("shape_", vid)]]
      shapes[i] <- if (!is.null(si)) as.numeric(si)
      else def_shapes[(i - 1) %% length(def_shapes) + 1]
      
      sel <- input[[paste0("color_selector_", vid)]]
      if (!is.null(sel) && sel == "custom")
        colors[i] <- normalize_hex_color(input[[paste0("color_", vid)]])
      else if (!is.null(sel))
        colors[i] <- sel
      else if (!is.null(rv$palette_colors)) {
        idx       <- match(vn, rv$od_vars)
        colors[i] <- if (!is.na(idx)) rv$palette_colors[idx]
        else def_colors[(i - 1) %% length(def_colors) + 1]
      } else
        colors[i] <- def_colors[(i - 1) %% length(def_colors) + 1]
      
      lt <- input[[paste0("line_type_", vid)]]
      line_types[i] <- if (!is.null(lt)) lt else "solid"
      
      ll      <- input[[paste0("legend_label_", vid)]]
      base_ll <- if (!is.null(ll) && nchar(ll) > 0) ll else vn
      wrap_w  <- if (!is.null(input$legend_wrap_width)) input$legend_wrap_width else 20L
      leg_labels[i] <- stringr::str_wrap(base_ll, width = wrap_w)
      
      fi <- input[[paste0("shape_filled_", vid)]]
      filled_map[i] <- if (!is.null(fi)) fi else TRUE
    }
    list(shapes = shapes, colors = colors, line_types = line_types,
         leg_labels = leg_labels, filled_map = filled_map)
  }
  
  # ── Core plot builder ────────────────────────────────────────────────────────
  build_plot <- function(plot_data, samples, aes_vals, highlight_samples = NULL) {
    shapes     <- aes_vals$shapes;     colors     <- aes_vals$colors
    line_types <- aes_vals$line_types; leg_labels <- aes_vals$leg_labels
    filled_map <- aes_vals$filled_map; n_vars     <- length(samples)
    
    display_colors <- colors
    if (!is.null(highlight_samples))
      for (i in seq_len(n_vars))
        if (!samples[i] %in% highlight_samples) display_colors[i] <- "#DDDDDD"
    
    plot_data$variable <- factor(plot_data$variable, levels = samples)
    p <- ggplot(plot_data,
                aes(x = time, y = mean_value, group = variable,
                    color = variable, shape = variable, linetype = variable))
    
    # Region highlighting
    if (isTRUE(input$enable_highlighting) && !is.null(input$region_count)) {
      for (i in seq_len(input$region_count)) {
        x1 <- input[[paste0("region_x_min_", i)]]; x2 <- input[[paste0("region_x_max_", i)]]
        y1 <- input[[paste0("region_y_min_", i)]]; y2 <- input[[paste0("region_y_max_", i)]]
        rc  <- input[[paste0("region_color_", i)]]
        col <- if (!is.null(rc) && rc == "custom") input[[paste0("region_color_custom_", i)]] else rc
        alp <- input[[paste0("region_alpha_", i)]]
        if (!is.null(x1) && !is.null(col) && !is.null(alp))
          p <- p + annotate("rect", xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=col, alpha=alp)
      }
    }
    
    # Time markers
    if (isTRUE(input$enable_time_markers) && !is.null(input$marker_count)) {
      yr <- range(plot_data$mean_value, na.rm = TRUE)
      if (input$y_scale_type == "log") yr[1] <- max(1e-9, yr[1])
      for (i in seq_len(input$marker_count)) {
        tp <- input[[paste0("marker_time_", i)]]
        lt <- input[[paste0("marker_line_type_", i)]]
        mc <- input[[paste0("marker_color_", i)]]
        lc <- if (!is.null(mc) && mc == "custom") input[[paste0("marker_color_custom_", i)]] else mc
        ls <- input[[paste0("marker_size_", i)]]
        if (!is.null(tp) && !is.null(lc))
          p <- p + geom_vline(xintercept = tp, linetype = lt, color = lc, linewidth = ls)
        if (isTRUE(input[[paste0("marker_label_", i)]])) {
          lbl  <- input[[paste0("marker_text_",         i)]]
          lsz  <- input[[paste0("marker_label_size_",   i)]]
          lpos <- input[[paste0("marker_label_position_",i)]]
          lhj  <- input[[paste0("marker_label_hjust_",  i)]]; if (is.null(lhj)) lhj <- 0
          yp <- switch(lpos,
                       top    = if (input$y_scale_type == "log") 10^(log10(yr[2]) - diff(log10(yr)) * 0.08) else yr[1] + diff(yr) * 0.92,
                       bottom = if (input$y_scale_type == "log") 10^(log10(yr[1]) + diff(log10(yr)) * 0.08) else yr[1] + diff(yr) * 0.08,
                       middle = if (input$y_scale_type == "log") 10^(mean(log10(yr))) else mean(yr))
          p <- p + annotate("text", x = tp + lhj, y = yp, label = lbl, size = lsz, color = lc)
        }
      }
    }
    
    # ── Variability / error display ──────────────────────────────────────────────
    using_shadow <- FALSE
    dm <- if (!is.null(input$error_display_mode)) input$error_display_mode else "bars"

    # Helper: compute stat-based error bounds into plot_data
    compute_stat_bounds <- function() {
      et <- if (!is.null(input$error_type)) input$error_type else "sem"
      ec <- switch(et, sd = "sd_value", sem = "sem_value", ci95 = "ci95_value", "sem_value")
      em <- if (!is.null(input$error_multiplier)) input$error_multiplier else 1
      ev <- plot_data[[ec]]; ev[is.na(ev)] <- 0
      lo <- if (!is.null(input$y_scale_type) && input$y_scale_type == "log")
        pmax(plot_data$mean_value - ev * em, 1e-9)
      else if (isTRUE(input$asymmetric_error))
        pmax(plot_data$mean_value - ev * em, 0)
      else
        plot_data$mean_value - ev * em
      plot_data$err_lo <<- lo
      plot_data$err_hi <<- plot_data$mean_value + ev * em
    }

    if (dm == "shadow") {
      compute_stat_bounds()
      using_shadow <- TRUE
      p <- p + geom_ribbon(data = plot_data,
                           aes(x = time, ymin = err_lo, ymax = err_hi,
                               group = variable, fill = variable),
                           alpha = if (!is.null(input$shadow_alpha)) input$shadow_alpha else 0.2,
                           color = NA, inherit.aes = FALSE)

    } else if (dm == "bars") {
      compute_stat_bounds()
      pos <- if (!is.null(input$error_bar_position) && input$error_bar_position == "dodge")
        position_dodge(0.2) else position_identity()
      blt <- if (!is.null(input$error_bar_style) && input$error_bar_style == "dashed") "dashed" else "solid"
      if (!is.null(input$error_bar_style) && input$error_bar_style == "T") {
        p <- p + geom_errorbar(data = plot_data,
                               aes(x = time, ymin = err_lo, ymax = err_hi, color = variable),
                               width     = if (!is.null(input$error_bar_width))     input$error_bar_width     else 1,
                               linewidth = if (!is.null(input$error_bar_thickness)) input$error_bar_thickness else 0.8,
                               linetype  = blt, position = pos, inherit.aes = FALSE)
      } else {
        p <- p + geom_linerange(data = plot_data,
                                aes(x = time, ymin = err_lo, ymax = err_hi, color = variable),
                                linewidth = if (!is.null(input$error_bar_thickness)) input$error_bar_thickness else 0.8,
                                linetype  = blt, position = pos, inherit.aes = FALSE)
      }

    } else if (dm == "spaghetti") {
      # Individual replicate traces rendered behind the mean line
      rep_d <- tryCatch(prepare_replicate_data(samples), error = function(e) NULL)
      if (!is.null(rep_d) && nrow(rep_d) > 0 && "replicate" %in% names(rep_d)) {
        s_alpha <- if (!is.null(input$spaghetti_alpha)) input$spaghetti_alpha else 0.35
        p <- p + geom_line(data = rep_d,
                           aes(x = time, y = rep_value,
                               group = interaction(variable, replicate),
                               color = variable),
                           alpha = s_alpha, linewidth = input$line_thickness * 0.55,
                           inherit.aes = FALSE)
      }

    } else if (dm == "quantile_bands") {
      # Nested empirical quantile bands: IQR inner, 2.5–97.5% outer
      rep_d <- tryCatch(prepare_replicate_data(samples), error = function(e) NULL)
      if (!is.null(rep_d) && nrow(rep_d) > 0) {
        q_inner <- if (!is.null(input$qband_inner_alpha)) input$qband_inner_alpha else 0.40
        q_outer <- if (!is.null(input$qband_outer_alpha)) input$qband_outer_alpha else 0.15
        qd <- rep_d %>%
          group_by(time, variable) %>%
          summarise(q25  = quantile(rep_value, 0.25,  na.rm = TRUE),
                    q75  = quantile(rep_value, 0.75,  na.rm = TRUE),
                    q025 = quantile(rep_value, 0.025, na.rm = TRUE),
                    q975 = quantile(rep_value, 0.975, na.rm = TRUE),
                    .groups = "drop")
        using_shadow <- TRUE
        p <- p +
          geom_ribbon(data = qd,
                      aes(x = time, ymin = q025, ymax = q975,
                          group = variable, fill = variable),
                      alpha = q_outer, color = NA, inherit.aes = FALSE) +
          geom_ribbon(data = qd,
                      aes(x = time, ymin = q25, ymax = q75,
                          group = variable, fill = variable),
                      alpha = q_inner, color = NA, inherit.aes = FALSE)
      }

    } else if (dm == "jitter") {
      # Raw replicate observations as jittered points
      rep_d <- tryCatch(prepare_replicate_data(samples), error = function(e) NULL)
      if (!is.null(rep_d) && nrow(rep_d) > 0) {
        j_alpha <- if (!is.null(input$jitter_alpha)) input$jitter_alpha else 0.55
        j_size  <- if (!is.null(input$jitter_size))  input$jitter_size  else 1.5
        j_w     <- if (!is.null(input$jitter_width)) input$jitter_width else 0.0
        p <- p + geom_jitter(data = rep_d,
                             aes(x = time, y = rep_value, color = variable),
                             width = j_w, height = 0,
                             alpha = j_alpha, size = j_size, inherit.aes = FALSE)
      }

    } else if (dm == "combo") {
      # Semi-transparent replicate traces + light CI ribbon behind the mean line
      rep_d <- tryCatch(prepare_replicate_data(samples), error = function(e) NULL)
      if (!is.null(rep_d) && nrow(rep_d) > 0 && "replicate" %in% names(rep_d)) {
        s_alpha <- if (!is.null(input$spaghetti_alpha)) input$spaghetti_alpha else 0.25
        p <- p + geom_line(data = rep_d,
                           aes(x = time, y = rep_value,
                               group = interaction(variable, replicate),
                               color = variable),
                           alpha = s_alpha, linewidth = input$line_thickness * 0.45,
                           inherit.aes = FALSE)
      }
      compute_stat_bounds()
      using_shadow <- TRUE
      p <- p + geom_ribbon(data = plot_data,
                           aes(x = time, ymin = err_lo, ymax = err_hi,
                               group = variable, fill = variable),
                           alpha = if (!is.null(input$shadow_alpha)) input$shadow_alpha * 0.6 else 0.12,
                           color = NA, inherit.aes = FALSE)
    }
    
    p <- p + geom_line(aes(group = variable), linewidth = input$line_thickness)
    
    if (input$show_points) {
      shape_map <- integer(n_vars)
      fill_map  <- character(n_vars)
      for (i in seq_len(n_vars)) {
        pt_s <- as.integer(shapes[i])
        if (!filled_map[i]) {
          pt_s         <- switch(as.character(pt_s),
                                 "16" = 21L, "15" = 22L, "17" = 24L, "18" = 23L, pt_s)
          fill_map[i]  <- "white"
        } else {
          fill_map[i]  <- display_colors[i]
        }
        shape_map[i] <- pt_s
      }
      names(shape_map) <- samples
      names(fill_map)  <- samples
      
      if (using_shadow) {
        p <- p + geom_point(
          data        = plot_data,
          aes(x = time, y = mean_value, color = variable, shape = variable,
              fill = variable),
          size        = input$shape_size,
          stroke      = input$point_stroke,
          inherit.aes = FALSE
        ) +
          scale_shape_manual(values = shape_map, labels = leg_labels,
                             name = NULL, guide = "none") +
          scale_fill_manual(values = setNames(display_colors, samples),
                            labels = leg_labels, guide = "none")
      } else {
        pt_data          <- plot_data
        pt_data$pt_fill  <- fill_map[pt_data$variable]
        
        p <- p + geom_point(
          data        = pt_data,
          aes(x = time, y = mean_value, color = variable, shape = variable,
              fill = pt_fill),
          size        = input$shape_size,
          stroke      = input$point_stroke,
          inherit.aes = FALSE
        ) +
          scale_shape_manual(values = shape_map, labels = leg_labels,
                             name = NULL, guide = "none") +
          scale_fill_identity()
      }
    } else if (using_shadow) {
      p <- p + scale_fill_manual(values = setNames(display_colors, samples),
                                 labels = leg_labels, guide = "none")
    }
    
    if (input$show_end_labels) {
      mt  <- max(plot_data$time, na.rm = TRUE)
      off <- mt * (input$label_offset / 100)
      ep  <- plot_data %>% group_by(variable) %>% filter(time == max(time)) %>% ungroup()
      p   <- p + geom_text_repel(data = ep,
                                 aes(label = variable, color = variable, x = Inf, y = mean_value),
                                 direction = "y", xlim = c(mt + off, Inf),
                                 min.segment.length = Inf, hjust = 0,
                                 size = input$label_font_size / 2.835,
                                 fontface = if (input$label_bold) "bold" else "plain")
    }
    
    p <- p +
      scale_color_manual(   values = setNames(display_colors, samples), labels = leg_labels, name = NULL) +
      scale_linetype_manual(values = line_types,                         labels = leg_labels, name = NULL)
    
    if (!isTRUE(input$show_points))
      p <- p + scale_shape_manual(values = shapes, labels = leg_labels, name = NULL)
    
    p <- p +
      guides(shape = guide_legend(override.aes = list(alpha = 1)),
             color = guide_legend(override.aes = list(alpha = 1)))
    
    x_exp  <- expansion(mult = c(input$x_expand_left,   input$x_expand_right))
    y_exp  <- expansion(mult = c(input$y_expand_bottom, input$y_expand_top))
    x_lims <- if (isTRUE(input$custom_x_limits)) c(input$x_min, input$x_max) else NULL
    y_lims <- if (isTRUE(input$custom_y_limits)) c(input$y_min, input$y_max) else NULL
    
    if (input$x_scale_type == "log") {
      p <- p + scale_x_log10(limits = x_lims, expand = x_exp)
    } else if (input$x_scale_type == "sqrt") {
      p <- p + scale_x_sqrt(limits = x_lims, expand = x_exp)
    } else if (input$x_scale_type == "reverse") {
      p <- p + scale_x_reverse(limits = x_lims, expand = x_exp)
    } else if (isTRUE(input$use_advanced_ticks)) {
      mt2 <- if (isTRUE(input$custom_x_limits)) input$x_max
      else max(plot_data$time, na.rm = TRUE)
      x0  <- if (isTRUE(input$custom_x_limits)) input$x_min else 0
      
      extra_ticks <- tryCatch({
        raw <- trimws(input$x_extra_ticks)
        if (nchar(raw) == 0) numeric(0)
        else {
          vals <- suppressWarnings(as.numeric(strsplit(raw, "[,;[:space:]]+")[[1]]))
          vals[is.finite(vals)]
        }
      }, error = function(e) numeric(0))
      
      manual_interval <- input$x_tick_interval
      use_manual <- !is.null(manual_interval) &&
        !is.na(manual_interval) &&
        is.numeric(manual_interval) &&
        manual_interval > 0
      
      if (use_manual) {
        iv  <- manual_interval
        mjb <- seq(x0, ceiling((mt2 - x0) / iv) * iv + x0, by = iv)
        if (!isTRUE(input$custom_x_limits) && !mt2 %in% mjb)
          mjb <- sort(unique(c(mjb, mt2)))
        mnb <- seq(x0, max(mjb), by = iv / 2)
      } else {
        if (mt2 <= 60)       { iv <- 10;  mjb <- seq(0, ceiling(mt2/10)*10,  10);  mnb <- seq(0, ceiling(mt2/10)*10,  5) }
        else if (mt2 <= 120) { iv <- 30;  mjb <- seq(0, ceiling(mt2/30)*30,  30);  mnb <- seq(0, ceiling(mt2/30)*30, 10) }
        else                 { iv <- 60;  mjb <- seq(0, ceiling(mt2/60)*60,  60);  mnb <- seq(0, ceiling(mt2/60)*60, 30) }
        if (!isTRUE(input$custom_x_limits) && !mt2 %in% mjb)
          mjb <- sort(unique(c(mjb, mt2)))
      }
      
      if (length(extra_ticks) > 0)
        mjb <- sort(unique(c(mjb, extra_ticks)))
      
      p <- p + scale_x_continuous(limits = x_lims, expand = x_exp, breaks = mjb,
                                  minor_breaks = mnb, guide = guide_prism_minor())
    } else {
      p <- p + scale_x_continuous(limits = x_lims, expand = x_exp)
    }
    
    if (input$y_scale_type == "log") {
      if (isTRUE(input$use_advanced_ticks)) {
        minE <- input$y_log_min_exponent; maxE <- input$y_log_max_exponent
        if (is.null(y_lims)) y_lims <- c(10^minE, if (maxE >= 0) 10^(maxE+1) else 10^maxE)
        ymn <- c(rep(1:9, maxE-minE+1) * 10^rep(minE:maxE, each=9))
        if (maxE >= 0) ymn <- c(ymn, 10^(maxE+1))
        p <- p + scale_y_log10(limits = y_lims, expand = y_exp, minor_breaks = ymn,
                               guide = guide_prism_minor())
      } else {
        p <- p + scale_y_log10(limits = y_lims, expand = y_exp)
      }
    } else if (input$y_scale_type == "sqrt") {
      p <- p + scale_y_sqrt(limits = y_lims, expand = y_exp)
    } else if (input$y_scale_type == "reverse") {
      p <- p + scale_y_reverse(limits = y_lims, expand = y_exp)
    } else {
      p <- p + scale_y_continuous(limits = y_lims, expand = y_exp)
    }
    
    tf      <- if (isTRUE(input$bold_title))        "bold"   else "plain"
    af      <- if (isTRUE(input$italic_axis_labels)) "italic" else "plain"
    base_t  <- if (isTRUE(input$use_advanced_ticks)) theme_prism(border=TRUE) else theme_pubr()
    maj_col <- if (!is.null(input$major_gridline_color) && input$major_gridline_color == "custom")
      input$major_gridline_color_custom else input$major_gridline_color
    min_col <- if (!is.null(input$minor_gridline_color) && input$minor_gridline_color == "custom")
      input$minor_gridline_color_custom else input$minor_gridline_color
    x_ang   <- as.numeric(input$axis_text_angle)
    
    p <- p + base_t + theme(
      text             = element_text(family = input$font_family),
      plot.title       = element_text(size = input$title_font_size,       face = tf, hjust = 0.5),
      axis.title       = element_text(size = input$axis_label_font_size,  face = af),
      axis.text        = element_text(size = input$axis_text_font_size),
      axis.text.x      = element_text(angle = x_ang,
                                      hjust = if (x_ang > 0) 1 else 0.5,
                                      vjust = if (x_ang > 0) 1 else 0.5),
      panel.grid.major = if (isTRUE(input$show_major_gridlines))
        element_line(color = maj_col, linewidth = input$major_gridline_size) else element_blank(),
      panel.grid.minor = if (isTRUE(input$show_minor_gridlines))
        element_line(color = min_col, linewidth = input$minor_gridline_size) else element_blank(),
      axis.ticks        = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position   = {
        lp <- if (!is.null(input$legend_position)) input$legend_position else "right"
        if (isTRUE(input$show_end_labels)) "none"
        else if (lp == "inside") c(
          if (!is.null(input$legend_x)) input$legend_x else 0.85,
          if (!is.null(input$legend_y)) input$legend_y else 0.95
        )
        else lp
      },
      legend.justification = if (!is.null(input$legend_position) && input$legend_position == "inside" && !isTRUE(input$show_end_labels)) c("right", "top") else "center",
      legend.background = if (isTRUE(input$legend_no_box)) element_blank()
                          else element_rect(fill = "white", color = "gray80"),
      legend.box.background = if (isTRUE(input$legend_no_box)) element_blank()
                              else element_rect(fill = NA),
      legend.key        = element_rect(fill = NA),
      legend.text       = element_text(size  = if (!is.null(input$axis_text_font_size))  input$axis_text_font_size  else 12,
                                       face  = if (!is.null(input$legend_text_face))  input$legend_text_face  else "plain"),
      legend.title      = element_text(size  = if (!is.null(input$axis_label_font_size)) input$axis_label_font_size else 14,
                                       face  = if (!is.null(input$legend_text_face))  input$legend_text_face  else "plain"),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

    # Reserve a fixed right margin when the legend isn't consuming space externally,
    # so the panel stays the same size whether legend is visible or not.
    eff_leg <- {
      lp <- if (!is.null(input$legend_position)) input$legend_position else "right"
      if (isTRUE(input$show_end_labels)) "none" else lp
    }
    if (eff_leg %in% c("none", "inside")) {
      rsv <- if (!is.null(input$legend_reserve_space)) input$legend_reserve_space else 130
      p <- p + theme(plot.margin = margin(5, rsv, 5, 5, "pt"))
    }

    if (isTRUE(input$custom_aspect_ratio))
      p <- p + theme(aspect.ratio = 1 / input$aspect_ratio)
    
    # Notes caption
    notes_cap <- if (isTRUE(input$notes_show_caption)) {
      tryCatch(build_notes_caption(), error = function(e) "")
    } else ""
    if (nchar(notes_cap) > 0) {
      cap_sz <- if (!is.null(input$notes_caption_size)) input$notes_caption_size else 9
      p <- p + theme(plot.caption = element_text(size = cap_sz, hjust = 0,
                                                  color = "#555555", face = "plain"),
                     plot.caption.position = "plot")
    }

    p + coord_cartesian(clip = "off") +
      labs(x        = input$x_axis_label,
           y        = input$y_axis_label,
           title    = input$plot_title,
           subtitle = if (nchar(trimws(input$plot_subtitle)) > 0) input$plot_subtitle else NULL,
           caption  = if (nchar(notes_cap) > 0) notes_cap else NULL)
  }
  
  generate_plot <- function() {
    req(rv$data, input$selected_samples)
    tryCatch(
      build_plot(prepare_plot_data(), input$selected_samples,
                 resolve_aesthetics(input$selected_samples)),
      error = function(e) {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste("Plot error:", conditionMessage(e)),
                   color = "red", size = 4, hjust = 0.5) +
          theme_void()
      }
    )
  }
  
  plot_inputs_debounced <- debounce(reactive({
    list(
      selected_samples    = input$selected_samples,
      error_type          = input$error_type,
      error_display_mode  = input$error_display_mode,
      error_multiplier    = input$error_multiplier,
      line_thickness      = input$line_thickness,
      show_points         = input$show_points,
      shape_size          = input$shape_size,
      x_scale_type        = input$x_scale_type,
      y_scale_type        = input$y_scale_type,
      custom_x_limits     = input$custom_x_limits,
      custom_y_limits     = input$custom_y_limits,
      x_min = input$x_min, x_max = input$x_max,
      y_min = input$y_min, y_max = input$y_max,
      enable_time_filter  = input$enable_time_filter,
      time_filter_range   = input$time_filter_range,
      exclude_timepoints  = input$exclude_timepoints
    )
  }), millis = 400)
  
  # ── Fixed-panel helpers ───────────────────────────────────────────────────────
  # Forces the ggplot panel to constant pt dimensions so the curve-drawing area
  # never changes regardless of legend width or display choices.
  # px → pt conversion at 96 DPI: 1 px = 0.75 pt.
  # Panel occupies 68% of figure width and 78% of figure height (axes/margins
  # take the rest on a typical pubr/prism theme).
  fix_panel_size <- function(p, base_w_px, base_h_px) {
    panel_w_pt <- max(base_w_px * 0.75 * 0.68, 60)
    panel_h_pt <- max(base_h_px * 0.75 * 0.78, 50)
    g   <- ggplotGrob(p)
    pos <- g$layout[g$layout$name == "panel", , drop = FALSE]
    g$widths[pos$l]  <- unit(panel_w_pt, "pt")
    g$heights[pos$t] <- unit(panel_h_pt, "pt")
    g
  }

  # Estimates how many extra pixels the legend will add to the right of the
  # panel when legend.position == "right".  Accounts for the wrapped labels
  # already applied by resolve_aesthetics().
  legend_extra_px <- reactive({
    lp       <- if (!is.null(input$legend_position)) input$legend_position else "right"
    show_end <- isTRUE(input$show_end_labels)
    if (show_end || lp != "right") return(0L)
    req(input$selected_samples)
    ll <- tryCatch(resolve_aesthetics(input$selected_samples)$leg_labels,
                   error = function(e) character(0))
    if (length(ll) == 0) return(0L)
    # Labels may already contain newlines from str_wrap — measure longest line
    max_chars <- max(nchar(unlist(strsplit(ll, "\n"))), na.rm = TRUE)
    font_sz   <- if (!is.null(input$axis_text_font_size)) input$axis_text_font_size else 12
    # 1 pt = 1.333 px at 96 DPI.  Key ≈ 3 em wide; text ≈ 0.65 em/char; 40 px padding.
    em_px <- font_sz * 1.333
    as.integer(ceiling(3 * em_px + max_chars * 0.65 * em_px + 40))
  })

  output$plot_container <- renderUI({
    # Container expands rightward to fit the legend; panel size stays constant.
    extra_px <- legend_extra_px()
    plotOutput("od_plot",
               height = paste0(input$plot_height, "px"),
               width  = paste0(input$plot_width + extra_px, "px"),
               click  = "plot_click")
  })

  # Click-to-place legend
  observeEvent(input$plot_click, {
    if (isTRUE(input$legend_click_mode) &&
        !is.null(input$legend_position) && input$legend_position == "inside") {
      cx <- input$plot_click$coords_npc$x
      cy <- input$plot_click$coords_npc$y
      if (!is.null(cx) && !is.null(cy)) {
        updateNumericInput(session, "legend_x", value = round(cx, 2))
        updateNumericInput(session, "legend_y", value = round(cy, 2))
      }
    }
  })

  # Corner preset buttons
  observeEvent(input$leg_tl, { updateNumericInput(session, "legend_x", value = 0.05); updateNumericInput(session, "legend_y", value = 0.95) })
  observeEvent(input$leg_tr, { updateNumericInput(session, "legend_x", value = 0.85); updateNumericInput(session, "legend_y", value = 0.95) })
  observeEvent(input$leg_bl, { updateNumericInput(session, "legend_x", value = 0.05); updateNumericInput(session, "legend_y", value = 0.10) })
  observeEvent(input$leg_br, { updateNumericInput(session, "legend_x", value = 0.85); updateNumericInput(session, "legend_y", value = 0.10) })
  
  output$od_plot <- renderPlot({
    plot_inputs_debounced()
    p <- generate_plot()
    g <- tryCatch(
      fix_panel_size(p,
                     base_w_px = if (!is.null(input$plot_width))  input$plot_width  else 820,
                     base_h_px = if (!is.null(input$plot_height)) input$plot_height else 600),
      error = function(e) ggplotGrob(p)
    )
    grid::grid.draw(g)
  }, res = 96)
  
  output$downloadPlot <- downloadHandler(
    filename = function()
      paste0("od_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$export_format),
    content  = function(file) {
      p    <- generate_plot()
      w_in <- input$export_width
      h_in <- input$export_height
      dpi  <- input$export_dpi
      fmt  <- input$export_format
      # Extra inches for the legend (use screen px estimate → convert to inches)
      extra_in <- legend_extra_px() / 96
      g <- tryCatch(
        fix_panel_size(p, base_w_px = w_in * 96, base_h_px = h_in * 96),
        error = function(e) ggplotGrob(p)
      )
      total_w_in <- w_in + extra_in
      switch(fmt,
        pdf  = grDevices::pdf( file, width = total_w_in, height = h_in),
        svg  = grDevices::svg( file, width = total_w_in, height = h_in),
        png  = grDevices::png( file, width = total_w_in * dpi, height = h_in * dpi,
                               res = dpi, units = "px"),
        tiff = grDevices::tiff(file, width = total_w_in * dpi, height = h_in * dpi,
                               res = dpi, units = "px"),
        jpeg = grDevices::jpeg(file, width = total_w_in * dpi, height = h_in * dpi,
                               res = dpi, units = "px"),
        grDevices::pdf(file, width = total_w_in, height = h_in)
      )
      tryCatch(grid::grid.draw(g), finally = grDevices::dev.off())
    }
  )
  
  # ── PowerPoint Export ────────────────────────────────────────────────────────
  if (has_officer && has_rvg) {
    output$downloadPPTX <- downloadHandler(
      filename = function()
        paste0("od_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pptx"),
      content = function(file) {
        library(officer); library(rvg)
        req(rv$data, input$selected_samples)
        samps  <- input$selected_samples
        aes_v  <- resolve_aesthetics(samps)
        pd_all <- prepare_plot_data(samples = samps)
        sw     <- max(input$export_width,  4)
        sh     <- max(input$export_height, 3)
        
        title_h  <- 0.5
        plot_top <- title_h
        plot_h   <- sh - title_h
        
        title_pt <- max(12, round(sw * 2.2))
        
        prs <- read_pptx()
        
        for (k in seq_along(samps)) {
          visible_samps <- samps[1:k]
          new_samp      <- samps[k]
          
          aes_k <- list(
            shapes     = aes_v$shapes[    visible_samps],
            colors     = aes_v$colors[    visible_samps],
            line_types = aes_v$line_types[visible_samps],
            leg_labels = aes_v$leg_labels[visible_samps],
            filled_map = aes_v$filled_map[visible_samps]
          )
          
          pd_k <- pd_all[pd_all$variable %in% visible_samps, , drop = FALSE]
          
          p_k <- build_plot(pd_k, visible_samps, aes_k, highlight_samples = NULL)
          
          prs <- add_slide(prs, layout = "Blank", master = "Office Theme")
          
          prs <- ph_with(prs,
                         value    = new_samp,
                         location = ph_location(left = 0, top = 0, width = sw, height = title_h))
          
          prs <- ph_with(prs,
                         dml(ggobj = p_k, bg = "white"),
                         location = ph_location(left = 0, top = plot_top,
                                                width = sw, height = plot_h))
        }
        
        p_final <- build_plot(pd_all, samps, aes_v, highlight_samples = NULL)
        prs <- add_slide(prs, layout = "Blank", master = "Office Theme")
        prs <- ph_with(prs,
                       dml(ggobj = p_final, bg = "white"),
                       location = ph_location(left = 0, top = 0, width = sw, height = sh))

        # Optional experiment notes slide
        if (isTRUE(input$notes_pptx_slide)) {
          n <- get_notes()
          fmt_field <- function(label, val) {
            if (nchar(trimws(val)) > 0) paste0(label, ": ", trimws(val)) else NULL
          }
          lines <- c(
            fmt_field("Experiment ID",    n$exp_id),
            fmt_field("Date",             n$date),
            fmt_field("Experimenter",     n$experimenter),
            fmt_field("Project",          n$project),
            fmt_field("Institution",      n$institution),
            "",
            fmt_field("Experiment Type",  n$exp_type),
            fmt_field("Host Strain(s)",   n$host_strain),
            fmt_field("Phage / Plasmid",  n$phage_plasmid),
            fmt_field("MOI",              n$moi),
            fmt_field("Replicate / Batch",n$replicate),
            fmt_field("Passage #",        n$passage),
            "",
            fmt_field("Growth Medium",    n$media),
            fmt_field("Temperature (°C)", n$temperature),
            fmt_field("Inducer",          n$inducer),
            fmt_field("Inducer Conc.",    n$inducer_conc),
            fmt_field("Time of Inf/Ind",  n$time_inf),
            fmt_field("Other Condition",  n$extra_cond)
          )
          if (nchar(trimws(n$custom_key1)) > 0) lines <- c(lines, fmt_field(n$custom_key1, n$custom_val1))
          if (nchar(trimws(n$custom_key2)) > 0) lines <- c(lines, fmt_field(n$custom_key2, n$custom_val2))
          if (nchar(trimws(n$custom_key3)) > 0) lines <- c(lines, fmt_field(n$custom_key3, n$custom_val3))
          free_sections <- c(
            if (nchar(trimws(n$observations)) > 0) c("", paste0("Observations: ", trimws(n$observations))) else NULL,
            if (nchar(trimws(n$issues))       > 0) c("", paste0("Issues: ",       trimws(n$issues)))       else NULL,
            if (nchar(trimws(n$next_steps))   > 0) c("", paste0("Next Steps: ",   trimws(n$next_steps)))   else NULL,
            if (nchar(trimws(n$tags))         > 0) c("", paste0("Tags: ",         trimws(n$tags)))         else NULL
          )
          lines <- c(Filter(Negate(is.null), lines), free_sections)
          body_text <- paste(lines, collapse = "\n")
          prs <- add_slide(prs, layout = "Blank", master = "Office Theme")
          prs <- ph_with(prs,
                         value    = "Experiment Notes",
                         location = ph_location(left = 0.3, top = 0.2, width = sw - 0.6, height = 0.5))
          prs <- ph_with(prs,
                         value    = body_text,
                         location = ph_location(left = 0.3, top = 0.8, width = sw - 0.6, height = sh - 1.0))
        }

        print(prs, target = file)
      }
    )
  }
  
  # ── GIF Export ───────────────────────────────────────────────────────────────
  if (has_gifski) {
    output$downloadGIF <- downloadHandler(
      filename = function()
        paste0("od_anim_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".gif"),
      content = function(file) {
        library(gifski)
        req(rv$data, input$selected_samples)
        samps  <- input$selected_samples
        aes_v  <- resolve_aesthetics(samps)
        pd_all <- prepare_plot_data(samples = samps)
        w_px   <- max(as.integer(input$gif_width),  200L)
        h_px   <- max(as.integer(input$gif_height), 150L)
        fps    <- max(input$gif_fps, 0.2)
        
        tmp_dir <- tempfile(); dir.create(tmp_dir)
        on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
        
        frame_files <- character(length(samps))
        for (k in seq_along(samps)) {
          visible_samps <- samps[1:k]
          
          aes_k <- list(
            shapes     = aes_v$shapes[    visible_samps],
            colors     = aes_v$colors[    visible_samps],
            line_types = aes_v$line_types[visible_samps],
            leg_labels = aes_v$leg_labels[visible_samps],
            filled_map = aes_v$filled_map[visible_samps]
          )
          
          pd_k <- pd_all[pd_all$variable %in% visible_samps, , drop = FALSE]
          
          p_k   <- build_plot(pd_k, visible_samps, aes_k, highlight_samples = NULL)
          fpath <- file.path(tmp_dir, sprintf("frame_%03d.png", k))
          ggsave(fpath, plot = p_k, device = "png",
                 width = w_px / 96, height = h_px / 96,
                 units = "in", dpi = 96)
          frame_files[k] <- fpath
        }
        
        gifski(frame_files, gif_file = file,
               width = w_px, height = h_px,
               delay = 1 / fps, loop = TRUE)
      }
    )
  }
  
  # ══════════════════════════════════════════════════════════════════════════════
  # ── Analysis Tab Server Logic ───────────────────────────────────────────────
  # ══════════════════════════════════════════════════════════════════════════════
  
  # ── Experiment Notes helpers ─────────────────────────────────────────────────

  # Collects all note inputs into a named list
  get_notes <- function() {
    list(
      exp_id        = input$notes_exp_id,
      date          = as.character(input$notes_date),
      experimenter  = input$notes_experimenter,
      project       = input$notes_project,
      institution   = input$notes_institution,
      exp_type      = input$notes_exp_type,
      host_strain   = input$notes_host_strain,
      phage_plasmid = input$notes_phage_plasmid,
      moi           = input$notes_moi,
      replicate     = input$notes_replicate,
      passage       = input$notes_passage,
      media         = input$notes_media,
      temperature   = input$notes_temperature,
      time_inf      = input$notes_time_inf,
      inducer       = input$notes_inducer,
      inducer_conc  = input$notes_inducer_conc,
      extra_cond    = input$notes_extra_cond,
      observations  = input$notes_observations,
      issues        = input$notes_issues,
      next_steps    = input$notes_next_steps,
      tags          = input$notes_tags,
      custom_key1   = input$notes_custom_key1,
      custom_val1   = input$notes_custom_val1,
      custom_key2   = input$notes_custom_key2,
      custom_val2   = input$notes_custom_val2,
      custom_key3   = input$notes_custom_key3,
      custom_val3   = input$notes_custom_val3
    )
  }

  # Builds the caption string from selected fields
  build_notes_caption <- function(notes = NULL, fields = NULL) {
    if (is.null(notes))  notes  <- get_notes()
    if (is.null(fields)) fields <- input$notes_caption_fields
    if (is.null(fields) || length(fields) == 0) return("")
    field_map <- list(
      exp_id        = c("ID",          notes$exp_id),
      date          = c("Date",        notes$date),
      experimenter  = c("By",          notes$experimenter),
      exp_type      = c("Type",        notes$exp_type),
      host_strain   = c("Host",        notes$host_strain),
      phage_plasmid = c("Phage/Plasmid", notes$phage_plasmid),
      moi           = c("MOI",         notes$moi),
      replicate     = c("Rep",         notes$replicate),
      inducer       = c("Inducer",     notes$inducer),
      inducer_conc  = c("Conc",        notes$inducer_conc),
      media         = c("Media",       notes$media),
      temperature   = c("Temp",        notes$temperature)
    )
    parts <- vapply(fields, function(f) {
      m <- field_map[[f]]
      if (!is.null(m) && nchar(trimws(m[2])) > 0) paste0(m[1], ": ", trimws(m[2])) else ""
    }, character(1))
    parts <- parts[nchar(parts) > 0]
    if (length(parts) == 0) return("")
    paste(parts, collapse = "  |  ")
  }

  # Caption preview
  output$notes_caption_preview <- renderText({
    cap <- build_notes_caption()
    if (nchar(cap) == 0) "(no fields selected or all fields empty)" else cap
  })

  # TXT download
  output$download_notes_txt <- downloadHandler(
    filename = function() {
      id  <- trimws(input$notes_exp_id)
      pfx <- if (nchar(id) > 0) paste0(gsub("[^A-Za-z0-9_-]", "_", id), "_") else ""
      paste0(pfx, "notes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      n   <- get_notes()
      hdr <- paste0(strrep("=", 60), "\n",
                    "EXPERIMENT NOTES\n",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
                    strrep("=", 60), "\n\n")
      section <- function(title, pairs) {
        lines <- vapply(pairs, function(p) {
          if (nchar(trimws(p[2])) > 0) sprintf("  %-24s %s", paste0(p[1], ":"), trimws(p[2])) else ""
        }, character(1))
        lines <- lines[nchar(lines) > 0]
        if (length(lines) == 0) return("")
        paste0(title, "\n", strrep("-", 40), "\n", paste(lines, collapse = "\n"), "\n\n")
      }
      identity_sec <- section("EXPERIMENT IDENTITY", list(
        c("Experiment ID",  n$exp_id),       c("Date",       n$date),
        c("Experimenter",   n$experimenter),  c("Project",    n$project),
        c("Institution",    n$institution)
      ))
      biology_sec <- section("BIOLOGICAL PARAMETERS", list(
        c("Experiment Type",  n$exp_type),   c("Host Strain(s)",  n$host_strain),
        c("Phage / Plasmid",  n$phage_plasmid), c("MOI",          n$moi),
        c("Replicate / Batch",n$replicate),  c("Passage #",       n$passage)
      ))
      cond_sec <- section("CONDITIONS", list(
        c("Growth Medium",   n$media),       c("Temperature (°C)",  n$temperature),
        c("Inducer",         n$inducer),     c("Inducer Conc.",      n$inducer_conc),
        c("Time of Inf/Ind", n$time_inf),    c("Other Condition",    n$extra_cond)
      ))
      custom_pairs <- list(
        if (nchar(trimws(n$custom_key1)) > 0) c(n$custom_key1, n$custom_val1) else NULL,
        if (nchar(trimws(n$custom_key2)) > 0) c(n$custom_key2, n$custom_val2) else NULL,
        if (nchar(trimws(n$custom_key3)) > 0) c(n$custom_key3, n$custom_val3) else NULL
      )
      custom_pairs <- Filter(Negate(is.null), custom_pairs)
      custom_sec <- if (length(custom_pairs) > 0) section("CUSTOM FIELDS", custom_pairs) else ""
      free_sec <- ""
      for (nm in c("observations", "issues", "next_steps")) {
        val <- trimws(n[[nm]])
        if (nchar(val) > 0) {
          lbl <- switch(nm, observations = "OBSERVATIONS / RESULTS",
                            issues       = "ISSUES / ANOMALIES",
                            next_steps   = "NEXT STEPS")
          free_sec <- paste0(free_sec, lbl, "\n", strrep("-", 40), "\n", val, "\n\n")
        }
      }
      tags_sec <- if (nchar(trimws(n$tags)) > 0)
        paste0("TAGS\n", strrep("-", 40), "\n  ", trimws(n$tags), "\n\n") else ""
      writeLines(paste0(hdr, identity_sec, biology_sec, cond_sec, custom_sec, free_sec, tags_sec), file)
    }
  )

  # CSV download (one-row log entry)
  output$download_notes_csv <- downloadHandler(
    filename = function() {
      id  <- trimws(input$notes_exp_id)
      pfx <- if (nchar(id) > 0) paste0(gsub("[^A-Za-z0-9_-]", "_", id), "_") else ""
      paste0(pfx, "notes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      n <- get_notes()
      custom_cols <- list()
      if (nchar(trimws(n$custom_key1)) > 0) custom_cols[[trimws(n$custom_key1)]] <- trimws(n$custom_val1)
      if (nchar(trimws(n$custom_key2)) > 0) custom_cols[[trimws(n$custom_key2)]] <- trimws(n$custom_val2)
      if (nchar(trimws(n$custom_key3)) > 0) custom_cols[[trimws(n$custom_key3)]] <- trimws(n$custom_val3)
      row <- data.frame(
        ExperimentID      = n$exp_id,
        Date              = n$date,
        Experimenter      = n$experimenter,
        Project           = n$project,
        Institution       = n$institution,
        ExperimentType    = n$exp_type,
        HostStrain        = n$host_strain,
        PhagePlasmid      = n$phage_plasmid,
        MOI               = n$moi,
        Replicate         = n$replicate,
        Passage           = n$passage,
        GrowthMedium      = n$media,
        Temperature_C     = n$temperature,
        TimeOfInfection   = n$time_inf,
        Inducer           = n$inducer,
        InducerConc       = n$inducer_conc,
        OtherCondition    = n$extra_cond,
        Observations      = n$observations,
        Issues            = n$issues,
        NextSteps         = n$next_steps,
        Tags              = n$tags,
        ExportedAt        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors  = FALSE
      )
      if (length(custom_cols) > 0) {
        for (k in names(custom_cols)) row[[k]] <- custom_cols[[k]]
      }
      write.csv(row, file, row.names = FALSE)
    }
  )

  # ── Notes: shared renderer (PDF + PNG) ───────────────────────────────────────
  # Draws a formatted lab-notes page onto the currently-open graphics device.
  draw_notes_page <- function(n) {
    grid::grid.newpage()
    # Outer margin viewport
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.93, height = 0.95))
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

    # ── Title bar ──────────────────────────────────────────────────────────────
    grid::grid.rect(x = 0, y = 1, width = 1, height = 0.065,
                    just = c("left", "top"),
                    gp = grid::gpar(fill = "#2c3e50", col = NA))
    title_text <- if (nchar(trimws(n$exp_id)) > 0)
      paste0("Experiment Notes  —  ", trimws(n$exp_id)) else "Experiment Notes"
    grid::grid.text(title_text, x = 0.015, y = 0.978,
                    just = c("left", "top"),
                    gp = grid::gpar(fontsize = 13, col = "white", fontface = "bold"))
    grid::grid.text(format(Sys.time(), "%Y-%m-%d %H:%M"),
                    x = 0.985, y = 0.978, just = c("right", "top"),
                    gp = grid::gpar(fontsize = 9, col = "#bdc3c7"))

    # ── Helpers ────────────────────────────────────────────────────────────────
    y <- 0.918   # current vertical position (0 = bottom, 1 = top)
    lh <- 0.030  # line height

    add_section <- function(label) {
      y <<- y - 0.008
      grid::grid.rect(x = 0, y = y, width = 1, height = 0.030,
                      just = c("left", "top"),
                      gp = grid::gpar(fill = "#ecf0f1", col = NA))
      grid::grid.text(toupper(label), x = 0.012, y = y - 0.005,
                      just = c("left", "top"),
                      gp = grid::gpar(fontsize = 8, fontface = "bold", col = "#2c3e50"))
      y <<- y - 0.035
    }

    add_field <- function(label, value, indent = 0.012, label_w = 0.30) {
      val <- trimws(value)
      if (nchar(val) == 0) return(invisible(NULL))
      grid::grid.text(paste0(label, ":"), x = indent, y = y,
                      just = c("left", "top"),
                      gp = grid::gpar(fontsize = 8, col = "#555555"))
      grid::grid.text(val, x = indent + label_w, y = y,
                      just = c("left", "top"),
                      gp = grid::gpar(fontsize = 8, col = "#111111"))
      y <<- y - lh
    }

    add_two <- function(l1, v1, l2, v2) {
      v1 <- trimws(v1); v2 <- trimws(v2)
      if (nchar(v1) == 0 && nchar(v2) == 0) return(invisible(NULL))
      if (nchar(v1) > 0) {
        grid::grid.text(paste0(l1, ":"), x = 0.012, y = y, just = c("left","top"),
                        gp = grid::gpar(fontsize = 8, col = "#555555"))
        grid::grid.text(v1, x = 0.22, y = y, just = c("left","top"),
                        gp = grid::gpar(fontsize = 8, col = "#111111"))
      }
      if (nchar(v2) > 0) {
        grid::grid.text(paste0(l2, ":"), x = 0.52, y = y, just = c("left","top"),
                        gp = grid::gpar(fontsize = 8, col = "#555555"))
        grid::grid.text(v2, x = 0.72, y = y, just = c("left","top"),
                        gp = grid::gpar(fontsize = 8, col = "#111111"))
      }
      y <<- y - lh
    }

    add_block <- function(label, value, indent = 0.012) {
      val <- trimws(value)
      if (nchar(val) == 0) return(invisible(NULL))
      grid::grid.text(paste0(label, ":"), x = indent, y = y,
                      just = c("left", "top"),
                      gp = grid::gpar(fontsize = 8, fontface = "bold", col = "#333333"))
      y <<- y - lh
      # Wrap text at ~100 chars and render each line
      wrapped <- strwrap(val, width = 100)
      for (ln in wrapped) {
        grid::grid.text(ln, x = indent + 0.015, y = y,
                        just = c("left", "top"),
                        gp = grid::gpar(fontsize = 8, col = "#111111"))
        y <<- y - lh
      }
      y <<- y - 0.006
    }

    # ── Section 1: Identity ────────────────────────────────────────────────────
    add_section("Experiment Identity")
    add_two("Date", n$date, "Experimenter", n$experimenter)
    add_two("Project", n$project, "Institution", n$institution)
    add_field("Experiment Type", n$exp_type)

    # ── Section 2: Biology ────────────────────────────────────────────────────
    add_section("Biological Parameters")
    add_two("Host Strain", n$host_strain, "Phage / Plasmid", n$phage_plasmid)
    add_two("MOI", n$moi, "Replicate / Batch", n$replicate)
    add_field("Passage #", n$passage)

    # ── Section 3: Conditions ─────────────────────────────────────────────────
    add_section("Conditions")
    add_two("Growth Medium", n$media, "Temperature (°C)", n$temperature)
    add_two("Time of Infection", n$time_inf, "Inducer", n$inducer)
    add_two("Inducer Conc.", n$inducer_conc, "Other", n$extra_cond)

    # ── Custom fields ─────────────────────────────────────────────────────────
    has_custom <- any(nchar(trimws(c(n$custom_key1, n$custom_key2, n$custom_key3))) > 0)
    if (has_custom) {
      add_section("Custom Fields")
      if (nchar(trimws(n$custom_key1)) > 0) add_field(trimws(n$custom_key1), n$custom_val1)
      if (nchar(trimws(n$custom_key2)) > 0) add_field(trimws(n$custom_key2), n$custom_val2)
      if (nchar(trimws(n$custom_key3)) > 0) add_field(trimws(n$custom_key3), n$custom_val3)
    }

    # ── Section 4: Free text ──────────────────────────────────────────────────
    add_section("Notes")
    add_block("Observations / Results", n$observations)
    add_block("Issues / Anomalies",     n$issues)
    add_block("Next Steps",             n$next_steps)
    if (nchar(trimws(n$tags)) > 0) add_field("Tags", n$tags)

    # ── Footer line ───────────────────────────────────────────────────────────
    grid::grid.lines(x = c(0, 1), y = c(0.012, 0.012),
                     gp = grid::gpar(col = "#bdc3c7", lwd = 0.5))
    grid::grid.text("Generated by Lysis Curve App", x = 0.5, y = 0.005,
                    just = c("center", "bottom"),
                    gp = grid::gpar(fontsize = 7, col = "#aaaaaa"))
    grid::popViewport()
  }

  # ── Notes: PDF ────────────────────────────────────────────────────────────────
  output$download_notes_pdf <- downloadHandler(
    filename = function() {
      id  <- trimws(input$notes_exp_id)
      pfx <- if (nchar(id) > 0) paste0(gsub("[^A-Za-z0-9_-]", "_", id), "_") else ""
      paste0(pfx, "notes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      n <- get_notes()
      grDevices::pdf(file, width = 8.5, height = 11, paper = "letter")
      tryCatch(draw_notes_page(n), finally = grDevices::dev.off())
    }
  )

  # ── Notes: PNG (image / photo) ────────────────────────────────────────────────
  output$download_notes_png <- downloadHandler(
    filename = function() {
      id  <- trimws(input$notes_exp_id)
      pfx <- if (nchar(id) > 0) paste0(gsub("[^A-Za-z0-9_-]", "_", id), "_") else ""
      paste0(pfx, "notes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      n <- get_notes()
      grDevices::png(file, width = 2550, height = 3300, res = 300, units = "px")
      tryCatch(draw_notes_page(n), finally = grDevices::dev.off())
    }
  )

  # ── Notes: HTML ──────────────────────────────────────────────────────────────
  output$download_notes_html <- downloadHandler(
    filename = function() {
      id  <- trimws(input$notes_exp_id)
      pfx <- if (nchar(id) > 0) paste0(gsub("[^A-Za-z0-9_-]", "_", id), "_") else ""
      paste0(pfx, "notes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
    },
    content = function(file) {
      n   <- get_notes()
      ts  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      ttl <- if (nchar(trimws(n$exp_id)) > 0) htmltools::htmlEscape(trimws(n$exp_id)) else "Experiment Notes"

      row_html <- function(label, value) {
        val <- trimws(value)
        if (nchar(val) == 0) return("")
        sprintf('<tr><td style="color:#555;padding:4px 10px 4px 0;width:35%%;font-weight:600;">%s</td><td style="padding:4px 0;">%s</td></tr>',
                htmltools::htmlEscape(label),
                htmltools::htmlEscape(val))
      }
      section_html <- function(title, rows_html) {
        rows <- paste(rows_html[nchar(rows_html) > 0], collapse = "\n")
        if (nchar(rows) == 0) return("")
        sprintf('<h3 style="background:#2c3e50;color:white;padding:6px 12px;margin:18px 0 6px;font-size:13px;border-radius:3px;">%s</h3><table style="width:100%%;border-collapse:collapse;font-size:13px;">%s</table>',
                htmltools::htmlEscape(title), rows)
      }
      block_html <- function(label, value) {
        val <- trimws(value)
        if (nchar(val) == 0) return("")
        sprintf('<div style="margin-bottom:12px;"><b style="color:#2c3e50;">%s</b><p style="margin:4px 0 0 12px;white-space:pre-wrap;font-size:13px;">%s</p></div>',
                htmltools::htmlEscape(label),
                htmltools::htmlEscape(val))
      }
      custom_rows <- c(
        if (nchar(trimws(n$custom_key1)) > 0) row_html(trimws(n$custom_key1), n$custom_val1) else "",
        if (nchar(trimws(n$custom_key2)) > 0) row_html(trimws(n$custom_key2), n$custom_val2) else "",
        if (nchar(trimws(n$custom_key3)) > 0) row_html(trimws(n$custom_key3), n$custom_val3) else ""
      )
      html <- paste0(
        '<!DOCTYPE html><html><head><meta charset="UTF-8">',
        '<title>', ttl, '</title>',
        '<style>body{font-family:Arial,sans-serif;max-width:800px;margin:30px auto;color:#222;line-height:1.5;}',
        'h1{background:#2c3e50;color:white;padding:14px 18px;margin:0 0 4px;border-radius:4px;}',
        '.meta{color:#888;font-size:12px;margin-bottom:20px;padding-left:4px;}</style></head><body>',
        '<h1>', ttl, '</h1><p class="meta">Generated: ', htmltools::htmlEscape(ts), '</p>',
        section_html("Experiment Identity", c(
          row_html("Date", n$date), row_html("Experimenter", n$experimenter),
          row_html("Project", n$project), row_html("Institution", n$institution),
          row_html("Experiment Type", n$exp_type))),
        section_html("Biological Parameters", c(
          row_html("Host Strain", n$host_strain), row_html("Phage / Plasmid", n$phage_plasmid),
          row_html("MOI", n$moi), row_html("Replicate / Batch", n$replicate),
          row_html("Passage #", n$passage))),
        section_html("Conditions", c(
          row_html("Growth Medium", n$media), row_html("Temperature (°C)", n$temperature),
          row_html("Time of Infection", n$time_inf), row_html("Inducer", n$inducer),
          row_html("Inducer Concentration", n$inducer_conc),
          row_html("Other Condition", n$extra_cond))),
        if (any(nchar(trimws(custom_rows)) > 0)) section_html("Custom Fields", custom_rows) else "",
        '<h3 style="background:#2c3e50;color:white;padding:6px 12px;margin:18px 0 6px;font-size:13px;border-radius:3px;">Notes</h3>',
        block_html("Observations / Results", n$observations),
        block_html("Issues / Anomalies", n$issues),
        block_html("Next Steps", n$next_steps),
        if (nchar(trimws(n$tags)) > 0) sprintf('<p style="font-size:12px;color:#666;"><b>Tags:</b> %s</p>', htmltools::htmlEscape(n$tags)) else "",
        '<hr style="margin-top:30px;border:none;border-top:1px solid #ddd;">',
        '<p style="font-size:11px;color:#aaa;text-align:center;">Generated by Lysis Curve App</p>',
        '</body></html>'
      )
      writeLines(html, file)
    }
  )

  # ── Notes: JSON ──────────────────────────────────────────────────────────────
  output$download_notes_json <- downloadHandler(
    filename = function() {
      id  <- trimws(input$notes_exp_id)
      pfx <- if (nchar(id) > 0) paste0(gsub("[^A-Za-z0-9_-]", "_", id), "_") else ""
      paste0(pfx, "notes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json")
    },
    content = function(file) {
      n  <- get_notes()
      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      # Build JSON manually — no extra package required
      esc <- function(x) gsub('"', '\\"', gsub("\\\\", "\\\\\\\\", trimws(x)))
      field <- function(key, val) sprintf('  "%s": "%s"', key, esc(val))
      custom_fields <- c(
        if (nchar(trimws(n$custom_key1)) > 0) field(esc(trimws(n$custom_key1)), n$custom_val1) else NULL,
        if (nchar(trimws(n$custom_key2)) > 0) field(esc(trimws(n$custom_key2)), n$custom_val2) else NULL,
        if (nchar(trimws(n$custom_key3)) > 0) field(esc(trimws(n$custom_key3)), n$custom_val3) else NULL
      )
      lines <- c(
        "{",
        field("exported_at",    ts),    ",",
        field("exp_id",         n$exp_id),        ",",
        field("date",           n$date),           ",",
        field("experimenter",   n$experimenter),   ",",
        field("project",        n$project),        ",",
        field("institution",    n$institution),    ",",
        field("exp_type",       n$exp_type),       ",",
        field("host_strain",    n$host_strain),    ",",
        field("phage_plasmid",  n$phage_plasmid),  ",",
        field("moi",            n$moi),            ",",
        field("replicate",      n$replicate),      ",",
        field("passage",        n$passage),        ",",
        field("media",          n$media),          ",",
        field("temperature",    n$temperature),    ",",
        field("time_inf",       n$time_inf),       ",",
        field("inducer",        n$inducer),        ",",
        field("inducer_conc",   n$inducer_conc),   ",",
        field("extra_cond",     n$extra_cond),     ",",
        field("observations",   n$observations),   ",",
        field("issues",         n$issues),         ",",
        field("next_steps",     n$next_steps),     ",",
        field("tags",           n$tags)
      )
      if (length(custom_fields) > 0) {
        lines[length(lines)] <- paste0(lines[length(lines)], ",")
        lines <- c(lines, paste(custom_fields, collapse = ",\n"))
      }
      lines <- c(lines, "}")
      writeLines(lines, file)
    }
  )

  rv_analysis <- reactiveValues(
    metrics        = NULL,
    metrics_raw    = NULL,
    metrics_inf    = NULL,
    metrics_inf_raw = NULL,
    metrics_long   = NULL,
    stats_result   = NULL,
    deriv_data     = NULL
  )
  
  observe({
    samps <- input$selected_samples
    if (!is.null(samps) && length(samps) > 0)
      updateSelectInput(session, "ref_sample", choices = samps, selected = samps[1])
  })
  
  # ── Calculate Metrics ────────────────────────────────────────────────────────
  observeEvent(input$calc_metrics, {
    req(rv$data, input$selected_samples)
    pd <- tryCatch(prepare_metrics_data(), error = function(e) NULL)
    if (is.null(pd) || nrow(pd) == 0) return()
    withProgress(message = "Calculating growth metrics...", value = 0.2, {
      m_raw <- calculate_growth_metrics(pd, smooth_window = input$metrics_smooth_window)
      setProgress(0.6, detail = "Computing derivatives...")
      rv_analysis$deriv_data <- calculate_derivative(prepare_plot_data(), smooth_window = input$metrics_smooth_window)
      setProgress(0.8, detail = "Summarising...")
      m_summary <- if (!is.null(m_raw) && "replicate" %in% names(m_raw)) {
        m_raw %>% group_by(sample) %>%
          summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
      } else m_raw
      rv_analysis$metrics        <- m_summary
      rv_analysis$metrics_raw    <- m_raw
      rv_analysis$metrics_inf    <- NULL
      rv_analysis$metrics_inf_raw <- NULL
      rv_analysis$metrics_long   <- summarise_metrics(m_raw)
      rv_analysis$stats_result   <- NULL
      updateSelectInput(session, "barplot_metric",
                        choices = c(core_metric_choices, lysis_metric_choices))
      updateSelectInput(session, "stats_metric",
                        choices = c(core_metric_choices, lysis_metric_choices))
      setProgress(1.0)
    })
  })
  
  # ── Infection Strength ───────────────────────────────────────────────────────
  observeEvent(input$calc_infection, {
    req(rv_analysis$metrics_raw, input$ref_sample)
    m_inf_raw <- calculate_infection_metrics(rv_analysis$metrics_raw, input$ref_sample)
    if (!is.null(m_inf_raw)) {
      m_inf_summary <- if ("replicate" %in% names(m_inf_raw)) {
        m_inf_raw %>% group_by(sample) %>%
          summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")
      } else m_inf_raw
      rv_analysis$metrics_inf     <- m_inf_summary
      rv_analysis$metrics_inf_raw <- m_inf_raw
      rv_analysis$metrics_long    <- summarise_metrics(m_inf_raw)
      rv_analysis$stats_result    <- NULL
      updateSelectInput(session, "barplot_metric",
                        choices = c(core_metric_choices, lysis_metric_choices, infection_metric_choices))
      updateSelectInput(session, "stats_metric",
                        choices = c(core_metric_choices, lysis_metric_choices, infection_metric_choices))
    }
  })
  
  # ── Metrics Table ────────────────────────────────────────────────────────────
  output$metrics_table <- DT::renderDT({
    req(rv_analysis$metrics)
    m <- rv_analysis$metrics
    num_cols <- names(m)[sapply(m, is.numeric)]
    for (col in num_cols) m[[col]] <- round(m[[col]], 4)
    display_names <- c(sample="Sample", initial_od="Initial OD", max_od="Max OD",
                       time_max_od="Time Max OD", final_od="Final OD", auc="AUC",
                       mu_max="μmax (h⁻¹)", doubling_time="Doubling Time",
                       lag_phase="Lag Phase", stat_phase_dur="Stationary Dur.",
                       lysis_time="Lysis Time", lysis_rate="Lysis Rate",
                       od_drop="OD Drop", residual_od="Residual OD", recovery_slope="Recovery Slope")
    nms <- names(m)
    new_names <- ifelse(nms %in% names(display_names), display_names[nms], nms)
    datatable(m, colnames = new_names, rownames = FALSE,
              options = list(scrollX = TRUE, pageLength = 50, dom = "Bfrtip"),
              class = "compact stripe hover")
  })
  
  # ── Infection Table ──────────────────────────────────────────────────────────
  output$infection_table <- DT::renderDT({
    req(rv_analysis$metrics_inf)
    m <- rv_analysis$metrics_inf
    show_cols <- intersect(c("sample","infection_strength","relative_growth",
                             "relative_mu_max","relative_max_od","lysis_onset_delta"), names(m))
    m <- m[, show_cols, drop = FALSE]
    num_cols <- names(m)[sapply(m, is.numeric)]
    for (col in num_cols) m[[col]] <- round(m[[col]], 4)
    datatable(m, rownames = FALSE,
              options = list(scrollX = TRUE, pageLength = 50, dom = "Bfrtip"),
              class = "compact stripe hover")
  })
  
  # ── Bar Plot Sample Selection ────────────────────────────────────────────────
  output$barplot_sample_selector <- renderUI({
    req(rv_analysis$metrics)
    samps <- unique(rv_analysis$metrics$sample)
    checkboxGroupInput("barplot_samples", "Show samples in bar plot:",
                       choices = samps, selected = samps, inline = TRUE)
  })
  observeEvent(input$barplot_select_all, {
    req(rv_analysis$metrics)
    updateCheckboxGroupInput(session, "barplot_samples",
                             selected = unique(rv_analysis$metrics$sample))
  })
  observeEvent(input$barplot_deselect_all, {
    updateCheckboxGroupInput(session, "barplot_samples", selected = character(0))
  })

  # ── Bar Plot Custom Order Widget ───────────────────────────────────────────
  output$barplot_custom_order_ui <- renderUI({
    req(rv_analysis$metrics)
    current_samps <- if (!is.null(input$barplot_samples) && length(input$barplot_samples) > 0)
      input$barplot_samples
    else
      as.character(unique(rv_analysis$metrics$sample))
    selectInput("barplot_custom_order", label = NULL,
                choices = current_samps, selected = current_samps,
                multiple = TRUE)
  })

  # ── Stats Sample Selector ──────────────────────────────────────────────────
  output$stats_sample_selector <- renderUI({
    req(rv_analysis$metrics)
    samps <- unique(rv_analysis$metrics$sample)
    checkboxGroupInput("stats_samples", "Samples to Compare:",
                       choices = samps, selected = samps, inline = TRUE)
  })
  observeEvent(input$stats_select_all, {
    req(rv_analysis$metrics)
    updateCheckboxGroupInput(session, "stats_samples",
                             selected = unique(rv_analysis$metrics$sample))
  })
  observeEvent(input$stats_deselect_all, {
    updateCheckboxGroupInput(session, "stats_samples", selected = character(0))
  })

  # ── Bar plot color pickers ───────────────────────────────────────────────────
  output$barplot_color_pickers <- renderUI({
    req(rv_analysis$metrics)
    samps <- unique(rv_analysis$metrics$sample)
    def_colors <- if (!is.null(input$selected_samples)) {
      aes_v <- resolve_aesthetics(input$selected_samples); aes_v$colors
    } else setNames(rep("#333333", length(samps)), samps)
    tagList(
      p(style = "font-size:.82em;color:#555;margin-bottom:6px;", "Custom color per bar (HEX):"),
      fluidRow(lapply(seq_along(samps), function(i) {
        vn <- samps[i]; dc <- if (vn %in% names(def_colors)) def_colors[vn] else "#333333"
        column(3, textInput(paste0("barcol_", safe_id(vn)), label = vn, value = dc))
      }))
    )
  })
  
  get_bar_colors <- function() {
    req(rv_analysis$metrics)
    samps <- unique(rv_analysis$metrics$sample)
    if (isTRUE(input$barplot_custom_colors)) {
      cols <- sapply(samps, function(vn) {
        val <- input[[paste0("barcol_", safe_id(vn))]]
        if (!is.null(val) && nchar(trimws(val)) > 0) normalize_hex_color(val) else "#333333"
      })
      return(setNames(cols, samps))
    }
    if (isTRUE(input$barplot_use_sample_colors) && !is.null(input$selected_samples)) {
      aes_v <- resolve_aesthetics(input$selected_samples)
      return(aes_v$colors[names(aes_v$colors) %in% samps])
    }
    NULL
  }
  
  # ── Bar Plot ─────────────────────────────────────────────────────────────────
  output$barplot_container <- renderUI({
    req(rv_analysis$metrics)
    total_samps <- length(unique(rv_analysis$metrics$sample))
    active <- input$barplot_samples
    n_active <- if (!is.null(active) && length(active) > 0) length(active) else total_samps
    base_width  <- if (!is.null(input$barplot_width))  input$barplot_width  else 700
    base_height <- if (!is.null(input$barplot_height)) input$barplot_height else 500
    sizing_mode <- if (!is.null(input$barplot_sizing_mode)) input$barplot_sizing_mode else "shrink"
    is_horiz <- isTRUE(input$barplot_horizontal)
    if (sizing_mode == "shrink" && n_active < total_samps) {
      scaled <- max(300, round(base_width * n_active / total_samps))
      if (is_horiz) {
        plot_w <- base_width
        plot_h <- max(200, round(base_height * n_active / total_samps))
      } else {
        plot_w <- scaled
        plot_h <- base_height
      }
    } else {
      plot_w <- base_width
      plot_h <- base_height
    }
    plotOutput("metrics_barplot",
               height = paste0(plot_h, "px"),
               width  = paste0(plot_w, "px"))
  })

  build_barplot <- function() {
    req(rv_analysis$metrics, input$barplot_metric)
    met <- input$barplot_metric
    m_raw <- if (!is.null(rv_analysis$metrics_inf_raw) && met %in% names(rv_analysis$metrics_inf_raw))
      rv_analysis$metrics_inf_raw else rv_analysis$metrics_raw
    m_use <- if (!is.null(m_raw) && met %in% names(m_raw)) m_raw else rv_analysis$metrics
    if (!met %in% names(m_use)) return(NULL)

    plot_df <- m_use %>% select(sample, value = !!sym(met)) %>% filter(is.finite(value))
    if (nrow(plot_df) == 0) return(NULL)

    summary_df <- plot_df %>% group_by(sample) %>%
      summarise(mean_val = mean(value), sd_val = sd(value, na.rm = TRUE),
                n = n(), .groups = "drop") %>%
      mutate(sd_val = ifelse(is.na(sd_val), 0, sd_val),
             se_val = ifelse(n > 1, sd_val / sqrt(n), 0),
             ci_val = ifelse(n > 1, qt(0.975, df = n - 1) * se_val, 0))

    err_type <- if (!is.null(input$barplot_error_type)) input$barplot_error_type else "sem"
    summary_df$err <- switch(err_type, sd = summary_df$sd_val, sem = summary_df$se_val,
                             ci95 = summary_df$ci_val, none = 0, summary_df$se_val)

    bar_colors <- get_bar_colors()
    all_samps <- as.character(unique(rv_analysis$metrics$sample))
    active_samps <- if (!is.null(input$barplot_samples) && length(input$barplot_samples) > 0)
      input$barplot_samples else if (!is.null(input$selected_samples))
      input$selected_samples else as.character(summary_df$sample)
    samp_order <- intersect(active_samps, as.character(summary_df$sample))
    if (length(samp_order) == 0) return(NULL)

    # Apply sort order
    sort_mode <- if (!is.null(input$barplot_sort_order)) input$barplot_sort_order else "default"
    samp_order <- switch(sort_mode,
      "alpha"  = sort(samp_order),
      "asc"    = { ord <- summary_df %>% filter(as.character(sample) %in% samp_order) %>% arrange(mean_val); as.character(ord$sample) },
      "desc"   = { ord <- summary_df %>% filter(as.character(sample) %in% samp_order) %>% arrange(desc(mean_val)); as.character(ord$sample) },
      "custom" = {
        co <- input$barplot_custom_order
        if (!is.null(co) && length(co) > 0) {
          valid <- intersect(co, samp_order)
          c(valid, setdiff(samp_order, valid))
        } else {
          samp_order
        }
      },
      samp_order
    )

    plot_df    <- plot_df    %>% filter(as.character(sample) %in% samp_order)
    summary_df <- summary_df %>% filter(as.character(sample) %in% samp_order)

    # Fixed sizing mode: keep all sample levels so deselected bars leave gaps
    sizing_mode <- if (!is.null(input$barplot_sizing_mode)) input$barplot_sizing_mode else "shrink"
    if (sizing_mode == "fixed") {
      all_levels <- if (sort_mode == "alpha") {
        sort(all_samps)
      } else if (sort_mode == "custom" && !is.null(input$barplot_custom_order) && length(input$barplot_custom_order) > 0) {
        co <- input$barplot_custom_order
        c(intersect(co, all_samps), setdiff(all_samps, co))
      } else {
        all_samps
      }
      summary_df$sample <- factor(summary_df$sample, levels = all_levels)
      plot_df$sample    <- factor(plot_df$sample,    levels = all_levels)
    } else {
      summary_df$sample <- factor(summary_df$sample, levels = samp_order)
      plot_df$sample    <- factor(plot_df$sample,    levels = samp_order)
    }

    y_label <- if (met %in% names(metric_labels)) metric_labels[met] else met
    err_label <- switch(err_type, sd=" (\u00b1 SD)", sem=" (\u00b1 SEM)", ci95=" (\u00b1 95% CI)", "")

    bp_bar_width <- if (!is.null(input$barplot_bar_width_ctrl)) input$barplot_bar_width_ctrl else 0.7
    bp_bar_alpha <- if (!is.null(input$barplot_bar_alpha_ctrl)) input$barplot_bar_alpha_ctrl else 0.85
    bp_type <- if (!is.null(input$barplot_type)) input$barplot_type else "bar"

    # Build plot based on type
    if (bp_type == "bar") {
      p <- ggplot(summary_df, aes(x = sample, y = mean_val, fill = sample)) +
        geom_col(width = bp_bar_width, color = "black", linewidth = 0.3, alpha = bp_bar_alpha) +
        geom_jitter(data = plot_df, aes(x = sample, y = value), width = 0.15,
                    size = 2.5, shape = 21, fill = "white", color = "black",
                    stroke = 0.5, inherit.aes = FALSE, show.legend = FALSE)
      if (err_type != "none" && any(summary_df$err > 0))
        p <- p + geom_errorbar(aes(ymin = mean_val - err, ymax = mean_val + err),
                               width = 0.25, linewidth = 0.6, color = "black")

    } else if (bp_type == "dot") {
      p <- ggplot(summary_df, aes(x = sample, y = mean_val, color = sample)) +
        geom_jitter(data = plot_df, aes(x = sample, y = value), width = 0.15,
                    size = 2, shape = 16, alpha = 0.4,
                    inherit.aes = FALSE, show.legend = FALSE) +
        geom_point(size = 4, shape = 18, show.legend = FALSE) +
        geom_errorbar(aes(ymin = mean_val - err, ymax = mean_val + err),
                      width = 0.25, linewidth = 0.6)

    } else if (bp_type == "box") {
      p <- ggplot(plot_df, aes(x = sample, y = value, fill = sample)) +
        geom_boxplot(width = bp_bar_width, alpha = bp_bar_alpha, color = "black",
                     linewidth = 0.3, outlier.shape = NA) +
        geom_jitter(width = 0.15, size = 2, shape = 21, fill = "white",
                    color = "black", stroke = 0.5, show.legend = FALSE)

    } else if (bp_type == "violin") {
      p <- ggplot(plot_df, aes(x = sample, y = value, fill = sample)) +
        geom_violin(width = bp_bar_width, alpha = bp_bar_alpha, color = "black",
                    linewidth = 0.3, trim = FALSE) +
        geom_jitter(width = 0.15, size = 2, shape = 21, fill = "white",
                    color = "black", stroke = 0.5, show.legend = FALSE) +
        stat_summary(fun = mean, geom = "crossbar", width = 0.3,
                     linewidth = 0.4, color = "black", show.legend = FALSE)
    }

    # Fixed sizing: keep empty factor levels visible
    if (sizing_mode == "fixed")
      p <- p + scale_x_discrete(drop = FALSE)

    # Significance brackets (only for non-horizontal bar/dot types for clean rendering)
    sig_text_sz <- if (!is.null(input$barplot_sig_text_size)) input$barplot_sig_text_size else 4.5
    if (isTRUE(input$barplot_show_sig) && !is.null(rv_analysis$stats_result) && !isTRUE(input$barplot_horizontal)) {
      sr <- rv_analysis$stats_result
      sr <- sr[sr$metric == met & sr$significance %in% c("*","**","***") &
                 !grepl("Overall", sr$comparison), , drop = FALSE]
      if (nrow(sr) > 0) {
        sr <- sr %>% group_by(comparison) %>% slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>% ungroup()
        if (bp_type %in% c("bar", "dot")) {
          y_max  <- max(summary_df$mean_val + summary_df$err, na.rm = TRUE)
          y_min  <- min(summary_df$mean_val - summary_df$err, na.rm = TRUE)
        } else {
          y_max  <- max(plot_df$value, na.rm = TRUE)
          y_min  <- min(plot_df$value, na.rm = TRUE)
        }
        y_span <- max(y_max - y_min, abs(y_max), 1)
        y_step <- y_span * 0.08
        for (j in seq_len(min(nrow(sr), 5))) {
          pair <- strsplit(sr$comparison[j], " vs ")[[1]]
          if (length(pair) == 2 && all(pair %in% levels(summary_df$sample))) {
            x1 <- which(levels(summary_df$sample) == pair[1])
            x2 <- which(levels(summary_df$sample) == pair[2])
            yb <- y_max + y_step * j
            p <- p +
              annotate("segment", x=x1, xend=x2, y=yb, yend=yb, linewidth=0.4) +
              annotate("segment", x=x1, xend=x1, y=yb, yend=yb - y_step*0.3, linewidth=0.4) +
              annotate("segment", x=x2, xend=x2, y=yb, yend=yb - y_step*0.3, linewidth=0.4) +
              annotate("text", x=(x1+x2)/2, y=yb + y_step*0.15,
                       label=sr$significance[j], size=sig_text_sz, fontface="bold")
          }
        }
      }
    }

    bp_title    <- if (!is.null(input$barplot_title) && nchar(trimws(input$barplot_title)) > 0) input$barplot_title else y_label
    bp_x_angle  <- if (!is.null(input$barplot_x_angle)) as.numeric(input$barplot_x_angle) else 45
    bp_title_sz <- if (!is.null(input$barplot_title_size)) input$barplot_title_size else 16
    bp_ax_sz    <- if (!is.null(input$barplot_axis_text_size)) input$barplot_axis_text_size else 12
    bp_ax_lab   <- if (!is.null(input$barplot_axis_label_size)) input$barplot_axis_label_size else 14
    bp_leg      <- if (!is.null(input$barplot_legend_pos)) input$barplot_legend_pos else "none"
    bp_leg_sz   <- if (!is.null(input$barplot_legend_size)) input$barplot_legend_size else 12
    bp_leg_face <- if (!is.null(input$barplot_legend_face)) input$barplot_legend_face else "plain"

    p <- p + labs(x = NULL, y = paste0(y_label, err_label), title = bp_title, fill = NULL) +
      theme_pubr() +
      theme(text = element_text(family = input$font_family),
            plot.title = element_text(size = bp_title_sz, face = "bold", hjust = 0.5),
            axis.title = element_text(size = bp_ax_lab),
            axis.text  = element_text(size = bp_ax_sz),
            axis.text.x = element_text(angle = bp_x_angle, hjust = if (bp_x_angle > 0) 1 else 0.5),
            legend.position = bp_leg,
            legend.text = element_text(size = bp_leg_sz, face = bp_leg_face),
            legend.title = element_text(size = bp_leg_sz + 1, face = bp_leg_face),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

    # Apply colors
    if (!is.null(bar_colors) && length(bar_colors) > 0) {
      if (bp_type == "dot") {
        p <- p + scale_color_manual(values = bar_colors)
      } else {
        p <- p + scale_fill_manual(values = bar_colors)
      }
    }

    # Horizontal bars (coord_flip)
    if (isTRUE(input$barplot_horizontal))
      p <- p + coord_flip()

    p
  }
  
  output$metrics_barplot <- renderPlot({ build_barplot() }, res = 96)
  
  output$download_barplot <- downloadHandler(
    filename = function() paste0("barplot_", input$barplot_metric, "_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$barplot_export_fmt),
    content = function(file) {
      p <- build_barplot(); req(p)
      ggsave(file, plot = p, device = input$barplot_export_fmt,
             width = input$barplot_export_w, height = input$barplot_export_h,
             units = "in", dpi = input$export_dpi)
    }
  )
  
  # ── Derivative Plot ──────────────────────────────────────────────────────────
  output$deriv_plot_container <- renderUI({
    req(rv_analysis$deriv_data)
    plotOutput("deriv_plot", height = paste0(input$deriv_height, "px"),
               width = paste0(input$deriv_width, "px"))
  })
  
  build_deriv_plot <- function() {
    req(rv_analysis$deriv_data)
    dd <- rv_analysis$deriv_data
    samps <- input$selected_samples
    if (!is.null(samps)) dd <- dd %>% filter(variable %in% samps)
    if (nrow(dd) == 0) return(NULL)
    line_colors <- if (!is.null(samps)) resolve_aesthetics(samps)$colors else NULL
    p <- ggplot(dd, aes(x = time, y = dODdt, color = variable)) +
      geom_line(linewidth = 0.8) +
      labs(x = input$x_axis_label, y = "dOD/dt (rate of change)",
           title = "Derivative Plot: Rate of OD Change", color = NULL) +
      theme_pubr() +
      theme(text = element_text(family = input$font_family),
            plot.title  = element_text(size = if (!is.null(input$deriv_title_size)) input$deriv_title_size else 16,
                                       face = if (!is.null(input$deriv_title_face)) input$deriv_title_face else "bold",
                                       hjust = 0.5),
            axis.title  = element_text(size = if (!is.null(input$deriv_axis_label_size)) input$deriv_axis_label_size else 14,
                                       face = if (!is.null(input$deriv_axis_face)) input$deriv_axis_face else "plain"),
            axis.text   = element_text(size = if (!is.null(input$deriv_axis_text_size)) input$deriv_axis_text_size else 12),
            axis.text.x = element_text(angle = as.numeric(if (!is.null(input$deriv_x_angle)) input$deriv_x_angle else "0"),
                                       hjust = if (!is.null(input$deriv_x_angle) && input$deriv_x_angle != "0") 1 else 0.5),
            legend.position = if (!is.null(input$deriv_legend_pos)) input$deriv_legend_pos else "right",
            legend.text  = element_text(size = if (!is.null(input$deriv_legend_size)) input$deriv_legend_size else 12),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
    if (isTRUE(input$deriv_show_zero))
      p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5)
    if (!is.null(line_colors)) p <- p + scale_color_manual(values = line_colors)
    p
  }
  
  output$deriv_plot <- renderPlot({ build_deriv_plot() }, res = 96)
  output$download_deriv_plot <- downloadHandler(
    filename = function() paste0("derivative_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$deriv_export_fmt),
    content = function(file) {
      p <- build_deriv_plot(); req(p)
      ggsave(file, plot = p, device = input$deriv_export_fmt,
             width = input$deriv_export_w, height = input$deriv_export_h, units = "in", dpi = input$export_dpi)
    }
  )
  
  # ── Annotated Growth Curves ──────────────────────────────────────────────────
  output$annot_plot_container <- renderUI({
    req(rv_analysis$metrics)
    plotOutput("annot_plot", height = paste0(input$annot_height, "px"),
               width = paste0(input$annot_width, "px"))
  })
  
  build_annot_plot <- function() {
    req(rv$data, input$selected_samples, rv_analysis$metrics)
    pd <- tryCatch(prepare_plot_data(), error = function(e) NULL)
    if (is.null(pd) || nrow(pd) == 0) return(NULL)
    m <- rv_analysis$metrics; aes_v <- resolve_aesthetics(input$selected_samples)
    pd$variable <- factor(pd$variable, levels = input$selected_samples)
    p <- ggplot(pd, aes(x = time, y = mean_value, color = variable)) +
      geom_line(linewidth = 0.9) +
      labs(x = input$x_axis_label, y = input$y_axis_label,
           title = "Annotated Growth Curves", color = NULL) +
      scale_color_manual(values = aes_v$colors) +
      theme_pubr() +
      theme(text = element_text(family = input$font_family),
            plot.title  = element_text(size = if (!is.null(input$annot_title_size)) input$annot_title_size else 16,
                                       face = if (!is.null(input$annot_title_face)) input$annot_title_face else "bold",
                                       hjust = 0.5),
            axis.title  = element_text(size = if (!is.null(input$annot_axis_label_size)) input$annot_axis_label_size else 14,
                                       face = if (!is.null(input$annot_axis_face)) input$annot_axis_face else "plain"),
            axis.text   = element_text(size = if (!is.null(input$annot_axis_text_size)) input$annot_axis_text_size else 12),
            axis.text.x = element_text(angle = as.numeric(if (!is.null(input$annot_x_angle)) input$annot_x_angle else "0"),
                                       hjust = if (!is.null(input$annot_x_angle) && input$annot_x_angle != "0") 1 else 0.5),
            legend.position = if (!is.null(input$annot_legend_pos)) input$annot_legend_pos else "right",
            legend.text  = element_text(size = if (!is.null(input$annot_legend_size)) input$annot_legend_size else 12),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
    for (i in seq_len(nrow(m))) {
      s <- m$sample[i]; sc <- pd %>% filter(variable == s)
      if (nrow(sc) == 0) next
      if (isTRUE(input$annot_lag) && is.finite(m$lag_phase[i])) {
        lag_od <- approx(sc$time, sc$mean_value, xout = m$lag_phase[i])$y
        if (!is.null(lag_od) && is.finite(lag_od))
          p <- p + annotate("point", x = m$lag_phase[i], y = lag_od, shape = 17, size = 3.5, color = "#27AE60")
      }
      if (isTRUE(input$annot_mumax) && is.finite(m$time_max_od[i]) && is.finite(m$max_od[i]))
        p <- p + annotate("point", x = m$time_max_od[i], y = m$max_od[i], shape = 18, size = 4, color = "#2980B9")
      if (isTRUE(input$annot_lysis) && is.finite(m$lysis_time[i])) {
        lys_od <- approx(sc$time, sc$mean_value, xout = m$lysis_time[i])$y
        if (!is.null(lys_od) && is.finite(lys_od))
          p <- p + annotate("point", x = m$lysis_time[i], y = lys_od, shape = 25, size = 3.5, color = "#E74C3C", fill = "#E74C3C")
      }
    }
    p
  }
  
  output$annot_plot <- renderPlot({ build_annot_plot() }, res = 96)
  output$download_annot_plot <- downloadHandler(
    filename = function() paste0("annotated_curves_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$annot_export_fmt),
    content = function(file) {
      p <- build_annot_plot(); req(p)
      ggsave(file, plot = p, device = input$annot_export_fmt,
             width = input$annot_export_w, height = input$annot_export_h, units = "in", dpi = input$export_dpi)
    }
  )
  
  # ── Heatmap Sample Selectors ─────────────────────────────────────────────────
  output$heatmap_sample_selector <- renderUI({
    req(rv_analysis$metrics)
    samps <- unique(rv_analysis$metrics$sample)
    checkboxGroupInput("heatmap_samples", "Show samples:",
                       choices = samps, selected = samps, inline = TRUE)
  })
  observeEvent(input$heatmap_select_all, {
    req(rv_analysis$metrics)
    updateCheckboxGroupInput(session, "heatmap_samples",
                             selected = unique(rv_analysis$metrics$sample))
  })
  observeEvent(input$heatmap_deselect_all, {
    updateCheckboxGroupInput(session, "heatmap_samples", selected = character(0))
  })

  output$od_heatmap_sample_selector <- renderUI({
    samps <- input$selected_samples
    if (is.null(samps) && !is.null(rv$od_vars)) samps <- rv$od_vars
    req(samps)
    checkboxGroupInput("od_heatmap_samples", "Show samples:",
                       choices = samps, selected = samps, inline = TRUE)
  })
  observeEvent(input$od_heatmap_select_all, {
    samps <- if (!is.null(input$selected_samples)) input$selected_samples else isolate(rv$od_vars)
    req(samps)
    updateCheckboxGroupInput(session, "od_heatmap_samples", selected = samps)
  })
  observeEvent(input$od_heatmap_deselect_all, {
    updateCheckboxGroupInput(session, "od_heatmap_samples", selected = character(0))
  })

  # ── Phenotype Heatmap (with metric selection) ────────────────────────────────
  output$heatmap_metric_checkboxes <- renderUI({
    req(rv_analysis$metrics)
    m <- if (!is.null(rv_analysis$metrics_inf)) rv_analysis$metrics_inf else rv_analysis$metrics
    num_cols <- names(m)[sapply(m, function(x) is.numeric(x) && any(is.finite(x)))]
    num_cols <- setdiff(num_cols, "sample")
    labels <- ifelse(num_cols %in% names(metric_labels), metric_labels[num_cols], num_cols)
    checkboxGroupInput("heatmap_metrics", "Select Metrics:", choices = setNames(num_cols, labels),
                       selected = num_cols, inline = TRUE)
  })
  
  output$heatmap_container <- renderUI({
    req(rv_analysis$metrics)
    plotOutput("heatmap_plot", height = paste0(input$heatmap_height, "px"),
               width = paste0(input$heatmap_width, "px"))
  })
  
  build_heatmap <- function() {
    req(rv_analysis$metrics, input$heatmap_metrics)
    m <- if (!is.null(rv_analysis$metrics_inf)) rv_analysis$metrics_inf else rv_analysis$metrics
    sel_cols <- intersect(input$heatmap_metrics, names(m))
    if (length(sel_cols) < 1) return(NULL)
    heat_df <- m %>% select(sample, all_of(sel_cols)) %>%
      pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
      filter(is.finite(value)) %>%
      group_by(metric) %>%
      mutate(z = if (sd(value, na.rm = TRUE) > 0) (value - mean(value)) / sd(value) else 0) %>%
      ungroup()
    heat_df$metric_label <- ifelse(heat_df$metric %in% names(metric_labels),
                                   metric_labels[heat_df$metric], heat_df$metric)
    samp_order <- if (!is.null(input$heatmap_samples) && length(input$heatmap_samples) > 0)
      intersect(input$heatmap_samples, unique(heat_df$sample))
    else if (!is.null(input$selected_samples))
      intersect(input$selected_samples, unique(heat_df$sample))
    else unique(heat_df$sample)
    heat_df <- heat_df %>% filter(sample %in% samp_order)
    if (nrow(heat_df) == 0) return(NULL)
    heat_df$sample <- factor(heat_df$sample, levels = rev(samp_order))
    ggplot(heat_df, aes(x = metric_label, y = sample, fill = z)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = round(value, 2)), size = 3, color = "black") +
      scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#E74C3C", midpoint = 0, name = "Z-score") +
      labs(x = NULL, y = NULL, title = "Phenotype Heatmap") +
      theme_minimal() +
      theme(text = element_text(family = input$font_family),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 11), panel.grid = element_blank(), legend.position = "right")
  }
  
  output$heatmap_plot <- renderPlot({ build_heatmap() }, res = 96)
  output$download_heatmap <- downloadHandler(
    filename = function() paste0("heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$heatmap_export_fmt),
    content = function(file) {
      p <- build_heatmap(); req(p)
      ggsave(file, plot = p, device = input$heatmap_export_fmt,
             width = input$heatmap_export_w, height = input$heatmap_export_h, units = "in", dpi = input$export_dpi)
    }
  )
  
  # ── OD Over Time Heatmap ─────────────────────────────────────────────────────
  output$od_heatmap_container <- renderUI({
    req(rv$data, input$selected_samples)
    plotOutput("od_heatmap_plot", height = paste0(input$od_heatmap_height, "px"),
               width = paste0(input$od_heatmap_width, "px"))
  })
  
  build_od_heatmap <- function() {
    req(rv$data, input$selected_samples)
    od_samps <- if (!is.null(input$od_heatmap_samples) && length(input$od_heatmap_samples) > 0)
      input$od_heatmap_samples else input$selected_samples
    pd <- tryCatch(prepare_plot_data(samples = od_samps), error = function(e) NULL)
    if (is.null(pd) || nrow(pd) == 0) return(NULL)
    pd$variable <- factor(pd$variable, levels = rev(intersect(od_samps, unique(pd$variable))))
    pal <- if (!is.null(input$od_heatmap_palette)) input$od_heatmap_palette else "inferno"
    ggplot(pd, aes(x = time, y = variable, fill = mean_value)) +
      geom_tile() +
      scale_fill_viridis_c(option = pal, name = "OD") +
      labs(x = input$x_axis_label, y = NULL, title = "OD Over Time Heatmap") +
      theme_minimal() +
      theme(text = element_text(family = input$font_family),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 10),
            panel.grid = element_blank(), legend.position = "right")
  }
  
  output$od_heatmap_plot <- renderPlot({ build_od_heatmap() }, res = 96)
  output$download_od_heatmap <- downloadHandler(
    filename = function() paste0("od_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", input$od_heatmap_export_fmt),
    content = function(file) {
      p <- build_od_heatmap(); req(p)
      ggsave(file, plot = p, device = input$od_heatmap_export_fmt,
             width = input$od_heatmap_export_w, height = input$od_heatmap_export_h, units = "in", dpi = input$export_dpi)
    }
  )
  
  # ── Statistical Tests ────────────────────────────────────────────────────────
  observeEvent(input$run_stats, {
    req(rv_analysis$metrics_long, input$stats_metric)
    ml <- rv_analysis$metrics_long
    if (!is.null(input$stats_samples) && length(input$stats_samples) > 0)
      ml <- ml %>% filter(sample %in% input$stats_samples)
    res <- tryCatch(
      compare_metrics(ml, input$stats_metric),
      error = function(e) { showNotification(paste("Stats error:", e$message), type = "error"); NULL })
    if (is.null(res)) {
      showNotification("Not enough replicate values for this metric to run stats.", type = "warning")
      rv_analysis$stats_result <- NULL
    } else {
      rv_analysis$stats_result <- res
    }
  })
  
  output$stats_table <- DT::renderDT({
    req(rv_analysis$stats_result)
    s <- rv_analysis$stats_result; s$p_value <- signif(s$p_value, 4)
    datatable(s, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 25, dom = "Bfrtip"),
              class = "compact stripe hover") %>%
      formatStyle("significance", backgroundColor = styleEqual(
        c("***","**","*",".","ns"), c("#d4edda","#d4edda","#fff3cd","#f8f9fa","#f8f9fa")))
  })
  
  # ── Download handlers ────────────────────────────────────────────────────────
  output$download_metrics_csv <- downloadHandler(
    filename = function() paste0("metrics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      m <- if (!is.null(rv_analysis$metrics_inf)) rv_analysis$metrics_inf else rv_analysis$metrics
      req(m); write.csv(m, file, row.names = FALSE)
    }
  )
  output$download_infection_csv <- downloadHandler(
    filename = function() paste0("infection_metrics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) { req(rv_analysis$metrics_inf); write.csv(rv_analysis$metrics_inf, file, row.names = FALSE) }
  )
  output$download_stats_csv <- downloadHandler(
    filename = function() paste0("stats_", input$stats_metric, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) { req(rv_analysis$stats_result); write.csv(rv_analysis$stats_result, file, row.names = FALSE) }
  )
  output$download_all_stats <- downloadHandler(
    filename = function() paste0("all_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      req(rv_analysis$metrics_long)
      all_mets <- unique(rv_analysis$metrics_long$metric)
      all_res <- lapply(all_mets, function(met) tryCatch(compare_metrics(rv_analysis$metrics_long, met), error = function(e) NULL))
      combined <- bind_rows(all_res[!sapply(all_res, is.null)])
      if (nrow(combined) > 0) { combined$p_adj <- p.adjust(combined$p_value, method = "BH"); write.csv(combined, file, row.names = FALSE) }
      else write.csv(tibble(note = "No valid comparisons"), file, row.names = FALSE)
    }
  )
  
}
shinyApp(ui, server)
