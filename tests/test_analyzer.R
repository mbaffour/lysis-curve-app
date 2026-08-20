# ============================================================================
# Verification suite for the merged "OD Growth Curve Analyzer"
# (lysis-curve-app + killcurveplot v2.8 fixes). Drives the REAL server code
# via shiny::testServer.
# ============================================================================
suppressPackageStartupMessages({
  library(shiny); library(tidyverse); library(ggpubr); library(scales)
  library(ggrepel); library(ggprism); library(svglite); library(jsonlite)
  library(zoo); library(DT)
})
if (!file.exists("Lysis Curve Ap 26.03.13.R")) stop("Run from the repository root: Rscript tests/test_analyzer.R")

pass <- 0L; fail <- 0L
check <- function(desc, cond) {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat("PASS:", desc, "\n") }
  else              { fail <<- fail + 1L; cat("FAIL:", desc, "\n") }
}
near <- function(a, b, tol = 1e-9) is.finite(a) && abs(a - b) < tol

env <- new.env()
app <- source("Lysis Curve Ap 26.03.13.R", local = env)$value
check("app sources -> shiny.appobj", inherits(app, "shiny.appobj"))

cat("\n--- helpers ---\n")
check("rgb_to_hex no alpha suffix", toupper(env$rgb_to_hex(255, 0, 0)) == "#FF0000")
dd <- env$make_demo_data()
check("demo: 39 rows (3 stacked blocks x 13 times)", nrow(dd) == 39)
check("demo: 3 sample columns", ncol(dd) == 4)

cat("\n--- infection metrics fix (ground truth) ---\n")
mk <- function(t, y, name) data.frame(time = t, variable = name, mean_value = y)
multi <- rbind(
  mk(seq(0, 180, 15), c(0.05,0.08,0.13,0.22,0.35,0.55,0.80,1.00,1.15,1.22,1.26,1.28,1.29), "Uninfected"),
  mk(seq(0, 180, 15), c(0.05,0.08,0.13,0.22,0.35,0.30,0.15,0.08,0.05,0.04,0.04,0.05,0.06), "Infected"))
gm <- as.data.frame(env$calculate_growth_metrics(multi, smooth_window = 5))
im <- as.data.frame(env$calculate_infection_metrics(gm, "Uninfected"))
u  <- im[im$sample == "Uninfected", ]; inf <- im[im$sample == "Infected", ]
check("ref vs itself: relative_growth = 1",    near(u$relative_growth, 1))
check("ref vs itself: infection_strength = 0", near(u$infection_strength, 0))
check("ref vs itself: relative_max_od = 1",    near(u$relative_max_od, 1))
check("infected: relative_growth = auc ratio", near(inf$relative_growth, inf$auc / u$auc))

flat <- mk(seq(0, 60, 10), rep(0.3, 7), "F")
fm <- as.data.frame(env$calculate_growth_metrics(flat, smooth_window = 5))
check("flat curve: mu_max NA (epsilon guard)", is.na(fm$mu_max) && is.na(fm$doubling_time))

cat("\n--- server pipelines via testServer ---\n")
testServer(app, {
  # 1. stacked-block replicates (this app's documented convention)
  session$setInputs(metrics_smooth_window = 5)
  apply_data_to_rv(env$make_demo_data(), "demo")
  session$setInputs(selected_samples = rv$od_vars, error_type = "sem",
                    exclude_timepoints = "")
  pd <- prepare_plot_data()
  check("stacked: 3 samples x 13 times = 39 rows", nrow(pd) == 39)
  check("stacked: n = 3 everywhere", all(pd$n == 3))

  md <- prepare_metrics_data()
  check("stacked: metrics data has replicate column", "replicate" %in% names(md))
  check("stacked: 3 replicates per sample",
        length(unique(md$replicate[md$variable == md$variable[1]])) == 3)

  rd <- prepare_replicate_data()
  check("stacked: replicate data 3x39 rows", nrow(rd) == 117)

  # 2. duplicate-column replicates (killcurveplot convention)
  tmp <- tempfile(fileext = ".csv")
  writeLines(c("Time,WT,WT,WT,Mutant,Mutant,Mutant",
               "0,0.05,0.06,0.04,0.05,0.04,0.06",
               "30,0.12,0.15,0.10,0.08,0.09,0.07"), tmp)
  apply_data_to_rv(read.csv(tmp, stringsAsFactors = FALSE, check.names = FALSE), "dup.csv")
  session$setInputs(selected_samples = c("WT", "Mutant"))
  check("dup cols: col_map WT -> 3 columns", length(rv$col_map[["WT"]]) == 3)
  pd2 <- prepare_plot_data()
  wt0 <- pd2[pd2$variable == "WT" & pd2$time == 0, ]
  check("dup cols: WT t=0 n = 3",        wt0$n == 3)
  check("dup cols: WT t=0 mean = 0.05",  near(wt0$mean_value, 0.05))
  check("dup cols: WT t=0 sd = 0.01",    near(wt0$sd_value, 0.01, 1e-12))
  check("dup cols: sem = sd/sqrt(3)",    near(wt0$sem_value, 0.01/sqrt(3), 1e-12))
  rd2 <- prepare_replicate_data()
  check("dup cols: replicate data has 3 reps per sample",
        length(unique(rd2$replicate)) == 3)

  # 3. NA replicate: n counts only non-NA
  tmp3 <- tempfile(fileext = ".csv")
  writeLines(c("Time,A,A,A", "0,0.1,NA,0.3", "10,0.2,0.4,0.6"), tmp3)
  apply_data_to_rv(read.csv(tmp3, stringsAsFactors = FALSE, check.names = FALSE), "na.csv")
  session$setInputs(selected_samples = "A")
  a0 <- prepare_plot_data(); a0 <- a0[a0$time == 0, ]
  check("NA rep: n = 2", a0$n == 2)
  check("NA rep: sem uses n = 2", near(a0$sem_value, sd(c(0.1, 0.3))/sqrt(2), 1e-12))

  # 4. REGRESSION: time filter + stacked replicates must not empty the metrics
  apply_data_to_rv(env$make_demo_data(), "demo")
  session$setInputs(selected_samples = rv$od_vars,
                    enable_time_filter = TRUE,
                    time_filter_range = c(0, 120),
                    exclude_timepoints = "")
  mdf <- tryCatch(prepare_metrics_data(), error = function(e) NULL)
  check("time filter + stacked reps: metrics data non-empty (was silently empty)",
        !is.null(mdf) && nrow(mdf) > 0)
  check("time filter: range respected", !is.null(mdf) && max(mdf$time) <= 120)
  check("time filter: replicate structure survives",
        !is.null(mdf) && length(unique(mdf$replicate)) == 3)

  # 5. long format with explicit replicate column
  tmp5 <- tempfile(fileext = ".csv")
  writeLines(c("Time,Sample,replicate,OD",
               "0,WT,1,0.05", "0,WT,2,0.06", "0,WT,3,NA",
               "30,WT,1,0.12", "30,WT,2,0.15", "30,WT,3,0.10"), tmp5)
  apply_data_to_rv(read.csv(tmp5, stringsAsFactors = FALSE, check.names = FALSE), "long.csv")
  session$setInputs(selected_samples = "WT", enable_time_filter = FALSE)
  pl <- prepare_plot_data(); l0 <- pl[pl$time == 0, ]
  check("long + NA: n = 2, mean = 0.055", l0$n == 2 && near(l0$mean_value, 0.055))
})

cat("\n============ SUMMARY:", pass, "passed /", fail, "failed ============\n")
if (fail > 0) quit(status = 1)
