# ============================================================================
# STRESS TEST — OD Growth Curve Analyzer
# Pathological inputs, scale, and crash-safety across ingestion, stats,
# metrics, QC, and report generation. A "PASS" means either correct output
# or a *graceful* refusal (req()/validate/NULL) — never an uncaught crash.
# ============================================================================
suppressPackageStartupMessages({
  library(shiny); library(tidyverse); library(ggpubr); library(scales)
  library(ggrepel); library(ggprism); library(svglite); library(jsonlite)
  library(zoo); library(DT)
})
if (!file.exists("Lysis Curve Ap 26.03.13.R")) stop("Run from the repository root: Rscript tests/stress_test_analyzer.R")

pass <- 0L; fail <- 0L
check <- function(desc, cond) {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat("PASS:", desc, "\n") }
  else              { fail <<- fail + 1L; cat("FAIL:", desc, "\n") }
}
# graceful = TRUE result, or a shiny req() halt; hard error = fail
graceful <- function(expr) {
  tryCatch({ force(expr); TRUE },
           shiny.silent.error = function(e) TRUE,   # req() refusal is graceful
           validation         = function(e) TRUE,
           error = function(e) { cat("   HARD ERROR:", conditionMessage(e), "\n"); FALSE })
}
csv <- function(...) { f <- tempfile(fileext = ".csv"); writeLines(c(...), f); f }
rd  <- function(f) tryCatch(read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
                            error = function(e) read.csv(f, stringsAsFactors = FALSE))

env <- new.env()
app <- source("Lysis Curve Ap 26.03.13.R", local = env)$value
check("app sources", inherits(app, "shiny.appobj"))

# ── Pure-function stress (no server) ─────────────────────────────────────────
cat("\n--- metrics engine edge cases ---\n")
mk <- function(t, y, nm = "S") data.frame(time = t, variable = nm, mean_value = y)

check("1-point curve: metrics don't crash",
      graceful(env$calculate_growth_metrics(mk(0, 0.5), smooth_window = 5)))
check("2-point curve: metrics don't crash",
      graceful(env$calculate_growth_metrics(mk(c(0, 10), c(0.5, 0.4)), smooth_window = 5)))
m2 <- env$calculate_growth_metrics(mk(c(0, 10), c(0.5, 0.4)), smooth_window = 5)
check("2-point curve: mu_max NA (window needs >= 3)", is.na(m2$mu_max))
empty_pd <- data.frame(time = numeric(0), variable = character(0),
                       mean_value = numeric(0))
check("NULL / empty input returns NULL",
      is.null(env$calculate_growth_metrics(NULL)) &&
      is.null(env$calculate_growth_metrics(empty_pd)))
check("all-identical times: no crash",
      graceful(env$calculate_growth_metrics(mk(c(5, 5, 5), c(0.1, 0.2, 0.3)), 5)))
check("huge values (1e6): finite metrics",
      graceful({ m <- env$calculate_growth_metrics(mk(seq(0, 100, 10), seq(1, 11) * 1e5), 5)
                 stopifnot(is.finite(m$auc)) }))
check("tiny values (1e-9): no crash",
      graceful(env$calculate_growth_metrics(mk(seq(0, 100, 10), rep(1e-9, 11)), 5)))
check("negative times: no crash",
      graceful(env$calculate_growth_metrics(mk(seq(-60, 60, 15), runif(9, 0.1, 1), "N"), 5)))
check("window larger than n: no crash",
      graceful(env$calculate_growth_metrics(mk(seq(0, 30, 10), c(0.1, 0.2, 0.3, 0.4)), 15)))

cat("\n--- local virulence edge cases ---\n")
lc <- mk(seq(0, 120, 10), c(0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.2, 1.3, 1.33, 1.34, 1.34, 1.34, 1.34), "Ref")
tr <- mk(seq(0, 120, 10), lc$mean_value * 0.5, "Trt")
both <- rbind(lc, tr)
check("local virulence: reference missing -> graceful NULL/empty",
      graceful({ r <- env$calculate_local_virulence(both, "NotThere", 5)
                 stopifnot(is.null(r) || nrow(r) == 0 || all(is.na(r$local_virulence))) }))
check("local virulence: flat reference -> graceful (no div-by-zero)",
      graceful(env$calculate_local_virulence(
        rbind(mk(seq(0, 60, 10), rep(0, 7), "Z"), mk(seq(0, 60, 10), runif(7), "T")), "Z", 5)))
check("local virulence: t_stat override beyond run -> graceful",
      graceful(env$calculate_local_virulence(both, "Ref", 5, t_stat_override = 1e6)))
check("local virulence: single timepoint -> graceful",
      graceful(env$calculate_local_virulence(mk(0, 0.5, "Ref"), "Ref", 5)))

cat("\n--- qc_flags edge cases ---\n")
check("qc_flags: NULL input graceful", graceful(env$qc_flags(NULL)))
check("qc_flags: zero-row input graceful",
      graceful(env$qc_flags(empty_pd)))
check("qc_flags: single point graceful", graceful(env$qc_flags(mk(0, 2, "S"))))

cat("\n--- report builders under stress ---\n")
big_stats <- data.frame(time = rep(seq(0, 990, 10), 50),
                        sample = rep(paste0("S", 1:50), each = 100),
                        mean = runif(5000), sd = runif(5000, 0, 0.1), n = 3L,
                        sem = runif(5000, 0, 0.05), ci95_half_width = runif(5000, 0, 0.1))
big_metrics <- env$calculate_growth_metrics(
  do.call(rbind, lapply(1:50, function(i)
    mk(seq(0, 180, 15), pmax(cumsum(rnorm(13, 0.05, 0.1)), 0.01), paste0("S", i)))), 5)
meta <- list(title = "Stress <script>alert(1)</script> & Co", author = "A&B <i>",
             notes = paste(rep("long note line", 200), collapse = "\n"),
             generated = "now", data_file = "x.csv", format = "wide", time_col = "T",
             samples_txt = paste(paste0("S", 1:50), collapse = ", "),
             filter_txt = "none", error_txt = "SEM",
             settings_df = data.frame(setting = "X", value = "1"),
             notebook_df = data.frame(field = "Obs", value = strrep("verylongword", 50)),
             qc_df = data.frame(sample = "S1", flag = "CEILING", severity = "warn",
                                detail = "x > 1"),
             font_pt = 14)
html <- NULL
check("HTML report: 50 samples x 100 times builds",
      graceful({ html <<- env$build_html_report(meta, "AAAA", big_stats, big_metrics) }))
check("HTML report: script tag escaped (no XSS)",
      !grepl("<script>alert", html, fixed = TRUE) &&
       grepl("&lt;script&gt;", html, fixed = TRUE))
pdff <- tempfile(fileext = ".pdf")
check("PDF report: 5000-row stats table paginates without error",
      graceful({ p <- ggplot(mk(0:10, runif(11)), aes(time, mean_value)) + geom_line()
                 env$build_pdf_report(pdff, meta, p, big_stats, big_metrics) }))
check("PDF: >30 KB and valid magic",
      file.exists(pdff) && file.info(pdff)$size > 30000 &&
      identical(readBin(pdff, "raw", 4), charToRaw("%PDF")))

# ── Server-level stress via testServer ───────────────────────────────────────
cat("\n--- server ingestion stress ---\n")
testServer(app, {
  ok_ingest <- function(desc, df) {
    check(paste0("ingest: ", desc),
          graceful({ apply_data_to_rv(df, desc); TRUE }))
  }
  ok_ingest("header-only (0 rows)", rd(csv("Time,A,B")))
  ok_ingest("single timepoint",     rd(csv("Time,A,B", "0,0.1,0.2")))
  ok_ingest("single sample col",    rd(csv("Time,A", "0,0.1", "10,0.2")))
  ok_ingest("all-NA sample",        rd(csv("Time,A,B", "0,NA,0.2", "10,NA,0.3")))
  ok_ingest("junk text column",     rd(csv("Time,A,Note", "0,0.1,hello", "10,0.2,world")))
  ok_ingest("negative times",       rd(csv("Time,A", "-30,0.1", "0,0.2", "30,0.3")))
  ok_ingest("duplicate times",      rd(csv("Time,A", "0,0.1", "0,0.12", "10,0.2")))
  ok_ingest("unicode names",        rd(csv("Time,\u0394cI phage,WT\u00b5", "0,0.1,0.2", "10,0.2,0.3")))
  ok_ingest("blank column name",    rd(csv("Time,A,", "0,0.1,9", "10,0.2,9")))
  ok_ingest("unequal stacked blocks", rd(csv("Time,A", "0,0.1", "10,0.2", "20,0.3", "0,0.11", "10,0.21")))

  # correctness after the weird ones: reload clean demo and verify
  apply_data_to_rv(env$make_demo_data(), "demo")
  session$setInputs(selected_samples = rv$od_vars, error_type = "sem",
                    exclude_timepoints = "", metrics_smooth_window = 5,
                    extinction_threshold = 0.05)
  pd <- prepare_plot_data()
  check("recovery after stress: demo pipeline intact (39 rows, n=3)",
        nrow(pd) == 39 && all(pd$n == 3))

  # mixed convention: dup columns AND stacked blocks simultaneously
  mixed <- rd(csv("Time,A,A,B", "0,0.1,0.11,0.5", "10,0.2,0.21,0.4",
                  "0,0.12,0.13,0.52", "10,0.22,0.23,0.42"))
  apply_data_to_rv(mixed, "mixed")
  session$setInputs(selected_samples = c("A", "B"))
  pdm <- prepare_plot_data()
  a0 <- pdm[pdm$variable == "A" & pdm$time == 0, ]
  check("mixed conventions: A pools 2 cols x 2 blocks = n 4",
        nrow(a0) == 1 && a0$n == 4)
  check("mixed conventions: mean over all 4 replicates",
        abs(a0$mean_value - mean(c(0.1, 0.11, 0.12, 0.13))) < 1e-12)
  rdm <- prepare_replicate_data()
  check("mixed conventions: 4 distinct replicate ids for A",
        length(unique(rdm$replicate[rdm$variable == "A"])) == 4)

  # scale: 100 samples x 100 timepoints
  nT <- 100; nS <- 100
  bigw <- data.frame(Time = seq(0, 990, 10))
  for (i in seq_len(nS)) bigw[[paste0("S", i)]] <- pmax(cumsum(rnorm(nT, 0.02, 0.05)), 0.01)
  t0 <- Sys.time()
  apply_data_to_rv(bigw, "big")
  session$setInputs(selected_samples = rv$od_vars)
  pdb <- prepare_plot_data()
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  check("scale: 100x100 pipeline < 10 s", nrow(pdb) == nT * nS && el < 10)
  cat(sprintf("   (100 samples x 100 timepoints: %.2f s)\n", el))
  t1 <- Sys.time()
  mb <- env$calculate_growth_metrics(prepare_metrics_data(), smooth_window = 5)
  el2 <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  check("scale: metrics for 100 samples < 60 s", nrow(mb) >= nS && el2 < 60)
  cat(sprintf("   (metrics on 100 samples: %.2f s)\n", el2))

  # time filter interactions on the big data
  session$setInputs(enable_time_filter = TRUE, time_filter_range = c(100, 500),
                    exclude_timepoints = "200, 300")
  pdf2 <- prepare_plot_data()
  check("scale + filter + exclusions: correct window",
        min(pdf2$time) >= 100 && max(pdf2$time) <= 500 &&
        !any(pdf2$time %in% c(200, 300)))
  session$setInputs(enable_time_filter = FALSE)

  # qc flags fire correctly through the reactive path
  qcd <- rd(csv("Time,Hot,Cold", "0,0.5,0.1", "10,1.6,0.1", "20,1.8,0.1",
                "60,1.9,0.1", "70,2.0,0.1"))
  apply_data_to_rv(qcd, "qc")
  session$setInputs(selected_samples = c("Hot", "Cold"))
  fl <- env$qc_flags(prepare_plot_data(), ref_sample = "Cold", od_cap = 1.0)
  check("qc: CEILING fires for Hot", any(fl$flag == "CEILING" & fl$sample == "Hot"))
  check("qc: GAP fires (20->60)",    any(fl$flag == "GAP"))
  check("qc: FLAT_REFERENCE fires for Cold", any(fl$flag == "FLAT_REFERENCE"))
})

cat("\n============ STRESS SUMMARY:", pass, "passed /", fail, "failed ============\n")
if (fail > 0) quit(status = 1)
