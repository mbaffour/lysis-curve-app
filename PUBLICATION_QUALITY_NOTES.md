# Publication-Quality Review Notes

This review audited `Lysis Curve Ap 26.03.13.R` against a publication-quality
rubric for figures and statistical reporting. The app is already strongly
publication-oriented (vector PDF/SVG + high-DPI PNG/TIFF export with a
configurable DPI defaulting to 300, configurable figure width/height, a
`theme_pubr()` / `theme_prism()` clean theme, editable axis titles with units,
per-element font-size controls, selectable error statistic SD/SEM/95% CI, an
Okabe-Ito colorblind palette, viridis, significance brackets on the bar plot,
and named tests in the stats output). The changes below are conservative
additions and clarifications only.

> **IMPORTANT — UNVERIFIED:** These edits were made without an R runtime. The
> app was **not executed**. Please launch the app and smoke-test the Plot,
> Analysis (Statistics), and export paths before merging.

## (a) Changes applied

1. **Palette label now names the colorblind-safe scheme (Okabe-Ito).**
   The "Colorblind-friendly" dropdown entry already returned the Okabe-Ito
   8-colour palette; relabeled to **"Colorblind-safe (Okabe-Ito)"** so users can
   cite it correctly. Pure UI string change; the underlying value (`"colorblind"`)
   and colours are unchanged.

2. **Statistical results table now reports sample size (n).**
   `compare_metrics()` previously computed per-group replicate counts but
   discarded them. The results tibble now carries per-comparison `n1`/`n2`
   (2-group t-test / Wilcoxon and pairwise post-hoc) and `N_total` (ANOVA row).
   These flow automatically into the on-screen stats table (`output$stats_table`)
   and the `stats_*.csv` / `all_stats_*.csv` exports without any change to the
   test logic, the BH correction, or the significance thresholds.

3. **Optional error-bar definition + n line on the main lysis-curve figure.**
   New checkbox **"Add error-bar definition + n to caption"** (default OFF,
   independent of the existing experiment-info caption). When enabled, the figure
   caption gains a line stating what the variability display represents — e.g.
   `Error bars = SEM; n = 3` / `Shaded band = 95% CI; n = 3-4` /
   `Bands = IQR (inner) and 2.5-97.5% (outer)` — derived from the current
   `error_display_mode`, `error_type`, `error_multiplier`, and the replicate
   count already carried in `plot_data$n`. Purely additive, guarded with the same
   `!is.null(input$...)` pattern used throughout, and reuses the existing caption
   theming path.

All three changes preserve existing reactivity, defaults, and statistical
methodology. The settings save/restore path uses `reactiveValuesToList(input)`,
so the new checkbox round-trips automatically.

## (b) Recommendations deferred (need author's scientific decision or an R test-run)

These were intentionally **not** applied — they either change methodology (the
author's deliberate choice) or need an R run to verify safely:

1. **Effect sizes and confidence intervals in the stats output.** Consider
   reporting Cohen's d / Hedges' g (2-group) and the mean-difference 95% CI
   alongside p-values. `t.test()` already returns `$estimate` and `$conf.int`;
   surfacing them is high value but adds columns/logic that should be R-tested.

2. **Which multiple-comparison correction is the default.** Pairwise post-hoc is
   hard-coded to Benjamini-Hochberg (`p.adjust.method = "BH"`). Recommend exposing
   a correction selector (BH / Holm / Bonferroni / none) as an **option**, leaving
   BH as the default. Not changing the default silently.

3. **Parametric vs non-parametric default.** For 2 groups both `t.test` and
   `wilcox.test` are reported (good). For >=3 groups only ANOVA + pairwise t is
   offered. Consider adding Kruskal-Wallis + Dunn as a user-selectable
   non-parametric path — an addition, not a default change.

4. **Name the test in the bar-plot significance annotation.** The brackets show
   stars but the test name lives only in the stats table. Consider an optional
   sub-caption on the bar plot, e.g. "pairwise t-test, BH-adjusted; * p<0.05".
   Deferred because bar-plot caption layout should be visually checked in R.

5. **Reproducible methods/script export.** Consider a "copy methods sentence" /
   `sessionInfo()` + settings-JSON bundle so a figure's exact parameters and
   package versions travel with it. The settings JSON already captures inputs;
   adding `sessionInfo()` capture would complete reproducibility.

6. **Assumption checks.** Optional normality (Shapiro) / variance (Levene)
   diagnostics could guide the parametric-vs-non-parametric choice, surfaced as
   advisory output rather than gating the analysis.

7. **State n and error definition in the bar-plot too.** The bar plot already
   encodes the error type in its y-axis label and shows individual points;
   consider mirroring the new "n = ..." annotation there for consistency.
