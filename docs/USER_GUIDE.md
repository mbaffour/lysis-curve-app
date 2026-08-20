# User Guide — OD Growth Curve Analyzer

This is the complete walkthrough of the app. For data-format rules see
[DATA_FORMATS.md](DATA_FORMATS.md); for every metric formula see
[METRICS_REFERENCE.md](METRICS_REFERENCE.md); for exports and reports see
[REPORTS_AND_EXPORTS.md](REPORTS_AND_EXPORTS.md).

---

## 1. Launching

| How | Steps |
|---|---|
| **Double-click (recommended)** | Windows: `Run OD Growth Curve Analyzer.bat` · macOS: `Run OD Growth Curve Analyzer.command` · Linux: `./run.sh`. First run auto-installs missing R packages into a user-writable library; afterwards the app opens instantly in your browser. |
| **Terminal** | `Rscript run_app.R` |
| **RStudio** | Open `Lysis Curve Ap 26.03.13.R` → **Run App** |

Requirements: R ≥ 4.2. Everything else installs itself.

**First time?** Click **Load Demo Data** under the file picker: a built-in
phage-infection dataset (uninfected control + two MOIs, 3 replicates each,
stacked-block format) loads instantly so you can explore every feature
without a file.

## 2. The layout

- **Sidebar (left, sticky)** — all plot controls, grouped into pill tabs:
  **Axes | Style | Marks | Error | Export**. On desktop the sidebar scrolls
  in its own pane, so the plot never leaves your view. The file picker,
  **Load Demo Data**, the green **data summary box**, and Night Mode are
  always visible above the tabs.
- **Main panel (right)** — the working tabs: **Plot · Analysis · Experiment
  Notes · Data · Curve Fitting · Compare · Export**.

**Always check the data summary box after loading.** It states exactly what
the app detected: format (wide/long), time column, number of time points,
each sample, and its replicate structure (columns and/or stacked blocks).
If it does not match your intent, fix the file — silent misreads are how
wrong figures happen.

## 3. Styling the plot

Everything is live: change anything and the plot re-renders.

- **Click-to-edit** (Style tab toggle, on by default): click the plot's
  title, an axis label, or the legend region to edit exactly that element —
  text, font size, font family, and for legend entries the sample's label,
  color, shape, line type. Clicking the legend region preselects the sample
  nearest your click.
- **Sub/superscripts**: highlight characters in any text field, then click
  **x₂ Subscript** / **x² Superscript** (or press **Ctrl+=** / **Ctrl+Shift+=**,
  as in Word). `A550` → select `550` → click → `A₅₅₀`. Click again to undo.
  Power path: `A_{550}` / `10^{-7}` markup is also understood everywhere.
- **Fonts**: default is Microsoft Sans Serif (falls back to standard sans on
  machines without it). Separate size controls for title, axis labels, axis
  text, and legend.
- **Per-sample styles** (Style → Variable Styling): color (presets, HEX, or
  RGB sliders), point shape and fill, line type, and the legend label for
  each sample. Selection order = legend order.
- **Palettes**: Viridis, Plasma, Okabe-Ito colorblind-safe, publication,
  rainbow, grayscale, or fully custom.

## 4. Legend

Legend → Position offers:

| Mode | What it does |
|---|---|
| Right / Left / Top / Bottom | Classic legend box outside the panel |
| Inside (custom) | Box inside the plot: corner preset buttons, exact X/Y, or **click the plot to place it** |
| **Inside — auto** | The app counts data points in each corner and puts the legend in the *least crowded* one, so it never covers your curves |
| **Direct labels at line ends** | No legend box at all: each curve is named at its end, in its own color, using your custom legend labels; overlap is automatically avoided |
| None | No legend |

Long labels wrap at a width you control, so the panel size stays constant.

## 5. Error display and replicates

Error type: SD, SEM, or 95% CI (t-based), with a multiplier; displayed as
T-bars, line ranges, or a shaded ribbon. Seven replicate display modes,
including spaghetti (every replicate as a faint trace), jitter points, and
quantile bands. `n` counts only non-missing values — replicates with a
missing reading do not silently shrink your error bars.

## 6. Analysis tab

1. **QC banners** appear automatically when something is off: saturated OD
   (> the cap you set), negative values, acquisition gaps, or a flat
   reference culture (which invalidates every ratio metric).
2. **Calculate Metrics** — per sample (and per replicate): μmax, lag phase,
   doubling time, AUC, lysis onset/rate, OD drop, residual OD, recovery
   slope, centroid metrics, first-peak and extinction time, and more. See
   [METRICS_REFERENCE.md](METRICS_REFERENCE.md) for every formula and source.
3. Pick a **Reference sample** (your uninfected control) and **Compute
   Infection Metrics** for relative growth, infection strength, relative
   μmax/max OD, lysis-onset Δ, **Centroid Index**, and **Storms local
   virulence** (with the integration window shown, and overridable).
4. **Statistics** — t-test / Wilcoxon, ANOVA + Tukey, Kruskal–Wallis + Dunn,
   BH correction, significance brackets on bar plots.
5. **Plots / Heatmaps / Replicate QC** sub-tabs: metric bar/box/violin
   plots, phenotype Z-score heatmap, OD-time heatmap, derivative plot,
   CV% and outlier flagging.

## 7. Publication presets (Export tab)

One click sets fonts (real print points), figure dimensions, DPI, format,
and the colorblind palette together:

| Preset | Size | Output |
|---|---|---|
| Single column | 89 mm (3.5 in) | TIFF 600 dpi (LZW) |
| 1.5 column | 140 mm (5.51 in) | TIFF 600 dpi (LZW) |
| Double column | 183 mm (7.2 in) | TIFF 600 dpi (LZW) |
| PDF vector | 7 × 5 in | PDF (cairo, fonts embedded) |
| Presentation slide | 10 × 7.5 in | PNG 300 dpi |
| Prism-style look | visual only | — |

Presets use Arial/Helvetica deliberately — that is what journals require —
even though the app's everyday default is Microsoft Sans Serif.

## 8. Project themes (Style tab)

Save the current per-sample styles (color, shape, line type, label) as a
named theme. Themes are JSON files in a `themes/` folder next to the app —
commit that folder to your project repository to share with collaborators.
Set an **Active theme** with auto-apply on, and every dataset you load gets
matching sample names styled identically: "WT" is the same black solid
circle in every figure of the project, forever.

## 9. Experiment Notes and reports

Fill in the **Experiment Notes** tab (experiment ID, strains, phage, MOI,
conditions, observations, custom fields). Every non-empty field is embedded
automatically into the one-click **single-file HTML report** and **multi-page
PDF report** (Export tab), alongside the figure, summary statistics, all
computed metrics with their definitions, QC flags, data provenance, and the
exact R/package versions. One file = the complete, defensible record of the
analysis. Details: [REPORTS_AND_EXPORTS.md](REPORTS_AND_EXPORTS.md).

## 10. Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| Data summary shows 1 replicate when you have 3 | Check [DATA_FORMATS.md](DATA_FORMATS.md) — replicates must be duplicate column names or stacked time blocks (or a `replicate` column in long format) |
| Error bars absent | Only one replicate per sample/time — nothing to show |
| Metrics table empty | Click **Calculate Metrics**; check the QC banners |
| Infection metrics absent | Pick a reference sample and click **Compute Infection Metrics** |
| Sub/superscript button does nothing | Select (highlight) the characters first |
| A ratio metric looks absurd | Read the QC banner: a flat or saturated reference invalidates ratios; blank-subtract your data |
