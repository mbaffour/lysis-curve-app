# Lysis Curve OD Visualization App

**Interactive R Shiny platform for bacterial growth and lysis curve phenotyping**

> **Author:** Michael Baffour Awuah — Ramsey Lab  
> **Source file:** `Lysis Curve Ap 26.03.13.R`  
> **Blog post / overview:** [mbaffour.github.io/lysis-curve-app](https://mbaffour.github.io/lysis-curve-app) *(see `blogpost.html`)*

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://cran.r-project.org)
[![Shiny](https://img.shields.io/badge/Shiny-1.7%2B-brightgreen)](https://shiny.posit.co)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Data Format](#data-format)
- [Features](#features)
  - [Plot Tab](#plot-tab)
  - [Variability & Error Display](#variability--error-display)
  - [Fixed Graph Dimensions](#fixed-graph-dimensions)
  - [Experiment Notes](#experiment-notes)
  - [Analysis Tab](#analysis-tab)
  - [Data Tab](#data-tab)
  - [Curve Fitting Tab](#curve-fitting-tab)
  - [Compare Tab](#compare-tab)
  - [Export Tab](#export-tab)
  - [Save & Load Settings](#save--load-settings)
- [All 19 Growth Metrics](#all-19-growth-metrics)
- [Statistical Methods](#statistical-methods)
- [Changelog](#changelog)
- [Troubleshooting](#troubleshooting)
- [Contributing / Suggesting Changes](#contributing--suggesting-changes)
- [Citation and Contact](#citation-and-contact)

---

## Overview

The Lysis Curve OD Visualization App is a self-contained R Shiny application that turns raw plate-reader OD (optical density) CSV exports into publication-quality figures and quantitative phenotypic profiles — with no programming required beyond launching the app.

### What it does

1. Ingests CSV time-series data from any microplate reader (wide or long format, auto-detected).
2. Visualizes growth and lysis curves with fully customizable aesthetics and **7 variability display modes**.
3. Computes **19 phenotypic metrics** per sample: growth rate, lag phase, lysis onset, AUC, and five infection-context metrics.
4. Runs t-tests, Wilcoxon tests, ANOVA, and pairwise BH-corrected comparisons, with significance brackets overlaid on plots.
5. Exports figures, tables, and statistics in 6+ formats (PDF, PNG, SVG, TIFF, JPEG, PPTX, GIF).
6. Keeps a full **Experiment Notes** record exportable as TXT, CSV, PDF, HTML, JSON, or an image (PNG).
7. Saves and restores complete session state — plot settings, sample styles, and experiment notes — via a JSON settings file.
8. **Data tab** — view, inline-edit, and download the working dataset; enter data manually or paste directly from Excel.
9. **Curve Fitting** — fits Logistic and Gompertz growth models, overlays predictions, and reports A, k, t₀, and R².
10. **Compare tab** — overlay or facet two independent CSV experiments side-by-side.
11. **Batch Export** — one-click ZIP of all plots, heatmaps, and tables in your chosen format and DPI.
12. **Replicate QC** — CV% summary, outlier flagging, and per-replicate exclusion before analysis.
13. **Threshold annotation** — draw a labeled horizontal OD reference line on the main plot.

### Who it is for

Experimental microbiologists, phage biologists, and any researcher who uses a microplate reader to track bacterial cultures. No R experience needed. Supported experiment types:

- Bacterial growth curves
- Phage lysis assays
- Plasmid expression / toxicity experiments
- Any OD-based time series (antibiotics, chemical perturbations, competition assays)

---

## Installation

### Prerequisites

- **R ≥ 4.2** — [cran.r-project.org](https://cran.r-project.org)
- **RStudio (recommended)** — [posit.co/download/rstudio-desktop](https://posit.co/download/rstudio-desktop)

### Packages

```r
# Required
install.packages(c(
  "shiny", "tidyverse", "ggpubr", "scales", "ggrepel",
  "ggprism", "svglite", "jsonlite", "zoo", "DT"
))

# Optional — Excel-paste grid in Data Entry tab
install.packages("rhandsontable")

# Optional — PowerPoint export
install.packages(c("officer", "rvg"))

# Optional — animated GIF export
install.packages("gifski")
```

| Package | Role |
|---|---|
| `shiny` | Web application framework |
| `tidyverse` | Data manipulation (dplyr, tidyr) and plotting (ggplot2) |
| `ggpubr` | Statistical annotation, publication themes |
| `scales` | Axis formatting (log scales, labels) |
| `ggrepel` | Non-overlapping end-of-curve labels |
| `ggprism` | Publication-style ggplot2 theme |
| `svglite` | SVG vector graphics device |
| `jsonlite` | JSON for Save/Load Settings |
| `zoo` | Rolling window functions for metric calculation |
| `DT` | Interactive data tables |
| `rhandsontable` | Excel-paste spreadsheet grid in Data Entry *(optional)* |
| `officer` | PowerPoint generation *(optional)* |
| `rvg` | Vector graphics in PowerPoint *(optional)* |
| `gifski` | Animated GIF encoding *(optional)* |

---

## Quick Start

**Option 1 — RStudio:**
1. Open RStudio → `File > Open File` → select `Lysis Curve Ap 26.03.13.R`
2. Click **Run App** in the script editor toolbar.

**Option 2 — R console:**
```r
shiny::runApp("path/to/Lysis Curve Ap 26.03.13.R")
```

The app opens in your default browser. Keep the R console open while using it.

---

## Data Format

### Wide format *(recommended)*

Most plate readers export wide format — one column per sample, one row per timepoint.

```csv
time,SampleA,SampleB,Control
0,0.05,0.06,0.04
10,0.08,0.09,0.07
20,0.14,0.15,0.12
```

**Replicates in wide format:** Stack replicate blocks vertically. The app detects where the time counter resets and assigns replicate IDs automatically.

```csv
time,SampleA,Control
0,0.05,0.04
10,0.08,0.07
20,0.14,0.12
0,0.06,0.05   ← replicate 2 starts here
10,0.09,0.08
20,0.15,0.13
```

### Long format

A grouping column named `variable`, `condition`, `treatment`, `sample`, or `group` triggers long-format detection.

```csv
time,condition,OD
0,phage_treated,0.05
0,untreated,0.06
10,phage_treated,0.08
10,untreated,0.09
```

A column named `replicate`, `rep`, or `well` is used for replicate-aware statistics.

---

## Features

### Plot Tab

- **Sample selection** — choose any subset of uploaded samples
- **Axis scales** — linear or log (y), with manual tick control
- **Customizable aesthetics** per sample: color (palette, custom hex, or RGB sliders), point shape, line type, fill
- **Annotations** — custom text labels, vertical time markers, region highlighting, end-of-curve labels
- **Legend** — position (right, left, top, bottom, inside, none), auto label wrapping, click-to-place
- **Themes** — ggpubr clean or ggprism bordered with advanced tick marks
- **Grid lines** — major/minor, independently toggled
- **Aspect ratio** — lockable

### Variability & Error Display

Seven modes for showing data spread, selected from the **Variability / Error Display** panel:

| Mode | What is shown |
|---|---|
| **None** | Mean line only |
| **Error Bars** | ±SD, SEM, or 95 % CI as T-bars or line bars |
| **Shadow / Ribbon** | Filled band ±SD, SEM, or 95 % CI |
| **Spaghetti (Replicates)** | Every individual replicate as a semi-transparent trace |
| **Quantile Bands** | Nested shaded bands: IQR (25–75 %) inner + 2.5–97.5 % outer |
| **Jitter Points** | Raw replicate observations scattered at each timepoint |
| **Combo (Traces + CI Band)** | Replicate spaghetti traces + a light CI ribbon, mean line on top |

All modes use the actual replicate structure of your data (wide block detection or explicit replicate column).

### Fixed Graph Dimensions

The plotting area (panel) is locked to constant pixel dimensions regardless of:
- How long your sample names are
- Whether a legend is visible or not
- Which variability mode is selected

Long labels are automatically **word-wrapped** at a configurable character limit (default 20 chars, adjustable in the Legend panel). The figure container expands rightward to fit whatever the legend needs — the graph itself never shrinks.

### Experiment Notes

A dedicated **Experiment Notes** tab captures a full experimental record:

- **Identity** — Experiment ID, date, experimenter, project, institution
- **Biology** — Experiment type, host strain, phage/plasmid, MOI, replicate, passage
- **Conditions** — Growth medium, temperature, time of infection, inducer, concentration, other
- **Free text** — Observations, issues, next steps, tags
- **Custom fields** — 3 extensible key-value pairs for any extra metadata

Notes can be:
- Embedded as a caption on the main plot
- Added as a slide in PPTX exports
- Downloaded in **6 formats**: TXT, CSV (one-row log), PDF, HTML (open in browser and print), JSON (machine-readable), or **PNG image** (300 DPI, letter-size — shareable as a photo)

### Analysis Tab

After uploading data and computing metrics, the Analysis tab provides:

- **Metrics table** — 19 growth/lysis metrics per sample, sortable and downloadable
- **Bar / dot / box / violin plots** — compare any metric across samples with error bars
- **Derivative plot** — rate of OD change (dOD/dt) vs time, highlighting lysis dynamics
- **Annotated growth curves** — lag, peak OD, and lysis time marked on the main curve
- **Phenotype heatmap** — Z-score normalized matrix across all metrics and samples
- **OD-over-time heatmap** — time × sample color map for spotting temporal patterns
- **Replicate QC** — CV% summary, outlier flagging by SD threshold, per-replicate exclusion checkboxes

### Data Tab

- **View & Edit** — full interactive table of the loaded dataset; double-click any cell to edit inline; changes propagate immediately to plots and metrics
- **Revert** — restore the original uploaded file at any point
- **Download Current Data** — export the working dataset (including any edits) as CSV
- **Data Entry** — two ways to enter data without a CSV file:
  - *Option A:* Excel-paste grid (requires `rhandsontable`) — resize, paste with Ctrl+V, click Apply
  - *Option B:* Paste comma/tab/semicolon-separated text, auto-detect separator, preview, and apply
- **Rename Samples** — rename sample display labels and assign group tags without re-uploading; changes update legends and all downstream analysis

### Curve Fitting Tab

- Fit **Logistic** (3-parameter) and/or **Gompertz** growth models to each sample using nonlinear least squares (`nls`)
- Option to restrict fitting to the growth phase only (up to max OD)
- Interactive overlay plot: observed mean curve (solid) vs fitted prediction (dashed/dotted)
- Parameters table: A (carrying capacity), k (growth rate), t₀ (inflection point), R²
- Downloadable fitted plot and CSV of parameters

### Compare Tab

- Load a second CSV experiment alongside the primary dataset
- Choose **Overlay** (same axes, linetype differentiates datasets) or **Facet** (side-by-side panels)
- Customizable dataset labels
- Downloadable comparison plot

### Export Tab (Batch)

- Select any combination of outputs: main plot, derivative plot, annotated curves, both heatmaps, curve fit plot, metrics CSV, stats CSV, raw data CSV
- Set format (PNG/PDF/SVG/TIFF), DPI, and dimensions once — applied to all image outputs
- Download a single **ZIP archive** with all selected files numbered and named

### Threshold / Annotation Line

A collapsible sidebar panel lets you draw a labeled horizontal OD reference line on the main plot:
- Set the OD value, color, linetype, and line width
- Optionally add a text label anchored to the right edge of the curve
- Integrates with all export formats

### Per-Sample Time Filtering

Located inside the **Time Point Filtering** sidebar panel. Enable per-sample ranges to give each sample its own independent time window:
- One **From / To** row per currently selected sample, pre-filled with the global range
- Works on top of the existing global time filter (global runs first, per-sample narrows further)
- **Reset All to Full Range** button restores every sample at once
- Affects plots, metrics, statistics, heatmaps, and all replicate display modes consistently

### Save & Load Settings

Click **Save Settings** to download a JSON file that captures:
- All plot settings (axis scales, colors, shapes, line types, error mode, etc.)
- All per-sample aesthetic customizations
- The complete Experiment Notes record

Load the JSON on your next session to restore everything exactly. Sliders, checkboxes, dropdowns, and multi-select fields all round-trip correctly.

### Export Options

| Output | Formats |
|---|---|
| Main plot | PDF, PNG, SVG, TIFF, JPEG |
| Animated GIF | Frame-by-frame reveal of samples |
| PowerPoint | Vector plot + optional notes slide |
| Analysis plots | PDF, PNG, SVG per plot type |
| Curve fit plot | PDF, PNG, SVG, TIFF |
| Compare plot | PDF, PNG, SVG, TIFF |
| Metrics table | CSV |
| Statistics | CSV |
| Curve fit parameters | CSV |
| Notes | TXT, CSV, PDF, HTML, JSON, PNG |
| Settings | JSON (save/restore full session) |
| Batch ZIP | All of the above in one archive |

---

## All 19 Growth Metrics

| # | Metric | Description |
|---|---|---|
| 1 | **μmax** | Maximum specific growth rate (log-linear slope) |
| 2 | **Lag phase** | Time before exponential growth begins |
| 3 | **Max OD** | Peak optical density reached |
| 4 | **Time to max OD** | Time at which peak OD occurs |
| 5 | **AUC** | Area under the OD curve (trapezoidal) |
| 6 | **Final OD** | OD at the last recorded timepoint |
| 7 | **Doubling time** | ln(2) / μmax |
| 8 | **Time to half-max OD** | Time to reach 50% of peak OD |
| 9 | **Lysis onset time** | First timepoint where OD begins sustained decline |
| 10 | **Lysis rate** | Slope of OD decline during lysis |
| 11 | **Min OD post-lysis** | Lowest OD reached after lysis onset |
| 12 | **Time to min OD** | Time at which minimum post-lysis OD occurs |
| 13 | **OD drop** | Magnitude of OD decline from peak to minimum |
| 14 | **Lysis efficiency** | OD drop as a fraction of peak OD |
| 15 | **Infection strength** | Ratio of AUC (infected) to AUC (uninfected reference) |
| 16 | **Growth suppression** | 1 − (peak OD infected / peak OD uninfected) |
| 17 | **Lysis speed** | Inverse of time from infection to lysis onset |
| 18 | **Burst proxy** | Post-lysis recovery slope |
| 19 | **Recovery ratio** | Final OD / peak OD |

---

## Statistical Methods

Statistics are computed in the Analysis tab after selecting a reference sample:

- **Two-sample comparisons** — Student's t-test (parametric) or Wilcoxon rank-sum (non-parametric), based on normality
- **Multi-sample comparisons** — one-way ANOVA with Tukey HSD, or Kruskal-Wallis with Dunn's test
- **Pairwise corrections** — Benjamini-Hochberg (BH) false-discovery-rate correction
- **Significance brackets** — overlaid automatically on bar plots

---

## Changelog

### v2.1.0 — 2026-04-13

**Sidebar addition:**
- **Per-Sample Time Filtering** — independent From/To time window per sample inside the Time Point Filtering panel; stacks on top of the global filter; Reset All button; affects plots, metrics, stats, heatmaps, and all replicate modes

---

### v2.0.0 — 2026-04-13

**New tabs:**
- **Data tab** — inline cell editing, revert to upload, CSV download, manual data entry (Excel-paste grid + textarea parse), sample renaming with group tags
- **Curve Fitting tab** — Logistic and Gompertz `nls` fitting, overlay plot, parameters table (A, k, t₀, R²)
- **Compare tab** — overlay or facet two independent CSV experiments
- **Export tab** — batch ZIP of all plots and tables with unified format/DPI controls

**New Analysis sub-tab:**
- **Replicate QC** — CV% table, outlier flagging, per-replicate exclusion

**Sidebar additions:**
- **Threshold / Annotation Line** — labeled horizontal OD reference line on the main plot

**Optional new dependency:**
- `rhandsontable` — Excel-paste grid in Data Entry (app works without it)

---

### v1.0.0 — 2026-03-13

- Initial release

---

## Troubleshooting

| Problem | Likely cause | Fix |
|---|---|---|
| App does not launch | Missing packages | Re-run `install.packages(...)` from the Installation section |
| No samples appear after upload | File not parsed | Check that the file is CSV and the first column is numeric time |
| Error bars missing | Only one replicate | Replicate blocks must be present for SD/SEM to be computed |
| Metrics tab empty | Metrics not calculated | Click **Calculate Metrics** in the Analysis tab |
| Plot panel shrinks with long names | Fixed-panel feature | Increase **Wrap labels at** slider in the Legend panel to wrap long names |
| Settings file does not load | Wrong JSON structure | Only load files saved by this app |
| Data Entry grid not visible | `rhandsontable` not installed | Run `install.packages("rhandsontable")` and restart the app |
| Curve fit did not converge | Too few points or no growth phase | Try toggling **Fit growth phase only** or use a dataset with a clearer growth curve |
| Batch ZIP is empty | No data or metrics loaded | Load data, calculate metrics, then export |
| Per-sample filter has no effect | Checkbox not ticked | Enable **"Enable per-sample time ranges"** in the Time Point Filtering panel |

---

## Contributing / Suggesting Changes

Feature requests and bug reports are welcome via GitHub Issues:

**[github.com/mbaffour/lysis-curve-app/issues](https://github.com/mbaffour/lysis-curve-app/issues)**

Please include:
- A description of what you expected vs what happened
- Your R version (`R.version.string`) and OS
- A minimal CSV file that reproduces the issue (if applicable)

---

## Citation and Contact

If you use this app in published work, please cite it as:

> Baffour Awuah, M. (2026). *Lysis Curve OD Visualization App* [R Shiny application]. Ramsey Lab. https://github.com/mbaffour/lysis-curve-app

**Contact:** Open a GitHub Issue, or reach out via LinkedIn.

**Lab:** Ramsey Lab — PI: Dr. Jolene Ramsey
