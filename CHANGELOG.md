# Changelog

## 2026-08-20 — figure-craft wave

### Draggable floating legend labels (`7f3eae3`)
- New legend mode "Floating labels (drag to place)": each sample name is an
  individually draggable element on the plot; positions are drawn into the
  ggplot (panel-npc annotations via a new `annotation_npc()` geom, since
  `annotation_custom()` collapses on log10 axes), so all exports match the
  placement. Positions persist in settings and project themes. Panel-rect
  measurement accurate to 0.01 px; drag guarded against click-to-edit.
- SVG exports switched to svglite with `fix_text_size = FALSE`: real
  editable `<text>` elements, no `textLength` stretching in Inkscape.
- Settings-loader hardening for partial/hand-edited files.

### Exact graph-panel size + presentation parity (this commit)
- "Exact graph-panel size" mode (default): the user sets the panel itself
  (px on screen, inches for export); the canvas is measured — not
  estimated — to fit axis decorations and the actual legend. A 3-sample
  and a 30-sample legend produce exactly identical panels. Tall legends
  gain pad rows rather than clipping. "Reserved right space" is obsolete
  in this mode.
- PPTX/GIF cumulative builds: every slide/frame shares the final frame's
  canvas and legend-cell width, so the full-sized graph box is present
  from frame 1 and the panel never shifts as lines/legend entries appear.
  (PPTX additionally fits the canvas to the slide box, scaling uniformly.)
- Publication presets set the export panel dimensions too.

## 2026-08-19 — merge with killcurveplot + major feature wave

The app became the single merged base for the kill-curve tooling
(absorbing the killcurveplot v2.8 work), with an accuracy audit, new
science, reports, and UI overhaul. Chronological by commit:

### Accuracy fixes (`4f126c5`)
- **Infection metrics recycling bug**: `relative_growth`,
  `infection_strength`, `relative_mu_max`, `relative_max_od` were computed
  with a length-1 `ifelse()`, giving every sample the first sample's
  values. Each sample now gets its own ratio.
- **Wide-format replicate handling rebuilt** (`wide_to_long` + `col_map`):
  duplicate column headers are pooled by position (previously only the
  first column of each name was used — n=1, no error bars), and stacked-
  block replicate ids are recomputed on filtered rows (previously time
  filtering + stacked replicates silently emptied Calculate Metrics).
- `n` counts only non-missing values (SEM/CI no longer understated).
- Settings save no longer drops `shape_size` / `color_palette`.
- `rgb_to_hex` no longer appends a spurious alpha channel.
- `mu_max` requires slope > 1e-10 (flat curves report NA, not noise).
- Load Demo Data button + data summary box (detected format/replicates).

### Reports & launchers (`e9dfb09`)
- Single-file **HTML** and multi-page **PDF analysis reports**: figure,
  experiment notes, summary statistics, metrics with definitions,
  provenance, settings, session info. No pandoc/LaTeX.
- Double-click launchers (`.bat` / `.command` / `run.sh`) + auto-installing
  `run_app.R` (the repo previously had no launcher).
- TIFF exports use lossless LZW.

### Science quick wins (`e8e83eb`)
- **Centroid Index** (Hosseini 2024) with exact trapezoid centroids.
- **PNAS-2025 metrics**: `first_peak_od`/`t_first_peak` (prominence-based)
  and `extinction_time` (interpolated threshold crossing); caveat notes on
  `mu_max`/`lysis_rate`.
- **Storms local virulence** with auto-detected t_stat window + override;
  `infection_strength` relabeled as the full-run AUC ratio (≡ PhageScore/100).
- **QC flags** (CEILING / NEGATIVE / GAP / FLAT_REFERENCE) in-app and in
  reports.
- Demo growth rate tuned so the control reaches stationary phase.

### UI overhaul (`1ae0281`)
- Sticky sidebar (plot stays visible while scrolling controls).
- Sidebar reorganized into pill tabs: Axes | Style | Marks | Error | Export.
- Format-detection fix: numeric flat columns no longer flip wide files to
  long format.

### Click-to-edit, presets, themes (`f7fa5ee`)
- Click the plot title / axis labels / legend region to edit that element.
- Publication presets: 89/140/183 mm TIFF-600-LZW, PDF vector, slide,
  Prism-style.
- Project themes: named per-sample style sets in `themes/*.json`,
  auto-applied to matching sample names on load.

### Legend & fonts (`fd4d508`)
- Legend modes: **Direct labels at line ends** (colored, custom names,
  overlap-avoiding) and **Inside — auto** (least-crowded corner).
- Legend font size control; font family + sizes editable in every
  click-to-edit dialog.
- Microsoft Sans Serif + Segoe UI font choices; MS Sans Serif default;
  `cairo_pdf` export for font embedding.
- Legibility defaults: point size 4, line 1.2, legend 18 pt.

### Sub/superscripts (`3eb799c`, `ae173d1`)
- `_{...}` / `^{...}` plotmath markup in all plot text (titles, axes,
  legend, end labels, annotations), consistent across every export.
- **Highlight-and-click**: select characters → x₂/x² buttons or
  Ctrl+= / Ctrl+Shift+= convert to real Unicode sub/superscripts in place
  (toggle to undo). WYSIWYG; works in every launch mode.

### Documentation & tests (this commit)
- `docs/` (user guide, metrics reference, data formats, reports/exports,
  developer notes), this changelog, and the verification suites moved
  into `tests/` (26 integration + 41 stress checks).

## Earlier
See the header changelog inside `Lysis Curve Ap 26.03.13.R` and the git
history (v2.1.0 per-sample time filtering, etc.).
