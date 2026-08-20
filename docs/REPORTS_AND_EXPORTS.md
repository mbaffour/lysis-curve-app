# Reports & Exports

Everything the app can produce, and what ends up inside each artifact.
All outputs reflect the **current** state of the plot: sample selection,
time filtering, styling — what you see is what exports.

---

## Figure exports (sidebar → Export tab)

| Format | Notes |
|---|---|
| PDF | `cairo_pdf` — system fonts (Microsoft Sans Serif, Segoe UI, …) embed correctly |
| SVG | vector, editable in Illustrator/Inkscape |
| PNG / JPEG | raster at your chosen DPI |
| TIFF | **lossless LZW compression** (journal requirement), your DPI |
| PPTX | one slide per sample (cumulative build) + final all-samples slide, vector graphics |
| Animated GIF | cumulative build, one line added per frame |

Export width/height are in **inches at print scale** with a separate DPI —
so "3.5 in wide, 600 dpi TIFF" is exactly what a journal's figure pipeline
expects. Legend width is compensated so the *panel* keeps the stated size.

### Publication presets

One click sets size + DPI + format + fonts (real print points) + the
colorblind palette together: Single column (89 mm), 1.5 column (140 mm),
Double column (183 mm) — all TIFF 600 dpi LZW — plus PDF vector,
Presentation slide (PNG), and a Prism-style visual preset. Presets use
Arial/Helvetica on purpose (journal requirement); your interactive default
font is unaffected.

## Data & statistics exports

- **Summary statistics CSV** — per sample × time: mean, SD, n (non-missing
  only), SEM, 95% CI half-width.
- **Metrics CSV** — every metric in
  [METRICS_REFERENCE.md](METRICS_REFERENCE.md), per sample (and per
  replicate where structure exists).
- **Raw data CSV** — the working dataset as currently edited.
- **Batch ZIP** (main Export tab) — all selected plots (main, derivative,
  annotated, heatmaps, curve fit) and tables in one archive at your chosen
  format/DPI.
- **Experiment notes** — TXT / CSV log-row / PDF / HTML / JSON / PNG.

## Single-file analysis reports (main Export tab)

One click, one file, everything inside — built for "attach it to the lab
notebook / send it to your PI" moments.

**HTML report** (self-contained, opens anywhere, print-to-PDF friendly):
1. Title, author, your summary notes
2. **Quality control** — the QC flags table, or "all checks passed"
3. **Data provenance** — file name, detected format, time column, samples,
   filters applied, error-display settings
4. **Experiment notebook** — every non-empty Experiment Notes field,
   including custom key–value pairs
5. **Figure** — embedded at your export settings
6. **Summary statistics** table
7. **Curve metrics** table **with the definition of every metric present**
8. Key visual settings (including μmax smoothing window and reference
   sample), and the exact R + package versions

**PDF report** — the same content as a paginated multi-page PDF (title/
provenance page, QC page, notebook page, figure page, chunked tables).
No LaTeX or pandoc needed.

Report controls: title, author, free-text notes, and a **base font size**
that scales all report text and tables in both formats. Sub/superscripts
in the report title render as real `<sub>/<sup>` in the HTML.

## Settings & themes

- **Save/Load Settings** (Export sidebar tab) — a JSON snapshot of every
  visual setting *plus* per-sample aesthetics keyed by sample name and the
  experiment notes. Loading restores settings and applies sample styles to
  matching names; your data is never touched.
- **Project themes** (Style sidebar tab) — named per-sample style sets
  stored as `themes/*.json` next to the app. With an active theme set,
  matching sample names are styled automatically on every load. Commit the
  `themes/` folder to a project repo to share the look with collaborators.
