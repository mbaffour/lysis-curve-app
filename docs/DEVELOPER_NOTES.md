# Developer Notes

Architecture, key functions, invariants, and how to test. The app is a
single R file — `Lysis Curve Ap 26.03.13.R` (~8,000 lines) — deliberately,
so "install" is "copy one file".

---

## Dependency policy

Required: shiny, tidyverse, ggpubr, scales, ggrepel, ggprism, svglite,
jsonlite, zoo, DT. Optional (feature-gated with `has_*` flags checked at
startup): officer + rvg (PPTX), gifski (GIF), base64enc (HTML report),
gridExtra (PDF report), rhandsontable (Excel-paste grid). `run_app.R`
installs anything missing into a user-writable library on first run.
New features should be pure base-R/tidyverse where possible; anything else
must be optional-gated.

## File layout (top to bottom)

1. **Header changelog**, `library()` calls, `has_*` guards
2. **Top-level pure helpers** — safe to unit-test by `source()`ing into an
   env: `safe_id`, `normalize_hex_color`, `rgb_to_hex`,
   `parse_excluded_timepoints`, `make_demo_data`, `infer_wide_replicates`,
   `pick_legend_corner`, `pub_preset_settings`, markup converters
   (`markup_has` / `markup_to_plotmath` / `markup_label` / `markup_to_html`),
   report builders (`fmt_sig`, `html_escape`, `html_table`,
   `build_html_report`, `build_pdf_report`), `METRIC_DEFINITIONS`,
   the metrics engine (`calculate_growth_metrics`,
   `calculate_infection_metrics`, `calculate_local_virulence`,
   `calculate_derivative`, `summarise_metrics`), `qc_flags`
3. **UI** — `tags$head` (JS: color preview, sub/superscript converter
   `csConvert`/`csConvertLast` + Ctrl+= shortcuts; CSS: sticky sidebar,
   pill styling, dark `.well` theme), then `sidebarLayout`: sidebar =
   always-visible header + `tabsetPanel(id="sidebar_tabs", type="pills")`
   (Axes | Style | Marks | Error | Export); main = `tabsetPanel(id=
   "main_tabs")` (Plot / Analysis / Experiment Notes / Data / Curve
   Fitting / Compare / Export)
4. **Server** — ingestion, settings/themes, dynamic per-sample UI,
   data prep, plot builder, click handlers/modals, analysis observers,
   download handlers

## The data pipeline (the part you must not break)

```
apply_data_to_rv(data, source_name)     <- ALL ingestion goes through here
  ├─ sanitises rows/columns, detects format (is_long_format_detect)
  ├─ wide: detect_od_column_idx -> rv$col_map (name -> column indices;
  │        duplicates preserved BY POSITION), infer_wide_replicates
  ├─ long: group/value/replicate column detection
  └─ theme auto-apply hook

wide_to_long(d, meas, value_name)       <- THE single wide->long conversion
  - recomputes stacked-block ids on the rows it is given (never reuse a
    cached block vector on filtered rows — that bug silently emptied
    metrics under time filtering)
  - expands duplicate columns via rv$col_map; replicate id = block or
    block_c<k> when both conventions are present

prepare_plot_data()     -> time/variable/mean/sd/n/sem/ci95 (n = non-NA only)
prepare_metrics_data()  -> per-replicate means for the metrics engine
prepare_replicate_data()-> raw replicate traces for spaghetti/jitter/bands
```

**Invariants** (all enforced by tests):
- `n` counts only non-missing values; `ci95 = qt(.975, n-1) * sem`.
- Duplicate-name columns all contribute (never subset a data.frame by a
  duplicated column *name* — use `rv$col_map` positions).
- Numeric columns are never long-format grouping candidates.
- `calculate_infection_metrics` must use `if/else`, never `ifelse()` with a
  length-1 condition (recycling bug, fixed 2026-08-19).
- `mu_max` requires rolling slope > 1e-10 (flat curve = NA, not 1e-18).
- Plain text labels take the identical old code path; plotmath conversion
  only activates when `markup_has()` is TRUE.

## Metrics engine

`calculate_growth_metrics(pd, smooth_window, extinction_threshold)` groups
by sample (+ replicate) and computes everything in
[METRICS_REFERENCE.md](METRICS_REFERENCE.md). The rolling regression is
plain base R replicating `zoo::rollapply(..., align="center", fill=NA)`
semantics. `calculate_local_virulence` finds t_stat from the reference's
rolling ln-slope and integrates trapezoids up to it (interpolating a point
at t_stat when it falls between samples).

## Testing

```bash
Rscript tests/test_analyzer.R          # 26 integration checks (shiny::testServer)
Rscript tests/stress_test_analyzer.R   # 41 stress checks (pathological inputs, scale, XSS)
```

Patterns used: `source()` the app into an env for pure helpers;
`shiny::testServer(app, { ... })` to drive the real server (call
`apply_data_to_rv`, `session$setInputs`, then the `prepare_*` functions);
render-diff the UI against `git show HEAD:` to prove no input id was lost
after UI refactors (strip the random `tab-NNNN-*` tabset keys). A "PASS"
for pathological input means correct output **or a graceful `req()`
refusal — never an uncaught crash**.

The sibling repo `killcurveplot` carries `tests/test_metrics_parity.R`,
which verifies its port of this metrics engine reproduces this app's
numbers exactly (it looks for this repo as a sibling directory).

## Gotchas

- **File encoding**: UTF-8 **without BOM** (R's `source()` rejects a BOM).
  PowerShell's `Set-Content -Encoding UTF8` adds one — use
  `[IO.File]::WriteAllText` with `UTF8Encoding($false)` instead.
- The last top-level expression must remain `shinyApp(ui, server)` —
  `run_app.R` does `source(APP)$value`.
- Shiny download handlers bind their URLs when their tab first becomes
  visible; the Export tab must be opened before its buttons have live hrefs
  (relevant for scripted testing only).
- `tabsetPanel` mints random ids per render — ignore them in UI diffs.
- The plot is drawn as a fixed-panel grob (`fix_panel_size`), so plot
  clicks carry pixel/CSS coordinates rather than a full ggplot coordmap;
  the click-to-edit zone logic works in image fractions.
- `legend_extra_px()` reserves horizontal space so panel width is constant
  across legend modes.

## Release checklist

1. `Rscript tests/test_analyzer.R` and `tests/stress_test_analyzer.R` green.
2. Source check prints `OK TRUE`.
3. If UI structure changed: render-diff against HEAD (no lost input ids).
4. Load demo data in a browser; export one HTML report; eyeball it.
5. Update `CHANGELOG.md` and the header changelog in the app file.
