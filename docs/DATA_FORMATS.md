# Data Formats

The app ingests CSV. Format, time column, samples, and replicate structure
are auto-detected — and whatever was detected is displayed in the green
**data summary box** after every load. Always read it.

---

## Wide format (recommended; what plate readers export)

One row per time point; first (or `time`-named) column is time; every other
column is a sample.

```csv
time,SampleA,SampleB,Control
0,0.05,0.06,0.04
10,0.08,0.09,0.07
```

### Replicates in wide format — two conventions, combinable

**1. Stacked blocks** — repeat the whole time series vertically; the app
detects each place the time counter resets and assigns block replicate IDs:

```csv
time,SampleA,Control
0,0.05,0.04
10,0.08,0.07
0,0.06,0.05    <- replicate 2 starts here
10,0.09,0.08
```

**2. Duplicate column headers** — replicate columns share one name; all of
them are pooled per sample (resolved internally by column position):

```csv
time,SampleA,SampleA,SampleA,Control,Control,Control
0,0.05,0.06,0.04,0.04,0.05,0.04
10,0.08,0.09,0.07,0.07,0.08,0.06
```

Both at once also works: a file with 2 duplicate columns × 2 stacked blocks
yields n = 4 per sample/time, with distinct replicate IDs (`block_c<k>`)
for the replicate-level display modes.

> **Accuracy history (fixed 2026-08-19):** duplicate-name columns were
> previously collapsed to the *first* column of each name (means from one
> replicate, no error bars), and enabling time filtering with stacked
> blocks silently emptied Calculate Metrics. Both fixed and regression-
> tested; the data summary box exists so you can always confirm what was
> detected.

## Long format

Detected when a grouping column named `variable`, `condition`, `treatment`,
`sample`, or `group` (any case) exists — or, for small files, when a
**non-numeric** low-cardinality column looks like a grouping factor.
(Numeric columns are never treated as grouping candidates: a flat/constant
OD column must not flip a wide file into long format — this was a real bug,
fixed 2026-08-19.)

```csv
time,condition,OD
0,phage_treated,0.05
0,untreated,0.06
10,phage_treated,0.08
```

A column named `replicate`, `rep`, `rep_id`, `well`, `well_id`,
`technical_replicate`, `biological_replicate`, or `repeat` is used for
replicate-aware statistics and display modes.

## Detection rules (what the app actually does)

1. Rows whose time value is not a finite number are dropped up front.
2. Blank-named columns are dropped.
3. Time column: exact `time`/`t` name match first, then partial `time`
   match, then the first monotonically non-decreasing numeric column, else
   column 1.
4. Wide OD columns: name-matched (`od`, `absorbance`, `value`) if any,
   else all numeric columns — **by position**, so duplicates survive.
5. Missing values: fine anywhere. Statistics count only non-missing
   replicates (`n`, SEM, CI are not diluted by NAs).

## Edge cases that are tested and safe

Header-only files, single time point, single sample, all-NA samples,
non-numeric junk columns, negative times, duplicate time values, Unicode
sample names (Δ, µ, …), unequal stacked blocks, 100 samples × 100 time
points (≈4 s). See `tests/stress_test_analyzer.R`.

## Practical tips

- **Blank-subtract before importing.** Every ratio metric (virulence,
  Centroid Index, infection strength) is biased toward "no effect" by an
  unsubtracted media blank.
- Keep sample names consistent across experiments so **project themes**
  can restyle them automatically.
- The Data tab can view/edit cells inline, paste from Excel, and rename
  samples after import.
