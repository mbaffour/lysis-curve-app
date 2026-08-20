# Metrics Reference

Every metric the Analyzer computes, with the exact formula, the literature
source where one exists, and the pitfalls. All metrics are computed from the
per-sample **mean curve as plotted** — i.e. after time filtering and
point exclusion — and additionally **per replicate** when replicate
structure exists (the table then also shows per-sample means of the
per-replicate values).

Notation: `t_i` time points, `y_i` mean OD at `t_i`, `n` number of points.

---

## Growth metrics

| Metric | Definition |
|---|---|
| `initial_od` | `y_1` (first plotted time point) |
| `max_od`, `time_max_od` | Global maximum of the mean curve and its time |
| `min_od`, `time_min_od` | Global minimum and its time |
| `final_od` | `y_n` |
| `auc` | Trapezoidal area: `Σ (t_{i+1}−t_i)·(y_i+y_{i+1})/2` |
| `mu_max` | Maximum specific growth rate: the largest **centered rolling-window regression slope of ln(y) vs t** (window = the μmax smoothing setting, default 5 points; window is clamped to n). Slopes must exceed 1e-10 — an exactly flat curve reports `NA`, not floating-point noise. |
| `doubling_time` | `ln(2) / mu_max`; `NA` when there is no growth phase |
| `lag_phase` | First time at which the rolling ln-slope reaches **10% of mu_max** |
| `stat_phase_dur` | Time spanned by all points ≥ **95% of the peak** |

## Lysis metrics

| Metric | Definition |
|---|---|
| `lysis_time` | Lysis onset: first post-peak time at which the decline from the peak exceeds **5% of the peak**; `NA` if never |
| `lysis_rate` | Linear-regression slope of the mean from the **peak to the post-peak trough** (negative = decline) |
| `od_drop` | Peak minus post-peak minimum |
| `residual_od` | Post-peak minimum (what lysis left behind) |
| `recovery_slope` | Positive regression slope from the post-peak trough to the end (regrowth); `NA` if no recovery |
| `lysis_efficiency`* | `od_drop / max_od` |
| `recovery_ratio`* | `final_od / max_od` |
| `t_half_max_decline`* | First time after the peak at which the mean falls below **50% of the peak**, linearly interpolated; `NA` if never |
| `steepest_decline_rate/interval`* | Most negative slope between consecutive mean points, and where |

\* companion metrics carried over from the killcurveplot v2.8 engine.

**Caveat (Blazanin et al., PNAS 122(37):e2513377122, 2025).** In simulation,
exponential-phase growth rate and the rate of population decline do **not**
track phage infectivity; **time of peak density** and **extinction time**
do. The metric labels carry this note.

## PNAS-2025 recommended metrics

| Metric | Definition |
|---|---|
| `first_peak_od`, `t_first_peak` | The **first** local maximum of the lightly smoothed curve (centered rolling mean, window ≤ 3) with topographic prominence ≥ 5% of the curve's range. A regrowing culture's *global* maximum may be its second peak — this metric is not fooled. Falls back to the global max when no interior peak exists. |
| `extinction_time` | First time the curve falls **below the extinction threshold** (Analysis tab input, default 0.05), linearly interpolated between the bracketing points; `NA` when never crossed. Note: OD cannot see the 10⁴ CFU/mL threshold used in the paper's simulations — plate-reader floors are ~10⁶–10⁷ CFU/mL; the threshold here is OD-native. |

## Centroid metrics — Hosseini et al., *Commun Biol* 7:673 (2024), doi:10.1038/s42003-024-06379-z

Center of the area under the OD–time curve, computed with **exact
trapezoid-strip centroids** (this choice matters at coarse sampling; the
app's methods text states it):

```
per strip:  h = t_{i+1}−t_i,  a = y_i,  b = y_{i+1}
A_i  = h(a+b)/2
x̄_i  = t_i + h(a+2b) / (3(a+b))
ȳ_i  = (a² + ab + b²) / (3(a+b))
centroid_x = Σ x̄_i A_i / Σ A_i      centroid_y = Σ ȳ_i A_i / Σ A_i
```

| Metric | Definition |
|---|---|
| `centroid_x`, `centroid_y` | As above; `NA` when total area ≤ 0 or n < 2 |
| `centroid_index` (needs a reference) | `CI = 1 − (x̄·ȳ) / mean(x̄_ref·ȳ_ref)` — the denominator is the mean of the reference replicates' **products**, so the reference vs itself is exactly 0 under replicate averaging. CI = 1: maximal efficacy · CI = 0: none · **CI < 0: resistance/overgrowth**. CI is sensitive to late regrowth that AUC ratios miss — when CI and the AUC-based metrics disagree, *that disagreement is the finding*. Doubly sensitive to baseline offset: blank your data. |

## Reference-normalized infection metrics (pick an uninfected reference)

| Metric | Definition |
|---|---|
| `relative_growth` | `auc / auc_ref` |
| `infection_strength` | `1 − auc/auc_ref` over the **full run**. This is algebraically **PhageScore/100** (Vandersteegen lineage), *not* the Storms virulence index — do not conflate them. |
| `relative_mu_max` | `mu_max / mu_max_ref` |
| `relative_max_od` | `max_od / max_od_ref` |
| `lysis_onset_delta` | `lysis_time − lag_phase` |
| `local_virulence` | **Storms et al., *PHAGE* 1(1):27–36 (2020), doi:10.1089/phage.2019.0001.** `v = 1 − A_i/A_0` with **both** areas integrated only from 0 to `t_stat` — the onset of stationary phase in the phage-free reference (auto-detected as the first time after μmax at which the rolling ln-slope drops below 10% of μmax; a manual override input exists, and the t_stat actually used is displayed). Differs from `infection_strength` exactly when the treated culture's behavior after the control plateaus matters. Requires blanked data; not comparable across hosts/media/temperature. Negative values are legal and are not clamped. |

> **History.** Before 2026-08-19 the reference implementation computed the
> four relative metrics with a length-1 `ifelse()`, silently recycling the
> first sample's value to every sample. Fixed; a regression test asserts the
> reference sample's self-relative values are exactly 1 / 0 / 1 / 1.

## QC flags

| Flag | Trigger | Why it matters |
|---|---|---|
| `CEILING` (warn) | max OD > cap (default 1.0) | Plate-reader linearity is lost above ~0.5–1.0; ratios biased |
| `NEGATIVE` (error) | any OD < 0 | Over-subtracted blank |
| `GAP` (warn) | any time step > 3× the modal interval | Derivative/timing metrics unreliable across the gap |
| `FLAT_REFERENCE` (error) | reference range < 0.1 | Every ratio metric (virulence, CI, infection strength) is invalid |

Flags appear as banners in the Analysis tab and as a **Quality control**
section in the HTML/PDF reports. They never auto-exclude data — you decide.

## Statistical methods (Analysis tab)

- Summary statistics: mean, SD, `n` (**non-missing values only** — NA
  replicates do not inflate n and deflate SEM), `SEM = SD/√n`,
  `CI95 half-width = t(0.975, n−1) × SEM`.
- Two samples: Student's t-test or Wilcoxon rank-sum (by normality).
- Multi-sample: one-way ANOVA + Tukey HSD, or Kruskal–Wallis + Dunn.
- Multiplicity: Benjamini–Hochberg FDR correction.
- Significance brackets are drawn on the metric bar plots.

## Key sources

- Storms ZJ et al. *PHAGE* 1(1):27–36, 2020 — virulence index. doi:10.1089/phage.2019.0001
- Hosseini N et al. *Commun Biol* 7:673, 2024 — Centroid Index. doi:10.1038/s42003-024-06379-z
- Blazanin M et al. *PNAS* 122(37):e2513377122, 2025 — which metrics track infectivity. doi:10.1073/pnas.2513377122
- Classification of phage–host growth dynamics: *Microorganisms* 9(12):2470, 2021. doi:10.3390/microorganisms9122470
