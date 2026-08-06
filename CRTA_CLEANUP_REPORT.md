# `compare_retention_time_alignment.py` — Cleanup Report

First cleanup pass on `bin/compare_retention_time_alignment.py` (clusterless RT-alignment
step). Branch `alignment_clean`. **Goal of this pass: remove dead code / improve readability
with ZERO change to output.** All output-affecting issues below are *flagged only* — nothing
in section C was changed.

File shrank **1984 → 1617 lines** (−367).

---

## A0. Readability pass (verified output-neutral)

- Removed the dev-only `--show-plots` flag and its interactive `fig.show()` path (never used
  by the pipeline; needs a browser).
- Added a module-docstring table of contents and **7 numbered section banners**
  (formatting/outlier helpers · curve-fitting models · feature loading · matching/comparison ·
  HTML report · corrections export · entry point), sub-banners for the picklable wrappers and
  the fitting dispatcher, and **Phase 1-7 markers through `main()`**.
- Added block comments across the denser code (LOESSInterpolator, dispatcher, main phases).
- All comment/structure only — verified byte-identical JSON + HTML for polynomial and loess.

## A. What was changed (verified output-neutral)

All removals below are either never-called code or a provably-dead branch. Verified by A/B
diff (see section D): identical HTML and identical JSON for every fitting method.

| # | Removed | Why safe |
|---|---------|----------|
| A1 | `create_comparison_plot()` (~115 lines) | Never called. main() only builds the RT-vs-RT-diff plot, never the M/Z-vs-RT-diff plot. |
| A2 | `create_rt_transformation_plots()` (~141 lines) | Never called anywhere. |
| A3 | `serialize_model_to_json()` (~48 lines) | Never called. main() uses the `serialization_data` returned by `fit_retention_time_model()` instead. |
| A4 | `find_decimal_precision()` (~13 lines) | Never called. |
| A5 | `import make_subplots` | Unused. |
| A6 | Whole `sklearn` import block (`GaussianProcessRegressor`, `SKLearnRBF`, `ConstantKernel`, `HAS_SKLEARN`) | `HAS_SKLEARN` is set but never read; none of the symbols are used. RBF fitting uses `scipy.interpolate.Rbf`, not sklearn. |
| A7 | Dead pre-fit branch inside `create_rt_vs_rtdiff_plot()` (~26 lines) + its two vestigial params (`corrections_dict`, `equations_plain_dict`) | The branch tested `file1_name in corrections_dict`, but main() passes the **flat** `{"f1_vs_f2": …}` dict, so a bare filename is never a key → branch never ran. It also referenced non-existent `'model'`/`'pikle:'` keys. The report always silently re-fits via the fallback path (which is now the only path). Confirmed HTML unchanged. |

External importers were checked: `compare_retention_time_alignment_2d.py` imports only
`load_feature_json, get_coordinate, resolve_mz_collisions,
match_features_by_closest_mz_and_rt, round_with_precision` (all kept); pairing imports only
`LOESSInterpolator` / `PiecewisePolynomial` (both kept, unchanged).

---

## C. Flagged problems — **would affect output, NOT changed** (need your call)

### C1. ✅ FIXED — `--loess_iterations` was silently ignored by the actual fit  *(highest impact)*
`fit_loess_equation()` hardcoded `it=3` in the `lowess(...)` call. The `loess_iterations` value
was threaded through `fit_retention_time_model` but **never passed to `lowess`**, and the
equation label also hardcoded `"it=3"`. (Note: the live JSON path stores only `model_pickle`
for loess — the earlier claim that `iterations` was written to JSON metadata was via the now
-deleted `serialize_model_to_json`, not the live path.)
- **Fix applied (option a, approved):** added an `iterations` parameter to
  `fit_loess_equation`, passed it as `lowess(it=iterations)`, and threaded
  `loess_iterations` into all four call sites. The label is now dynamic (`it={iterations}`).
- **Verified behaviour change:** with the pipeline default (`--loess_iterations 100`), all
  20/20 corrections' fitted models now differ from the old `it=3` output — i.e. the knob now
  genuinely takes effect. Output is deterministic run-to-run under `PYTHONHASHSEED=0`.
- **⚠ This changes `fitted_corrections.json` for the loess method** (the pipeline default),
  which is the intended effect. Other methods (polynomial/spline/rbf/piecewise) are untouched.

### C2. ✅ FIXED — `--spline_degree` / `--spline_smoothing` defaults (was a real bug, not just a doc mismatch)
- `--spline_degree` default was **8**, help said "default 3"; `--spline_smoothing` default was
  **1000**, help said "default 100".
- **It was broken, not just inconsistent:** `UnivariateSpline` requires `1 ≤ k ≤ 5`, so a
  flagless spline run (default degree 8) fitted **0 corrections** (`k should be 1 <= k <= 5` on
  every pair). Verified against the original script.
- **Fix:** aligned the argparse defaults to the documented/valid values (`spline_degree=3`,
  `spline_smoothing=100`); flagless spline now fits all pairs. Pipeline output unchanged (the
  `.nf` always passes both explicitly, 3 / 50).

### C3. ✅ FIXED — removed the inert `'direction': 'forward'` field
Each `corrections[pair_key]` carried a hardcoded `direction: 'forward'` that never varied and
carried no information. Verified **no consumer reads it** (pairing/orchestrator and the older
pair_features variants use only `method` + `serialization` + the `f1_vs_f2` key; the
`direction` hits elsewhere in the repo are unrelated edge-direction logic).
- **Fix:** dropped the `direction` field. The genuinely-useful inspection fields
  (`file1`, `file2`, `n_points`, `n_points_fitted`, `n_points_filtered`, `equation`) were
  **kept** — they aid debugging and are not misleading.
- **Verified:** the new JSON equals the old JSON with `direction` stripped (nothing else
  changed); pairing still loads it since those keys were never read.

### C4. ✅ FIXED — `generate_curve_approximation()` used an x-*value* as a point *count* (HTML bloat)
`n_points = int(max_x)` rendered one curve point per second of RT — ~9,442 points per fitted
curve on real data. Combined with a **dead `customdata`** on each fit trace (9,442 copies of an
identical string that nothing read — the toggle JS matches by trace `name`), the fitted-curve
traces were **~12 MB of a 13 MB report**.
- **Fix:** fixed render resolution `CURVE_RENDER_POINTS = 200` (a report line needs no more),
  and removed the dead `customdata`.
- **Result:** on a **real 15-file run** the report went **133 MB → 17 MB (~8×)**;
  `fitted_corrections.json` byte-identical (HTML-only change). (In the 133 MB original, the dead
  `customdata` alone was 48 MB and the oversized fit curves ~89 MB.)
- **Remaining 17 MB** is the raw scatter: ~8 MB per-point hover strings + ~4 MB full-precision
  coordinate arrays. Two further levers were identified — round coordinates to ~3 decimals
  (~2 MB) and move hover to a Plotly `hovertemplate` + `customdata` (~6 MB), together ~7 MB —
  but **per user decision, 17 MB is acceptable and these were not applied.**

### C5. (Known, from handoff) domain mismatch, `rt_source=start` vs pairing's `y_center`
Not a code defect in this file — a config choice. The `.nf` fits on `rt_start`/`mz_start`
(`--rt_source start`) while pairing evaluates the curve at `y_center` (median ~8.9 s field
gap). Already documented; the direction half was fixed on the pairing side. No action here. Further tests required, at this point the start field seems to provide better results, but that was without applying the new unweighted center calculations

---

## E. Output-neutral nits

- **E1 ✅ FIXED** — `calculate_percentile_threshold()` comments/docstring said "IQR method";
  it's plain symmetric percentile trimming. Reworded (code comments/docstring only — the one
  "IQR" string that lives *inside the HTML report* at the config info-box was left untouched,
  since editing it would change the HTML bytes; see E4-adjacent note).
- **E2 ✅ FIXED** — `fit_polynomial_equation()` / `fit_polynomial_equation_plain()` built the
  equation with self-replacing no-ops (`.replace("+ ", "+ ")`, `.replace("- ", "- ")`).
  Removed both; kept the defensive `.replace("+-", "-")`. Verified byte-identical polynomial
  JSON + HTML vs baseline.
- **E3 — not an issue, no change.** `decimal_points` defaults to `35` consistently across the
  helper functions *and* argparse; the pipeline overrides to `40` explicitly. There is no
  internal inconsistency, and changing the default would alter output for CLI-without-flag
  users. Left as-is.
- **E4 — re-scoped: NOT output-neutral. Decision: keep current code (won't fix).** The HTML
  report loads Plotly from `https://cdn.plot.ly/plotly-latest.min.js`, which is **frozen at
  plotly.js 1.58.5**, while the installed **plotly.py is 6.8.0** (generates figure JSON for
  plotly.js 3.x) — a real version mismatch. Any fix (e.g. `include_plotlyjs='cdn'`, or
  pinning/inlining a version) would **change the rendered HTML report**, so per user decision it
  is left untouched. Report-only; `fitted_corrections.json` is unaffected either way.
- **E5 ✅ RESOLVED by the C1 fix** — `create_rt_vs_rtdiff_plot()` now passes `loess_iterations`
  into its `fit_loess_equation()` call, so the param is no longer vestigial.
- **E6 ✅ FIXED** — vestigial `decimal_points` parameter. It was accepted by
  `fit_spline_equation`, `fit_loess_equation`, `fit_rbf_equation`, and
  `fit_piecewise_polynomial` (and dutifully forwarded by every caller) but **never used** in
  any of those four bodies — only the two polynomial builders format coefficients with it.
  Removed the param from those four signatures/docstrings and stripped it from all 12
  spline/loess/rbf/piecewise call sites (`fit_retention_time_model`, `create_rt_vs_rtdiff_plot`,
  and the main equation-gen loop). `decimal_points` is retained everywhere polynomial. Verified
  output-neutral across all five methods (polynomial/piecewise/loess/rbf byte-identical to
  baseline; spline identical modulo its pre-existing pickle nondeterminism; HTML identical).

---

## D. Verification methodology & results

- Inputs: 5 real `*_feature_data.json` from `work/5a/bdc155…/`, staged as `features_*.json`.
- Ran the **backup original** vs the **cleaned** script with the pipeline's real args
  (`--mz_source start --rt_source start --outlier_threshold 0.95 --one_to_one_only
  --polynomial_degree 8 --spline_degree 3 --spline_smoothing 50 --loess_fraction 0.04
  --loess_iterations 100 --decimal_points 40`), for **every** fitting method.
- **`PYTHONHASHSEED=0` is required** for comparison: the code iterates over `set(common_peps)`,
  whose order is per-process randomized; without a fixed seed even the *original* script
  differs run-to-run (this is itself a latent reproducibility note, not caused by cleanup).
- Compared JSON (minus `timestamp`) and HTML (minus Plotly's random UUID div-ids).

| method | JSON | HTML |
|--------|------|------|
| polynomial | identical | identical |
| loess | identical | identical |
| spline | identical except the base64 pickle blob | identical |
| rbf | identical | identical |
| piecewise | identical | identical |

The spline pickle blob differs **run-to-run in the original script too** (old-vs-old also
differs) — `UnivariateSpline` pickling is nondeterministic. Not caused by this cleanup; noted
as a separate reproducibility concern.

---

# FINAL REVIEW (full-file pass)

A last read of the whole file (2026-07-22), focused on possible errors, then readability &
maintainability. Nothing below is changed yet — this is the report. Items are grouped by
severity. Line numbers are approximate (current file).

## F. Possible errors / latent bugs

- **F1 (low-risk, real) — `LOESSInterpolator._eval_scalar` can divide by zero.** [`~L328`]
  `y0 + (y1-y0)*(x-x0)/(x1-x0)` divides by `x1-x0`. `statsmodels.lowess` does not deduplicate
  its output x's, so two adjacent fitted points with an equal x would give `x1-x0 == 0` → a
  `nan`/`inf` correction that then flows into `fitted_corrections.json` and pairing. Rare (needs
  duplicate fitted RTs) but unguarded. Fix: guard `x1==x0` (return `y0`).
- **F2 (robustness) — `get_coordinate` silently returns `0.0` for a missing key.** [`~L546`]
  `feature.get(key, 0.0)`. If a feature lacks e.g. `rt_start`, it contributes RT=0 → a huge bogus
  `rt_diff` that (if not trimmed) skews the fit, with no warning. Consider raising / logging
  instead of defaulting.
- **F3 (edge) — spline degree is not clamped to ≤5.** [`~L214`] `adjusted_degree =
  min(degree, len-1)` but `UnivariateSpline` requires `k ≤ 5`. Default is now 3 and the .nf
  passes 3, so it's safe in practice, but `--spline_degree 6` would raise (same failure class as
  the just-fixed C2). Fix: `min(degree, 5, len-1)`.
- **F4 (config smell, not a crash) — `match_features_by_closest_mz_and_rt` default
  `rt_weight=0.001` contradicts its own docstring.** [`~L605`] The docstring says it penalizes
  large RT differences, but 0.001 makes the distance essentially m/z-only. Harmless for the
  pipeline (one_to_one_only=True ⇒ one feature per pep per file ⇒ weights never matter), but
  misleading, and `compare_files` never forwards `mz_weight`/`rt_weight` anyway (so they're
  effectively fixed at the defaults). Decide the intended weight and document, or drop the knobs.

## G. Dead code / maintainability (output-neutral to fix)

- **G1 — `fit_retention_time_model` builds `fit_params` but never returns/uses it.** [`~L448`,
  and every branch: 454-455, 464-465, 474-475, 482-483, 490-491]. Also `x_array`/`y_array`
  [`L445-446`] are computed and unused. All dead.
- **G2 — `all_equations` dict** [`L1299`, `L1325`] is populated and never read.
- **G3 — `file2_first_per_pep`** [`L726`] is built and never used.
- **G4 — `file_mapping` is functionally a no-op.** [built `L1098`; looked up `L1239-40/1321/1531-32`]
  It maps `{basename: real_name}`, but every lookup is `file_mapping.get(<real_name>, <real_name>)`
  — and `file_data` is already keyed by real names, so the key is never present and it always
  returns the default. Output is correct (real names), but the mapping layer is confusing dead
  weight; can be removed.
- **G5 — redundant `HAS_STATSMODELS` check** in the loess branch of the dispatcher [`L471`];
  `fit_loess_equation` already checks it [`L239`].
- **G6 — dead defensiveness** at [`L936`]: `equation_text` is always a `str` in every branch,
  so the `isinstance(..., bytes)`/`.decode()` fallback can never fire.
- **G7 — `create_rt_vs_rtdiff_plot` recomputes `sorted(comparisons.keys()).index(file2_name)`
  per iteration** [`L945`] (O(n²)); use `enumerate` over the sorted items. `percentile_used`
  [`L893`] is also unused here.
- **G8 — the report re-fits each model independently** from the JSON pass (main pass 5b vs
  `create_rt_vs_rtdiff_plot`). Same inputs ⇒ same result for deterministic methods, but it's
  duplicated fitting logic that can drift; consider reusing the already-fitted models.

## H. Readability / documentation inaccuracies

- **H1 — HTML still says "(IQR method)".** [`~L1477`, inside `html_content`] It's percentile
  trimming, not IQR. Fixing changes the HTML report bytes (output-affecting) — flagged.
- **H2 — final summary prints "Outlier threshold: {value} seconds".** [`~L1688`] The value is a
  percentile fraction (e.g. 0.95), not seconds. stdout only.
- **H3 — `fit_loess_equation` docstring** [`~L233`] advises `frac` "best set between 0.95-0.99",
  which contradicts the pipeline default (0.04) and LOESS locality (≈1.0 ⇒ near-global). Verify
  and correct the guidance.
- **H4 — spline label hardcodes "Cubic Spline"** [`~L220-221`] even when `k != 3`.
- **H5 — stale type hints / docstrings:** `compare_files` is annotated `-> Tuple[Dict[float,
  float], Dict]` and its docstring says "keyed by m/z", but it returns `(List[Dict], Dict)`
  [`~L670-685`]. `export_fitted_corrections_json` uses `Dict[str, Dict[str, any]]` (builtin
  `any`, should be `typing.Any`) and its docstring still says `[direction]` even though the
  `direction` key was removed (C3) [`~L994-1000`].

## Overall assessment

No high-severity/output-corrupting bug in the pipeline's actual configuration (loess, 1:1,
source=start). The one genuine latent correctness risk is **F1** (div-by-zero → NaN correction);
**F2/F3** are robustness gaps that only bite on malformed data or non-default CLI use. Everything
in **G** is safe, output-neutral cleanup; **H** is docs/labels (H1 touches HTML output).

---

# FINAL REVIEW — RESOLUTIONS APPLIED (2026-07-22)

Per user decisions, the F/G/H items were resolved. Verified: polynomial & loess JSON identical
to the prior baseline (the only diff is the already-approved C3 `direction` removal); HTML
changes to exactly the intended H1 line; F1/F2/F3 behaviours unit-checked; all methods run.

- **F1 ✅ fixed** — `LOESSInterpolator._eval_scalar` now guards `x1 == x0` (returns `y0`) so
  duplicate fitted x's can't produce a NaN/inf correction.
- **F2 ✅ fixed** — `get_coordinate` now **raises `ValueError`** on a missing key or a `0.0`
  value (a real m/z or RT is always > 0), instead of silently defaulting to 0.0. Fails the run
  loudly on malformed input rather than corrupting the fit.
- **F3 ✅ fixed** — spline degree is clamped `min(degree, 5, len-1)` (UnivariateSpline's hard
  limit is 5); `--spline_degree 8` now yields k=5 instead of raising.
- **F4 — kept as intended** — RT is down-weighted deliberately (mirrors
  `params.mmf_mz_rt_weight_ratio = 0.001`: the RT axis spans a much larger, more variable
  range). Docstring corrected to explain this rather than contradict it.
- **G1-G7 ✅ removed** (all verified truly unused): dead `fit_params` + `x_array`/`y_array` in
  the dispatcher; `all_equations`; `file2_first_per_pep`; the no-op `file_mapping` layer (dict,
  its diagnostic print, and all identity `.get()` lookups — replaced with the real names
  directly); the redundant loess `HAS_STATSMODELS` re-check; the dead `bytes/.decode()`
  fallback; the O(n^2) `.index()` (now `enumerate`); and the unused `percentile_used` (plot +
  both main passes). Comments updated to describe only the new state.
- **H1 ✅** removed the false "(IQR method)" from the HTML report. **H2 ✅** summary now prints
  "Outlier threshold: X (fraction of central data kept)" instead of "seconds". **H3 ✅** the
  "best 0.95-0.99 / lower RT shift -> higher" guidance moved from the loess `frac` docstring to
  the `--outlier_threshold` help, where it belongs. **H4 ✅** spline label is now "Spline (k=…)"
  (accurate for any k), not "Cubic Spline". **H5 ✅** fixed stale hints/docstrings
  (`compare_files -> Tuple[List[Dict], Dict]`; `export_...` uses `typing.Any` and drops the
  `[direction]` mention).

**Output impact of this round:** loess/polynomial `fitted_corrections.json` unchanged; the HTML
report loses the "(IQR method)" phrase; spline runs now write "Spline" (not "Cubic Spline") in
the inspection-only `equation` field. F2 will now abort a run on a missing/zero coordinate
(intended).
