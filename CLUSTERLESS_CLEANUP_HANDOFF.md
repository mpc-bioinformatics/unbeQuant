# Clusterless Workflow Cleanup — Handoff Notes

Working context for the ongoing cleanup of the `map_mzml_features_v2_clusterless.nf`
workflow and its helper scripts. Written to survive a context reset. Branch: `alignment_clean`.

## Project goal & scope

Clean up (readability, dead code, correctness) the scripts behind
`map_mzml_features_v2_clusterless.nf`, one at a time:

1. `bin/pair_features_1_clusterless.py` — **DONE** (this file documents it)
2. `bin/compare_retention_time_alignment.py` — **DONE** (see `CRTA_CLEANUP_REPORT.md`)
3. `bin/extract_feature_data.py` — **DONE** (see section below)
4. `bin/build_network_graph_clusterless.py` — not started — **NEXT UP** (biggest file, ~3600 lines)
5. `bin/generate_pairing_report_v2_clusterless.py` — not started
6. `bin/orchestrate_clusterless_cycles.py` — touched (parameter threading) but not fully reviewed
7. `bin/clusters_json_to_consensusxml_clusterless.py` — not started (final quant artifact writer)

**Hard rule from the user:** preserve output. If a change would alter results, flag it and
get explicit approval first. Behaviour-changing fixes are added behind a flag (default =
legacy) OR hardcoded only after A/B validation.

## The workflow call chain

```
map_mzml_features_v2_clusterless.nf
  ├─ extract_feature_data.py            (per file → *_feature_data.pkl/.json)
  ├─ compare_retention_time_alignment.py (→ fitted_corrections.json)   ← NEXT TARGET
  └─ orchestrate_clusterless_cycles.py   (iterative cycles)
        ├─ pair_features_1_clusterless.py        (pairing → edges)
        └─ build_network_graph_clusterless.py    (components + filtering)
     → generate_clusterless_report / clusters_json_to_consensusxml_clusterless.py
```

Cycles: pair → build components → accept clean components → re-pool `unpaired + deleted`
vertices → repeat. Re-pooled vertices carry `openms_fid` (stable identity across cycles).

## Conventions established (follow these)

- **Four-layer parameter threading** for any new knob:
  1. pairing/helper script: `--flag` argparse + function kwarg + metadata record
  2. `orchestrate_clusterless_cycles.py`: argparse + `config` dict + conditional forward
     (`if config.get('x', default) != default:`)
  3. `.nf`: `params.mmf_*` declaration + `def x_arg = params.mmf_x != <default> ? "--flag ..." : ""` + pass-through line
  4. grid script: add as a dimension
- **Groovy truthiness gotcha:** the string `"false"` is truthy in Groovy. Only pass a
  boolean-ish flag to the `.nf` when enabling it; omit otherwise and rely on the `.nf` default.
- **Docker/stale-code caveat:** `.nf` invokes helpers by absolute path
  (`python3 /workspaces/unbeQuant/bin/...`) inside `luxii/unbequant:latest` (workspace
  bind-mounted). Nextflow hashes the `.nf` + params but **NOT** the `bin/*.py` contents, so
  `-resume` can silently reuse results computed with old code. Grid scripts keep
  `NF_RESUME=""` by default. A long-running orchestrator process also holds old code in
  memory while spawning freshly-read pairing subprocesses — caused a real mid-run version-skew
  failure once.
- **Testing pattern:** verified changes on real `*_feature_data.pkl` from `work/` dirs (find
  the newest per sample: `find work -name '*_feature_data.pkl' -printf '%f\t%p\n' | sort -u -k1,1`).
  Scratchpad-based equivalence harness compared old vs new output byte-for-byte.
  `--rt_correction_json` models are pickled `LOESSInterpolator`/`PiecewisePolynomial` that
  unpickle against `__main__`, so test scripts must import those classes at module top.

## `pair_features_1_clusterless.py` — changes made (all verified, byte-equivalent unless noted)

**Removed no-op / vestigial flags** (from pairing only):
- `--rt_source` (pairing side — always used `y_center` anyway; alignment `rt_source` is a
  DIFFERENT, live thing — see below), `--skip-matchfinder`, `--analyze_pep_idents`,
  `--disable_multiprocessing`. `disable_multiprocessing` is still live in the orchestrator/
  `.nf`/`build_network_graph` (real filtering-parallelism control) — only removed from pairing.

**Removed dead features:** `normalize_edge_distances` (fn+flag), `apply_postpair_coordinate_normalization`
(+ `_lookup_edge_features` helper + `--postpair_normalize_coordinates` + its 2 metadata keys).
Neither was ever enabled in the clusterless chain.

**Removed redundant outputs** (edges were serialized 4×; only `output_pkl` is consumed):
- dropped the full `output_json` dump, `edges_edges.pkl`, `edges_edges.json`, and the
  `--output_edges_pkl` + `--skip_json_output` args.
- **Kept (the only consumed files):** `paired_features_edges.pkl` (→ `build_network_graph
  --input_pkl`: reads `network_edges`, `file_index_map`, `cutoff_params`),
  `..._ident_lookup.json` (→ `--ident-lookup`), `paired_features_unpaired.json`
  (→ `--unpaired-json`), and `..._metadata.json` (kept for inspection, unconsumed).

**Correctness fixes (behaviour-changing, done deliberately):**
- **RT-correction curve direction:** hardcoded to prefer the **reverse key** `(query, ref)`,
  which is evaluated on its own fitted domain → EXACT. The forward key `(ref, query)` is an
  off-domain first-order approximation. (Was briefly a toggle, then hardcoded ON after A/B.)
  Impact was real: at `mz_rt_weight_ratio=0.001` this + the k=N fix changed edge counts a lot.
- **k=N "best-passing" matching** (`--max_match_candidates`, nf `mmf_pairing_max_match_candidates`,
  default 1 = legacy): the old code took only the single KD-tree nearest neighbour and dropped
  the pair if it failed the cutoff, even when a valid in-cutoff neighbour existed just behind it
  (severe with small weight ratios, where "nearest" is chosen almost purely on m/z). Now scans
  up to N nearest ascending and keeps the first that passes. `N=10` gave +325% edges on 3 files
  at weight 0.001 — monotonic (only adds). Partially mitigated by the cycle re-pooling, so it's
  more a convergence-speed/quality knob than a pure bug. **Left default N=1 (off).**

**Readability cleanups (byte-identical output):** dropped the always-`[]` `matches` return
(`feature_pairing_optimized` now returns `(edges, file_index_map)`); simplified single-file
branch; direct edge-key access in `extract_and_save_unpaired_vertices`; removed dead
`skip_matchfinder` metadata key (kept `best_match_only`); fixed docstrings/typos;
**`--normalize_coordinates default=False → None`** so the documented "AUTO if omitted" path
(OFF for coord cutoffs, ON otherwise) actually fires — was dead code. Pipeline unaffected
(`.nf` passes it explicitly).

**New `feature_mode` values** (in `bin/extract_feature_data.py`, documented in the `.nf`):
- `first_mean_unweighted_mean` — first isotope, geometric midrange `(min+max)/2`, NO intensity
  weighting (`x/y_center == x/y_center_geo`).
- `first_mean_unweighted_first` — first isotope, earliest-eluting peak (min RT).
- Existing modes (`CoM`, `first_mean`, `single_peak`) verified byte-identical after the change.

## Key domain findings (RT alignment) — directly relevant to the NEXT script

`compare_retention_time_alignment.py` fits the RT-correction curves consumed by pairing.

- **What it fits:** for matched shared-pep_ident feature pairs, `rt_diff = f2_rt - f1_rt`;
  fits `model(x = rt_file1) = rt_diff`. So `model_(A,B)(rt_A) = rt_B - rt_A`. x-axis is FILE 1.
- **Coordinate source:** `get_coordinate(feature, 'rt', source)` maps `'start'→rt_start`,
  `'center'→y_center`, `'geo_center'→y_center_geo` (m/z analogous). The `.nf` passes
  `--mz_source start --rt_source start` (`mmf_rt_alignment_mz_source/rt_source = 'start'`), so
  curves are fitted on **`rt_start`/`mz_start`** (raw elution-window boundaries).
- **JSON keys** `"file1_vs_file2"` → parsed to `(file1, file2)` tuples in pairing. A
  `'direction'` field exists but is **hardcoded `'forward'`** (carries no info). All 210 = 15×14
  ordered pairs are emitted (both directions present).
- **Two domain mismatches** (pairing evaluates the curve at `y_center`, an intensity-weighted
  centroid computed from `feature_ls`; `rt_start` is `min(raw rt boundaries)`, computed OUTSIDE
  `feature_ls`, so `feature_mode` never affects `rt_start`):
  1. **Field mismatch** — fitted on `rt_start`, evaluated at `y_center` (median **8.9 s** gap on
     real data). Fixable by setting `mmf_rt_alignment_rt_source='center'` (fits on `y_center`,
     same field pairing uses). NOTE: `center` also changes `mz_source`→`x_center`, which changes
     *which* pairs are matched during fitting, so it's not a pure domain fix.
  2. **Direction mismatch** — forward-key off-domain evaluation (fixed on the pairing side, above).
- Both are ~1% of the drift (drifts 23–87 s, errors 0.1–4.4 s), well inside `rt_cutoff=100`.

## Grid script: `run_parameter_grid_feature_modes.sh`

- 20 runs = 5 `FEATURE_MODES` × 2 `COORD_SOURCES` (`start`/`center`, sets mz+rt source in
  lockstep) × 2 `MAX_MATCH_CANDIDATES` (`1`=nofix / `10`=fix). Run names encode all three.
- `RUN_STEP2`/`RUN_STEP3` toggles (default **off** — only step 1, the pairing/mapping, runs;
  steps 2/3 are consensusXML + ProteoBench and were failing at step 2 independently).
- Disk pressure is real: `work/` + published results filled the volume once. `cleanup=false`,
  no `-resume`, so every run makes fresh work dirs. Prune `work/` between runs.

## `compare_retention_time_alignment.py` — DONE

See `CRTA_CLEANUP_REPORT.md` for the full write-up. Summary of output-affecting changes now
LIVE (the next pipeline run's `fitted_corrections.json` will differ accordingly):
- **C1**: `--loess_iterations` now actually reaches `lowess(it=)` (was silently ignored;
  `lowess` hard-coded `it=3`) → changes loess corrections. NOTE the `.nf` default was later
  set to `mmf_rt_alignment_loess_iterations = 10`.
- **C2**: spline defaults 8/1000 → 3/100 (old degree 8 broke `UnivariateSpline`, which needs k≤5).
- **C3**: dropped the inert `direction` field from the JSON (only `method`+`serialization` are consumed).
- **F1/F2/F3**: LOESSInterpolator div-by-zero guard; `get_coordinate` now RAISES on missing/0.0
  coord (aborts on malformed input) instead of silently defaulting; spline degree clamped ≤5.
- Removed dead code (fit_params, all_equations, the no-op file_mapping layer, `--show-plots`,
  `--decimal_points`), report size 133MB→17MB. File ~1668 lines.
- Verification gotcha: output only deterministic under `PYTHONHASHSEED=0` (iterates
  `set(common_peps)`); spline's base64 pickle blob is nondeterministic even old-vs-old.

### PERF FIX (later chat): every RT model was fitted THREE times per file pair

Found by `/code-review`. Each ordered pair was fitted three separate times on the SAME points:
Pass 5a (to build an equation string, model discarded), Pass 5b (for the JSON), and again inside
`create_rt_vs_rtdiff_plot` (for the report curve). With 25 files that is 600 pairs × 3 = **1800
lowess fits**, at `it=10` since C1 made the flag live — the reason this process needed 4G → 24G.

Now ONE fit per pair, reused for all three consumers. Changes:
- **`format_polynomial_equation(coeffs, degree, decimals, html=)`** extracted — rendering split
  from fitting. `fit_polynomial_equation_plain` was a byte-for-byte copy of
  `fit_polynomial_equation` differing only in `x<sup>2</sup>` vs `*x^2`; it refit the data to
  change a string. Deleted.
- **`fit_retention_time_model`'s polynomial branch** now builds `plain_eq` from the coefficients
  it just fitted. It previously did `plain_eq = html_eq`, i.e. put HTML markup in the JSON's
  "plain" field — reachable whenever Pass 5a had no entry for the pair.
- **Passes 5a and 5b merged.** Verified they computed IDENTICAL filtered inputs (same threshold
  from the same `abs(d['y'])` list, same predicate, same iteration order), so 5a existed only to
  recompute 5b's fit and keep a string.
- **`create_rt_vs_rtdiff_plot` no longer fits** — takes `fitted_models`
  `{file2: (model, html_eq, degree)}`. Six now-unused params and the dead `equations_plain`
  return went with it (the caller already discarded it as `_`).

**Verified on the real 25-file input at production params, `PYTHONHASHSEED=0`:**
- `fitted_corrections.json` — all **600 correction records identical in every field**, including
  the base64 pickled models and equation strings. Only the top-level `timestamp` differs, which
  `datetime.now()` makes differ between ANY two runs. ⚠️ **md5 of this file is therefore useless
  as an equivalence check** — compare `json.load(...)['corrections']` instead.
- Also verified with `--fitting_method polynomial` (6-file subset): all records identical.
- Speed, 6-file subset (30 pairs): **loess 85s → 37s (2.3×)**, polynomial 16s → 14s. Since
  loading/matching is unchanged, the fitting portion improved by considerably more.

⚠️ **The HTML report's fitted curve changes by ~1e-11 s** (max relative 1e-9, on values of scale
~7 s). Cause, proven empirically: `fit_loess_equation` sorts with `np.argsort`, whose default
quicksort is **NOT stable**. Real RT values contain ties, so sorting an already-sorted array (what
the plot used to pass) vs an unsorted one orders tied points differently and lowess sees a
different `y` order. A no-ties control is bit-identical. This is a convergence, not a regression:
the report previously drew a SEPARATE fit from the one pairing actually applies, and now draws the
same model that is in the JSON.

Unrelated observation while verifying: **`rt_alignment_report.html` is 1.2 GB** on 25 files (it
embeds every data point). The earlier "133MB→17MB" note was measured on a much smaller input.

## `extract_feature_data.py` — DONE (this chat)

Per-file feature TSV → `*_feature_data.pkl/.json`. First stage of the workflow. ~440 lines now.

**Output-affecting changes (LIVE — re-extract before next run):**
- 🔴 **CRITICAL intensity bug FIXED.** A loop variable `intensity` shadowed the feature's TSV
  intensity, so the emitted `intensity` field held the *last peak's* intensity for the weighted
  modes (`CoM`, `first_mean`, `single_peak`) — wrong for 1200/1200 features, median ~1/25th of
  true, some 0 (→ downstream div-by-zero). Fixed by renaming the loop var to `peak_intensity`.
  `x_center`/`y_center` were always correct, so pairing/edges/components are UNCHANGED — only the
  `intensity` field (→ `ident_lookup` → consensusXML `it=` → ProteoBench) changes. The two
  `*_unweighted_*` modes were already correct.
- Removed the `isotope_dist_mean` output field (computed every row, read by nobody). Only field
  removed; every other field byte-identical across all 5 modes.

**Loud-failure changes (behaviour on BAD input only; valid input unaffected):**
- All 4 parsers (`parse_list`/`parse_list_flexible`/`parse_nested_list`/`parse_scalar`) now ABORT
  with column+row+value on an unparseable cell (was: silent `[]`/`None`). Fixes the NaN-charge
  crash (empty `charge` cell → clean error, not a raw traceback).
- `--feature_mode` now `required` with `choices=[...]`; unknown/missing mode errors instead of
  silently running `CoM`. (Removed the `else: # CoM` fallback in `calculate_feature_center`.)
- Column validation (`_validate_feature_tsv`): all 12 read columns mandatory + two either/or
  pairs (`id`|`openms_fid`, `ms2_scans`|`l_ms2_scans`). All 6 sample TSVs carry all 14 cols, so
  no valid input breaks. On a full 56k-row file: **0 rows skipped** — the skip paths are purely
  defensive.
- `skipped_count` split into 3 reported reasons (no bounding box / no usable peaks / outside RT window).

**Readability:** merged `reduce_feature_to_point` + the row-loop centre math into ONE function
`calculate_feature_center` with a single mode-dispatch block (mode → peaks → weighted? in one
place); deleted the `UNWEIGHTED_FEATURE_MODES` / `FIRST_ISOTOPE_MODES` module sets (user dislikes
them); removed dead numpy import, `--mzml`/`--spectrum_pkl`/`--round_up_to`/`--generate_diagnostic`
args. **Removed `--round_up_to` + `--generate_diagnostic` from ALL THREE calling workflows**
(`map_mzml_features_v2_clusterless.nf`, `map_mzml_features.nf`, `map_mzml_features_v2.nf`) so they
don't pass now-unknown args. (The `--round_up_to` still in v1/v2 belongs to `process_mzml_file.py`,
a different script — leave it.)

**Left as-is / flagged, not changed:** `feature_count` output field is misnamed (it's the PEAK
count, not a feature count) and is misused as an index fallback in the LEGACY `pair_features.py`/
`pair_features_old.py` (not the clusterless chain). Kept because those legacy pipelines read the
same pkl; rename to `n_peaks` or remove is a pending user decision.

Weighted-mean formula was CHALLENGED and PROVEN correct on live data (identical to textbook
`ΣcI/ΣI` to float epsilon; the `geo` seed cancels). The in-code "weighted mean anders berechnen"
TODO was unfounded and removed.

## `build_network_graph_clusterless.py` — IN PROGRESS

### Fix 1 (DONE, output-affecting): non-reproducible component order

`export_accepted_components` iterated `component_info` in dict-insertion order.
`filtered_cluster_info` is filled in two phases: no-dup components in sorted `comp_idx`
order (PHASE 1), then filtered components in **`imap_unordered` worker-completion order**.
So the `components` array order was race-dependent, and the orchestrator assigns global
`component_N` IDs **by position in that array** → global component IDs shuffled run to run.

- **Fix:** iterate `sorted(component_info.items(), key=lambda item: item[0])`. One line, at the
  serialization boundary so it covers the `--filter-duplicate-file-vertices`-disabled branch too.
  Keys are int `comp_idx` in both branches.
- **Verified:** cycle 20 unmodified → **2 distinct outputs across 6 runs** (4×/2× split); post-fix
  → 6/6 byte-identical, and equal to the pre-computed normalized baseline hash
  `fec4ed47dd8ab52ba540f1b5b33f9374`. Content proven identical pre/post (same 15 components;
  cycle 5: same 689). `deleted_vertices.json` byte-identical pre/post (already sorted at the
  `deleted_vertices_list.sort(...)` line). Parallel and serial (`--disable-multiprocessing`) now
  also agree byte-for-byte (cycle 5 → `7451f124...`).
- **Race is load-dependent** — cycle 5 (479 filtered components, big worker chunks) was stable
  across 6 runs while cycle 20 (12 components, chunks finishing near-simultaneously) flipped.
  "Looked reproducible" was luck. Always test on a SMALL/late cycle, not a big early one.
- Downstream bytes will differ from previous runs (component IDs are now canonical, not
  race-dependent); component *contents* are unchanged.

### Fix 2 (DONE, byte-identical): dead-code removal — 3658 → 2360 lines (−1298, 35.5%)

Removed, all verified byte-identical on the full harness matrix (36 output files):
- **17 never-called functions (1136 lines)** — the leftover recovery/singleton/cluster-export
  workflow from the clustered version (`recover_deleted_vertices_edge_based`,
  `save_deleted_vertices` (209 lines; signature said 2 return values, returned 5),
  `create_singleton_clusters_for_un{recovered_deleted,paired}`, `load_unpaired_vertices_json`,
  `export_{clusters,filtered_clusters}_to_json`, `save_deleted_vertices_report`,
  `calculate_file_linkage_scores`, `detect_*`, `create_cluster_size_distribution`,
  `_compute_cluster_edge_weights`, both `_batch_components_by_complexity*`, and
  `aggregate_cluster_{pep,prot}_idents` which were reachable only from the above).
- **The orphaned code block** (docstring + body, no `def`, referenced undefined
  `components_batch`) — only parsed because it sat after `return results` inside
  `_lazy_worker_with_self_batching`.
- **`_get_actual_available_cpus` (66 lines)** — unreachable cgroup/container CPU detection; main
  always passes `os.cpu_count() or 1`. The dead `num_cores is None` branch was collapsed to a
  plain `os.cpu_count() or 1` fallback (kept, so the function stays safe if called with None).
- **The PHASE 0 `edge_index`** and its module global — built O(E), forked to every worker,
  read by nobody (workers use `graph.neighbors()`).
- Unused imports: `shutil`, `math`, `random`, `get_context`, `minimize_scalar`, `pdb ... as bp`.

`pyflakes`: **0 undefined names** after removal (remaining warnings are pre-existing noise —
f-strings without placeholders, unused locals).

**Measured perf win** (all from the `edge_index` removal; the rest was never executed).
Cycle 3 = 694,977 edges, 3 runs each, non-overlapping ranges:
filtering phase 2.48–2.57s → **1.53–1.64s** (−37%); IPC overhead 1.67–1.71s → 1.35–1.45s
(the dict was being forked to every worker); total 6.59–6.96s → **5.70–5.79s** (−15%).
Scales linearly at **~1.23 us/edge**, so a big early cycle (cycle 1 pairs 1.37M features over
25 files) saves tens of seconds plus a multi-GB structure. NOTE: on a small late cycle
(cycle 20, 50k edges) the difference is pure noise — measure perf on cycle 3+, not cycle 20.

### Fix 3 (DONE, byte-identical): small correctness oddities

All verified byte-identical on the full harness matrix (36 files); the two unreachable-path fixes
were unit-tested directly since the harness can't reach them.
- **`start_time` moved above the `try`** in `_filter_component_worker`. The `except` handler
  reports elapsed time, so an early failure raised `UnboundLocalError` *inside the handler* and
  masked the real error. Now returns a clean error record (tested: reports the real `KeyError`).
- **`scoring_path_counts` rewritten as an accumulating dict.** The old fixed key list could never
  count the fallback paths (they embed `type(e).__name__` in the name) and had two phantom keys
  (`fallback_from_exception`, `modularity_fallback`) plus one mismatched
  (`modularity_graph_neighbors_fallback` vs `modularity_no_graph_fallback`).
- **`float(best_modularity) if best_modularity else None` → `... is not None else None`.** A real
  score of 0 was being written as `None`. Field is unconsumed downstream, so output is unaffected.
- **Removed the verbatim copy-pasted complexity-estimation preamble** in
  `_generate_filter_candidates` (5 lines duplicated), and the never-read `num_duplicate_clusters`.
- **`mega_complexity_threshold` docs corrected** — said "default: 1e15"; the real production value
  is **5000**, set by `params.mmf_mega_complexity_threshold` (`map_mzml_features_v2_clusterless.nf:73`)
  and threaded through the orchestrator. Both the docstring and the argparse help now say so.

### Fix 4 (DONE, byte-identical): one exact complexity formula, `count_filter_candidates()`

There were TWO formulas. PHASE 1.5 predicted the heuristic/exhaustive choice with an
**approximation** (`int(mean duplicate-set size) ** number_of_dup_sets`, capped 1e15); the real
decision in `_generate_filter_candidates` used the **exact product** of duplicate-set sizes.
They drifted: at threshold 5000, cycle 3 **44/1846 (2.4%)**, cycle 5 **6/479 (1.3%)**, cycle 20
0/12 — always one-directional (the approximation *underestimates*, predicting exhaustive when
reality goes heuristic, because `int()` truncates the mean). Real example: dup-set sizes
`[4,3,3,3,3,2,2,2,2]` → mean 2.67 truncates to 2 → `2**9 = 512`, vs the true `5,184`.

**Both sites now call one module-level helper** `count_filter_candidates(clusters_with_duplicates)`
returning the exact product, so the prediction is a guarantee rather than a guess.
- The old running-product check is equivalent to a single final comparison (every factor is >= 2,
  so the product only grows). Proven by randomized equivalence test: 3000 shapes x 4 thresholds.
- PHASE 1.5 compares the **UNCLAMPED int** against the threshold, then clamps with
  `float(min(count, 10**15))` for the load-balancing weight. Clamping before comparing would
  misclassify huge components whenever the threshold is set above the clamp; converting an
  unclamped product to float can raise OverflowError (products of 1e50+ are real here).
- Verified on cycle 3: predicted split was 328 heur / 1518 exh vs an ACTUAL 370 / 1476 (off by
  42); it is now **370 / 1476, an exact match**. Output byte-identical (36/36) — Fix 1 made the
  output invariant to work distribution, so only weights and log text moved.

### Fix 5 (DONE, byte-identical): B — removed the no-op cutoff machinery

Edge filtering here could never drop an edge (see the finding below), so the whole layer went:
- Graph script: `--edge_cutoff`, `--mz_cutoff`, `--rt_cutoff`, the `cutoff_params` reconciliation
  block in main, and `add_cutoff_suffix_to_filename` (with the cutoffs gone it could only ever
  return the filename unchanged — `fraction` was always None). `build_network_graph()` no longer
  takes cutoffs; `load_edges_and_cutoff_params_from_paired_features` → `load_edges_from_paired_features`
  returning `(edges, file_index_map)`.
- Orchestrator: dropped the three forwards to the graph step. `config['edge_cutoff']` was never
  set by anything, so that branch could never fire either.
- **Untouched:** the pairing cutoffs (`--ppm` / `--rt_cutoff`, `.nf` lines ~229-231). They are the
  only place cutoffs are enforced and they stay.
- Verified: dropping `--rt_cutoff` from the invocation was ALREADY byte-identical before any code
  changed (isolation check), then 36/36 byte-identical after.

### Fix 6 (DONE): C1 + C2 — removed `filter_method` from the clusterless chain entirely

`--filter-method` selected between two modularity variants that are unreachable here (main builds
one cluster per component, so `len(clusters_with_duplicates)` is always 1 — measured: cycle 3
`{1: 1846}`, cycle 5 `{1: 479}`, cycle 20 `{1: 12}`). Production had it set to `'modularity'`, so
it had been requesting a path that never ran.
- Graph script: ~100 lines of modularity scoring, the arg, the `filter_method` module global and
  its whole threading chain (signature → global → work_item → worker), and `scoring_path_counts`
  (which only ever tracked the removed branches). Scoring collapsed to the one live line:
  `score = sum(graph.degree(v) for v in test_vertices)`.
- Orchestrator: argparse, config entry, forward, log line. `.nf`: `params.mmf_filter_method` and
  both its uses.
- consensusXML writer: the `--mmf_filter_method` arg and the `<UserParam name="mmf_filter_method">`
  line. **This is an output change** — verified as EXACTLY one line removed from a 225MB
  consensusXML (159,621 components / 1,372,578 vertices otherwise identical).
- Graph script output: 36/36 byte-identical.
- **Legacy pipelines keep `filter_method`** (`map_mzml_features.nf`, `map_mzml_features_v2.nf`,
  `bin/build_network_graph.py`, `bin/clusters_json_to_consensusxml.py`) — they do real clustering,
  so it genuinely applies there. Do NOT propagate this removal to them.

### Fix 7 (DONE, byte-identical): the two selection methods are documented in-place

`SELECTION METHOD 1 of 2: EXHAUSTIVE (JOINT OPTIMUM)` above the candidate loop in
`_filter_component_worker`, and `SELECTION METHOD 2 of 2: HEURISTIC (GREEDY PER DUPLICATE SET)`
above the `if use_heuristic:` block in `_generate_filter_candidates`. They sit ~1000 lines apart
so each cross-references the other. Key facts captured there:
- Exhaustive scores by TOTAL DEGREE of the retained set and compares all combinations → JOINT
  optimum. Ties keep the first candidate (strict `>`); order is deterministic (sorted iteration).
- Heuristic picks per duplicate set by `internal_degree - external_degree`, but **external_degree
  is structurally always 0** — main sets `'vertices'` and `'clusters' {0: ...}` to the SAME list,
  so `other_cluster_vertices` is always empty. It reduces to "highest-degree vertex per file",
  decided greedily per set.
- This is why `mega_complexity_threshold` is OUTPUT-AFFECTING: raising it moves components from
  the greedy choice to the joint one. Measured on cycle 3 at 5000 vs 100000: heuristic/exhaustive
  370/1476 → 293/1553, deleted count identical (31,717) but a DIFFERENT hash — same number of
  vertices removed, different ones.

### Fix 8 (DONE, byte-identical): D and F marked as optional diagnostic output

User decision: keep both, they are useful for manual analysis. Measured footprint ~283 of 2260
lines (13%), so keeping them is cheap.
- `analyze_graph`, `analyze_graph_composition`, `print_analysis` docstrings now say DIAGNOSTIC
  OUTPUT ONLY / not read by any downstream step, plus a COST NOTE that `diameter()` on the largest
  component is the expensive part and the reason `--skip-analysis` exists. Same marking on the
  main() block and on the `--output_analysis` / `--output_composition` argparse help.
- Deleted-vertex record marked as diagnostic bookkeeping, listing the only 3 consumed keys.
  **CORRECTION to the earlier assumption:** those diagnostic fields never reach disk.
  `export_deleted_vertices_for_next_cycle` takes only `vertex_id` and re-reads every attribute
  from `vertex_attrs`, so `reason` and `best_modularity` — the two fields that exist ONLY in that
  record — are currently discarded. Noted in the code. Surface them in the export if wanted.

### Fix 9 (DONE, byte-identical): G — removed the log lines that could never fire

`mega_components` / `large_components` compared `combo_estimate` against
`COMBINATIONS_BATCH_LIMIT * 10` and `* 1`, but a component's weight is capped by
`mega_complexity_threshold`, so both always printed 0 (confirmed across all 11 filtering runs) and
their ">10M"/">1M" labels matched neither threshold. The Min/Max/Avg complexity stats are kept.

### Fix 10 (DONE, byte-identical): E — one `--num_workers` knob, threaded to the .nf

Four independent things used to decide parallelism: `os.cpu_count()`, a hard cap of 10 "cores",
`num_workers = min(max(4, 4), cores_to_use)` (i.e. `min(4, …)`, with a comment claiming 8), and a
CPU-pinning modulus taken from yet another variable. All replaced by one function:

`resolve_worker_count(requested) -> (num_workers, cpu_list)` — the single source of truth.
- `cpu_list` comes from `sched_getaffinity`, so it respects cgroup/container limits (the old
  `_init_worker` pinned to `worker_index % num_cpus`, a modulus of a plain count that could point
  at a CPU the container does not own; it now indexes into the real affinity list).
- Default `DEFAULT_NUM_WORKERS = 1`, and a count of 1 SKIPS the Pool entirely and runs in-process
  (same path as `--disable-multiprocessing`). Parallelism is opt-in. This is deliberate: a
  1-worker Pool still pays fork + IPC with nothing to overlap and measured ~23% SLOWER than
  in-process (3.24s vs 2.64s on cycle 3). After the change, default / `--num_workers 1` /
  `--disable-multiprocessing` all land on the same path (2.61 / 2.63 / 2.51s = noise).
- Requests above the available CPU count are clamped, with a printed note.
- Threaded 4 layers per the project convention: `--num_workers` (graph script) →
  `--num_workers` (orchestrator argparse + config + conditional forward) → `params.mmf_num_workers`
  (`.nf`, null = script default, `workers_arg` builder + pass-through).
- `--disable-multiprocessing` still wins over `--num_workers` (orchestrator uses `elif`).

Verified: 36/36 byte-identical; `nextflow inspect` exit 0; and default / 1 / 8 / 999 workers all
produce the SAME md5 on cycle 3 (999 clamps to 16 with the message). Output really is invariant
to parallelism.

**Measured on cycle 3 (694,977 edges / 3,499 components), filtering phase, best of 2:**
| workers | 1 (default, in-process) | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| time | 2.61s | 1.90s | 1.58s | **1.34s** | 1.50s |
Diminishing past 4, best at 8, regression at 16 (pinning contention + IPC). So ~2x is available
on this cycle by setting `mmf_num_workers = 8` AND `mmf_disable_multiprocessing = false`.

### ⚠ NOTE: parallel filtering is currently OFF in production

`params.mmf_disable_multiprocessing = true` in the .nf means the flag IS passed, so filtering runs
single-process and **`mmf_num_workers` has no effect until that is set to false**. The inline
comment there reads "Enable parallel processing for filtering", which is inverted - the param is
named after the flag, not the feature. Left as-is deliberately (user declined the comment change).

### Fix 11 (DONE, byte-identical): diagnostic output compressed, minor lint cleared

Logic was cleaned through Fix 10 but the INSTRUMENTATION never got the same treatment.
Stdout on cycle 3 (1846 filtered components): **13,852 -> 267 lines (-98%)**. File 2335 -> 2141.
- **Five overlapping report blocks collapsed into one** (120 -> 33 lines): FILTER TIMING SUMMARY,
  PROCESS POOL UTILIZATION DIAGNOSTICS, PER-WORKER BATCH TIMING ANALYSIS, IPC OVERHEAD ANALYSIS
  and Task timing all re-derived the same handful of timestamps. Now one block: wall clock,
  fork/workers/collect split, per-process counts, avg+max per component, imbalance warning.
- **Per-worker summary 34 -> 13 lines** (one line per worker instead of an 11-line table each).
- **Per-component logging sampled**, `LOG_EVERY_N_COMPONENTS = 250` with a per-worker-process
  counter `_components_seen`. Was ~5 lines x N components. Anomalies (AFFINITY CHANGED, graph is
  None) are ALWAYS logged regardless of cadence — only the "affinity held" confirmation is sampled.
- Removed per-batch chatter, phase banners, "Freed X to recover memory" lines, the duplicate
  per-worker chunk announcements, and a stale "Details saved in output JSON report" (that report
  was deleted back in Fix 2).
- Minor lint cleared: unused locals (`components_to_process`, `comp_vertex_set`, `batch_size`),
  unused `traceback` import, redundant `global graph`, and 27 f-strings without placeholders.
  **pyflakes is now completely clean (0 findings).**

⚠️ **Gotcha for whoever does this next:** stripping the `f` prefix with a regex is easy to get
wrong. My first attempt used nested capture groups (`r'\bf(("...")))'`) and emitted group1+group2,
which are the SAME text — so it silently DUPLICATED 13 string literals. Python concatenates
adjacent literals, so those lines printed twice and no test caught it (stdout only). Detect with
`re.findall(r'("(?:[^"\\\n]|\\.)*")\1', src)`; fix by substituting `\1`.

### Fix 12 (DONE, ⚠ ONE OUTPUT CHANGE): `_generate_filter_candidates` flattened

The function took `{cluster_id: {filename: [vertex_ids]}}` and cross-multiplied combinations
across clusters. The caller only ever built ONE synthetic cluster `{0: ...}` — it is clusterless —
so the whole cluster dimension was ceremony. Removing it: **177 → 81 lines, nesting depth 10 → 2**.
Both callers now build `file_vertex_map` directly from `comp_vertices`. Renamed
`_generate_filter_candidates` → `generate_filter_candidates` (it is not private to a worker).

Two collapses that followed from "there is only one cluster":
- Heuristic score was `internal_degree − external_degree`. With one cluster per component every
  neighbour of a member IS a member, so `external_degree` is structurally 0 → plain `degree`.
- `heuristic_complexity` multiplied by a cluster count that was always 1.

**Equivalence proven three ways**, after the first harness run FAILED and forced a real diagnosis:
1. **Algorithm level.** Reimplemented the OLD generator verbatim (synthetic cluster, per-cluster
   combos, cross-multiply, old heuristic) and compared the ENTIRE ORDERED candidate sequence plus
   the winner, per component: cycle 20 @5000 (12 comps), cycle 20 @1 (12, all heuristic), cycle 30
   @5000 (5), cycle 5 @5000 (479), cycle 3 @5000 (1846 = 370 heuristic + 1476 exhaustive).
   **0 sequence mismatches, 0 selected-set mismatches**, `(combo_estimate, use_heuristic)` equal
   on every one. Both branches covered, 2354 components total.
2. **Output level.** Full-payload comparison keyed on `(filename, feature_idx, openms_fid)` —
   NOT on labels — across 5680 components / 60174 vertices: every field identical, no label
   collisions, `num_vertices` identical.
3. `deleted_vertices.json` byte-identical to baseline in all 13 runs (accepted+deleted partition
   the vertex set, so this is independent confirmation).

⚠️ **THE OUTPUT CHANGE — needs a decision if not already taken.** The first harness run after
flattening failed on `accepted_components.json`. Membership was identical; the VERTEX ORDER inside
each component had changed. Cause: `list(best_candidate)` on a `set` — CPython set iteration order
depends on how the set was built, so the emitted order silently tracked the candidate-construction
sequence. Fixed with `sorted(best_candidate)`, making it canonical (ascending vertex id) and
independent of construction — same class of fix as Fix 1, and it makes the file MORE stable, not
less. But `accepted_components.json` bytes now differ from pre-Fix-12 runs, and the consensusXML's
`groupedElementList` ordering follows. Options: keep the sort (recommended — membership is provably
unchanged so quantification is unaffected), or revert both the sort and the flattening.

**Set-equality caveat, so nobody over-reads test 1:** Python `set == set` is content-based, so the
sequence test proves MEMBERSHIP equality, not iteration order. The ordering difference is a
separate consequence, diagnosed and pinned by `sorted()` above.

### Fix 13 (DONE, byte-identical): restructured into titled sections + naming pass

Pure readability. 2014 → 2053 lines (net +39: comments added, dead assignments removed).
**All 36 harness outputs byte-identical; stdout identical too** (diffed with timing lines filtered
— only the run-directory name in `✓ Saved analysis:` paths differs). pyflakes still 0.

**Eight numbered sections**, with a header block each explaining what it is and how it fits:
1 tuning constants · 2 shared state for worker processes (the fork/COW contract) · 3 input ·
4 graph construction · 5 diagnostic analysis · 6 duplicate-vertex filtering (sub-numbered 6.1–6.5,
ordered bottom-up: parallelism setup → candidate generation → per-component → chunk driver →
orchestrator) · 7 output · 8 CLI entry point. Module docstring now states the pipeline position,
the determinism contract (the two sorts), and the vestigial cluster layer.

**Renames (nothing outside this file referenced them — `build_network_graph_multithread_v2.py` is
the separate non-clusterless variant with its own copies):**

| old | new | why |
|---|---|---|
| `best_modularity` | `best_score` | it is total degree, not modularity — C1 removed modularity |
| `'best_modularity'` (dict key) | `'best_score'` | same; in-memory only, never reaches disk |
| `metadata['combo_estimate']` | `metadata['load_weight']` | held the load-balancing WEIGHT, not a combination count — the two diverge for heuristic components, and `'actual_combo_estimate'` sat next to it holding the real one |
| `_filter_component_worker` | `filter_one_component` | runs in-process too; "worker" was misleading |
| `_lazy_worker_with_self_batching` | `process_component_chunk` | says what it does |
| `_generate_filter_candidates` | `generate_filter_candidates` | not worker-private (Fix 12) |
| `test_vertices` | `kept_vertices` | they are the survivors, nothing is being tested |
| `filtered_vertices_all` | `kept_vertices_all` | consistency |
| `components_to_estimate` | `components_needing_filter` | estimation is not what the list is for |
| `worker_pids` / `task_times` | `components_per_pid` / `component_times_ms` | contents + units |
| `chunk_complexities` | `chunk_weights` | matches `load_weight` |
| `save_graph(format=)` | `save_graph(graph_format=)` | shadowed the builtin |
| `pass_idx` | `scanned` | it is a progress counter |
| `v_dict` | `deleted` | says what it holds |

**Dead assignments removed** (verified unread, no output effect): `metadata['actual_combo_estimate']`
(set "for logging", the log used the other key); `will_use_heuristic`, `combo_estimate` and
`__estimated_size_bytes` on the WORK ITEM (the worker recomputes both flags; the size estimate is
still computed for the LARGE COMPONENT warning, just not stored); `feature_ids`/`feature_id` in
`analyze_graph_composition` (collected then discarded, and `'feature_id'` is not even a vertex attr
— it is `'feature_idx'`, so the `.get` always hit its default); the duplicated init block and the
half-used `worker_timing` dict in the chunk driver; `filtered_edges = edges` and `edges = edges_full`
aliases; redundant function-local `import os/sys/time`; the `try/except NameError` around the memory
`del`s (`ident_lookup` is always bound).

Also refactored the two batch bodies in `process_component_chunk` into one `run_batch()` closure —
they were copy-pasted, ~15 lines each.

**Two things deliberately NOT changed, flagged instead:**
- **The two edge-parsing passes in `build_network_graph` are asymmetric.** Pass 1 requires BOTH
  file indices in `file_index_map` before using it (else synthesises `file_N` names for both) and
  `continue`s on an unrecognised edge; pass 2 does a per-index `.get` with fallback and treats
  "neither `ref_file_idx` nor `file1`" as the legacy nested format, so it would `KeyError` where
  pass 1 skipped. On well-formed input they agree exactly. Unifying them is a behaviour change on
  malformed input, so it is documented in a comment and left alone.
- **The `cluster_id` layer** (`vertex_to_cluster_map`, `component_cluster_info[c]['clusters']`,
  the `cluster_id` in every deleted-vertex record). Always `(comp_id, 0)`. The names correctly
  describe the data as it stands, so renaming without removing would make things worse. Removing it
  is a real structural refactor — it reaches the deleted-list sort key — and belongs on its own.
  Documented in the module docstring and the SECTION 6 header.
- Component-level `'pep_ident'`/`'prot_ident'` in `accepted_components.json` are never set by this
  script, so they always serialise as `null`. Kept (the consumer expects the keys) and commented.

### Fix 14 (DONE, byte-identical): collapsed four representations of the partition into one

The script derived `graph.components()` → `component_info` → `component_cluster_info` →
`vertex_to_cluster_map`, then a fifth copy `filtered_vertex_map`. All five described the same
partition. 2053 → 1963 lines, and a real memory win: two of them held **one entry per vertex**,
plus a per-component slice copied into every work item.

What the tracing showed (all verified by grep before removing anything):
- `components[i]` is **already a plain Python list** — picklable, enumerable. Nothing needed
  "extracting" from igraph at all; `VertexClustering` is not an opaque handle.
- `component_info[i]['component_id']` and `['num_vertices']` were **never read**. The export takes
  `component_id` from the dict KEY and recomputes the count.
- The `component_cluster_info` module GLOBAL was **never read** — published for fork-sharing,
  `global`-declared and `del`'d in the worker, but nothing ever consumed it. Pure waste in every
  worker's address space.
- `vertex_to_cluster_map[v]` returned `(comp_idx, 0)` — and `comp_idx` is already a parameter of
  the function doing the lookup. A per-vertex dict answering a question the caller had answered.
- `filtered_vertex_map` was a full per-vertex dict copy whose ONLY consumer was
  `len(filtered_vertex_map)`, i.e. `graph.vcount() - len(deleted_vertices_list)`.

Root cause: the cluster dimension is inherited from the CLUSTERED variant
(`build_network_graph_multithread_v2.py`), where Leiden makes `cluster_id` vary and the map is
genuinely needed. Clusterless pinned it to 0 and the whole layer became derivable.

Now: `filter_duplicate_file_vertices(graph, components, vertex_attrs, ...)` returns
`(accepted_components, deleted_vertices_list)` where `accepted_components` is
`{comp_idx: {'vertices': [...]}}` covering every component. `cluster_id` is gone from the deleted
record and from the sort key (it was a constant wedged in the middle of a 3-tuple).

**Verification: all 36 harness outputs byte-identical.** stdout identical across all 13 runs except
the `⚠ load imbalance` line — which is a wall-clock ratio crossing a 5x threshold and **varies
between logically identical builds** (15x vs 16x on cycle 3 for two same-logic runs; on cycle 5 it
appeared in one run and not the other). Not a behaviour change.

⚠️ The `LARGE COMPONENT DETECTED` log lost its `clusters=` field. Confirmed it never fires on any
harness run (it needs ~460k vertices in one component), so stdout is unaffected in practice.

### DOCUMENTATION RULE (user instruction, applies to every script)

**Never reference removed or previous code in the script itself** — not in comments, not in
docstrings. "This used to be…", "the old formula…", "X was removed because…" are all out. Users
read the script to learn how it works NOW; history is noise there. Rationale that survives without
the history (e.g. *why* sorting is required) stays, phrased in the present tense. The
before/after record belongs in THIS handoff, not in the source.

A sweep applied this across the whole script: the `resolve_worker_count` docstring, both selection
methods, the per-worker summary, the load-weight note, the report block, the sort key, the module
docstring, and the component-extraction comment.

### OPEN: `unique_id` in the consensusXML is randomised every run (pre-existing)

Found while verifying Fix 6, NOT caused by it. `generate_unique_id()` in
`clusters_json_to_consensusxml_clusterless.py` uses `hash(bytes)`, which Python randomises per
process — three runs gave three different sets of IDs; stable only under `PYTHONHASHSEED=0`.
Same class of defect as Fix 1, but in the FINAL artifact consumed by OpenMS/ProteoBench.
Blast radius looks contained: called once per map, and `unique_id` appears only in the 25
`<map ...>` lines with no internal cross-references — self-consistent within a run. Risk is for
anything correlating IDs across files/runs. One-line fix (stable digest, e.g.
`int.from_bytes(hashlib.md5(combined).digest()[:8], 'big')`), but it changes emitted IDs → needs
a decision.

### Assessment findings — status

**CLOSED by the fixes above:** edge filtering no-op → Fix 5 (removed). `--filter-method`
unreachable → Fix 6 (removed). D/F diagnostic outputs → Fix 8 (kept, marked). Worker-count
tangle (E) → Fix 10 (`--num_workers`). Impossible log lines (G) → Fix 9 (removed).

**STILL OPEN:**
- `--unpaired-json` is parsed and passed by the orchestrator but never read (the orchestrator
  collects unpaired vertices itself). Removing it means touching the orchestrator + `.nf`.
- `graph_analysis.json` / `graph_composition.json` are written every cycle, read by nobody in
  this chain. Decision taken in Fix 8: KEEP, marked diagnostic-only — useful for manual analysis.
- `unique_id` randomisation (section above) — needs a decision.
- The asymmetric edge-parsing passes in `build_network_graph` (flagged in Fix 13). The
  `cluster_id` layer, the other Fix 13 item, was removed in Fix 14.
- `filter_one_component` still logs `graph is None! Heuristic and induced_subgraph will fall back`,
  but nothing in the filtering path builds an induced subgraph — the only `induced_subgraph` call
  left is in `analyze_graph`, which never runs during filtering. Half that message is stale.
- The stdout banner still reads `Component Building (Clusterless)`; the code finds components
  rather than building them.
- Scoring uses `graph.degree(v)` — degree in the FULL graph — so a candidate is credited for edges
  to vertices that same candidate is deleting. With more than one duplicate set the "joint optimum"
  therefore optimises a proxy, not the connectivity of the surviving subgraph. Pre-existing;
  changing it would alter which vertices survive.

**Historical detail on the two closed findings**, kept because it is the evidence:

- **Edge filtering here was a no-op under EVERY config** (verified). ppm mode: `mz_cutoff` is None
  so the coordinate branch is skipped, and the load-cutoffs-from-file fallback requires
  `rt_cutoff is None` too (orchestrator always passes `--rt_cutoff`) → `edge_cutoff` stays `inf`
  → log shows `49930 total edges, 49930 pass distance cutoff inf`. `mz_cutoff` mode: filters on
  `mz_distance`/`rt_distance`, which pairing edges never carry → `.get(...,0)` → all pass.
  **Harmless** (pairing already enforces cutoffs before emitting an edge), but 3 CLI args and the
  `cutoff_params` block are decorative. This CLOSES the open question from the pairing review.
- **`--filter-method` is unreachable in this chain.** main builds exactly one cluster per
  component, so `clusters_with_duplicates` always has len 1 → scoring always takes the
  `total_degree` branch. Verified: all `scoring_path_counts` are 0 in a real run. ~100 lines
  of modularity scoring dead, plus the arg the orchestrator threads through.
- **E:** the old `num_workers = min(max(4,4), cores_to_use)` — i.e. `min(4, …)`, the comment
  claimed 8, and `num_cores` was separately capped at 10. Safe to tune because Fix 1 made output
  invariant to worker count (proven — parallel and `--disable-multiprocessing` produce identical
  bytes). Replaced by `resolve_worker_count()` in Fix 10.
- **G:** `mega_components` / `large_components` compared `combo_estimate` against
  `COMBINATIONS_BATCH_LIMIT * 10` and `* 1`, but a component's weight is capped by
  `mega_complexity_threshold`. With the constants in use (limit 500,000, threshold 5,000) both
  always printed 0 — confirmed across all 11 filtering runs. Their labels said ">10M"/">1M",
  matching neither. Removed in Fix 9.

### `COMBINATIONS_BATCH_LIMIT` — what it actually is (comment was a fossil)

The old comment said "ORCHESTRATOR-LEVEL", "Increased to 5 million", and "reduce IPC overhead".
All three were wrong: the orchestrator-level batching functions were dead code removed in Fix 2;
the value was 50,000 not 5 million (both added in the SAME commit, so it never matched); and it
cannot affect IPC because the metadata is already inside the worker by then. It is purely a
memory-pacing knob — how often a worker flushes a batch and calls `gc.collect()`. Rewritten to
say so (user then set it to 500,000).
Measured on cycle 3, varying ONLY this constant: components-per-worker identical (461/462/462/461)
at every value — distribution across workers is separate greedy balancing, independent of this.
Batches per worker went 25 → 3 → 1 → 1 for limits 5k → 50k → 500k → 100M, output md5 identical at
all four. Each worker's chunk totals ~110,000, so **at 500,000 batching is effectively OFF (one
batch per worker) and raising it further does nothing**. Peak RSS barely moves (631.9 → 632.3 MB)
because work_items hold REFERENCES into vertex_attrs, not copies.

### Verification setup for this script — USE THE HARNESS

`harness.sh` in the session scratchpad runs a 13-run matrix over REAL data and md5s every output
(36 files, ~30s). Re-create it if the scratchpad is gone; it is the reason Fixes 1 and 2 could be
proven byte-identical.
- `SCRIPT=<path> ./harness.sh <label>` — `SCRIPT` defaults to the live script, so an untouched
  copy (e.g. `pre_phase1.py`) can be baselined against the same matrix. Compare `runs/*/HASHES.txt`.
- Matrix: cycles **20, 30, 45, 5, 3** production-like (analysis ENABLED — the orchestrator never
  passes `--skip-analysis`), plus serial (`--disable-multiprocessing`), forced-heuristic
  (`--mega_complexity_threshold 1`), the `modularity` method, filtering disabled, and a repeat run
  to prove within-build determinism.
- Scale covered: 13,735 → 694,977 edges; 3 → 3,499 components. Cycle 20 exercises BOTH filter
  paths (7 heuristic + 5 exhaustive); cycle 45 is heuristic-only.
- **Never force exhaustive globally** (`--mega_complexity_threshold 1e18`): some components have
  ~1e31 combinations and it hangs forever. To widen exhaustive coverage use a bounded value (50000).
- Cycle 1 is a 1GB pkl — too heavy. Use cycle 3 as the "large" case.

## Next step

`build_network_graph_clusterless.py` is **DONE** through Fix 13 (3658 → 2053 lines, pyflakes clean,
sectioned and commented). What remains for it are the decisions listed under "STILL OPEN" above.

**Scripts not yet cleaned, in the order they run:**
1. `bin/orchestrate_clusterless_cycles.py` — touched repeatedly for parameter threading (Fixes 5,
   6, 10) but never reviewed end to end. It is the loop driver, so it is the right next target.
2. `bin/generate_pairing_report_v2_clusterless.py` (2143 lines) — see
   `MEMORY_OPTIMIZATION_REPORT_generate_pairing_report_v2.md`.
3. `bin/clusters_json_to_consensusxml_clusterless.py` (335 lines) — small, and carries the OPEN
   `unique_id` bug. Its output is the FINAL artifact (OpenMS/ProteoBench), so verify with a
   consensusXML diff, not just JSON hashes.

**Approach that has worked and should be repeated:** read fully → map consumed-vs-dead by tracing
actual consumers → build an equivalence harness over REAL data before touching anything → make
byte-identical changes first, flag output-affecting ones separately for a decision.

**Operational reminders before the next production run:**
- Re-extract the feature pkls — `extract_feature_data.py` fixed a real intensity bug and dropped a
  field, so old pkls are wrong.
- `luxii/unbequant:latest` is stale (missing py deps, old `bin/`); patch via `Dockerfile.unbequant_pydeps`.
- The next run will NOT be byte-comparable to old runs (Fix 1 and Fix 12 both changed ordering).
