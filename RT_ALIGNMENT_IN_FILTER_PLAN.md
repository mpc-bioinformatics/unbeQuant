# Plan (ON HOLD) — RT alignment inside the anchor/BFS filter

**Status: ON HOLD. Nothing implemented, nothing in `bin_anchor/` touched by this plan.**
Companion docs: `ANCHOR_FILTER_PLAN.md` (the filter itself), `FILTER_SCORING_HANDOFF.md`
(the duplicate-filter scoring decision underneath both).

## Why this is worth doing

`build_network_graph_clusterless.py` reads `y_center` out of `ident_lookup`, and that value is
**pre-alignment** — the loess corrections are applied to `f.copy()` objects inside pairing and never
reach the lookup. Measured on cycle 3: **1.031% of 589962 edges** exceed the RT cutoff in exported
coordinates, worst overshoot **29.5 s**. m/z is exact (0 edges violate, max lands exactly on the
cutoff).

Against the 150 s anchor tolerance that is noise. Against the **50 s BFS budget** now configured in
`map_mzml_features_v3_anchor.nf`, up to 30 s of a 50 s budget can be pure alignment artifact — 60%
of it. For the current configuration this is a correctness fix, not a refinement.

## The model to mirror

`pair_features_1_clusterless.py` already does this properly. `_resolve_rt_correction(query_file,
ref_file, corrections)` returns `(model, sign, status)`; the shift is `sign * model(rt)`:

| key found | mapping | status |
|---|---|---|
| reverse `(query, ref)` | `rt_on_ref = rt_query + model(rt_query)` | **exact** — evaluated on its own fitted domain |
| forward `(ref, query)` | implicit; one fixed-point iteration | **approximate**, ~0.1–4.4 s error |
| neither | no shift | **missing** — compared unaligned |

Pairing then builds a *reference frame* per reference file (`_build_frame`) with every other file
aligned onto it, and matches inside that frame. The filter should mirror exactly this.

Every RT comparison in the filter has a natural reference: the **anchor**, or the **BFS start**. So
"align the other vertex onto the reference's file scale" is well defined at every call site.

## Design

### One accessor, four call sites

```python
aligned_rt(v, ref_file) -> float      # v's RT expressed on ref_file's scale
```

replaces the raw `y_center` at:

| site | today | becomes |
|---|---|---|
| `_in_anchor_range(anchor, v)` | `abs(y[a] - y[v])` | `abs(y[a] - aligned_rt(v, file[a]))` |
| `_pair_conforms(u, v)` (prune mode) | `abs(y[u] - y[v])` | `abs(y[u] - aligned_rt(v, file[u]))` |
| `_anchor_scaled(v)` (ranking distance) | `y[v]` | `aligned_rt(v, ref_file)` |
| `grow_bfs` budget accumulation | `abs(y[v] - y[u])` | see "distance from the origin" below |

`_anchor_scaled` is the only signature change — it currently takes a vertex alone and must also take
the reference.

### ALREADY HANDLED — do not double-solve

The tree/candidate **scores** sum the stored edge `distance` values, and pairing computed those
*inside its reference frame*, so they are already alignment-corrected. Neither `score_kept_set` nor
the BFS total-distance ranking needs to change. **Only the coordinate-based admission tests are
wrong today.**

Verified in `pair_features_1_clusterless.py`: `_build_frame` shifts `rt` by `sign * model(rt)`
(`:362-363`) and derives `scaled` from the shifted array (`:365`); the stored `distance` is computed
from `query_group['scaled'][idx]` against `ref_point`. The corrected RTs deliberately live in the
frame and are never written back onto the feature dicts (`:348-350`) — which is exactly why
`ident_lookup` ends up holding raw values while edge distances do not.

### Scale check — the weight ratio does NOT suppress RT in the cutoff window

Worth recording because it is counter-intuitive and was initially got wrong. With
`mz_rt_weight_ratio = 0.001`, `range_x = 1049` Da and `range_y = 7196` s, contributions to the
SCALED distance are:

| quantity | scaled contribution |
|---|---|
| m/z, 5 ppm at m/z 1000 | 4.77e-06 |
| RT, 100 s | **1.39e-05** (2.9x larger) |
| RT, a 30 s misalignment | 4.17e-06 |

The 0.001 weight squashes RT's *full range*, but 5 ppm is only 4.8e-6 of the m/z range whereas 100 s
is 1.4% of the RT range. Net: inside the cutoff window the metric is RT-dominated, not m/z-dominated.

Two consequences: the alignment correction already in `distance` is doing real work rather than being
washed out by the weighting; and an uncorrected `y_center` perturbs the DOMINANT axis of the geometry
the filter reasons about, which makes this fix more valuable than the RT-budget arithmetic alone
suggests.

### Performance

Naive per-pair evaluation is O(n^2) model calls per component, ~59M in cycle 1. Too slow.

**Precompute one global table, vectorized per file pair.** `LOESSInterpolator.__call__` already
accepts arrays:

```
for each (query_file, ref_file):                  # 25 x 25 = 625 calls for a WHOLE cycle
    model, sign, status = _resolve_rt_correction(query_file, ref_file, corrections)
    aligned[ref_file][vertices of query_file] = rt + sign * model(rt)      # one array call
```

- **625 vectorized evaluations per cycle**, against 59M scalar ones
- memory `V x F` float32: 1.34M x 25 x 4 B = **~134 MB** in cycle 1
- built BEFORE the worker Pool forks, so workers inherit it copy-on-write exactly like `graph` and
  `vertex_attrs` (SECTION 2) — no per-worker copy

Fallback if that memory is unwelcome: compute lazily per component, ~19M scalar calls, roughly
60–95 s per cycle. The table is preferred.

### aligntree and loess are NOT symmetric

- **`aligntree`** — one global transformed RT per feature. Overwrite `y_center` at load, exactly as
  pairing does via `apply_aligntree_rt_transformations`. No table, no frames, trivial.
- **`loess` / `polynomial` / `spline` / `rbf` / `piecewise`** — pair-specific. No single "corrected
  RT" exists; the table is required.
- **`noalign`** — nothing to do.

## Decisions

### D1 — export: **INTERNAL ONLY** (settled)

The filter uses aligned RT internally; `final_components_all_cycles.json` keeps the raw `y_center`,
so the consensusXML centroid and everything downstream are untouched.

**Consequence that must not be forgotten:** every conformance number measured so far
(`summarize_filter_metrics.py`, `star_test.py`) reads the exported raw coordinates and would keep
reporting *unaligned* drift, understating the improvement. Those measurements have to be re-pointed
at aligned space as part of this work, or the results will look like a regression when they are not.

### D2 — how "distance from the origin" is measured: **OPEN**

Today the BFS budget is a **path sum**, not a direct offset:

```python
nrt = path_rt[u] + abs(rtof[v] - rtof[u])      # u = path parent, v = new vertex
```

Two designs follow, and they differ in what they guarantee:

**Option A — per-hop pairwise.** Each hop `u -> v` uses the directly-fitted `(file(v), file(u))`
model. Most faithful per individual comparison. But the accumulated sum then mixes frames — hop 1 in
file A's frame, hop 2 in file B's — and loess models are fitted per pair independently, so
`align(w, A) != align(align(w, B), A)`. The triangle-inequality guarantee that produced **zero
violating pairs** assumes a single metric, so under Option A the bound degrades from guaranteed to
approximate, by roughly the forward-key error (~0.1–4.4 s) per hop, growing with depth.

**Option B — direct offset in the start's frame (RECOMMENDED).**

```
admit v  iff  |rt_S - aligned_rt(v, file(S))| <= rt_budget
```

Still **fully pairwise**: `aligned_rt(v, file(S))` uses the directly-fitted `(file(v), file(S))`
model. No composition, no chaining. Every vertex is aligned to the start by its own pair-specific
curve. What is given up is path-length semantics; what is gained is an exact triangle bound, because
everything lives in one frame.

Since the budget's entire job is delivering that guarantee, trading it away to gain hop-by-hop
faithfulness is the wrong direction.

Consequence to accept under B: the budget becomes "radius from the start" rather than "path length".
At depth 1 the two are identical. At depth >= 2 B is more permissive — a vertex two hops out but
geometrically close to the start now qualifies. Depth moved average linkage by 0.002 in testing, so
this is expected to be negligible, but it is unmeasured.

**Still to decide:** go straight to B, or measure A vs B divergence in Phase 1 first.

## Threading

`--rt_correction_json`, `--rt_alignment_method`, `--trafoxmls_dir` -> the graph script, forwarded by
`orchestrate_clusterless_cycles.py` from config it already holds. **Nextflow needs no new params** —
`params.mmf_rt_alignment_method` exists and the correction JSON is already an orchestrator input.

Reuse rather than duplicate:

```python
from pair_features_1_clusterless import (load_rt_corrections_json,
                                         _resolve_rt_correction,
                                         build_rt_correction_model)
```

Same directory, and that module is `__main__`-guarded, so importing it is side-effect free.

## Phases

1. **Scratchpad prototype.** Build the table; re-measure cycle 3 for `bfs 2.5:50` and
   `anchor half`, aligned vs unaligned; if D2 is still open, measure A vs B divergence.
   Deliverable: how much the budgets actually move. **Nothing in `bin_anchor/` until reviewed.**
2. **Implement** in `build_network_graph_clusterless.py`, gated on `--rt_alignment_method` with
   default `noalign` — today's behaviour, so it is opt-in and gives a free regression test.
3. **Thread** through the orchestrator.
4. **Verify**: `duplicate_only` byte-identical; determinism across worker counts; `noalign`
   identical to current output; re-point the metrics scripts; then a full run.

## Risks

1. **Memory** — ~134 MB pre-fork in cycle 1. Small against the 40 G allocation, but new.
2. **Cycle 2+** — vertices arrive from `deleted_vertices.json`. Filenames persist so the models
   still resolve, but this needs an explicit check that the file-name keys still match.
3. **Approximate direction** — pairs served only by a forward key carry 0.1–4.4 s error. Log the
   exact/approximate/missing counts per cycle, as pairing already does.
4. **Output changes** — components will differ, and it cascades across cycles. Per-cycle deltas
   understate the end state.
5. **Do not edit `bin_anchor/*.py` while a run is in progress** — Nextflow invokes the scripts by
   absolute path, so an edit lands between cycles and different cycles execute different code.
