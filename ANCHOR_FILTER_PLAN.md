# Anchor-based component trimming — plan and implementation record

**Status: IMPLEMENTED on a separate file line. `bin/` and the original workflow are untouched.**
Companion doc: `FILTER_SCORING_HANDOFF.md` (the duplicate-filter scoring decision this builds on).

## The new file line

| new | copied from | changed |
|---|---|---|
| `bin_anchor/pair_features_1_clusterless.py` | `bin/` | persists `scaling_constants` |
| `bin_anchor/build_network_graph_clusterless.py` | `bin/` | scoring fix + anchor filter |
| `bin_anchor/orchestrate_clusterless_cycles.py` | `bin/` | forwards the new args |
| `bin_anchor/{clusters_json_to_consensusxml,extract_feature_data,compare_retention_time_alignment,generate_pairing_report_v2}_clusterless.py` | `bin/` | unchanged copies |
| `map_mzml_features_v3_anchor.nf` | `map_mzml_features_v2_clusterless.nf` | new params, paths rewired to `bin_anchor/`, outdir `feature_analysis_clusterless_anchor` |

**Baseline parity was verified before any change**: the untouched copies reproduce
`accepted_components.json`, `deleted_vertices.json`, `graph_analysis.json` and
`graph_composition.json` byte-for-byte on cycle 3.

## Measured results (cycle 3, ppm=5 / rt=150, 2769 components)

| run | avg linkage | connections | orphans | ppm conform | pep conflicts | minority |
|---|---|---|---|---|---|---|
| original `bin/` | 6.932 | 72589 | 141 | 94.2% | 18.1% | 3.68% |
| scoring fix only | 6.932 | 73492 (+903) | 100 | 94.3% | 18.0% | 3.65% |
| **anchor_connectivity** | 6.787 | **82179 (+13%)** | **37** | **94.8%** | **16.7%** | **3.39%** |
| anchor_nearest | 6.787 | 79410 | 58 | 94.6% | 16.8% | 3.65% |

The scoring fix delivers exactly the +903 connections predicted in `FILTER_SCORING_HANDOFF.md`,
which is the cross-check that it is doing what was designed.

### Radius modes, measured (prototype, all components)

| mode | avg linkage | file coverage | evicted non-dups | ppm conform | pep conflicts |
|---|---|---|---|---|---|
| `full` (chosen) | 6.787 | 97.9% | 232 | 94.8% | 16.9% |
| `prune` | 6.224 | 89.8% | 1107 | **100%** | 14.0% |
| `half` | 5.889 | 84.9% | 1842 | 95.6% | 12.6% |

`full` costs only 2.1% of linkage; `prune` buys total pairwise conformance for 10%; `half` is
dominated - it loses more linkage than `prune` AND conforms less, because halving the radius
splits pairs that `prune` would have kept.

## FULL-RUN RESULT (25 files, ppm=5 / rt=100, anchor rt=150, radius `full`)

Two complete pipeline runs, identical inputs and parameters, differing only in `--filter_method`.
The baseline reused the anchor run's cached extraction and RT alignment via `-resume`, so the
comparison isolates the filter exactly.

| | `duplicate_only` | `anchor_connectivity` |
|---|---|---|
| cycles to converge | 75 | 82 |
| wall time | 1h 46m | 1h 47m |
| components | 161287 | 161062 |
| vertices | 1372578 | 1372578 |
| **avg linkage size** | **8.510** | **8.522** |
| avg linkage (>=2) | 11.720 | 11.719 |
| pep_ident conflicts | 13.4% | **13.1%** |
| minority vertices | 1.47% | **1.34%** |
| PAIRWISE conforming | 94.33% | **95.54%** |
| STAR conforming | 97.80% | **99.67%** |
| worst m/z span | 18.54 ppm (**3.71x**) | 9.97 ppm (**1.99x**) |
| worst RT span | 7181 s (**47.9x**) | 299 s (**1.99x**) |
| components beyond 2x cutoff | 1612 (1.43%) | **0** |

**Linkage did not drop.** 8.510 -> 8.522, i.e. unchanged within noise, despite the anchor filter
evicting 145309 vertices in cycle 1 alone. The cycle loop recycles them and they find groups they
fit; the cost is 7 extra cycles, not lost linkage. The concern that linkage might fall below ~7 did
not materialise.

**The unambiguous win is the bound on transitive drift.** Connected components let membership chain:
A-B and B-C inside the cutoff does not keep A-C inside it, and nothing downstream re-checked. The
baseline shows how far that goes - one component spanning 7181 s of retention time, 47.9x the
cutoff, emitted as a single consensus feature whose centroid averages over two hours. Under the
anchor filter no component exceeds 2x on either axis, and 2x is exactly the theoretical maximum for
`full` radius (two members on opposite sides of the anchor). Confirmed empirically: worst observed
1.99x, zero components beyond 2x.

**The identification-quality gains are real but small**: conflicts 13.4% -> 13.1%, minority vertices
1.47% -> 1.34% (-8.8%). Consistent with everything else measured in this project - geometry-only
scoring is only weakly coupled to identification evidence.

**Known imperfection**: 0.33% of components are not STAR-conforming, i.e. have no member within
tolerance of all others. Cause: the duplicate-resolution step picks one vertex per file from the
pool, and when the anchor's own file has several pool members it can pick a different one, evicting
the anchor. The 2x bound survives regardless (every member was within 1x of the anchor), but the
star property does not. Fix if wanted: pin the anchor in `resolve_pool_*`.

## Two implementation findings worth keeping

1. **The greedy fallback had to be replaced.** When a pool has too many combinations to enumerate,
   the old path fell back to `generate_filter_candidates`' heuristic, which ranks by FULL-GRAPH
   degree - the separable measure this whole change exists to remove. It made `anchor_connectivity`
   score *worse* than `anchor_nearest` on 152/2769 components (150 of them on that path).
   Replaced with hill-climbing on the real score, seeded from nearest-to-anchor.
2. **Incremental scoring is not optional.** Rescanning the component's edge table per swap trial
   took 301 s for cycle 3. Scoring only the edges incident to the swapped pair brought it to 17 s.

## Why

Today the filter's only lever is choosing among duplicates — it keeps at most one vertex per file
and cannot remove anything else. Measured consequences (cycle 3, ppm=5/rt=100 run,
`work/c1/ed55a3f9ae28b013320e845be68e61`):

- **9.59%** of finished components violate the pairwise cutoff (`final_components_all_cycles.json`,
  ppm=10/rt=100 run): RT span > cutoff in 8.46%, m/z span > cutoff in 2.04%. Worst RT span
  **7174 s (71.7x)**, worst m/z span **37.88 ppm (3.8x)**. Components are connected components, so
  membership is transitive while the cutoff is pairwise — A–B and B–C inside the cutoff does not
  keep A–C inside it.
- A pure veto layer in front of the filter **cannot fix this**: of the 88 exhaustive components with
  no conforming candidate, **81 are infeasible because of forced vertices** — files that
  contributed exactly one vertex and so appear in every candidate. Adding the constraint changed
  **3 of 1267** components (ppm) and **0** (ppm+RT).
- Allowing eviction of non-duplicates does fix it: **141 non-duplicate vertices** evicted across
  84 components (6.6%) buys 100% pairwise ppm conformance, at **-1.84% file slots** and
  **-1.9% connections**. Zero components had no solution.

Evicted vertices flow into `deleted_vertices.json` and re-pair next cycle, so the cost is a
deferral rather than a loss — the steady-state effect is unmeasured.

## What it does

1. **Anchor** — for every vertex in the component, count how many *distinct files* have at least one
   vertex within cutoff range of it. Highest count wins.
2. **Pool** — every vertex within cutoff range of the anchor. Everything else is evicted to
   `deleted_vertices.json`.
3. **Resolve duplicates** inside the pool, by a selectable method:
   - `anchor_connectivity` — the scoring from `FILTER_SCORING_HANDOFF.md`: maximise distinct
     retained connections, then minimise total retained distance.
   - `anchor_nearest` — per file, keep the vertex closest to the anchor.

This replaces the one-per-file filter rather than layering on it: the output still satisfies
one-per-file, and additionally trims outliers.

## Decisions (settled 2026-07-28)

### D1 — anchor radius semantics: **`full`** (the cutoff), REVISABLE

"Within cutoff of the anchor" is weaker than "pairwise within cutoff": two vertices each 10 ppm from
the anchor can be 20 ppm apart. Modes considered:

| mode | radius | pairwise guaranteed | note |
|---|---|---|---|
| **`full`** | **the cutoff** | **no** | **CHOSEN for now** |
| `half` | cutoff / 2 | yes (triangle inequality) | strictest; evicts most |
| `prune` | cutoff, then drop violators | yes | middle |

**Decision: `full`.** Keep `--anchor_radius_mode` as a parameter with the other two implemented, so
this can be revisited without a code change — it is explicitly marked as changeable later. Phase 1
should still report all three so the trade is on record. For reference, the strict pairwise variant
costs 141 evicted vertices / 1.9% of connections.

### D2 — distance metric: **the pairing weight from the Nextflow workflow**

Use the same weighted space the pairing step uses, driven by `params.mmf_mz_rt_weight_ratio`
(currently `0.0001` in `map_mzml_features_v2_clusterless.nf:46`), rather than inventing a
cutoff-normalised metric:

```
x_scaled = (mz - min_mz) / range_mz
y_scaled = (rt - min_rt) / range_rt * mz_rt_weight_ratio
ranking  = euclidean distance in that space
in range = box test on the anchor tolerances (see D3), not on the scaled distance
```

**Neither the weight ratio nor the normalisation constants exist in the graph script today.**
Decision: **pairing writes them into `paired_features_edges.pkl`** — they are pairing-derived facts,
so the pkl is where they belong, and this avoids threading anything through
nextflow -> orchestrator. Add to `pair_features_1_clusterless.py`'s `pairing_parameters` block:

```
mz_rt_weight_ratio, min_x, range_x, min_y, range_y
```

(the four constants are already computed at `pair_features_1_clusterless.py:462-473`; they are
simply not persisted).

This makes Phase 2 touch `pair_features_1_clusterless.py` as well — additive only, no behaviour
change. **Fallback required:** pkls written before this change lack the fields, so the graph script
must recompute from its own vertices and warn. Every scratchpad measurement to date uses such pkls,
so Phase 1 runs on the fallback path by definition.

### D3 — RT coordinates: keep as-is, widen the RT tolerance to **150**

`ident_lookup` stores **pre-alignment** RT, so any coordinate-based RT test is wrong by up to ~30 s
(measured: 1.031% of 589962 edges exceed the RT cutoff in exported coordinates, max 129.5 s; m/z is
exact — 0 edges violate, max exactly at the cutoff).

**Decision: do not fix the coordinates now. Widen the anchor's RT tolerance from 100 to 150** to
absorb the distortion.

**The anchor gets its OWN tolerances, independent of the pairing cutoffs** — `--anchor_ppm` and
`--anchor_rt_cutoff` (default 150). Pairing behaviour is unchanged; `params.mmf_rt_cutoff` stays at
100 and continues to decide which edges exist.

**The m/z tolerance is computed in ppm exactly as pairing does it**
(`pair_features_1_clusterless.py:571`) — from the reference feature's own m/z, which here is the
anchor:

```
in range(v, anchor) :=
    |mz_v - mz_anchor| <= mz_anchor * (anchor_ppm / 1e6)
    AND |rt_v - rt_anchor| <= anchor_rt_cutoff
```

Asymmetric by construction, and correct for an anchor method — the anchor *is* the reference. Note
the scratchpad measurements so far used a symmetric ppm form (`2·Δmz/(mz_a+mz_b)`); the difference
is negligible at these tolerances but Phase 1 should switch to the pairing form for exactness.

**Later work (not now):** load the loess corrections for each file *combination* and apply them at
every distance calculation, rather than relying on one global corrected coordinate.

Recorded for whoever picks that up — the concern that per-pair alignment makes a single corrected RT
impossible does **not** apply to the current implementation.
`_apply_json_rt_corrections` loops per *file*, not per pair (`pair_features_1_clusterless.py:354`):
file 0 is the reference, each file `i>=1` gets the single `(file_i, file_0)` model, and every
feature receives one corrected `y_center` on file 0's scale (`:379`). It must collapse this way —
the KD-tree is built once per file, not once per compared pair. So a per-combination scheme would be
a genuine *improvement* in fidelity, not a repair of something broken.

Per RT method today:
- `aligntree` — already correct; `apply_aligntree_rt_transformations` mutates the feature lists
  directly, so `ident_lookup` holds corrected RT.
- `loess` and the other fitted methods — corrections applied to `f.copy()` objects, never reach
  `ident_lookup`. A `y_center_corrected` field would be additive and small.
- `noalign` — nothing to correct.

## Anchor selection — full algorithm

Ties are resolved by *evaluating* the tied anchors, not by an arbitrary key. Specified 2026-07-28.

```
1. For every vertex v in the component:
       files(v) = number of DISTINCT files having >=1 vertex within range of v
                  (box test on the raw cutoffs: ppm, and anchor_rt_cutoff)

2. A = { v : files(v) is maximal }          # all tied anchors are candidates

3. If |A| == 1: anchor = the single member.
   Otherwise, score every a in A by the method in use and take the best:

   anchor_nearest      : pool(a) = one vertex per file, each the nearest to a.
                         score = -SUM(ranking distance from a to each kept vertex)
                         i.e. prefer the anchor whose in-range unique vertices are closest.

   anchor_connectivity : pool(a) = all vertices within range of a.
                         Evaluate EVERY one-per-file combination inside pool(a).
                         Best combination = ( retained connections , -total distance )
                         where a connection counts ONLY if both endpoints are inside
                         pool(a) and kept.
                         score = that best combination's key.
                         i.e. prefer the anchor whose best achievable candidate is best.

4. Final tie-break, if scores are still equal: lowest vertex id.
   REQUIRED - the module's determinism contract (docstring :24-30) admits no
   dependence on set or dict iteration order.
```

**Cost note.** Step 3 runs the full combination enumeration once per tied anchor, so worst-case work
is `|A| x` the current per-component cost. `mega_complexity_threshold` bounds one enumeration, not
the product. Phase 1 must measure the distribution of `|A|` and, if it has a long tail, add either a
cap on tied anchors evaluated or a threshold divided by `|A|`.

## Phases

### Phase 1 — prototype in scratchpad, no `bin/` changes

Extend the existing harness to score side by side: today's filter, the `FILTER_SCORING_HANDOFF`
recommendation, `anchor_connectivity` and `anchor_nearest`, each across the three radius modes.
Metrics already implemented in scratchpad scripts: retained connections, file coverage, evicted
non-duplicates, orphan count and connectivity rate, pairwise conformance, pooled pep_ident accuracy.

Deliverable: a comparison table. **No pipeline changes until it is reviewed.**

### Phase 2 — `build_network_graph_clusterless.py`

- New args `--filter_method {duplicate_only,anchor_connectivity,anchor_nearest}` (default
  `duplicate_only`, bit-identical to today unless opted in) and `--anchor_radius_mode
  {full,half,prune}`.
- Read `cutoff_params` from the input pkl — already present and deliberately unread (`:132`), so no
  new plumbing is needed for the cutoffs themselves.
- New functions alongside section 6: `select_anchor()`, `anchor_pool()`, `resolve_nearest()`.
  `filter_one_component()` branches on the method. `generate_filter_candidates()` and the scoring
  become the duplicate-resolution step of `anchor_connectivity` — reused, not replaced.
- **Determinism is a hard contract** (module docstring `:24-30`). Implement the anchor tie-break
  exactly as specified in "Anchor selection", including the final vertex-id fallback.
- New args `--anchor_ppm` and `--anchor_rt_cutoff` (default 150). The weight ratio and
  normalisation constants come from the pkl (D2), not the CLI.
- `pair_features_1_clusterless.py`: persist `mz_rt_weight_ratio` + the four normalisation constants
  into `pairing_parameters` (additive; no behaviour change). Graph script falls back to recomputing
  them, with a warning, for pkls written before this.
- Evicted non-duplicates need a distinct `reason`. Note the export drops diagnostic fields today
  (`:1030-1034`), so surfacing eviction reasons on disk is a small extra change if wanted.

### Phase 3 — threading

- `orchestrate_clusterless_cycles.py`: two `add_argument`s -> `config` dict -> forward to the graph
  script at `:224`, exactly the shape `mega_complexity_threshold` already uses.
- `map_mzml_features_v2_clusterless.nf`: `params.mmf_filter_method` and
  `params.mmf_anchor_radius_mode`, appended to `filter_args` (~`:239`). Defaults preserve current
  behaviour.

### Phase 4 — verification

Recreate `harness.sh`, snapshot the current script, A/B against it. **The harness will fail on md5 —
that is expected; compare structurally.** Then full-pipeline runs across the method matrix,
comparing final components on: count, size distribution, pairwise conformance, peptide-conflict
rate, and cycles to convergence.

## Risks

1. **Convergence.** Eviction recycles more vertices per cycle, so the loop may need more cycles or
   settle on a different fixed point. The reference run used 74 cycles. Must be measured.
2. **Cascade.** Output changes compound across cycles; per-cycle deltas understate the end state.
3. **Anchor instability.** A different anchor in an early cycle reshapes everything downstream.
   Addressed by the evaluate-the-ties algorithm above plus the vertex-id fallback; the residual
   risk is cost, not correctness (see the cost note there).
4. **Not measurable by pep_ident.** Every scoring variant tested so far was statistically
   indistinguishable on identification ground truth (all p > 0.2, n=95 strong cases). Expect the
   same here; judge on conformance and coverage, not accuracy.
5. **Never edit `bin/*.py` while a pipeline run is in progress** — Nextflow invokes scripts by
   absolute path into the mounted workspace, so a running job picks up edits mid-run.

## Still to decide before Phase 1

All design questions are settled (D1, D2, D3 above). Remaining: **go-ahead for Phase 1** — scratchpad
prototype, nothing in `bin/`.

## Relationship to the scoring decision

Independent and composable. The scoring fix repairs machinery that currently does nothing
(exhaustive returns the greedy answer on 1267/1267 components); this adds a new capability.
`anchor_connectivity` *is* that scoring applied to a trimmed pool, so landing the scoring fix first
means it is exercised on its own before this is layered on top.
