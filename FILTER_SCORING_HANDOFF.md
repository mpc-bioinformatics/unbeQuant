# Handoff — exhaustive filter scoring in `build_network_graph_clusterless.py`

Open design decision. **Nothing has been changed yet** — the file is still byte-identical to the
verified state. This doc is the context needed to resume cold.

## The problem in one paragraph

The duplicate-vertex filter keeps at most one feature per input file per component. Two paths pick
the survivors: a greedy HEURISTIC (large components) and an EXHAUSTIVE search (small ones). Both
currently score with `sum(graph.degree(v) for v in kept)` — the degree in the FULL graph, which
counts edges to vertices that same candidate is deleting. That is not "connections retained", and
because it is a sum of per-vertex constants it is **separable**, so maximising it is identical to
choosing each duplicate set independently. Measured: the exhaustive path returns the same answer
as the greedy heuristic on **1267/1267** components. The enumeration buys nothing.

## Intended semantics (confirmed with the user)

- **HEURISTIC path — CORRECT AS IS. Do not touch.** "Keep the vertex with the most connections",
  decided per duplicate set. `max(vertex_ids, key=degree)` is exactly that.
- **EXHAUSTIVE path — WRONG.** Intent is: for each combination, measure how many connections
  *remain*, i.e. edges with BOTH endpoints in the kept set. That is non-separable, which is what
  makes enumerating combinations meaningful in the first place.

Note: every edge lies inside exactly one component and components are filtered independently, so
counting surviving edges within the component == that component's contribution to the whole graph.

### Why today's score is separable — the exact identity

For a kept set K inside a component (every neighbour of a component member is also a member):

```
degsum(K) = 2 * |edges with both endpoints in K| + |edges from K to deleted vertices|
```

each edge contributing once per endpoint in K. The intended score is the first term. The second
rewards keeping a vertex *because* it had many connections to features this same candidate deletes,
and it is also what makes the total decompose per-vertex — i.e. the sole cause of the 1267/1267
greedy agreement. Removing it is the whole fix.

### The joint choice does pay off under the corrected score

Scoring on the candidate's induced subgraph (`joint_vs_greedy.py`, cycle 3):

| | |
|---|---|
| joint optimum strictly beats greedy | **314 / 1267 (24.8%)** — today: 0/1267 |
| extra connections retained | **+903** |
| per-component gain | 119×+1, 74×+2, 39×+3, 21×+4, 17×+5, 16×+6, 11×+7, … one at +16 |

Worked example, component 355 (7 files, one duplicate set): greedy keeps `v25261` (degree 9) and
retains 16 connections; the joint optimum keeps `v25262` (degree 8) and retains 17. `v25262` holds
the strongest link to `v13639` — a vertex whose only other tie into the component is a single edge —
plus links to `v43127`/`v34694`/`v994` that `v25261` lacks. Locally worse, jointly better, and
invisible to any per-vertex rule.

## Measurements (cycle 3, threshold 5000, 1267 exhaustive components)

| question | answer |
|---|---|
| exhaustive winner == greedy pick, today | **1267/1267 (100%)** |
| today's pick differs from max-retained-edges | 314 (24.8%) |
| switching to retained-edge count: components changed | 693 (54.7%) [no secondary key — with a `degsum` secondary it is 314; the difference is tie winners flipping to an equally-scoring candidate] |
| vertices swapped kept↔deleted | 1066 |
| edges retained | 30060 → 30963 (**+903, +3.0%**) |
| exhaustive would then beat greedy in | 693 (54.7%) |
| scoring cost | 205 ms → 1283 ms (**6.3×**) |
| candidates tied at best edge count | **880 (69.5%)** — none tied at zero |
| degree-sum secondary key changes the tie winner | 452 (35.7%) |

**Ties are the common case, not an edge case.** A plain edge count would decide 69.5% of components
by "whichever candidate the generator emitted first" (an artifact of filename sort order). A
secondary key is mandatory, not optional.

## Candidate formulas

**Option 1 — lexicographic**

```
score(K) = ( |edges with both endpoints in K| ,  sum(deg_G(v) for v in K) )
```
maximise lexicographically. Primary = intent. Secondary = today's score, consulted only on ties
(69.5% of the time), and it is also the heuristic's criterion, so it is principled rather than
arbitrary.

**Option 2 — distance-weighted**

Normalise distances within the component, weight each retained edge by `1 - d_norm`, sum.

⚠️ **Known flaw, measured:** with min-max normalisation the LONGEST edge in every one of the 1267
components gets weight exactly **0** — retaining the weakest connection scores identically to
severing it. Worse, the zero point is set by whichever edge happens to be longest *in that
component*, so a component whose distances span 1.00e-5..1.01e-5 is stretched to the full 0..1
range exactly like one spanning 1e-9..1e-5. The score cannot distinguish "all edges nearly equally
good" from "one much worse".

Variants that avoid the zero:
- `w = d_min / d` — strongest = 1, weakest = `d_min/d_max` > 0, scale-free, never 0.
- `w = 2 - d_norm` ∈ [1,2] — every edge worth ≥1 so "amount" is preserved, strength adds ≤1 more.
- lexicographic `(count, strength)` — count dominates, strength only breaks ties.

## Variant comparison — RESULTS

Harness: `variant_compare2.py` (session scratchpad; read-only, imports the module and re-runs
candidate enumeration offline on the cycle-3 pkl, threshold 5000, 1267 exhaustive components).
Validated by reproducing two independent prior numbers exactly: today's score agrees with greedy
on **1267/1267**, and today's retained-distinct-edge total is **30060**.

| variant | primary → secondary | raw edges | distinct | strength `Σ dmin/d` | ties left | ≠ greedy |
|---|---|---|---|---|---|---|
| V0 today | degsum | 55598 | 30060 | 4047.6 | 247 | 0 |
| V1 | raw count → degsum | **56842** | 30557 | 4217.1 | 100 | 266 |
| V6 | raw count → `Σ dmin/d` | **56842** | 30513 | **4294.9** | 6 | 348 |
| V1d | distinct → degsum | 55981 | **30963** | 3978.7 | 134 | 314 |
| V5d | distinct → `Σ 1-d_norm` | 55992 | **30963** | 4014.8 | **4** | 418 |
| V6d | distinct → `Σ dmin/d` | 55718 | **30963** | **4113.4** | **4** | 484 |
| V7d | distinct → mutual → `Σ dmin/d` | 56077 | **30963** | 4048.4 | **4** | 409 |

Pairwise disagreement (kept set differs, of 1267): V1↔V6 119, V1d↔V5d 163, V5d↔V6d 191,
V1d↔V7d 142, V1↔V1d 232, V0↔V6d 484. Full matrix in `variant_results2.json`.

### The finding that reframes the choice: "one edge" is ambiguous

Pairing emits **directed** records, and both directions are attempted independently, so a mutual
nearest-neighbour pair lands in the graph **twice** — two parallel edges, always with identical
`distance` (verified: 0 mismatches). Measured on cycle 3: 589962 records over 399281 distinct
vertex pairs; **47.8%** of pairs are mutual. Within the denser exhaustive components the raw:distinct
ratio is ~1.85.

So "count the edges with both endpoints kept" has two readings, and they are **not** the same
objective — maximising one does not maximise the other:

- **raw records** (V1/V6): a mutual match counts 2, a one-way match counts 1. Reaches 56842 raw but
  only 30557 distinct — 406 fewer real connections than the distinct-optimum.
- **distinct connections** (V1d/V5d/V6d/V7d): every connected feature pair counts once. Reaches
  30963 distinct, the ceiling; this is the number Option 1 was written to mean.

Raw counting silently weights mutual agreement double. That is defensible evidence-wise, but it is
an undocumented weighting smuggled into a count. V7d makes it explicit instead — distinct first,
mutual agreement as an ordered secondary — and lands within 4 components of V5d/V6d.

### Recommendation

**V6d: `(distinct retained edges, Σ d_min/d over those edges)`.**

- Primary is the literal intent, and hits the retained-connection ceiling (30963, +3.0% over today).
- Secondary is scale-free, never zero, so it discriminates on all 880 count-ties instead of
  falling back to generator order. It leaves only **4/1267** components genuinely tied — and it
  fixes that at the same cost as the flawed min-max form (V5d, also 4).
- Highest retained strength of any distinct-optimal variant (4113.4).
- Trade: `degsum` drops to 137037 (lowest). That is expected and not a regression — degsum is the
  separable quantity we are deliberately abandoning.

If the mutual-match evidence should count, take **V7d** instead — same primary, same 4 ties, keeps
25114 mutual pairs vs V6d's 24755. Do **not** take V1/V6: they optimise a count that double-counts
by accident.

⚠️ **The pep_ident evaluation below could not confirm any of this.** All variants, including today's
score, are statistically indistinguishable on identification ground truth. The recommendation above
rests on the *argument* — which objective is correct, and which tie-break is principled rather than
sort-order — not on measured accuracy. Read that section before treating V6d as validated.

**NEXT STEP:** user picks V6d / V7d / V1d, then implement (see constraints below). The
implementation must dedupe parallel edges per component before counting — `graph.induced_subgraph`
returns both records.

## pep_ident ground truth — the variants are NOT distinguishable

Peptide identifications are independent of the geometry the score optimises, so they can referee
the choice. Within a component the majority peptide (vote over labelled vertices) is taken as the
true analyte; a duplicate set is then judged in one of two tiers:

- **STRONG** — the set holds a majority-carrier AND a member carrying a *different* peptide.
  Picking the wrong one is a real error.
- **WEAK** — majority-carrier vs an *unlabelled* rival. Only suggestive: the rival may be the same
  peptide, merely unidentified.

Pooled over cycles 2–73 (`ident_paired.py`), exact McNemar against today's score:

| variant | STRONG acc (n=95) | p vs V0 | WEAK acc (n=3599) | p vs V0 |
|---|---|---|---|---|
| V0 today | 62.6% (57/91) | — | 57.9% | — |
| V1 raw→deg | 66.7% | 0.219 | 58.6% | 0.206 |
| V6 raw→ratio | 65.6% | 0.453 | 58.5% | 0.341 |
| V1d dist→deg | 66.7% | 0.508 | 57.5% | 0.568 |
| V5d dist→minmax | 66.3% | 0.508 | 57.5% | 0.512 |
| V6d dist→ratio | 66.7% | 0.508 | 57.9% | **1.000** |
| V7d dist→mutual→ratio | 66.7% | 0.508 | 57.3% | 0.381 |

**Conclusions, stated carefully:**

1. **No variant beats today's score significantly** (every p > 0.2), and no variant beats another.
   The four distinct-primary variants differ on **0–1 of 95** STRONG cases — on this evidence they
   are the same classifier. The formula choice must be settled on which objective is *correct*, not
   on measured accuracy.
2. **This is low power, not proof of equivalence.** n=95 STRONG cases cannot detect a 4-point
   difference. Absence of significance is not evidence of absence.
3. **V6d changes 625 WEAK decisions and nets exactly zero** (313 gained, 312 lost, p=1.000). The
   churn on unlabelled duplicates is noise, not improvement.
4. **Every variant gets ~1/3 of checkable cases wrong.** Geometry-only scoring is only weakly
   aligned with identification evidence. That ceiling is a property of the approach, not of the
   tie-break.
5. Caveat: cycles are not independent samples — each re-pairs the previous cycle's deletions.

## pep_ident conflicts — measured at the WRONG stage first, corrected here

⚠️ **Correction.** An earlier version of this section reported "20.4% of identified components merge
two or more peptides" and called identification-aware filtering a larger lever than the scoring
change. That number came from the **pre-filter connected components inside a single cycle** — the
deliberately over-merged intermediate state that the cycle loop exists to resolve. It is not an
output-quality figure and should not be read as one. The claim is withdrawn.

**Measured on the actual output** (`final_conflicts.py`, `final_components_all_cycles.json`,
159727 components over 74 cycles):

| | |
|---|---|
| components with no identifications | 106775 (66.8%) |
| identified, unanimous | 45032 (**85.0%**) |
| identified, conflicting | 7920 (15.0%) |

Component-level conflict rate overstates the damage, because conflicting components are large
(median 23 vertices vs 15 unanimous) and the minority is one or two vertices. Per vertex:

```
labelled vertices in final components : 474933
  in their component's majority peptide: 465317  (97.98%)
  in a minority peptide                :   9616  ( 2.02%)
```

Conflict rate by component size — consistent with stray intruders, not genuine merges:
4.3% at size 2, 4.9% at 3-5, 9.8% at 6-10, 15.0% at 11-20, 22.2% at 21+; median minority
fraction 0.056.

**So the cycle loop does most of this work already.** Residual contamination is ~2% of identified
vertices, part of which will be MS2 identification error rather than pairing error — the two are not
separable from this data. There is no demonstrated lever here bigger than the scoring change.

Two supporting facts from the per-cycle analysis, still valid on their own terms (cycle 3, ppm=5
run, pre-filter): removing the minority leaves the majority connected in 95.6% of conflicting
components, and cross-peptide edges are 2.88x longer at the median (7.33e-06 vs 2.55e-06) but the
distributions overlap heavily (72.6% / 81.0%) — so a distance cut cannot separate them cleanly.

⚠️ The two work dirs used here are **different runs** (`work/c1/ed55...` is ppm=5,
`work/5a/c8cc...` is ppm=10). Do not read any before/after difference between them as a cycle-loop
effect.

## Facts about the data (cycle 3)

- edge `distance` range across components: 4.497e-09 .. 1.448e-05; median within-component span
  1.174e-05.
- **0** components where all edge distances are equal (so no 0/0 in practice — still guard it).
- **0** components with a single edge.

## Constraints before touching anything

0. **A pipeline run was ACTIVE when this comparison was made** (nextflow, `mode_first_mean_...
   _ppm10_rt100`, started 11:53). Nothing in `bin/` was touched — the comparison is a standalone
   scratchpad script that only imports the module. Confirm the run has finished before editing.
1. **This change alters output.** 1066 vertices move between kept and deleted on cycle 3 alone, and
   deleted vertices feed the next cycle's pairing, so it cascades across the 500-cycle loop. The
   per-cycle delta understates the end effect — measure across several cycles.
2. **Do not edit `bin/*.py` while a pipeline run is in progress.** Nextflow invokes scripts by
   absolute path into the mounted workspace, so a running job picks up edits mid-run and different
   cycles would execute different code.
3. `best_score` becomes a tuple under Option 1; the deleted-vertex record does `float(best_score)`
   and would raise. Store the parts as separate fields (`retained_edges`, `degree_sum`). Those
   records never reach disk, so this is free. **The user wants those records KEPT** — they are
   intentional for debugging.
4. Comments must never reference removed/previous code. History goes in handoffs, not source.

## Verification setup

- `harness.sh` (session scratchpad) — 13-run matrix, md5s 36 JSON outputs. **The scratchpad gets
  wiped between sessions; recreate it.** Full description in `CLUSTERLESS_CLEANUP_HANDOFF.md`.
- ⚠️ The old data dir is gone. Current one with cycles 20/30/45/5/3:
  `work/c1/ed55a3f9ae28b013320e845be68e61`
- Snapshot the script before editing so you can A/B against it (`SCRIPT=<path> ./harness.sh <label>`).
- Expect the harness to FAIL on this change — that is the point. Compare structurally (which
  vertices moved, how many edges retained), not by md5.

## Also still open, unrelated to scoring

The `/simplify` pass produced a triage that was never applied — see the end of
`CLUSTERLESS_CLEANUP_HANDOFF.md`. The one worth doing is hoisting degree lookups out of the
exhaustive candidate loop (33,600 igraph calls to retrieve 27 constants on one measured component).
**Note it becomes moot under a retained-edge score** — that objective genuinely varies per
candidate, so the work is real rather than redundant. Settle scoring first.
