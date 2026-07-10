#!/usr/bin/env python3
"""
Generate summary report for clusterless iterative component building workflow.

This simplified version reports on:
- Cycle convergence (iterations until no more deletions)
- Component size distribution (unified across all cycles)
- Peptide identification analysis and conflict detection
- File linkage analysis
- Basic statistics

Removed sections (not applicable to clusterless workflow):
- Cluster size distributions
- Resolution optimization results
- Clustering-specific metrics
"""

import argparse
import ast
import csv
import json
import os
import random
from collections import Counter
from itertools import combinations
import statistics
from typing import Dict, List, Optional, Tuple
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import base64
from io import BytesIO

# ---------------------------------------------------------------------------
# Monoisotopic residue masses (Da) — standard Unimod/NIST values
# ---------------------------------------------------------------------------
_AA_MASSES: Dict[str, float] = {
    'G': 57.021464, 'A': 71.037114, 'V': 99.068414, 'L': 113.084064,
    'I': 113.084064, 'P': 97.052764, 'F': 147.068414, 'W': 186.079313,
    'M': 131.040485, 'S': 87.032028, 'T': 101.047679, 'C': 103.009185,
    'Y': 163.063329, 'H': 137.058912, 'D': 115.026943, 'E': 129.042593,
    'N': 114.042927, 'Q': 128.058578, 'K': 128.094963, 'R': 156.101111,
}

_NTERM_MOD_MASSES: Dict[str, float] = {
    'none':        0.0,
    'acetyl':      42.010565,   # N-acetylation (+CH2CO)
    'biotin':      226.077598,  # N-biotinyl label
    'fluorescein': 358.047876,  # 5-FAM N-terminal label
}

_CYS_MOD_MASSES: Dict[str, float] = {
    'none':            0.0,
    'oxidized':        15.994915,  # Cys sulphoxide (+O per Cys residue)
    'carbamidomethyl': 57.021464,  # IAA alkylation (+CH2CONH2 per Cys)
}

_CTERM_MOD_MASSES: Dict[str, float] = {
    'none':  0.0,
    'amide': -0.984016,  # C-terminal NH2 instead of OH (16.018724 - 17.002740)
}

_H2O_MASS = 18.010565  # H (1.007825) + OH (17.002740)


def compute_monoisotopic_mass(sequence: str,
                              nterm_mod: str = 'none',
                              cys_mod: str = 'none',
                              cterm_mod: str = 'none') -> float:
    """Return the neutral monoisotopic mass for a peptide sequence.

    Handles inline bracket modifications (e.g. M[15.9949]) from raw_pep_ident.
    When a residue carries an explicit [delta], that delta is used directly and
    cys_mod is NOT applied to that residue (it already has a fixed modification).

    nterm_mod : none | acetyl | biotin | fluorescein
    cys_mod   : none | oxidized | carbamidomethyl  (fallback for unmodified Cys)
    cterm_mod : none | amide  (C-terminal -CONH2 instead of -COOH)
    """
    nterm_delta = _NTERM_MOD_MASSES.get(nterm_mod.lower(), 0.0)
    cys_delta   = _CYS_MOD_MASSES.get(cys_mod.lower(), 0.0)
    cterm_delta = _CTERM_MOD_MASSES.get(cterm_mod.lower(), 0.0)
    mass = _H2O_MASS + nterm_delta + cterm_delta
    seq = sequence.upper()
    i = 0
    while i < len(seq):
        aa = seq[i]
        i += 1
        residue = _AA_MASSES.get(aa)
        if residue is None:
            # skip non-AA chars; if we land on '[' from a previous skipped char, consume it
            if aa == '[':
                close = seq.find(']', i)
                if close != -1:
                    i = close + 1
            continue
        # Check for inline modification bracket immediately after residue
        if i < len(seq) and seq[i] == '[':
            close = seq.find(']', i + 1)
            if close != -1:
                try:
                    residue += float(seq[i + 1:close])
                except ValueError:
                    pass
                i = close + 1
        else:
            # No inline mod — apply global cys_mod to unmodified Cys
            if aa == 'C':
                residue += cys_delta
        mass += residue
    return mass


def build_raw_pep_ident_lookup(tsv_paths: List[str]) -> Dict[Tuple[str, str], List[str]]:
    """Scan TSV files once and build a compact lookup for raw pep_idents.

    Returns {(filename_stem, openms_fid): [raw_pep_ident, ...]}
    Only entries with at least one non-empty raw_pep_ident are stored.
    TSV file contents are not retained after this function returns.
    """
    csv.field_size_limit(2 ** 31 - 1)
    lookup: Dict[Tuple[str, str], List[str]] = {}
    for path in tsv_paths:
        stem = os.path.splitext(os.path.basename(path))[0]
        try:
            with open(path, newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                for row in reader:
                    raw_str = row.get('l_raw_pep_ident', '[]') or '[]'
                    if raw_str == '[]':
                        continue
                    try:
                        seqs = ast.literal_eval(raw_str)
                    except (ValueError, SyntaxError):
                        continue
                    seqs = [s for s in seqs if s]
                    if seqs:
                        lookup[(stem, row['openms_fid'])] = seqs
        except (OSError, KeyError):
            pass
    return lookup


def detect_conflicting_pep_idents_in_component(vertices: List[Dict]) -> bool:
    """
    Detect if a component has conflicting pep_idents.
    
    A conflict exists if vertices with pep_idents don't share at least one common pep_ident.
    Single vertices or vertices without pep_idents don't have conflicts.
    
    Logic:
    - Collect pep_idents from all vertices (handle both list and string types)
    - Find the intersection of all pep_ident sets
    - If intersection is non-empty: NOT conflicting (all share at least one)
    - If intersection is empty: CONFLICTING (no common pep_ident)
    
    Args:
        vertices: List of vertex dictionaries with 'pep_ident' field
        
    Returns:
        True if conflicts exist, False otherwise
    """
    # Collect pep_ident sets for each vertex that has pep_idents
    pep_ident_sets = []
    
    for vertex in vertices:
        pep_ident = vertex.get('pep_ident')
        if pep_ident:
            if isinstance(pep_ident, list):
                # Convert list to set
                pep_set = set(pep_ident) if pep_ident else set()
            else:
                # Single string pep_ident
                pep_set = {pep_ident}
            
            if pep_set:  # Only add non-empty sets
                pep_ident_sets.append(pep_set)
    
    # If 0 or 1 vertices have pep_idents, no conflict possible
    if len(pep_ident_sets) <= 1:
        return False
    
    # Find intersection of all pep_ident sets (what they all share)
    common_pep_idents = pep_ident_sets[0]
    for pep_set in pep_ident_sets[1:]:
        common_pep_idents = common_pep_idents.intersection(pep_set)
        if not common_pep_idents:
            # Early exit if no common pep_idents found
            break
    
    # If there's any common pep_ident, no conflict
    # If intersection is empty, there's a conflict
    return not common_pep_idents


def analyze_pep_idents(components: List[Dict],
                       nterm_mod: str = 'none',
                       cys_mod: str = 'none',
                       cterm_mod: str = 'none',
                       raw_pep_ident_lookup: Optional[Dict[Tuple[str, str], List[str]]] = None) -> Dict:
    """
    Analyze pep_ident assignments and conflicts within components.

    Returns:
        Dictionary with pep_ident statistics:
        - total_components, components_with_pep_idents, components_with_conflicts,
          conflict_percentage, conflict_examples
        - mass_diff_stats: {mean, median, n_components} or None — pairwise monoisotopic
          mass differences between conflicting sequences across all conflicted components
    """
    _empty: Dict = {
        'total_components': 0,
        'components_with_pep_idents': 0,
        'components_with_conflicts': 0,
        'conflict_percentage': 0.0,
        'conflict_examples': [],
        'mass_diff_stats': None,
        'all_pairwise_diffs': [],
        'multi_component_stats': None,
        'vertex_conflict_stats': None,
    }
    if not components:
        return _empty

    components_with_pep_idents = 0
    components_with_conflicts = 0
    conflict_examples: List[Dict] = []
    all_component_mean_diffs: List[float] = []
    all_pairwise_diffs: List[float] = []
    pep_ident_to_comp_ids: Dict[str, set] = {}

    # Vertex-level conflict counters
    total_vertices_all = 0
    total_vertices_with_pid = 0
    minimal_conflict_vertices = 0   # non-dominant pep_ident in a conflicted component
    moderate_conflict_vertices = 0  # any pep_ident vertex in a conflicted component
    maximal_conflict_vertices = 0   # every vertex (incl. unidentified) in a conflicted component

    for comp in components:
        comp_id = comp.get('component_id', 'unknown')
        vertices = comp.get('vertices', [])
        total_vertices_all += len(vertices)

        # Per-vertex: has any non-empty pep_ident?
        v_with_pid = [v for v in vertices if v.get('pep_ident')]
        total_vertices_with_pid += len(v_with_pid)

        if not v_with_pid:
            continue

        components_with_pep_idents += 1

        # Collect unique non-empty pep_ident strings for the whole component
        raw: List[str] = []
        for v in v_with_pid:
            pid = v.get('pep_ident')
            if isinstance(pid, list):
                raw.extend(pid)
            else:
                raw.append(pid)
        unique_seqs = list({s for s in raw if s})

        # Track cross-component pep_ident appearances
        for pid in unique_seqs:
            pep_ident_to_comp_ids.setdefault(pid, set()).add(comp_id)

        # Build clean→raw mapping for this component using the TSV lookup.
        # For each vertex, match its raw pep_idents back to the clean pep_idents
        # we already collected so mass calculation can use the modified form.
        clean_to_raw: Dict[str, str] = {}
        if raw_pep_ident_lookup:
            for v in v_with_pid:
                fname = v.get('filename', '')
                fid   = v.get('openms_fid', '')
                raw_seqs = raw_pep_ident_lookup.get((fname, fid), [])
                clean_pids = v.get('pep_ident', [])
                if isinstance(clean_pids, str):
                    clean_pids = [clean_pids]
                for clean in clean_pids:
                    if clean and clean not in clean_to_raw:
                        # Find the matching raw seq (first one whose stripped form equals clean)
                        for raw_seq in raw_seqs:
                            stripped = ''.join(
                                c for c in raw_seq if c.isalpha()
                            ) if '[' in raw_seq else raw_seq
                            if stripped.upper() == clean.upper():
                                clean_to_raw[clean] = raw_seq
                                break

        # Conflict analysis
        if detect_conflicting_pep_idents_in_component(vertices):
            components_with_conflicts += 1

            seq_masses = {
                seq: compute_monoisotopic_mass(
                    clean_to_raw.get(seq, seq), nterm_mod, cys_mod, cterm_mod
                )
                for seq in unique_seqs
            }

            if len(unique_seqs) >= 2:
                pairwise = [
                    abs(seq_masses[a] - seq_masses[b])
                    for a, b in combinations(unique_seqs, 2)
                ]
                all_component_mean_diffs.append(statistics.mean(pairwise))
                all_pairwise_diffs.extend(pairwise)

            if len(conflict_examples) < 10:
                conflict_examples.append({
                    'component_id': comp_id,
                    'num_vertices': comp.get('num_vertices', 0),
                    'pep_idents': unique_seqs,
                    'pep_ident_masses': seq_masses,
                    'pep_ident_raw': {seq: clean_to_raw.get(seq, seq) for seq in unique_seqs},
                    'cycle_num': comp.get('cycle'),
                })

            # ── Vertex-level estimates ────────────────────────────────────────
            maximal_conflict_vertices += len(vertices)
            moderate_conflict_vertices += len(v_with_pid)

            # Dominant pep_ident(s): the one(s) appearing in the most vertices.
            # A vertex that has a pep_ident list is counted once per matching pid.
            pid_vertex_count: Counter = Counter()
            for v in v_with_pid:
                pid = v.get('pep_ident')
                pids = pid if isinstance(pid, list) else [pid]
                for p in pids:
                    if p:
                        pid_vertex_count[p] += 1

            if pid_vertex_count:
                max_count = max(pid_vertex_count.values())
                dominant_pids = {p for p, c in pid_vertex_count.items() if c == max_count}

                for v in v_with_pid:
                    pid = v.get('pep_ident')
                    v_pids = set(pid) if isinstance(pid, list) else ({pid} if pid else set())
                    v_pids = {p for p in v_pids if p}
                    if v_pids and not (v_pids & dominant_pids):
                        minimal_conflict_vertices += 1

    conflict_percentage = 0.0
    if components_with_pep_idents > 0:
        conflict_percentage = (components_with_conflicts / components_with_pep_idents) * 100

    mass_diff_stats = None
    if all_component_mean_diffs:
        mass_diff_stats = {
            'mean': statistics.mean(all_component_mean_diffs),
            'median': statistics.median(all_component_mean_diffs),
            'n_components': len(all_component_mean_diffs),
        }

    # Multi-component redundancy
    multi_component_stats = None
    if pep_ident_to_comp_ids:
        distribution: Dict[int, int] = {}
        for pid, comp_ids in pep_ident_to_comp_ids.items():
            n = len(comp_ids)
            if n >= 2:
                distribution[n] = distribution.get(n, 0) + 1
        # Top-20 most redundant pep_idents (highest component count first)
        multi_sorted = sorted(
            ((pid, len(cids)) for pid, cids in pep_ident_to_comp_ids.items() if len(cids) >= 2),
            key=lambda x: x[1], reverse=True,
        )
        multi_component_stats = {
            'total_unique_pep_idents': len(pep_ident_to_comp_ids),
            'total_multi_component': sum(distribution.values()),
            'distribution': distribution,
            'top_redundant': multi_sorted[:20],  # list of (pep_ident, n_components)
        }

    # Vertex-level conflict stats
    ident_denom = total_vertices_with_pid or 1  # avoid div-by-zero
    all_denom   = total_vertices_all or 1
    vertex_conflict_stats = {
        'total_vertices': total_vertices_all,
        'total_with_pep_ident': total_vertices_with_pid,
        'pct_with_pep_ident': total_vertices_with_pid / all_denom * 100,
        'minimal': minimal_conflict_vertices,
        'minimal_pct': minimal_conflict_vertices / ident_denom * 100,
        'moderate': moderate_conflict_vertices,
        'moderate_pct': moderate_conflict_vertices / ident_denom * 100,
        'maximal': maximal_conflict_vertices,
        'maximal_pct': maximal_conflict_vertices / ident_denom * 100,
    }

    return {
        'total_components': len(components),
        'components_with_pep_idents': components_with_pep_idents,
        'components_with_conflicts': components_with_conflicts,
        'conflict_percentage': conflict_percentage,
        'conflict_examples': conflict_examples,
        'mass_diff_stats': mass_diff_stats,
        'all_pairwise_diffs': all_pairwise_diffs,
        'multi_component_stats': multi_component_stats,
        'vertex_conflict_stats': vertex_conflict_stats,
    }


def analyze_feature_masses(components: List[Dict]) -> Dict:
    """Compute mass-based statistics across all components.

    Feature mass = x_center * charge  (as directed; no proton subtraction).
    Returns three datasets:
      all_pairwise_diffs  – absolute mass diffs between every vertex pair within a component
      component_mean_diffs – per-component mean of those pairwise diffs
      nn_diffs            – distance from each component's mean mass to its nearest neighbour
    """
    all_pairwise_diffs: List[float] = []
    component_mean_diffs: List[float] = []
    component_mean_masses: List[float] = []

    for comp in components:
        vertices = comp.get('vertices', [])
        masses = []
        for v in vertices:
            x = v.get('x_center')
            z = v.get('charge')
            if x is not None and z is not None and z > 0:
                masses.append(float(x) * float(z))

        if not masses:
            continue

        component_mean_masses.append(sum(masses) / len(masses))

        if len(masses) >= 2:
            pairs = [abs(masses[i] - masses[j])
                     for i in range(len(masses))
                     for j in range(i + 1, len(masses))]
            all_pairwise_diffs.extend(pairs)
            component_mean_diffs.append(sum(pairs) / len(pairs))

    # Nearest neighbour by mean mass (sort → adjacent = exact NN in 1D)
    sorted_masses = sorted(component_mean_masses)
    nn_diffs: List[float] = []
    for i, m in enumerate(sorted_masses):
        cands = []
        if i > 0:
            cands.append(abs(m - sorted_masses[i - 1]))
        if i < len(sorted_masses) - 1:
            cands.append(abs(m - sorted_masses[i + 1]))
        if cands:
            nn_diffs.append(min(cands))

    def _stats(vals):
        if not vals:
            return None
        return {'mean': statistics.mean(vals), 'median': statistics.median(vals), 'n': len(vals)}

    return {
        'all_pairwise_diffs':  all_pairwise_diffs,
        'pairwise_stats':      _stats(all_pairwise_diffs),
        'component_mean_diffs': component_mean_diffs,
        'mean_diff_stats':     _stats(component_mean_diffs),
        'nn_diffs':            nn_diffs,
        'nn_stats':            _stats(nn_diffs),
        'n_components_with_masses': len(component_mean_masses),
    }


# ---------------------------------------------------------------------------
# Shared canvas-based histogram JS (emitted once per report)
# ---------------------------------------------------------------------------
# Shared canvas histogram renderer — emitted once per report.
# Works on pre-computed bin counts+edges (no raw data in JS).
_HISTOGRAM_JS = r"""
<script>
function renderHistogram(cfg) {
    var canvas = document.getElementById(cfg.cid);
    if (!canvas) return;
    var ctx = canvas.getContext('2d');

    var nBins   = parseInt(document.getElementById(cfg.bid).value);
    var xmaxPct = parseFloat(document.getElementById(cfg.xid).value);

    /* Resolve xMax from pre-computed percentile table */
    var pctIdx = Math.min(Math.round(xmaxPct * 2), cfg.pctEdges.length - 1);
    var xMin = cfg.edges[0];
    var xMax = cfg.pctEdges[pctIdx];
    if (!xMax || xMax <= xMin) xMax = cfg.edges[cfg.edges.length - 1];
    var xRange = xMax - xMin;

    /* Adaptive tick formatter — picks decimal places from the axis range */
    function fmtX(v) {
        if (xRange === 0) return v.toFixed(2);
        var mag = Math.floor(Math.log10(xRange));
        if (mag >= 3)  return v.toFixed(0);
        if (mag >= 1)  return v.toFixed(1);
        if (mag >= 0)  return v.toFixed(2);
        if (mag >= -1) return v.toFixed(3);
        if (mag >= -2) return v.toFixed(4);
        return v.toExponential(2);
    }

    document.getElementById(cfg.bvid).textContent = nBins;
    document.getElementById(cfg.xvid).textContent = fmtX(xMax);

    /* Find last pre-bin fully within [xMin, xMax] */
    var nPre = cfg.counts.length;
    var lastPre = nPre - 1;
    for (var i = 0; i < nPre; i++) {
        if (cfg.edges[i + 1] > xMax) { lastPre = i - 1; break; }
    }
    lastPre = Math.max(0, lastPre);
    var nVisible = lastPre + 1;

    /* Count clipped values */
    var nOver = 0;
    for (var i = nVisible; i < nPre; i++) nOver += cfg.counts[i];

    /* Re-bin: group nVisible pre-bins into nBins display bars */
    var grp = Math.max(1, Math.ceil(nVisible / nBins));
    var bars = [];
    for (var i = 0; i < nVisible; i += grp) {
        var ct = 0;
        for (var j = i; j < Math.min(i + grp, nVisible); j++) ct += cfg.counts[j];
        bars.push(ct);
    }
    var nBars  = bars.length;
    var maxCt  = Math.max.apply(null, bars) || 1;

    /* Resize canvas to current CSS width */
    var W = canvas.offsetWidth || 840;
    canvas.width  = W;
    canvas.height = 380;
    var H   = canvas.height;
    var pad = {l:72, r:24, t:40, b:54};
    var pw  = W - pad.l - pad.r;
    var ph  = H - pad.t - pad.b;

    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, W, H);

    /* Y grid */
    ctx.strokeStyle = '#e5e5e5'; ctx.lineWidth = 1;
    for (var g = 0; g <= 5; g++) {
        var gy = pad.t + ph * (1 - g / 5);
        ctx.beginPath(); ctx.moveTo(pad.l, gy); ctx.lineTo(pad.l + pw, gy); ctx.stroke();
        ctx.fillStyle = '#666'; ctx.font = '11px Arial'; ctx.textAlign = 'right';
        ctx.fillText(Math.round(maxCt * g / 5).toLocaleString(), pad.l - 5, gy + 4);
    }

    /* Bars */
    var bpx = pw / nBars;
    for (var i = 0; i < nBars; i++) {
        if (!bars[i]) continue;
        var bh = (bars[i] / maxCt) * ph;
        ctx.fillStyle = '#3498db';
        ctx.fillRect(pad.l + i * bpx + 0.5, pad.t + ph - bh, Math.max(bpx - 1, 0.5), bh);
    }

    /* Axes box */
    ctx.strokeStyle = '#999'; ctx.lineWidth = 1;
    ctx.strokeRect(pad.l, pad.t, pw, ph);

    /* Mean / Median */
    function vline(val, color, lbl) {
        if (isNaN(val)) return;
        var x = pad.l + ((val - xMin) / (xMax - xMin)) * pw;
        if (x < pad.l || x > pad.l + pw) return;
        ctx.save();
        ctx.strokeStyle = color; ctx.lineWidth = 1.5;
        ctx.setLineDash([5, 3]);
        ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = color; ctx.font = 'bold 10px Arial'; ctx.textAlign = 'center';
        ctx.fillText(lbl + ': ' + fmtX(val), x, pad.t - 8);
        ctx.restore();
    }
    vline(cfg.mean,   '#e74c3c', 'mean');
    vline(cfg.median, '#27ae60', 'median');

    /* X ticks */
    ctx.fillStyle = '#555'; ctx.font = '11px Arial'; ctx.textAlign = 'center';
    for (var t = 0; t <= 7; t++) {
        var tv = xMin + xRange * t / 7;
        ctx.fillText(fmtX(tv), pad.l + pw * t / 7, pad.t + ph + 16);
    }

    /* Axis labels */
    ctx.fillStyle = '#333'; ctx.font = '13px Arial'; ctx.textAlign = 'center';
    ctx.fillText(cfg.xlabel, pad.l + pw / 2, pad.t + ph + 40);
    ctx.save();
    ctx.translate(14, pad.t + ph / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('Count', 0, 0);
    ctx.restore();

    /* Clipped note */
    if (nOver > 0) {
        ctx.fillStyle = '#aaa'; ctx.font = '10px Arial'; ctx.textAlign = 'right';
        ctx.fillText(nOver.toLocaleString() + ' values clipped', pad.l + pw - 4, pad.t + ph + 40);
    }
}
</script>
"""


def _interactive_histogram_html(
    data: List[float],
    hist_id: str,
    title: str,
    xlabel: str,
    n_total: int,
) -> str:
    """Return an HTML block for one interactive canvas histogram.

    All data is represented: numpy pre-computes a 2000-bin histogram in Python;
    only the counts + edges (~4000 numbers) are embedded in the HTML.
    The JS re-bins those on the fly — no raw values, no sampling.
    """
    if not data:
        return f'<p><em>No data available for: {title}</em></p>'

    arr = np.asarray(data, dtype=np.float64)
    mean_val   = float(np.mean(arr))
    median_val = float(np.median(arr))

    # Pre-compute high-resolution histogram
    counts, edges = np.histogram(arr, bins=2000)
    counts_js = json.dumps(counts.tolist())
    edges_js  = json.dumps([round(float(e), 6) for e in edges])

    # Pre-compute percentile lookup (201 values: 0 % to 100 % in 0.5 % steps)
    pct_vals = np.percentile(arr, np.linspace(0, 100, 201))
    pct_js   = json.dumps([round(float(v), 6) for v in pct_vals])

    # Default max-X slider at 99th percentile = index 198 (198/2 = 99 %)
    default_pct = 99.0

    cid  = f'{hist_id}-c'
    bid  = f'{hist_id}-b'
    bvid = f'{hist_id}-bv'
    xid  = f'{hist_id}-x'
    xvid = f'{hist_id}-xv'

    return f"""
<div style="margin:20px 0; padding:12px; background:#fafafa; border:1px solid #ddd; border-radius:6px;">
  <h4 style="margin:0 0 6px">{title}</h4>
  <p style="margin:0 0 10px;font-size:.95em">
    Mean:&nbsp;<strong>{mean_val:.3f}</strong>&nbsp;&nbsp;
    Median:&nbsp;<strong>{median_val:.3f}</strong>&nbsp;&nbsp;
    N:&nbsp;<strong>{n_total:,}</strong>
  </p>
  <div style="display:flex;gap:28px;flex-wrap:wrap;margin-bottom:8px;font-size:.9em;">
    <label>Bins:&nbsp;<strong id="{bvid}">80</strong>
      <input type="range" id="{bid}" min="5" max="500" value="80"
             style="width:140px;vertical-align:middle" oninput="_rh_{hist_id}()">
    </label>
    <label>Max&nbsp;X&nbsp;(percentile):&nbsp;<strong id="{xvid}">–</strong>
      <input type="range" id="{xid}" min="0" max="100" step="0.5" value="{default_pct}"
             style="width:140px;vertical-align:middle" oninput="_rh_{hist_id}()">
    </label>
  </div>
  <canvas id="{cid}" style="width:100%;height:380px;display:block;"></canvas>
</div>
<script>
(function(){{
  var _cfg = {{
    cid:'{cid}', bid:'{bid}', bvid:'{bvid}', xid:'{xid}', xvid:'{xvid}',
    counts:{counts_js},
    edges:{edges_js},
    pctEdges:{pct_js},
    mean:{mean_val}, median:{median_val},
    xlabel:{json.dumps(xlabel)}
  }};
  window._rh_{hist_id} = function(){{ renderHistogram(_cfg); }};
  setTimeout(window._rh_{hist_id}, 0);
  window.addEventListener('resize', window._rh_{hist_id});
}})();
</script>
"""


def generate_feature_mass_html_section(feature_mass_stats: Dict) -> str:
    """Return the full HTML section for feature mass analysis (3 interactive histograms)."""
    if not feature_mass_stats:
        return ''

    parts = []
    parts.append('<h2>⚖️ Feature Mass Analysis</h2>')
    parts.append(
        f'<p>Components with mass data: '
        f'<strong>{feature_mass_stats["n_components_with_masses"]:,}</strong> &nbsp;|&nbsp; '
        f'Mass&nbsp;=&nbsp;x_center&nbsp;×&nbsp;charge</p>'
    )
    parts.append(_HISTOGRAM_JS)

    # --- Histogram 1: all pairwise diffs ---
    pd = feature_mass_stats['all_pairwise_diffs']
    ps = feature_mass_stats['pairwise_stats']
    parts.append('<h3>1 · All Pairwise Feature Mass Differences Within Components</h3>')
    if ps:
        parts.append(_interactive_histogram_html(
            pd, 'fmh1',
            'Absolute mass difference between each feature pair in the same component',
            'Mass difference (Da)',
            ps['n'],
        ))
    else:
        parts.append('<p><em>No multi-vertex components found.</em></p>')

    # --- Histogram 2: per-component mean diff ---
    md = feature_mass_stats['component_mean_diffs']
    ms = feature_mass_stats['mean_diff_stats']
    parts.append('<h3>2 · Mean Mass Difference per Component</h3>')
    if ms:
        parts.append(_interactive_histogram_html(
            md, 'fmh2',
            'Mean of all pairwise mass differences within each component',
            'Mean mass difference (Da)',
            ms['n'],
        ))
    else:
        parts.append('<p><em>No multi-vertex components found.</em></p>')

    # --- Histogram 3: nearest-neighbour component mass diff ---
    nd = feature_mass_stats['nn_diffs']
    ns = feature_mass_stats['nn_stats']
    parts.append('<h3>3 · Nearest-Neighbour Component Mass Difference</h3>')
    parts.append('<p style="font-size:.9em;color:#555">Distance from each component\'s mean mass '
                 'to its closest neighbouring component (sorted by mean mass, adjacent = exact NN in 1D).</p>')
    if ns:
        parts.append(_interactive_histogram_html(
            nd, 'fmh3',
            'Mass difference from each component to its nearest neighbour (by mean mass)',
            'Mass difference (Da)',
            ns['n'],
        ))
    else:
        parts.append('<p><em>Not enough components.</em></p>')

    return ''.join(parts)


def load_final_components(components_json: str) -> Dict:
    """Load final components JSON from orchestrator."""
    if not os.path.exists(components_json):
        raise FileNotFoundError(f"Components JSON not found: {components_json}")
    
    with open(components_json, 'r') as f:
        data = json.load(f)
    
    print(f"✓ Loaded components from: {os.path.basename(components_json)}")
    return data


def load_cycle_history(history_json: str):
    """Load cycle history JSON from orchestrator. Returns (cycles, metadata)."""
    if not os.path.exists(history_json):
        print(f"⚠ Cycle history not found: {history_json}")
        return [], {}

    with open(history_json, 'r') as f:
        data = json.load(f)

    cycles = data.get('cycles', [])
    metadata = data.get('metadata', {})
    print(f"✓ Loaded cycle history: {len(cycles)} cycles")
    return cycles, metadata


def analyze_component_sizes(components: List[Dict]) -> Dict:
    """
    Analyze component size distribution.
    
    Returns:
        Dictionary with statistics and size distribution
    """
    if not components:
        return {
            'total_components': 0,
            'total_vertices': 0,
            'avg_size': 0.0,
            'median_size': 0.0,
            'min_size': 0,
            'max_size': 0,
            'size_distribution': {},
            'singletons': 0
        }
    
    sizes = [comp['num_vertices'] for comp in components]
    sizes.sort()
    
    # Build size distribution histogram
    size_dist = {}
    for size in sizes:
        size_dist[size] = size_dist.get(size, 0) + 1
    
    # Calculate statistics
    total_components = len(components)
    total_vertices = sum(sizes)
    avg_size = total_vertices / total_components if total_components > 0 else 0.0
    median_size = sizes[len(sizes) // 2] if sizes else 0.0
    min_size = min(sizes) if sizes else 0
    max_size = max(sizes) if sizes else 0
    singletons = size_dist.get(1, 0)
    
    return {
        'total_components': total_components,
        'total_vertices': total_vertices,
        'avg_size': avg_size,
        'median_size': median_size,
        'min_size': min_size,
        'max_size': max_size,
        'size_distribution': size_dist,
        'singletons': singletons
    }


def analyze_file_linkage(components: List[Dict]) -> Dict:
    """
    Analyze how well files are linked through components.
    
    Returns:
        Dictionary with file linkage statistics
    """
    # Collect all unique files
    all_files = set()
    for comp in components:
        for vertex in comp.get('vertices', []):
            filename = vertex.get('filename')
            if filename:
                all_files.add(filename)
    
    total_files = len(all_files)
    if total_files == 0:
        return {
            'total_files': 0,
            'file_linkage_scores': {},
            'avg_linkage_score': 0.0
        }
    
    # Calculate linkage score for each file
    # Linkage score = proportion of components that contain this file
    file_linkage_scores = {}
    
    for target_file in all_files:
        components_with_file = 0
        for comp in components:
            files_in_comp = set(v.get('filename') for v in comp.get('vertices', []))
            if target_file in files_in_comp:
                components_with_file += 1
        
        # Score is percentage of components containing this file
        score = components_with_file / len(components) if components else 0.0
        file_linkage_scores[target_file] = score
    
    avg_linkage_score = sum(file_linkage_scores.values()) / len(file_linkage_scores) if file_linkage_scores else 0.0
    
    return {
        'total_files': total_files,
        'file_linkage_scores': file_linkage_scores,
        'avg_linkage_score': avg_linkage_score
    }


def plot_component_size_distribution(size_dist: Dict[int, int]) -> str:
    """
    Generate bar plot of component size distribution.
    
    Returns:
        Base64-encoded PNG image
    """
    if not size_dist:
        return ""
    
    # Prepare data
    sizes = sorted(size_dist.keys())
    counts = [size_dist[s] for s in sizes]
    vertices = [s * size_dist[s] for s in sizes]  # total vertices at each size

    fontsize = max(6, min(9, 120 // max(len(sizes), 1)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # --- Top: number of components per size ---
    bars1 = ax1.bar(sizes, counts, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.set_ylabel('Number of Components')
    ax1.set_title('Component Size Distribution — Components per Size')
    ax1.grid(axis='y', alpha=0.3)
    for bar, count in zip(bars1, counts):
        ax1.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height(),
            f'{count:,}',
            ha='center', va='bottom', fontsize=fontsize, clip_on=True
        )

    # --- Bottom: total vertices per size ---
    bars2 = ax2.bar(sizes, vertices, color='darkorange', edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Component Size (number of vertices)')
    ax2.set_ylabel('Total Vertices')
    ax2.set_title('Component Size Distribution — Vertices per Size')
    ax2.grid(axis='y', alpha=0.3)
    for bar, v in zip(bars2, vertices):
        ax2.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height(),
            f'{v:,}',
            ha='center', va='bottom', fontsize=fontsize, clip_on=True
        )

    if len(sizes) > 50:
        ax1.set_xlim(0, max(sizes) + 1)

    plt.tight_layout()

    # Convert to base64
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=100)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)

    return img_base64


def plot_cycle_convergence(cycles: List[Dict]) -> str:
    """
    Generate line plot showing convergence across cycles.
    
    Returns:
        Base64-encoded PNG image
    """
    if not cycles:
        return ""
    
    cycle_nums = [c['cycle'] for c in cycles]
    edges_counts = [c['edges_generated'] for c in cycles]
    deleted_counts = [c['vertices_deleted'] for c in cycles]
    accepted_counts = [c['components_accepted'] for c in cycles]
    
    # Create plot with dual y-axes
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # Plot edges and deleted vertices on left axis
    color1 = 'steelblue'
    ax1.set_xlabel('Cycle Number')
    ax1.set_ylabel('Count', color=color1)
    ax1.plot(cycle_nums, edges_counts, marker='o', label='Edges', color=color1, linewidth=2)
    ax1.plot(cycle_nums, deleted_counts, marker='s', label='Deleted Vertices', color='orangered', linewidth=2)
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(alpha=0.3)
    ax1.legend(loc='upper left')
    
    # Plot accepted components on right axis
    ax2 = ax1.twinx()
    color2 = 'forestgreen'
    ax2.set_ylabel('Accepted Components', color=color2)
    ax2.plot(cycle_nums, accepted_counts, marker='^', label='Accepted Components', color=color2, linewidth=2, linestyle='--')
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.legend(loc='upper right')
    
    plt.title('Cycle Convergence')
    plt.tight_layout()
    
    # Convert to base64
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=100)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    
    return img_base64


def plot_multi_component_histogram(distribution: Dict[int, int]) -> str:
    """
    Bar chart: for each N ≥ 2, how many distinct pep_idents appear in exactly N components.

    Returns:
        Base64-encoded PNG, or empty string if no data.
    """
    if not distribution:
        return ""

    ns     = sorted(distribution.keys())
    counts = [distribution[n] for n in ns]

    fig, ax = plt.subplots(figsize=(max(6, len(ns) * 0.6 + 2), 5))
    bars = ax.bar(ns, counts, color='mediumpurple', edgecolor='black', alpha=0.75)

    for bar, count in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height(),
            f'{count:,}',
            ha='center', va='bottom', fontsize=9
        )

    ax.set_xlabel('Number of Components the Pep_ident Appears In')
    ax.set_ylabel('Count of Distinct Pep_idents')
    ax.set_title('Pep_ident Redundancy — How Many Pep_idents Span Multiple Components')
    ax.set_xticks(ns)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()

    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=100)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return img_base64


def plot_mass_diff_histogram(all_pairwise_diffs: List[float]) -> str:
    """
    Generate a histogram of pairwise monoisotopic mass differences between
    conflicting pep_ident sequences across all conflicted components.

    Returns:
        Base64-encoded PNG image, or empty string if no data.
    """
    if not all_pairwise_diffs:
        return ""

    import numpy as np

    diffs = np.array(all_pairwise_diffs)
    mean_val   = float(np.mean(diffs))
    median_val = float(np.median(diffs))

    # Freedman-Diaconis bin edges, capped at 500 for very large datasets
    fd_edges = np.histogram_bin_edges(diffs, bins='fd')
    n_bins = min(500, max(20, len(fd_edges) - 1))

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.hist(diffs, bins=n_bins, color='steelblue', edgecolor='black', alpha=0.75)
    ax.axvline(mean_val,   color='orangered',  linestyle='--', linewidth=1.8,
               label=f'Mean {mean_val:.3f} Da')
    ax.axvline(median_val, color='forestgreen', linestyle=':',  linewidth=1.8,
               label=f'Median {median_val:.3f} Da')

    ax.set_xlabel('Pairwise Mass Difference (Da)')
    ax.set_ylabel('Count')
    ax.set_title(f'Distribution of Pairwise Mass Differences Between Conflicting Sequences  (n={len(diffs):,})')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()

    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=100)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)

    return img_base64


def generate_html_report(
    components_data: Dict,
    cycles: List[Dict],
    cycle_metadata: Dict,
    component_stats: Dict,
    file_linkage: Dict,
    pep_ident_stats: Dict,
    output_html: str,
    parameters: Dict = None,
    feature_mass_stats: Dict = None,
):
    """Generate HTML report."""
    
    # Generate plots
    size_dist_img = plot_component_size_distribution(component_stats.get('size_distribution', {}))
    convergence_img = plot_cycle_convergence(cycles) if cycles else ""
    
    # Build HTML
    html_parts = []
    
    # Header
    html_parts.append("""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Clusterless Component Building Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background-color: white; padding: 30px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; border-bottom: 2px solid #95a5a6; padding-bottom: 5px; }
        h3 { color: #7f8c8d; }
        .metric { display: inline-block; margin: 10px 20px 10px 0; padding: 10px 15px; background-color: #ecf0f1; border-radius: 5px; }
        .metric-label { font-weight: bold; color: #7f8c8d; }
        .metric-value { font-size: 1.2em; color: #2c3e50; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #bdc3c7; padding: 10px; text-align: left; }
        th { background-color: #34495e; color: white; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .plot { margin: 20px 0; text-align: center; }
        .plot img { max-width: 100%; height: auto; border: 1px solid #ddd; box-shadow: 0 0 5px rgba(0,0,0,0.1); }
        .success { color: #27ae60; font-weight: bold; }
        .warning { color: #f39c12; font-weight: bold; }
        .info { color: #3498db; font-weight: bold; }
    </style>
</head>
<body>
<div class="container">
""")
    
    # Title
    html_parts.append("<h1>🔄 Clusterless Component Building Report</h1>")
    
    initial_vertices = cycle_metadata.get('initial_vertices')
    unpaired_vertices_final = component_stats['singletons']
    vertices_after_cutoff = (initial_vertices - unpaired_vertices_final) if initial_vertices is not None else None

    # Overview metrics
    html_parts.append("<h2>📊 Overview</h2>")
    html_parts.append(f"""
<div class="metric">
    <span class="metric-label">Total Cycles:</span>
    <span class="metric-value">{len(cycles)}</span>
</div>""")
    if initial_vertices is not None:
        pct = 100 * vertices_after_cutoff / initial_vertices if initial_vertices else 0
        html_parts.append(f"""
<div class="metric">
    <span class="metric-label">Vertices Before Cutoff (input):</span>
    <span class="metric-value">{initial_vertices:,}</span>
</div>
<div class="metric">
    <span class="metric-label">Unpaired Vertices (final):</span>
    <span class="metric-value">{unpaired_vertices_final:,}</span>
</div>
<div class="metric">
    <span class="metric-label">Vertices After Cutoff:</span>
    <span class="metric-value">{vertices_after_cutoff:,} ({pct:.1f}%)</span>
</div>""")
    html_parts.append(f"""
<div class="metric">
    <span class="metric-label">Final Components:</span>
    <span class="metric-value">{component_stats['total_components']}</span>
</div>
<div class="metric">
    <span class="metric-label">Total Vertices in Components:</span>
    <span class="metric-value">{component_stats['total_vertices']}</span>
</div>
<div class="metric">
    <span class="metric-label">Avg Component Size:</span>
    <span class="metric-value">{component_stats['avg_size']:.2f}</span>
</div>
<div class="metric">
    <span class="metric-label">Singletons:</span>
    <span class="metric-value">{component_stats['singletons']}</span>
</div>
""")

    # Cycle convergence table
    if cycles:
        html_parts.append("<h2>🔁 Cycle Convergence</h2>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>Cycle</th><th>Edges</th><th>Components Accepted</th><th>Unpaired Vertices</th><th>Vertices Deleted</th><th>Status</th></tr>")
        
        for cycle in cycles:
            cycle_num = cycle['cycle']
            edges = cycle['edges_generated']
            accepted = cycle['components_accepted']
            unpaired = cycle.get('unpaired_vertices', '-')
            deleted = cycle['vertices_deleted']

            if edges == 0 or deleted == 0:
                status = '<span class="success">✓ Converged</span>'
            else:
                status = '<span class="info">→ Continue</span>'

            html_parts.append(f"<tr><td>{cycle_num}</td><td>{edges}</td><td>{accepted}</td><td>{unpaired}</td><td>{deleted}</td><td>{status}</td></tr>")
        
        html_parts.append("</table>")
        
        # Convergence plot
        if convergence_img:
            html_parts.append('<div class="plot">')
            html_parts.append(f'<img src="data:image/png;base64,{convergence_img}" alt="Cycle Convergence"/>')
            html_parts.append('</div>')
    
    # Component size distribution
    html_parts.append("<h2>📈 Component Size Distribution</h2>")
    total_c = component_stats['total_components']
    total_v = component_stats['total_vertices']
    avg_c_size = total_v / total_c if total_c > 0 else 0.0
    html_parts.append(f"<p>Component sizes range from {component_stats['min_size']} to {component_stats['max_size']} vertices.</p>")
    html_parts.append(f"<p>Median component size: <strong>{component_stats['median_size']:.1f}</strong> &nbsp;|&nbsp; "
                      f"Average component size: <strong>{avg_c_size:.2f}</strong> "
                      f"({total_v:,} vertices / {total_c:,} components)</p>")

    if size_dist_img:
        html_parts.append('<div class="plot">')
        html_parts.append(f'<img src="data:image/png;base64,{size_dist_img}" alt="Component Size Distribution"/>')
        html_parts.append('</div>')
    
    # Pep_ident conflict analysis
    if pep_ident_stats['total_components'] > 0:
        html_parts.append("<h2>🔬 Peptide Identification Analysis</h2>")
        html_parts.append(f"<p>Total components analyzed: <strong>{pep_ident_stats['total_components']}</strong></p>")
        html_parts.append(f"<p>Components with identifications: <strong>{pep_ident_stats['components_with_pep_idents']}</strong></p>")
        html_parts.append(f"<p>Components with conflicting identifications: <strong>{pep_ident_stats['components_with_conflicts']}</strong></p>")
        
        if pep_ident_stats['components_with_pep_idents'] > 0:
            conflict_pct = pep_ident_stats['conflict_percentage']
            if conflict_pct > 20:
                status_class = 'warning'
                status_msg = f"{conflict_pct:.1f}% of identified components have conflicts"
            elif conflict_pct > 0:
                status_class = 'info'
                status_msg = f"{conflict_pct:.1f}% of identified components have conflicts"
            else:
                status_class = 'success'
                status_msg = "No conflicts detected!"
            
            html_parts.append(f'<p class="{status_class}">Conflict rate: {status_msg}</p>')

            # Mass difference statistics + histogram
            mds = pep_ident_stats.get('mass_diff_stats')
            if mds:
                html_parts.append(
                    f'<p>Mean pairwise mass difference between conflicting sequences: '
                    f'<strong>{mds["mean"]:.3f} Da</strong> &nbsp;|&nbsp; '
                    f'Median: <strong>{mds["median"]:.3f} Da</strong> '
                    f'(across {mds["n_components"]:,} conflicted components)</p>'
                )
                mass_hist_img = plot_mass_diff_histogram(
                    pep_ident_stats.get('all_pairwise_diffs', [])
                )
                if mass_hist_img:
                    html_parts.append('<div class="plot">')
                    html_parts.append(
                        f'<img src="data:image/png;base64,{mass_hist_img}" '
                        f'alt="Pairwise Mass Difference Distribution"/>'
                    )
                    html_parts.append('</div>')

            # Show conflict examples if any
            if pep_ident_stats['conflict_examples']:
                html_parts.append("<h3>Example Conflicts</h3>")
                html_parts.append("<p><em>Components where vertices have different peptide identifications:</em></p>")
                html_parts.append("<table>")
                html_parts.append("<tr><th>Component ID</th><th>Size</th><th>Cycle</th><th>Conflicting Pep_idents (raw form · mass)</th></tr>")

                for example in pep_ident_stats['conflict_examples']:
                    comp_id = example['component_id']
                    size = example['num_vertices']
                    cycle = example.get('cycle_num', 'N/A')
                    seq_masses = example.get('pep_ident_masses', {})
                    raw_map   = example.get('pep_ident_raw', {})
                    shown = example['pep_idents'][:5]
                    parts = []
                    for seq in shown:
                        raw = raw_map.get(seq, seq)
                        mass_str = f' [{seq_masses[seq]:.3f}&nbsp;Da]' if seq in seq_masses else ''
                        # Show raw form in monospace; if different from clean, also show clean in grey
                        if raw != seq:
                            parts.append(f'<code>{raw}</code><span style="color:#888;font-size:0.85em"> ({seq})</span>{mass_str}')
                        else:
                            parts.append(f'<code>{seq}</code>{mass_str}')
                    pep_idents = '<br>'.join(parts)
                    if len(example['pep_idents']) > 5:
                        pep_idents += f'<br><em>+{len(example["pep_idents"]) - 5} more</em>'

                    html_parts.append(f"<tr><td>{comp_id}</td><td>{size}</td><td>{cycle}</td><td>{pep_idents}</td></tr>")

                html_parts.append("</table>")

            # ── Multi-component pep_ident redundancy ──────────────────────────
            mcs = pep_ident_stats.get('multi_component_stats')
            if mcs:
                html_parts.append("<h3>Pep_ident Redundancy Across Components</h3>")
                html_parts.append(
                    f'<p>Unique pep_idents with identifications: '
                    f'<strong>{mcs["total_unique_pep_idents"]:,}</strong> &nbsp;|&nbsp; '
                    f'Appearing in ≥2 components: '
                    f'<strong>{mcs["total_multi_component"]:,}</strong></p>'
                )
                multi_hist_img = plot_multi_component_histogram(mcs.get('distribution', {}))
                if multi_hist_img:
                    html_parts.append('<div class="plot">')
                    html_parts.append(
                        f'<img src="data:image/png;base64,{multi_hist_img}" '
                        f'alt="Pep_ident Redundancy Histogram"/>'
                    )
                    html_parts.append('</div>')

                top_redundant = mcs.get('top_redundant', [])
                if top_redundant:
                    html_parts.append(
                        f'<h4>Top {len(top_redundant)} Most Redundant Pep_idents</h4>'
                    )
                    html_parts.append(
                        '<table>'
                        '<tr><th>#</th><th>Pep_ident sequence</th><th>Components</th></tr>'
                    )
                    for rank, (pid, n_comp) in enumerate(top_redundant, 1):
                        html_parts.append(
                            f'<tr><td>{rank}</td><td><code>{pid}</code></td><td>{n_comp:,}</td></tr>'
                        )
                    html_parts.append('</table>')

            vcs = pep_ident_stats.get('vertex_conflict_stats')
            if vcs:
                html_parts.append("<h3>Vertex-Level Pep_ident Conflict Analysis</h3>")
                html_parts.append(
                    f'<p>Total vertices: <strong>{vcs["total_vertices"]:,}</strong> &nbsp;|&nbsp; '
                    f'Vertices with pep_ident: <strong>{vcs["total_with_pep_ident"]:,}</strong> '
                    f'({vcs["pct_with_pep_ident"]:.1f}% of all vertices)</p>'
                )
                html_parts.append(
                    '<table>'
                    '<tr><th>Estimate</th><th>Vertices</th><th>% of identified vertices</th><th>Description</th></tr>'
                    f'<tr><td><strong>Minimal</strong></td>'
                    f'<td>{vcs["minimal"]:,}</td>'
                    f'<td>{vcs["minimal_pct"]:.1f}%</td>'
                    f'<td>Vertices whose pep_ident is not the dominant pep_ident of their component</td></tr>'
                    f'<tr><td><strong>Moderate</strong></td>'
                    f'<td>{vcs["moderate"]:,}</td>'
                    f'<td>{vcs["moderate_pct"]:.1f}%</td>'
                    f'<td>All vertices with a pep_ident in any conflicted component</td></tr>'
                    f'<tr><td><strong>Maximal</strong></td>'
                    f'<td>{vcs["maximal"]:,}</td>'
                    f'<td>{vcs["maximal_pct"]:.1f}%</td>'
                    f'<td>All vertices (including unidentified) in any conflicted component</td></tr>'
                    '</table>'
                )
        else:
            html_parts.append("<p><em>No components with peptide identifications found.</em></p>")

    # Feature mass analysis
    if feature_mass_stats:
        html_parts.append(generate_feature_mass_html_section(feature_mass_stats))

    # File linkage analysis
    if file_linkage['total_files'] > 0:
        html_parts.append("<h2>🔗 File Linkage Analysis</h2>")
        html_parts.append(f"<p>Total files: <strong>{file_linkage['total_files']}</strong></p>")
        html_parts.append(f"<p>Average linkage score: <strong>{file_linkage['avg_linkage_score']:.4f}</strong></p>")
        html_parts.append("<p><em>Linkage score = proportion of components containing each file</em></p>")
        
        # Show top 20 files by linkage score
        sorted_files = sorted(file_linkage['file_linkage_scores'].items(), key=lambda x: x[1], reverse=True)
        html_parts.append("<h3>Top Files by Linkage Score</h3>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>File</th><th>Linkage Score</th><th>Components</th></tr>")
        
        for filename, score in sorted_files[:20]:
            comp_count = int(score * component_stats['total_components'])
            html_parts.append(f"<tr><td>{filename}</td><td>{score:.4f}</td><td>{comp_count}</td></tr>")
        
        html_parts.append("</table>")
    
    # Parameters (if provided)
    if parameters:
        html_parts.append("<h2>⚙️ Parameters</h2>")
        html_parts.append("<table>")
        html_parts.append("<tr><th>Parameter</th><th>Value</th></tr>")
        for key, value in sorted(parameters.items()):
            html_parts.append(f"<tr><td>{key}</td><td>{value}</td></tr>")
        html_parts.append("</table>")
    
    # Footer
    html_parts.append("""
</div>
</body>
</html>
""")
    
    # Write HTML file
    with open(output_html, 'w') as f:
        f.write(''.join(html_parts))
    
    print(f"✓ Generated HTML report: {output_html}")


def generate_text_summary(
    component_stats: Dict,
    cycles: List[Dict],
    cycle_metadata: Dict,
    file_linkage: Dict,
    pep_ident_stats: Dict,
    output_summary: str,
    parameters: Dict = None,
    feature_mass_stats: Dict = None,
):
    """Generate plain text summary."""

    initial_vertices = cycle_metadata.get('initial_vertices')
    unpaired_vertices_final = component_stats['singletons']
    vertices_after_cutoff = (initial_vertices - unpaired_vertices_final) if initial_vertices is not None else None

    lines = []
    lines.append("=" * 70)
    lines.append("CLUSTERLESS COMPONENT BUILDING SUMMARY")
    lines.append("=" * 70)
    lines.append("")

    # Run parameters (shown first so they're easy to find in the file)
    if parameters:
        lines.append("RUN PARAMETERS")
        lines.append("-" * 70)
        for key, value in parameters.items():
            lines.append(f"  {key:<28} {value}")
        lines.append("")

    # Overview
    lines.append("OVERVIEW")
    lines.append("-" * 70)
    lines.append(f"Total cycles: {len(cycles)}")
    if initial_vertices is not None:
        lines.append(f"Vertices before cutoff (input): {initial_vertices:,}")
        lines.append(f"Unpaired vertices (final):      {unpaired_vertices_final:,}")
        if vertices_after_cutoff is not None:
            pct = 100 * vertices_after_cutoff / initial_vertices if initial_vertices else 0
            lines.append(f"Vertices after cutoff:          {vertices_after_cutoff:,}  ({pct:.1f}% of input)")
    lines.append(f"Final components: {component_stats['total_components']}")
    lines.append(f"Total vertices in components: {component_stats['total_vertices']}")
    lines.append(f"Average component size: {component_stats['avg_size']:.2f}")
    lines.append(f"Median component size: {component_stats['median_size']:.1f}")
    lines.append(f"Size range: {component_stats['min_size']} - {component_stats['max_size']}")
    lines.append(f"Singleton components: {component_stats['singletons']}")
    lines.append("")

    # Cycle convergence
    if cycles:
        lines.append("CYCLE CONVERGENCE")
        lines.append("-" * 70)
        lines.append(f"{'Cycle':<8} {'Edges':<12} {'Accepted':<12} {'Unpaired':<12} {'Deleted':<12} {'Status':<15}")
        lines.append("-" * 70)
        
        for cycle in cycles:
            cycle_num = cycle['cycle']
            edges = cycle['edges_generated']
            accepted = cycle['components_accepted']
            unpaired = cycle.get('unpaired_vertices', '-')
            deleted = cycle['vertices_deleted']

            if edges == 0 or deleted == 0:
                status = "CONVERGED"
            else:
                status = "CONTINUE"

            lines.append(f"{cycle_num:<8} {edges:<12} {accepted:<12} {unpaired:<12} {deleted:<12} {status:<15}")

        lines.append("")
    
    # Pep_ident analysis
    if pep_ident_stats['total_components'] > 0:
        lines.append("PEPTIDE IDENTIFICATION ANALYSIS")
        lines.append("-" * 70)
        lines.append(f"Total components: {pep_ident_stats['total_components']}")
        lines.append(f"Components with identifications: {pep_ident_stats['components_with_pep_idents']}")
        lines.append(f"Components with conflicts: {pep_ident_stats['components_with_conflicts']}")
        
        if pep_ident_stats['components_with_pep_idents'] > 0:
            lines.append(f"Conflict rate: {pep_ident_stats['conflict_percentage']:.1f}%")

            mds = pep_ident_stats.get('mass_diff_stats')
            if mds:
                lines.append(
                    f"Mass diff (mean):   {mds['mean']:.3f} Da  "
                    f"(median: {mds['median']:.3f} Da, N={mds['n_components']} conflicted components)"
                )

            mcs = pep_ident_stats.get('multi_component_stats')
            if mcs:
                dist = mcs.get('distribution', {})
                dist_str = '  '.join(f"{n} comps: {c:,}" for n, c in sorted(dist.items()))
                lines.append(
                    f"Pep_ident redundancy: {mcs['total_multi_component']:,} of "
                    f"{mcs['total_unique_pep_idents']:,} unique pep_idents appear in ≥2 components"
                )
                if dist_str:
                    lines.append(f"  Distribution: {dist_str}")
                top_redundant = mcs.get('top_redundant', [])
                if top_redundant:
                    lines.append(f"  Top {len(top_redundant)} most redundant pep_idents:")
                    for rank, (pid, n_comp) in enumerate(top_redundant, 1):
                        lines.append(f"    {rank:2d}. {pid}  ({n_comp} components)")

            vcs = pep_ident_stats.get('vertex_conflict_stats')
            if vcs:
                lines.append("")
                lines.append("VERTEX-LEVEL PEP_IDENT CONFLICT ANALYSIS")
                lines.append(f"  Total vertices:               {vcs['total_vertices']:,}")
                lines.append(
                    f"  Vertices with pep_ident:      {vcs['total_with_pep_ident']:,}"
                    f"  ({vcs['pct_with_pep_ident']:.1f}% of all)"
                )
                lines.append(f"  Minimal estimate:             {vcs['minimal']:,}"
                             f"  ({vcs['minimal_pct']:.1f}% of identified)"
                             f"  [non-dominant pep_ident in conflicted component]")
                lines.append(f"  Moderate estimate:            {vcs['moderate']:,}"
                             f"  ({vcs['moderate_pct']:.1f}% of identified)"
                             f"  [any pep_ident vertex in conflicted component]")
                lines.append(f"  Maximal estimate:             {vcs['maximal']:,}"
                             f"  ({vcs['maximal_pct']:.1f}% of identified)"
                             f"  [all vertices in conflicted component]")

            if pep_ident_stats['conflict_examples']:
                lines.append("")
                lines.append("Example conflicts (first 5):")
                for i, example in enumerate(pep_ident_stats['conflict_examples'][:5], 1):
                    comp_id = example['component_id']
                    size = example['num_vertices']
                    seq_masses = example.get('pep_ident_masses', {})
                    raw_map   = example.get('pep_ident_raw', {})
                    shown = example['pep_idents'][:3]
                    parts = []
                    for seq in shown:
                        raw = raw_map.get(seq, seq)
                        mass_str = f' [{seq_masses[seq]:.3f} Da]' if seq in seq_masses else ''
                        if raw != seq:
                            parts.append(f'{raw} (clean: {seq}){mass_str}')
                        else:
                            parts.append(f'{seq}{mass_str}')
                    pep_idents = ', '.join(parts)
                    if len(example['pep_idents']) > 3:
                        pep_idents += '...'
                    lines.append(f"  {i}. {comp_id} (size {size}): {pep_idents}")
        
        lines.append("")
    
    # Feature mass analysis
    if feature_mass_stats:
        lines.append("FEATURE MASS ANALYSIS")
        lines.append("-" * 70)
        lines.append(f"Components with mass data: {feature_mass_stats['n_components_with_masses']:,}")
        ps = feature_mass_stats.get('pairwise_stats')
        if ps:
            lines.append(f"All pairwise diffs  — mean: {ps['mean']:.3f} Da  "
                         f"median: {ps['median']:.3f} Da  N: {ps['n']:,}")
        ms = feature_mass_stats.get('mean_diff_stats')
        if ms:
            lines.append(f"Component mean diffs — mean: {ms['mean']:.3f} Da  "
                         f"median: {ms['median']:.3f} Da  N: {ms['n']:,}")
        ns = feature_mass_stats.get('nn_stats')
        if ns:
            lines.append(f"Nearest-neighbour   — mean: {ns['mean']:.3f} Da  "
                         f"median: {ns['median']:.3f} Da  N: {ns['n']:,}")
        lines.append("")

    # File linkage
    if file_linkage['total_files'] > 0:
        lines.append("FILE LINKAGE")
        lines.append("-" * 70)
        lines.append(f"Total files: {file_linkage['total_files']}")
        lines.append(f"Average linkage score: {file_linkage['avg_linkage_score']:.4f}")
        lines.append("")
    
    # Component size distribution (summary)
    lines.append("COMPONENT SIZE DISTRIBUTION")
    lines.append("-" * 70)
    total_c = component_stats['total_components']
    total_v = component_stats['total_vertices']
    avg_c_size = total_v / total_c if total_c > 0 else 0.0
    lines.append(f"Average component size: {avg_c_size:.2f}  ({total_v:,} vertices / {total_c:,} components)")
    lines.append("")
    size_dist = component_stats.get('size_distribution', {})
    for size in sorted(size_dist.keys())[:20]:  # Show first 20 sizes
        count = size_dist[size]
        lines.append(f"  Size {size:3d}: {count:6d} components")

    if len(size_dist) > 20:
        lines.append(f"  ... and {len(size_dist) - 20} more sizes")
    
    lines.append("")
    lines.append("=" * 70)
    
    # Write summary file
    with open(output_summary, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"✓ Generated text summary: {output_summary}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate summary report for clusterless component building workflow"
    )
    parser.add_argument("--components_json", required=True,
                       help="Input: final_components_all_cycles.json from orchestrator")
    parser.add_argument("--cycle_history_json", default=None,
                       help="Input: cycle_history.json from orchestrator (optional)")
    parser.add_argument("--output_summary", required=True,
                       help="Output: plain text summary file")
    parser.add_argument("--output_html", default=None,
                       help="Output: HTML report file (optional)")
    parser.add_argument("--parameters_json", default=None,
                       help="Input: JSON file with parameters to include in report (optional)")
    parser.add_argument("--nterm_mod", default="none",
                       choices=["none", "acetyl", "biotin", "fluorescein"],
                       help="N-terminal modification for monoisotopic mass calculation (default: none)")
    parser.add_argument("--cys_mod", default="none",
                       choices=["none", "oxidized", "carbamidomethyl"],
                       help="Cysteine modification for monoisotopic mass calculation (default: none)")
    parser.add_argument("--cterm_mod", default="none",
                       choices=["none", "amide"],
                       help="C-terminal modification: none (free -COOH) | amide (-CONH2, -0.984 Da) (default: none)")
    parser.add_argument("--tsv_files", nargs='*', default=None,
                       help="TSV files from map_mzml_features (for raw pep_ident with inline mods)")

    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("Generating Clusterless Component Building Report")
    print("="*70 + "\n")
    
    # Load data
    components_data = load_final_components(args.components_json)
    components = list(components_data.get('components', {}).values())
    
    cycles = []
    cycle_metadata = {}
    if args.cycle_history_json:
        cycles, cycle_metadata = load_cycle_history(args.cycle_history_json)
    
    parameters = None
    if args.parameters_json and os.path.exists(args.parameters_json):
        with open(args.parameters_json, 'r') as f:
            parameters = json.load(f)
        print(f"✓ Loaded parameters from: {os.path.basename(args.parameters_json)}")
    
    # Analyze components
    print("\nAnalyzing components...")
    component_stats = analyze_component_sizes(components)
    print(f"  Total components: {component_stats['total_components']}")
    print(f"  Total vertices: {component_stats['total_vertices']}")
    print(f"  Average size: {component_stats['avg_size']:.2f}")
    
    # Analyze file linkage
    print("\nAnalyzing file linkage...")
    file_linkage = analyze_file_linkage(components)
    print(f"  Total files: {file_linkage['total_files']}")
    print(f"  Average linkage score: {file_linkage['avg_linkage_score']:.4f}")
    
    # Analyze pep_idents
    print("\nAnalyzing peptide identifications...")
    print(f"  N-terminal modification: {args.nterm_mod}")
    print(f"  Cysteine modification:   {args.cys_mod}")
    print(f"  C-terminal modification: {args.cterm_mod}")
    raw_lookup = None
    if args.tsv_files:
        print(f"  Loading raw pep_idents from {len(args.tsv_files)} TSV file(s)...")
        raw_lookup = build_raw_pep_ident_lookup(args.tsv_files)
        print(f"  Found {len(raw_lookup):,} features with raw pep_idents")
    pep_ident_stats = analyze_pep_idents(
        components, args.nterm_mod, args.cys_mod, args.cterm_mod, raw_lookup
    )
    print(f"  Components with identifications: {pep_ident_stats['components_with_pep_idents']}")
    print(f"  Components with conflicts: {pep_ident_stats['components_with_conflicts']}")
    if pep_ident_stats['components_with_pep_idents'] > 0:
        print(f"  Conflict rate: {pep_ident_stats['conflict_percentage']:.1f}%")
    
    # Analyze feature masses
    print("\nAnalyzing feature masses...")
    feature_mass_stats = analyze_feature_masses(components)
    ps = feature_mass_stats.get('pairwise_stats')
    ns = feature_mass_stats.get('nn_stats')
    if ps:
        print(f"  Pairwise diffs: mean={ps['mean']:.3f} Da  median={ps['median']:.3f} Da  N={ps['n']:,}")
    if ns:
        print(f"  Nearest-neighbour: mean={ns['mean']:.3f} Da  median={ns['median']:.3f} Da")

    # Generate outputs
    print("\nGenerating reports...")
    generate_text_summary(
        component_stats, cycles, cycle_metadata, file_linkage, pep_ident_stats,
        args.output_summary, parameters, feature_mass_stats,
    )

    if args.output_html:
        generate_html_report(
            components_data,
            cycles,
            cycle_metadata,
            component_stats,
            file_linkage,
            pep_ident_stats,
            args.output_html,
            parameters,
            feature_mass_stats,
        )
    
    print("\n" + "="*70)
    print("✅ Report generation complete!")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
