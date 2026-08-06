#!/usr/bin/env python3
"""
Spike-in intensity QC report.

WHAT IT DOES
------------
Within each experimental condition, compares the (normalized) intensity
distribution of every spike-in protein across that condition's measurement
files, and summarizes how the spike-in level relates between conditions.

INPUT
-----
The normalized quantification table produced by merge_normalized_intensities.py.
It is schema-identical to the reduced quantification table:
  * ``<FILE>_____intensity``      -> log2-normalized intensity for that file
  * ``<FILE>_____l_prot_ident``   -> per-file protein identifications (a list)
The condition is encoded in the column-name prefix (``C1_FLI04006`` -> ``C1``),
and can alternatively be taken from the group-columns JSON.

FEATURE SELECTION (important)
-----------------------------
A feature counts towards a spike-in only if its identification is EXACTLY and
UNIQUELY that spike-in: the union of all protein accessions across every file
for that feature must be a single accession equal to the spike-in. Features that
also map to any other protein (e.g. a mouse background protein) are excluded.

OUTPUT
------
A self-contained HTML report with:
  * one section per condition, each holding a grouped boxplot (x = spike-ins,
    one box per file) of the log2-normalized intensities; and
  * a relative-mean table (rows = spike-ins, cols = conditions) where, per
    spike-in row, the lowest condition mean is anchored to 0.1 and the others
    scaled proportionally. Means are computed in LINEAR space (2**log2) so a
    condition 20x more intense than the lowest becomes 2.0.

CODE LAYOUT
-----------
The code is split into clearly separated phases so the analysis is easy to read:
  SECTION 1  Configuration and constants
  SECTION 2  Argument parsing
  SECTION 3  Low-level parsers (FASTA, intensity cells, protein-ident cells)
  SECTION 4  Data loading + feature selection  -> the raw `measurements`
  SECTION 5  Calculations (pure numbers, NO plotting/HTML)
  SECTION 6  Plot building (plotly figures from computed data)
  SECTION 7  HTML rendering (assemble the document)
  SECTION 8  main() orchestration
"""

import argparse
import csv
import json
import sys
from ast import literal_eval
from collections import namedtuple

import numpy as np
import plotly.graph_objects as go

# One measured feature in one file: its (log2) intensity plus the peptide
# sequence and charge that identified it - used to enrich outlier tooltips.
Point = namedtuple("Point", ["intensity", "pep", "charge"])


# ===========================================================================
# SECTION 1 - CONFIGURATION AND CONSTANTS
# ===========================================================================

# Column-name suffixes in the normalized table.
INTENSITY_SUFFIX = "_____intensity"
PROT_IDENT_SUFFIX = "_____l_prot_ident"
PEP_IDENT_SUFFIX = "_____l_pep_ident"
CHARGE_SUFFIX = "_____charge"

# Cell contents that mean "no value".
NONE_VALUES = {"nan", "inf", "-inf", "#n/a", "n/a", "na", ""}

# The lowest condition mean in each spike-in row is anchored to this value.
RELATIVE_ANCHOR = 0.1

# Qualitative palette for the per-file boxes (cycled if a condition has more
# files than colours).
FILE_COLORS = [
    "#4f7832", "#365575", "#076ddb", "#59cc06", "#cc6b06",
    "#8b2fb0", "#159e8c", "#b0402f", "#7d7d1f", "#2f6db0",
]

csv.field_size_limit(sys.maxsize)


# ===========================================================================
# SECTION 2 - ARGUMENT PARSING
# ===========================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Spike-in intensity QC report.")
    parser.add_argument(
        "--input_tsv", "-in", required=True,
        help="Normalized quantification table (clusters_quantification_normalized.tsv).")
    parser.add_argument(
        "--spikeins_fasta", "-fasta", required=True,
        help="FASTA of the spike-in proteins (files/spikeins.txt).")
    parser.add_argument(
        "--group_columns", "-g", default=None,
        help="Optional group-columns JSON (condition -> [intensity column names]). "
             "If omitted, the condition is the column-name prefix before the first "
             "underscore (e.g. 'C1').")
    parser.add_argument(
        "--output_html", "-out", required=True,
        help="Output self-contained HTML report.")
    parser.add_argument(
        "--output_table", "-t", default=None,
        help="Optional output TSV of the relative-mean table.")
    parser.add_argument(
        "--log2_input", "-lip", type=int, default=1,
        help="Are intensities log2 transformed? 1 (default) = yes; the table "
             "un-logs (2**x) before averaging. 0 = treat values as linear.")
    return parser.parse_args()


# ===========================================================================
# SECTION 3 - LOW-LEVEL PARSERS
# ===========================================================================

def parse_spikeins(fasta_path):
    """Return an ordered dict {accession: label} from a spike-in FASTA.

    Handles the well-formed ``>sp|ACC|ENTRY Description OS=...`` headers as well
    as the odd ``>sp|ACC|ENTRY|Description`` variant (myoglobin) by taking the
    accession as the 2nd pipe field and the label as the first token of the 3rd.
    """
    spikeins = {}
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split("|")
            if len(parts) < 3:
                continue
            acc = parts[1].strip()
            entry = parts[2].strip().split()[0] if parts[2].strip() else acc
            spikeins[acc] = entry
    return spikeins


def parse_intensity(cell):
    """Parse an intensity cell to float, or None if empty/unparsable."""
    if cell is None:
        return None
    cell = cell.strip()
    if cell.lower() in NONE_VALUES:
        return None
    try:
        val = float(cell)
    except ValueError:
        return None
    if np.isnan(val) or np.isinf(val):
        return None
    return val


def parse_prot_set(cell):
    """Parse an ``l_prot_ident`` cell into a set of protein accessions.

    A cell looks like ``['P02754']`` or ``['O55042,Q91ZZ3,P37840']`` - a protein
    group is comma-joined inside a single string element, so every element is
    additionally split on commas.
    """
    if not cell or cell.strip().lower() in NONE_VALUES:
        return set()
    try:
        elements = literal_eval(cell)
    except (ValueError, SyntaxError):
        return set()
    accs = set()
    for element in elements:
        for acc in str(element).split(","):
            acc = acc.strip()
            if acc:
                accs.add(acc)
    return accs


def parse_pep_ident(cell):
    """Parse an ``l_pep_ident`` cell (e.g. ``['LLTSLR']``) into a display string.

    Multiple peptides (across list elements or comma-joined) are collapsed to a
    unique, order-preserving ``;``-separated string.
    """
    if not cell or cell.strip().lower() in NONE_VALUES:
        return ""
    try:
        elements = literal_eval(cell)
    except (ValueError, SyntaxError):
        return ""
    peps = []
    for element in elements:
        for pep in str(element).split(","):
            pep = pep.strip()
            if pep and pep not in peps:
                peps.append(pep)
    return ";".join(peps)


def parse_charge(cell):
    """Parse a ``charge`` cell (a plain int like ``2``) into a string, or ''."""
    if cell is None:
        return ""
    cell = cell.strip()
    return "" if cell.lower() in NONE_VALUES else cell


# ===========================================================================
# SECTION 4 - DATA LOADING + FEATURE SELECTION
# ===========================================================================
#
# Produces the single raw data structure everything downstream consumes:
#
#   Measurements:
#     .file_ids        sorted list of file identifiers (e.g. 'C1_FLI04006')
#     .conditions      sorted list of conditions (e.g. ['C1', ... 'C5'])
#     .files_by_cond   {condition: [file_id, ...]}
#     .values          {spike_acc: {file_id: [Point(intensity, pep, charge), ...]}}
#     .feature_counts  {spike_acc: number of qualifying unique features}
#
# No calculations or presentation happen here - only reading and selection.
# ===========================================================================

class Measurements:
    def __init__(self, file_ids, conditions, files_by_cond, values, feature_counts):
        self.file_ids = file_ids
        self.conditions = conditions
        self.files_by_cond = files_by_cond
        self.values = values
        self.feature_counts = feature_counts


def _condition_from_prefix(file_id):
    """Condition = the column-name prefix before the first underscore (e.g. C1)."""
    return file_id.split("_", 1)[0]


def _build_condition_map(intensity_columns, group_columns_path):
    """Map file_id -> condition, preferring the group-columns JSON when given."""
    file_to_cond = {}
    if group_columns_path:
        with open(group_columns_path) as fh:
            group_columns = json.load(fh)
        col_to_cond = {col: cond
                       for cond, cols in group_columns.items()
                       for col in cols}
        for col, file_id in intensity_columns.items():
            if col in col_to_cond:
                file_to_cond[file_id] = col_to_cond[col]
    # Fall back to the prefix for anything the JSON did not cover.
    for file_id in intensity_columns.values():
        file_to_cond.setdefault(file_id, _condition_from_prefix(file_id))
    return file_to_cond


def load_measurements(input_tsv, spikeins, group_columns_path):
    """Read the normalized table and collect intensities of unique spike-in features."""
    spike_set = set(spikeins)

    with open(input_tsv, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)

        # Map each relevant column to its file_id and its row index.
        intensity_columns = {c: c[:-len(INTENSITY_SUFFIX)]
                             for c in header if c.endswith(INTENSITY_SUFFIX)}
        protident_columns = {c: c[:-len(PROT_IDENT_SUFFIX)]
                             for c in header if c.endswith(PROT_IDENT_SUFFIX)}
        pepident_columns = {c: c[:-len(PEP_IDENT_SUFFIX)]
                            for c in header if c.endswith(PEP_IDENT_SUFFIX)}
        charge_columns = {c: c[:-len(CHARGE_SUFFIX)]
                          for c in header if c.endswith(CHARGE_SUFFIX)}
        intensity_idx = {intensity_columns[c]: header.index(c) for c in intensity_columns}
        protident_idx = {protident_columns[c]: header.index(c) for c in protident_columns}
        pepident_idx = {pepident_columns[c]: header.index(c) for c in pepident_columns}
        charge_idx = {charge_columns[c]: header.index(c) for c in charge_columns}

        file_ids = sorted(intensity_idx)
        file_to_cond = _build_condition_map(intensity_columns, group_columns_path)
        conditions = sorted({file_to_cond[f] for f in file_ids})
        files_by_cond = {cond: [f for f in file_ids if file_to_cond[f] == cond]
                         for cond in conditions}

        values = {acc: {f: [] for f in file_ids} for acc in spikeins}
        feature_counts = {acc: 0 for acc in spikeins}

        n_rows = 0
        for row in reader:
            if not row:
                continue
            n_rows += 1

            # Union of protein identifications across every file for this feature.
            prot_set = set()
            for idx in protident_idx.values():
                if idx < len(row):
                    prot_set |= parse_prot_set(row[idx])

            # SELECTION: keep only features that are exactly one spike-in.
            if len(prot_set) != 1:
                continue
            acc = next(iter(prot_set))
            if acc not in spike_set:
                continue

            feature_counts[acc] += 1

            # Feature-level peptide/charge, pooled across ALL files. Used as a
            # fallback for files that measured the feature (intensity present)
            # but carry no local identification (e.g. match-between-runs), which
            # would otherwise show "(no peptide)".
            feat_peps = []
            for p_idx in pepident_idx.values():
                if p_idx < len(row):
                    for pep in parse_pep_ident(row[p_idx]).split(";"):
                        if pep and pep not in feat_peps:
                            feat_peps.append(pep)
            feature_pep = ";".join(feat_peps)
            feat_charges = []
            for c_idx in charge_idx.values():
                if c_idx < len(row):
                    ch = parse_charge(row[c_idx])
                    if ch and ch not in feat_charges:
                        feat_charges.append(ch)
            feature_charge = "/".join(feat_charges)

            for file_id, idx in intensity_idx.items():
                if idx >= len(row):
                    continue
                val = parse_intensity(row[idx])
                if val is None:
                    continue
                # Prefer this file's own identification; fall back to the
                # feature's identification from other files when absent.
                p_idx = pepident_idx.get(file_id)
                c_idx = charge_idx.get(file_id)
                local_pep = parse_pep_ident(row[p_idx]) if p_idx is not None and p_idx < len(row) else ""
                local_charge = parse_charge(row[c_idx]) if c_idx is not None and c_idx < len(row) else ""
                if local_pep:
                    pep = local_pep
                elif feature_pep:
                    pep = f"{feature_pep} (filled from other run)"
                else:
                    pep = "(no peptide)"
                charge = local_charge or feature_charge
                values[acc][file_id].append(Point(val, pep, charge))

    print(f"Scanned {n_rows:,} feature row(s).")
    return Measurements(file_ids, conditions, files_by_cond, values, feature_counts)


# ===========================================================================
# SECTION 5 - CALCULATIONS  (pure numbers - no plotting, no HTML)
# ===========================================================================

def linear_mean(values, log2_input):
    """Arithmetic mean in LINEAR space (un-logging first if input is log2).

    We average in linear space so that ratios between conditions are meaningful
    (a 20x more intense condition really is 20x), which the relative table needs.
    """
    if not values:
        return None
    arr = np.asarray(values, dtype=float)
    if log2_input:
        arr = np.power(2.0, arr)
    return float(np.mean(arr))


def compute_relative_table(m: Measurements, log2_input):
    """Compute the relative-mean table (spike-in x condition).

    Returns two dicts keyed [spike_acc][condition]:
      abs_means -> the linear-space mean of ALL pooled measurements of that
                   (spike-in, condition), or None if there were none.
      rel_means -> abs_means rescaled per spike-in row so the lowest condition
                   mean equals RELATIVE_ANCHOR (0.1) and the others scale up
                   proportionally. None where abs mean is missing.
    """
    abs_means = {}
    rel_means = {}
    for acc in m.values:
        # 1) Pool every measurement of this spike-in within each condition and
        #    take the linear-space mean.
        row_abs = {}
        for cond in m.conditions:
            pooled = []
            for file_id in m.files_by_cond[cond]:
                pooled.extend(pt.intensity for pt in m.values[acc][file_id])
            row_abs[cond] = linear_mean(pooled, log2_input)

        # 2) Anchor the lowest (positive) condition mean to 0.1 and scale up.
        present = [v for v in row_abs.values() if v is not None and v > 0]
        row_min = min(present) if present else None
        row_rel = {
            cond: (None if (v is None or row_min is None)
                   else RELATIVE_ANCHOR * v / row_min)
            for cond, v in row_abs.items()
        }

        abs_means[acc] = row_abs
        rel_means[acc] = row_rel
    return abs_means, rel_means


# ===========================================================================
# SECTION 6 - PLOT BUILDING  (plotly figures from computed data)
# ===========================================================================

def build_condition_boxplot(m: Measurements, spikeins, condition, y_title):
    """One grouped boxplot for a condition: x = spike-ins, one box per file."""
    spike_labels = [spikeins[acc] for acc in spikeins]
    fig = go.Figure()
    for i, file_id in enumerate(m.files_by_cond[condition]):
        # Flatten (spike-in, points) into parallel x / y arrays; plotly groups
        # boxes by identical x labels, and the trace name (file) sets the group.
        # `hovertext` is aligned with y so the drawn outlier points carry the
        # peptide sequence and charge (non-outliers are not drawn, so they never
        # show a tooltip). Using hovertext (not hovertemplate) keeps the default
        # box-statistics hover intact.
        xs, ys, htext = [], [], []
        for acc, label in zip(spikeins, spike_labels):
            for pt in m.values[acc][file_id]:
                xs.append(label)
                ys.append(pt.intensity)
                htext.append(
                    f"{pt.pep or '(no peptide)'} &middot; charge {pt.charge or '?'}"
                    f"<br>log2 intensity: {pt.intensity:.2f}"
                )
        fig.add_trace(go.Box(
            x=xs, y=ys, name=file_id, hovertext=htext,
            marker_color=FILE_COLORS[i % len(FILE_COLORS)],
            boxpoints="outliers", line_width=1,
        ))
    fig.update_layout(
        boxmode="group",
        xaxis_title="Spike-in protein",
        yaxis_title=y_title,
        template="plotly_white",
        height=460,
        margin=dict(l=60, r=20, t=20, b=120),
        legend_title="File",
    )
    fig.update_xaxes(categoryorder="array", categoryarray=spike_labels, tickangle=-40)
    return fig


# ===========================================================================
# SECTION 7 - HTML RENDERING  (assemble the document; no calculations here)
# ===========================================================================

def _fmt(value, digits=2):
    return "" if value is None else f"{value:.{digits}f}"


def render_relative_table(spikeins, feature_counts, conditions, rel_means, abs_means):
    """Render the relative-mean table as an HTML <table>."""
    head = "".join(f"<th>{c}</th>" for c in conditions)
    body = []
    for acc in spikeins:
        cells = []
        for cond in conditions:
            abs_m = abs_means[acc][cond]
            title = "" if abs_m is None else f' title="linear mean={abs_m:.3g}"'
            cells.append(f"<td{title}>{_fmt(rel_means[acc][cond])}</td>")
        body.append(
            f"<tr><th class='rowh'>{spikeins[acc]} "
            f"<span class='acc'>({acc}, n={feature_counts[acc]})</span></th>"
            + "".join(cells) + "</tr>"
        )
    return ("<table class='rel'><thead><tr><th>Spike-in</th>" + head
            + "</tr></thead><tbody>" + "".join(body) + "</tbody></table>")


def render_count_table(m: Measurements, spikeins, condition):
    """Number of features in each box (spike-in x file) for one condition.

    Rows = spike-ins, columns = the condition's files, cell = features in that
    box, plus a per-spike-in total (sum across files = total measurements; a
    feature seen in several files is counted once per file).
    """
    files = m.files_by_cond[condition]
    # Drop the shared condition prefix from the column headers for brevity.
    short = [f.split("_", 1)[1] if "_" in f else f for f in files]
    head = "".join(f"<th>{s}</th>" for s in short) + "<th>total</th>"
    body = []
    for acc in spikeins:
        counts = [len(m.values[acc][f]) for f in files]
        cells = "".join(f"<td>{c}</td>" for c in counts)
        body.append(f"<tr><th class='rowh'>{spikeins[acc]}</th>{cells}"
                    f"<td><strong>{sum(counts)}</strong></td></tr>")
    return ("<table class='rel'><thead><tr><th>Spike-in</th>" + head
            + "</tr></thead><tbody>" + "".join(body) + "</tbody></table>")


def render_report(input_tsv, m: Measurements, spikeins, figures,
                  rel_means, abs_means, log2_input):
    """Assemble the full HTML document from pre-built figures and the table."""
    y_title = "Normalized intensity (log2)" if log2_input else "Normalized intensity"
    parts = [_HTML_HEAD, "<h1>Spike-in intensity QC</h1>",
             "<p class='sub'>Per-condition intensity distributions of spike-in "
             "proteins. Only features whose identification is <em>exactly and "
             "uniquely</em> one spike-in protein are included. "
             f"Source: <code>{input_tsv}</code>.</p>"]

    # One section per condition, with its pre-built figure.
    for idx, (cond, fig) in enumerate(zip(m.conditions, figures)):
        cond_files = m.files_by_cond[cond]
        parts.append(f"<h2>Condition {cond}</h2>")
        parts.append(f"<p class='sub'>{len(cond_files)} file(s): "
                     + ", ".join(cond_files) + f" &mdash; y-axis: {y_title}.</p>")
        # Embed plotly.js only once (with the first figure) to stay self-contained.
        parts.append(fig.to_html(full_html=False,
                                 include_plotlyjs=(idx == 0)))
        # Per-box feature counts underneath the plot.
        parts.append("<p class='sub'>Features per box (columns = files, "
                     "total = sum across files):</p>")
        parts.append(render_count_table(m, spikeins, cond))

    # Relative-mean table.
    parts.append("<h2>Relative mean intensity (across conditions)</h2>")
    parts.append(
        "<p class='sub'>For each spike-in the lowest condition mean is set to "
        f"<strong>{RELATIVE_ANCHOR}</strong> and the others scaled proportionally "
        "(a condition 20&times; the lowest becomes 2.0). Means are computed in "
        "linear space. Empty cells = no qualifying measurements.</p>")
    parts.append(render_relative_table(
        spikeins, m.feature_counts, m.conditions, rel_means, abs_means))

    parts.append("</body></html>")
    return "\n".join(parts)


def write_relative_table_tsv(path, spikeins, conditions, rel_means):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["spike_in", "accession"] + list(conditions))
        for acc in spikeins:
            writer.writerow([spikeins[acc], acc]
                            + [_fmt(rel_means[acc][cond]) for cond in conditions])


_HTML_HEAD = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Spike-in intensity QC</title>
<style>
  body { font-family: -apple-system, Segoe UI, Roboto, sans-serif; margin: 24px;
         color: #1a1a1a; }
  h1 { font-size: 22px; }
  h2 { font-size: 18px; margin-top: 34px; border-bottom: 1px solid #ddd;
       padding-bottom: 4px; }
  p.sub { color: #555; font-size: 13px; max-width: 900px; }
  code { background: #f2f2f2; padding: 1px 4px; border-radius: 3px; }
  table.rel { border-collapse: collapse; margin-top: 12px; font-size: 13px; }
  table.rel th, table.rel td { border: 1px solid #ccc; padding: 5px 10px;
                               text-align: right; }
  table.rel thead th { background: #f2f2f2; text-align: center; }
  table.rel th.rowh { text-align: left; font-weight: 600; }
  table.rel .acc { color: #888; font-weight: 400; font-size: 11px; }
</style></head><body>
"""


# ===========================================================================
# SECTION 8 - MAIN ORCHESTRATION
# ===========================================================================

def main():
    args = parse_args()
    log2_input = bool(args.log2_input)

    # (3) Parse the spike-in list.
    spikeins = parse_spikeins(args.spikeins_fasta)
    if not spikeins:
        raise SystemExit(f"ERROR: no spike-ins parsed from {args.spikeins_fasta}")
    print(f"Loaded {len(spikeins)} spike-in(s): {', '.join(spikeins.values())}")

    # (4) Load data + select unique spike-in features.
    m = load_measurements(args.input_tsv, spikeins, args.group_columns)
    print(f"Found {len(m.file_ids)} file(s) across conditions: {', '.join(m.conditions)}")
    for acc in spikeins:
        print(f"    {spikeins[acc]} ({acc}): {m.feature_counts[acc]} unique feature(s)")

    # (5) Calculations.
    abs_means, rel_means = compute_relative_table(m, log2_input)

    # (6) Plot building.
    y_title = "Normalized intensity (log2)" if log2_input else "Normalized intensity"
    figures = [build_condition_boxplot(m, spikeins, cond, y_title)
               for cond in m.conditions]

    # (7) HTML rendering + optional TSV.
    html = render_report(args.input_tsv, m, spikeins, figures,
                         rel_means, abs_means, log2_input)
    with open(args.output_html, "w") as fh:
        fh.write(html)
    print(f"Wrote report: {args.output_html}")

    if args.output_table:
        write_relative_table_tsv(args.output_table, spikeins, m.conditions, rel_means)
        print(f"Wrote table: {args.output_table}")


if __name__ == "__main__":
    main()
