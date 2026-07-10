#!/usr/bin/env python3
"""
Generate summary report for the (non-cyclic) unbeQuant workflow, reading
directly from a raw_quantification_with_identifications_unbequant.tsv file
instead of the JSON component graph the clusterless orchestrator produces.

Each row of the TSV is one component ("openms_ceid"); each sample column
group present in that row (non-empty "<sample>_____openms_fid") is one
vertex of the component. There is no cycling in this workflow, so this
report has no "Cycle Convergence" section and only covers:
- Component Size Distribution
- Peptide Identification Analysis (conflict detection, mass differences,
  cross-component redundancy)

All analysis/plotting logic is imported unchanged from
generate_pairing_report_v2_clusterless.py so the numbers and charts are
produced by the exact same code as the clusterless report.
"""

import argparse
import ast
import csv
import os
import sys
from typing import Dict, List, Optional

csv.field_size_limit(sys.maxsize)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from generate_pairing_report_v2_clusterless import (
    analyze_component_sizes,
    analyze_pep_idents,
    plot_component_size_distribution,
    plot_mass_diff_histogram,
    plot_multi_component_histogram,
)

PER_FILE_FIELDS = [
    "intensity", "l_pep_ident", "l_raw_pep_ident", "l_prot_ident",
    "openms_fid", "charge", "l_ms2_scans", "l_mz_start", "l_mz_end",
    "l_rt_start", "l_rt_end", "l_retention_times", "l_mass_to_charges",
    "l_intensities",
]
FIELD_SUFFIX = "_____"


def _parse_list_field(raw: str) -> Optional[list]:
    """Parse a stringified Python list cell (e.g. "['MGSGIER']"). Empty/'[]' -> None."""
    if not raw or raw == "[]":
        return None
    try:
        parsed = ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return None
    if isinstance(parsed, list) and not parsed:
        return None
    return parsed


def _detect_sample_names(header: List[str]) -> List[str]:
    suffix = f"{FIELD_SUFFIX}openms_fid"
    samples = [col[: -len(suffix)] for col in header if col.endswith(suffix)]
    if not samples:
        raise ValueError("Could not detect any '<sample>_____openms_fid' columns in TSV header")
    return samples


def load_components_from_tsv(tsv_path: str) -> List[Dict]:
    """Load components from a raw_quantification_with_identifications_unbequant.tsv file.

    Returns a list of component dicts shaped like the clusterless orchestrator's
    JSON components, so they can be fed straight into analyze_component_sizes()
    and analyze_pep_idents():
        {'component_id': str, 'num_vertices': int, 'vertices': [ {...}, ... ]}
    """
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"TSV not found: {tsv_path}")

    components: List[Dict] = []

    with open(tsv_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        samples = _detect_sample_names(reader.fieldnames or [])

        for row in reader:
            vertices = []
            for sample in samples:
                fid = row.get(f"{sample}{FIELD_SUFFIX}openms_fid", "")
                if not fid:
                    continue

                intensity_raw = row.get(f"{sample}{FIELD_SUFFIX}intensity", "")
                charge_raw = row.get(f"{sample}{FIELD_SUFFIX}charge", "")

                vertices.append({
                    "filename": sample,
                    "openms_fid": fid,
                    "pep_ident": _parse_list_field(row.get(f"{sample}{FIELD_SUFFIX}l_pep_ident", "")),
                    "prot_ident": _parse_list_field(row.get(f"{sample}{FIELD_SUFFIX}l_prot_ident", "")),
                    "intensity": float(intensity_raw) if intensity_raw else None,
                    "charge": int(charge_raw) if charge_raw else None,
                })

            if not vertices:
                continue

            components.append({
                "component_id": row.get("openms_ceid", "unknown"),
                "num_vertices": len(vertices),
                "vertices": vertices,
            })

    print(f"✓ Loaded {len(components)} components from: {os.path.basename(tsv_path)}")
    return components


def generate_html_report(
    component_stats: Dict,
    pep_ident_stats: Dict,
    output_html: str,
    parameters: Dict = None,
):
    """Generate HTML report (no cycle/convergence section)."""

    size_dist_img = plot_component_size_distribution(component_stats.get("size_distribution", {}))

    html_parts = []

    html_parts.append("""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>UnbeQuant Component Report</title>
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

    html_parts.append("<h1>\U0001F9EC UnbeQuant Component Report</h1>")

    # Overview metrics
    html_parts.append("<h2>\U0001F4CA Overview</h2>")
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

    # Component size distribution
    html_parts.append("<h2>\U0001F4C8 Component Size Distribution</h2>")
    total_c = component_stats["total_components"]
    total_v = component_stats["total_vertices"]
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
    if pep_ident_stats["total_components"] > 0:
        html_parts.append("<h2>\U0001F52C Peptide Identification Analysis</h2>")
        html_parts.append(f"<p>Total components analyzed: <strong>{pep_ident_stats['total_components']}</strong></p>")
        html_parts.append(f"<p>Components with identifications: <strong>{pep_ident_stats['components_with_pep_idents']}</strong></p>")
        html_parts.append(f"<p>Components with conflicting identifications: <strong>{pep_ident_stats['components_with_conflicts']}</strong></p>")

        if pep_ident_stats["components_with_pep_idents"] > 0:
            conflict_pct = pep_ident_stats["conflict_percentage"]
            if conflict_pct > 20:
                status_class = "warning"
                status_msg = f"{conflict_pct:.1f}% of identified components have conflicts"
            elif conflict_pct > 0:
                status_class = "info"
                status_msg = f"{conflict_pct:.1f}% of identified components have conflicts"
            else:
                status_class = "success"
                status_msg = "No conflicts detected!"

            html_parts.append(f'<p class="{status_class}">Conflict rate: {status_msg}</p>')

            mds = pep_ident_stats.get("mass_diff_stats")
            if mds:
                html_parts.append(
                    f'<p>Mean pairwise mass difference between conflicting sequences: '
                    f'<strong>{mds["mean"]:.3f} Da</strong> &nbsp;|&nbsp; '
                    f'Median: <strong>{mds["median"]:.3f} Da</strong> '
                    f'(across {mds["n_components"]:,} conflicted components)</p>'
                )
                mass_hist_img = plot_mass_diff_histogram(pep_ident_stats.get("all_pairwise_diffs", []))
                if mass_hist_img:
                    html_parts.append('<div class="plot">')
                    html_parts.append(
                        f'<img src="data:image/png;base64,{mass_hist_img}" '
                        f'alt="Pairwise Mass Difference Distribution"/>'
                    )
                    html_parts.append('</div>')

            if pep_ident_stats["conflict_examples"]:
                html_parts.append("<h3>Example Conflicts</h3>")
                html_parts.append("<p><em>Components where vertices have different peptide identifications:</em></p>")
                html_parts.append("<table>")
                html_parts.append("<tr><th>Component ID</th><th>Size</th><th>Conflicting Pep_idents</th></tr>")

                for example in pep_ident_stats["conflict_examples"]:
                    comp_id = example["component_id"]
                    size = example["num_vertices"]
                    seq_masses = example.get("pep_ident_masses", {})
                    shown = example["pep_idents"][:5]
                    parts = [
                        f"{seq} [{seq_masses[seq]:.3f}&nbsp;Da]" if seq in seq_masses else seq
                        for seq in shown
                    ]
                    pep_idents = ", ".join(parts)
                    if len(example["pep_idents"]) > 5:
                        pep_idents += f" ... (+{len(example['pep_idents']) - 5} more)"

                    html_parts.append(f"<tr><td>{comp_id}</td><td>{size}</td><td>{pep_idents}</td></tr>")

                html_parts.append("</table>")

            mcs = pep_ident_stats.get("multi_component_stats")
            if mcs:
                html_parts.append("<h3>Pep_ident Redundancy Across Components</h3>")
                html_parts.append(
                    f'<p>Unique pep_idents with identifications: '
                    f'<strong>{mcs["total_unique_pep_idents"]:,}</strong> &nbsp;|&nbsp; '
                    f'Appearing in ≥2 components: '
                    f'<strong>{mcs["total_multi_component"]:,}</strong></p>'
                )
                multi_hist_img = plot_multi_component_histogram(mcs.get("distribution", {}))
                if multi_hist_img:
                    html_parts.append('<div class="plot">')
                    html_parts.append(
                        f'<img src="data:image/png;base64,{multi_hist_img}" '
                        f'alt="Pep_ident Redundancy Histogram"/>'
                    )
                    html_parts.append('</div>')

                top_redundant = mcs.get("top_redundant", [])
                if top_redundant:
                    html_parts.append(f'<h4>Top {len(top_redundant)} Most Redundant Pep_idents</h4>')
                    html_parts.append(
                        '<table>'
                        '<tr><th>#</th><th>Pep_ident sequence</th><th>Components</th></tr>'
                    )
                    for rank, (pid, n_comp) in enumerate(top_redundant, 1):
                        html_parts.append(
                            f'<tr><td>{rank}</td><td><code>{pid}</code></td><td>{n_comp:,}</td></tr>'
                        )
                    html_parts.append('</table>')

            vcs = pep_ident_stats.get("vertex_conflict_stats")
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

    html_parts.append("""
</div>
</body>
</html>
""")

    with open(output_html, "w") as f:
        f.write("".join(html_parts))

    print(f"✓ Generated HTML report: {output_html}")


def generate_text_summary(
    component_stats: Dict,
    pep_ident_stats: Dict,
    output_summary: str,
    parameters: Dict = None,
):
    """Generate plain text summary (no cycle/convergence section)."""

    lines = []
    lines.append("=" * 70)
    lines.append("UNBEQUANT COMPONENT ANALYSIS SUMMARY")
    lines.append("=" * 70)
    lines.append("")

    if parameters:
        lines.append("RUN PARAMETERS")
        lines.append("-" * 70)
        for key, value in parameters.items():
            lines.append(f"  {key:<28} {value}")
        lines.append("")

    lines.append("OVERVIEW")
    lines.append("-" * 70)
    lines.append(f"Final components: {component_stats['total_components']}")
    lines.append(f"Total vertices in components: {component_stats['total_vertices']}")
    lines.append(f"Average component size: {component_stats['avg_size']:.2f}")
    lines.append(f"Median component size: {component_stats['median_size']:.1f}")
    lines.append(f"Size range: {component_stats['min_size']} - {component_stats['max_size']}")
    lines.append(f"Singleton components: {component_stats['singletons']}")
    lines.append("")

    if pep_ident_stats["total_components"] > 0:
        lines.append("PEPTIDE IDENTIFICATION ANALYSIS")
        lines.append("-" * 70)
        lines.append(f"Total components: {pep_ident_stats['total_components']}")
        lines.append(f"Components with identifications: {pep_ident_stats['components_with_pep_idents']}")
        lines.append(f"Components with conflicts: {pep_ident_stats['components_with_conflicts']}")

        if pep_ident_stats["components_with_pep_idents"] > 0:
            lines.append(f"Conflict rate: {pep_ident_stats['conflict_percentage']:.1f}%")

            mds = pep_ident_stats.get("mass_diff_stats")
            if mds:
                lines.append(
                    f"Mass diff (mean):   {mds['mean']:.3f} Da  "
                    f"(median: {mds['median']:.3f} Da, N={mds['n_components']} conflicted components)"
                )

            mcs = pep_ident_stats.get("multi_component_stats")
            if mcs:
                dist = mcs.get("distribution", {})
                dist_str = "  ".join(f"{n} comps: {c:,}" for n, c in sorted(dist.items()))
                lines.append(
                    f"Pep_ident redundancy: {mcs['total_multi_component']:,} of "
                    f"{mcs['total_unique_pep_idents']:,} unique pep_idents appear in ≥2 components"
                )
                if dist_str:
                    lines.append(f"  Distribution: {dist_str}")
                top_redundant = mcs.get("top_redundant", [])
                if top_redundant:
                    lines.append(f"  Top {len(top_redundant)} most redundant pep_idents:")
                    for rank, (pid, n_comp) in enumerate(top_redundant, 1):
                        lines.append(f"    {rank:2d}. {pid}  ({n_comp} components)")

            vcs = pep_ident_stats.get("vertex_conflict_stats")
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

            if pep_ident_stats["conflict_examples"]:
                lines.append("")
                lines.append("Example conflicts (first 5):")
                for i, example in enumerate(pep_ident_stats["conflict_examples"][:5], 1):
                    comp_id = example["component_id"]
                    size = example["num_vertices"]
                    seq_masses = example.get("pep_ident_masses", {})
                    shown = example["pep_idents"][:3]
                    parts = [
                        f"{seq} [{seq_masses[seq]:.3f} Da]" if seq in seq_masses else seq
                        for seq in shown
                    ]
                    pep_idents = ", ".join(parts)
                    if len(example["pep_idents"]) > 3:
                        pep_idents += "..."
                    lines.append(f"  {i}. {comp_id} (size {size}): {pep_idents}")

        lines.append("")

    lines.append("COMPONENT SIZE DISTRIBUTION")
    lines.append("-" * 70)
    total_c = component_stats["total_components"]
    total_v = component_stats["total_vertices"]
    avg_c_size = total_v / total_c if total_c > 0 else 0.0
    lines.append(f"Average component size: {avg_c_size:.2f}  ({total_v:,} vertices / {total_c:,} components)")
    lines.append("")
    size_dist = component_stats.get("size_distribution", {})
    for size in sorted(size_dist.keys())[:20]:
        count = size_dist[size]
        lines.append(f"  Size {size:3d}: {count:6d} components")

    if len(size_dist) > 20:
        lines.append(f"  ... and {len(size_dist) - 20} more sizes")

    lines.append("")
    lines.append("=" * 70)

    with open(output_summary, "w") as f:
        f.write("\n".join(lines))

    print(f"✓ Generated text summary: {output_summary}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate component/peptide-ID report for the unbeQuant workflow from its raw TSV output"
    )
    parser.add_argument("--tsv", required=True,
                       help="Input: raw_quantification_with_identifications_unbequant.tsv")
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

    args = parser.parse_args()

    print("\n" + "=" * 70)
    print("Generating UnbeQuant Component Report")
    print("=" * 70 + "\n")

    components = load_components_from_tsv(args.tsv)

    parameters = None
    if args.parameters_json and os.path.exists(args.parameters_json):
        import json
        with open(args.parameters_json, "r") as f:
            parameters = json.load(f)
        print(f"✓ Loaded parameters from: {os.path.basename(args.parameters_json)}")

    print("\nAnalyzing components...")
    component_stats = analyze_component_sizes(components)
    print(f"  Total components: {component_stats['total_components']}")
    print(f"  Total vertices: {component_stats['total_vertices']}")
    print(f"  Average size: {component_stats['avg_size']:.2f}")

    print("\nAnalyzing peptide identifications...")
    print(f"  N-terminal modification: {args.nterm_mod}")
    print(f"  Cysteine modification:   {args.cys_mod}")
    print(f"  C-terminal modification: {args.cterm_mod}")
    pep_ident_stats = analyze_pep_idents(components, args.nterm_mod, args.cys_mod, args.cterm_mod)
    print(f"  Components with identifications: {pep_ident_stats['components_with_pep_idents']}")
    print(f"  Components with conflicts: {pep_ident_stats['components_with_conflicts']}")
    if pep_ident_stats["components_with_pep_idents"] > 0:
        print(f"  Conflict rate: {pep_ident_stats['conflict_percentage']:.1f}%")

    print("\nGenerating reports...")
    generate_text_summary(component_stats, pep_ident_stats, args.output_summary, parameters)

    if args.output_html:
        generate_html_report(component_stats, pep_ident_stats, args.output_html, parameters)

    print("\n" + "=" * 70)
    print("✅ Report generation complete!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
