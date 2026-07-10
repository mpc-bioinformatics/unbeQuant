#!/usr/bin/env python3
"""
Scan all parameter grid run directories and compile key metrics into a single
comparison table.  Reads clusterless_summary.txt from each run.

Usage:
    python3 bin/summarize_grid_results.py
    python3 bin/summarize_grid_results.py --results_dir results/parameter_grid
    python3 bin/summarize_grid_results.py --output comparison.txt
"""

import argparse
import os
import re
import sys


METRICS = [
    ("avg_component_size", r"^Average component size:\s*([\d.]+)"),
    ("median_component_size", r"^Median component size:\s*([\d.]+)"),
    ("final_components", r"^Final components:\s*([\d,]+)"),
    ("total_vertices", r"^Total vertices in components:\s*([\d,]+)"),
    ("conflict_rate", r"^Conflict rate:\s*([\d.]+)%"),
    ("components_with_conflicts", r"^Components with conflicts:\s*([\d,]+)"),
    ("components_with_pep_idents", r"^Components with identifications:\s*([\d,]+)"),
]


def parse_summary(path: str) -> dict:
    result = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.rstrip()
                for key, pattern in METRICS:
                    if key not in result:
                        m = re.match(pattern, line)
                        if m:
                            result[key] = m.group(1).replace(",", "")
    except OSError:
        pass
    return result


def main():
    parser = argparse.ArgumentParser(description="Summarise parameter grid results")
    parser.add_argument("--results_dir", default="results/parameter_grid",
                        help="Root directory containing run sub-directories")
    parser.add_argument("--output", default=None,
                        help="Output file path (default: print to stdout and write to <results_dir>/grid_comparison.txt)")
    args = parser.parse_args()

    results_dir = args.results_dir
    if not os.path.isdir(results_dir):
        print(f"ERROR: results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)

    # Collect all runs
    runs = []
    for entry in sorted(os.scandir(results_dir), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        summary_path = os.path.join(entry.path, "feature_analysis_clusterless", "clusterless_summary.txt")
        metrics = parse_summary(summary_path)
        runs.append((entry.name, summary_path, metrics))

    if not runs:
        print("No run directories found.", file=sys.stderr)
        sys.exit(1)

    # Build output lines
    lines = []
    lines.append("=" * 110)
    lines.append("PARAMETER GRID COMPARISON")
    lines.append("=" * 110)
    lines.append(f"Results directory: {os.path.abspath(results_dir)}")
    lines.append(f"Runs found: {len(runs)}")
    lines.append("")

    col_run   = 55
    col_avg   = 10
    col_med   = 10
    col_comp  = 10
    col_vert  = 12
    col_conf  = 12
    col_avail = 10

    header = (
        f"{'Run':<{col_run}} "
        f"{'AvgSize':>{col_avg}} "
        f"{'MedSize':>{col_med}} "
        f"{'#Comp':>{col_comp}} "
        f"{'Vertices':>{col_vert}} "
        f"{'ConflictRate':>{col_conf}} "
        f"{'w/Idents':>{col_avail}}"
    )
    lines.append(header)
    lines.append("-" * 110)

    for run_name, summary_path, m in runs:
        exists = os.path.isfile(summary_path)
        if not exists:
            status = "(no summary yet)"
        elif not m:
            status = "(parse error)"
        else:
            status = None

        if status:
            lines.append(f"{run_name:<{col_run}} {status}")
            continue

        avg   = m.get("avg_component_size", "—")
        med   = m.get("median_component_size", "—")
        comp  = m.get("final_components", "—")
        vert  = m.get("total_vertices", "—")
        conf  = m.get("conflict_rate", "—")
        if conf != "—":
            conf = conf + "%"
        avail = m.get("components_with_pep_idents", "—")

        lines.append(
            f"{run_name:<{col_run}} "
            f"{avg:>{col_avg}} "
            f"{med:>{col_med}} "
            f"{comp:>{col_comp}} "
            f"{vert:>{col_vert}} "
            f"{conf:>{col_conf}} "
            f"{avail:>{col_avail}}"
        )

    lines.append("=" * 110)
    lines.append("")

    text = "\n".join(lines)
    print(text)

    out_path = args.output or os.path.join(results_dir, "grid_comparison.txt")
    with open(out_path, "w") as f:
        f.write(text)
    print(f"Saved to: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
