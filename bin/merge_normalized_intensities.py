#!/usr/bin/env python3
"""
Merge ProtStatsWF-normalized intensities back into an unbeQuant quantification table.

ProtStatsWF::workflow_QC returns only an ID column plus the normalized intensity
columns (and it renames the ID column to 'prepared_data$ID'). The downstream
unbequant_ttest_and_plotting.py script, however, needs the full unbeQuant table:
'openms_ceid', the per-file identification columns
(_____l_pep_ident / _____l_raw_pep_ident / _____l_prot_ident / _____charge /
_____l_ms2_scans) and the eight global mz/rt columns.

This script rebuilds that table: it takes the original quantification table and
replaces the values of each '<file>_____intensity' column with the corresponding
normalized (log2) value, joining on the consensus element id. Column names and
column order are preserved, so the result is a drop-in replacement for the
original table with normalized intensities - exactly what
unbequant_ttest_and_plotting.py expects when run with
--log2_intensities_provided.

Also writes a small JSON mapping each condition to its intensity column names,
so that callers can expand a condition (e.g. "C1") into the -gA/-gB column
arguments that the t-test script takes.
"""

import argparse
import csv
import json
import os
import sys

csv.field_size_limit(2 ** 31 - 1)

INTENSITY_SUFFIX = "_____intensity"
ID_COLUMN = "openms_ceid"
# workflow_QC renames the ID column when it writes D_norm_ID
NORM_ID_CANDIDATES = ("prepared_data$ID", "ID", "openms_ceid")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge normalized intensities back into an unbeQuant quantification table"
    )
    parser.add_argument("--reduced_tsv", required=True,
                        help="Input: original quantification table (reduced) from "
                             "consensus_and_features_to_tsv.py")
    parser.add_argument("--normalized_tsv", required=True,
                        help="Input: D_norm_ID*.tsv written by workflow_QC")
    parser.add_argument("--column_map", required=True,
                        help="Input: *_column_map.json written by prepare_qc_input.py")
    parser.add_argument("--output_tsv", required=True,
                        help="Output: quantification table with normalized intensities")
    parser.add_argument("--output_group_columns", default=None,
                        help="Output: JSON mapping condition -> intensity column names "
                             "(default: alongside output_tsv)")
    parser.add_argument("--drop_unnormalized_rows", action="store_true",
                        help="Drop rows that have no normalized counterpart instead of "
                             "writing them with empty intensities")
    return parser.parse_args()


def load_normalized(path):
    """Read D_norm_ID into {id: {normalized_column: value}}."""
    values = {}
    columns = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        # workflow_QC quotes its output; csv handles that, but be defensive
        header = [h.strip().strip('"') for h in header]

        id_idx = 0
        for cand in NORM_ID_CANDIDATES:
            if cand in header:
                id_idx = header.index(cand)
                break
        else:
            print(f"  ! No known ID column in {os.path.basename(path)}; "
                  f"using first column '{header[0]}'")

        columns = [h for i, h in enumerate(header) if i != id_idx]
        for row in reader:
            if not row:
                continue
            key = row[id_idx].strip().strip('"')
            values[key] = {
                header[i]: row[i]
                for i in range(len(header))
                if i != id_idx and i < len(row)
            }
    return values, columns


def main():
    args = parse_args()

    print("=" * 70)
    print("Merging normalized intensities into the quantification table")
    print("=" * 70)

    with open(args.column_map) as fh:
        colmap = json.load(fh)
    # normalized column name (e.g. "C1_01") -> source column ("FLI04006_____intensity")
    norm_to_source = {
        name: meta["source_column"] for name, meta in colmap["columns"].items()
    }
    norm_to_condition = {
        name: meta["condition"] for name, meta in colmap["columns"].items()
    }
    print(f"✓ Column map: {args.column_map}")
    print(f"  {len(norm_to_source)} intensity column(s), "
          f"{len(colmap.get('conditions', []))} condition(s)")

    normalized, norm_columns = load_normalized(args.normalized_tsv)
    print(f"✓ Normalized table: {args.normalized_tsv}")
    print(f"  {len(normalized):,} row(s), {len(norm_columns)} value column(s)")

    unknown = [c for c in norm_columns if c not in norm_to_source]
    if unknown:
        print(f"  ! {len(unknown)} normalized column(s) not in the column map "
              f"(ignored): {', '.join(unknown[:5])}")

    with open(args.reduced_tsv, newline="") as in_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        header = next(reader)

        if ID_COLUMN not in header:
            raise SystemExit(
                f"ERROR: '{ID_COLUMN}' not found in {args.reduced_tsv}; "
                f"unbequant_ttest_and_plotting.py requires it."
            )
        id_idx = header.index(ID_COLUMN)

        # source intensity column -> its index in the reduced table
        source_to_idx = {
            col: idx for idx, col in enumerate(header)
            if col.endswith(INTENSITY_SUFFIX)
        }
        # index in reduced table -> normalized column name
        idx_to_norm = {}
        for norm_name, source_col in norm_to_source.items():
            if source_col in source_to_idx:
                idx_to_norm[source_to_idx[source_col]] = norm_name
            else:
                print(f"  ! Column map references '{source_col}' which is not in "
                      f"the quantification table (skipped)")

        print(f"✓ Quantification table: {args.reduced_tsv}")
        print(f"  {len(header)} column(s), {len(source_to_idx)} intensity column(s), "
              f"{len(idx_to_norm)} matched to normalized values")

        os.makedirs(os.path.dirname(os.path.abspath(args.output_tsv)), exist_ok=True)
        n_rows = n_matched = n_missing = n_written = 0
        n_values = 0
        with open(args.output_tsv, "w", newline="") as out_fh:
            writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            for row in reader:
                if not row:
                    continue
                n_rows += 1
                key = row[id_idx]
                norm_row = normalized.get(key)
                if norm_row is None:
                    n_missing += 1
                    if args.drop_unnormalized_rows:
                        continue
                    # No normalized counterpart: blank the intensities so the row
                    # cannot be mistaken for raw data mixed in with normalized data.
                    out = list(row) + [""] * max(0, len(header) - len(row))
                    for idx in idx_to_norm:
                        if idx < len(out):
                            out[idx] = ""
                else:
                    n_matched += 1
                    out = list(row) + [""] * max(0, len(header) - len(row))
                    for idx, norm_name in idx_to_norm.items():
                        value = norm_row.get(norm_name, "")
                        value = value.strip().strip('"')
                        if value.upper() in ("NA", "NAN", "INF", "-INF"):
                            value = ""
                        if idx < len(out):
                            out[idx] = value
                        if value:
                            n_values += 1
                writer.writerow(out[:len(header)])
                n_written += 1

    print(f"✓ Wrote: {args.output_tsv}")
    print(f"  {n_written:,} row(s) written  ({n_matched:,} normalized, "
          f"{n_missing:,} without a normalized counterpart)")
    print(f"  {n_values:,} normalized intensity value(s)")
    if n_rows:
        print(f"  coverage: {n_matched / n_rows * 100:.1f}% of rows")

    # Condition -> intensity column names, for building -gA/-gB arguments
    group_columns = {}
    for norm_name, condition in norm_to_condition.items():
        source_col = norm_to_source[norm_name]
        if source_col in source_to_idx:
            group_columns.setdefault(condition, []).append(source_col)
    for cond in group_columns:
        group_columns[cond].sort()

    group_path = args.output_group_columns
    if group_path is None:
        group_path = os.path.splitext(args.output_tsv)[0] + "_group_columns.json"
    with open(group_path, "w") as fh:
        json.dump(group_columns, fh, indent=2)
    print(f"✓ Wrote group columns: {group_path}")
    for cond in sorted(group_columns):
        print(f"    {cond}: {len(group_columns[cond])} column(s)")


if __name__ == "__main__":
    main()
