#!/usr/bin/env python3
"""
Prepare an unbeQuant quantification table for ProtStatsWF::workflow_QC.

The quantification tables produced by consensus_and_features_to_tsv.py name their
intensity columns "<filename>_____intensity". workflow_QC, however, derives the
sample grouping purely from the column name:

    group <- factor(limma::strsplit2(colnames(D), "_")[,1])

i.e. everything before the FIRST underscore becomes the group. Raw measurement
names therefore group incorrectly (e.g. "FLI04006" becomes its own one-sample
group, and "C2_FLI04008" would group as "C2" only by accident).

This script rewrites the table into a workflow_QC-ready form:
  * intensity columns renamed to <CONDITION>_<REPLICATE>  (e.g. C1_01)
  * ordered by condition, then replicate, so groups are contiguous
  * only 'openms_ceid' plus the intensity columns are kept

Only the ID + intensity columns are forwarded on purpose. The identification
columns hold Python-repr lists such as "['PEPTIDE']"; R's read.table treats the
single quote as a quoting character, which can corrupt field parsing. Those
columns are passengers for workflow_QC anyway, so they are left behind and can
be rejoined later on 'openms_ceid'.

Outputs:
  --output_tsv         table to hand to workflow_QC
  --output_column_map  JSON describing the renaming and the intensity column
                       range (1-based, inclusive) for workflow_QC
"""

import argparse
import csv
import json
import os
import sys

csv.field_size_limit(2 ** 31 - 1)

# Separator used by consensus_and_features_to_tsv.py between file and field name
FIELD_SEP = "_____"
INTENSITY_SUFFIX = FIELD_SEP + "intensity"
ID_COLUMN = "openms_ceid"
UNASSIGNED = "UNASSIGNED"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rename/reorder intensity columns of a quantification table "
                    "for ProtStatsWF::workflow_QC"
    )
    parser.add_argument("--input_tsv", required=True,
                        help="Input: quantification table from "
                             "consensus_and_features_to_tsv.py (reduced table)")
    parser.add_argument("--sample_conditions", required=True,
                        help="Input: TSV mapping filename -> condition "
                             "(see files/sample_conditions.tsv)")
    parser.add_argument("--output_tsv", required=True,
                        help="Output: workflow_QC-ready table")
    parser.add_argument("--output_column_map", default=None,
                        help="Output: JSON describing column renaming and the "
                             "intensity column range (default: alongside output_tsv)")
    parser.add_argument("--strict_names", action="store_true",
                        help="Require the table's file names to match the condition "
                             "file exactly. Without this, a name may also match after "
                             "its leading '<CONDITION>_' prefix is stripped, so tables "
                             "produced before the RAW files were renamed still work.")
    return parser.parse_args()


def load_sample_conditions(path):
    """Parse the condition mapping file.

    Returns (ordered_conditions, filename_to_condition, order_index) where
    order_index preserves the file order within each condition, which is what
    determines the replicate numbering.
    """
    filename_to_condition = {}
    ordered_conditions = []
    order_index = {}

    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 2:
                # Tolerate comma-separated files as well
                parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                raise SystemExit(
                    f"ERROR: {path}:{lineno}: expected 2 tab-separated columns, got: {line!r}"
                )
            fname, condition = parts[0], parts[1]
            if fname.lower() in ("filename", "file", "sample"):
                continue  # header row

            if condition == UNASSIGNED:
                raise SystemExit(
                    f"ERROR: {path}:{lineno}: file '{fname}' is still {UNASSIGNED}. "
                    f"Assign it a condition before running normalization."
                )
            if "_" in condition:
                raise SystemExit(
                    f"ERROR: {path}:{lineno}: condition '{condition}' contains an "
                    f"underscore. workflow_QC derives the group as the text before the "
                    f"first '_', so '{condition}' would silently merge with other "
                    f"conditions sharing the prefix '{condition.split('_')[0]}'."
                )
            if fname in filename_to_condition:
                raise SystemExit(
                    f"ERROR: {path}:{lineno}: duplicate filename '{fname}'"
                )

            filename_to_condition[fname] = condition
            if condition not in ordered_conditions:
                ordered_conditions.append(condition)
            order_index[fname] = len(order_index)

    if not filename_to_condition:
        raise SystemExit(f"ERROR: no usable rows found in {path}")

    return ordered_conditions, filename_to_condition, order_index


def build_name_resolver(filename_to_condition, strict):
    """Return a function mapping a table file name to (mapping_key, condition).

    Falls back to matching with the leading '<CONDITION>_' prefix stripped from
    either side, so tables generated before the RAW files were prefixed still
    resolve. Returns (None, None) when unresolvable.
    """
    aliases = {}
    if not strict:
        for fname, condition in filename_to_condition.items():
            prefix = condition + "_"
            if fname.startswith(prefix):
                aliases.setdefault(fname[len(prefix):], (fname, condition))

    def resolve(table_name):
        if table_name in filename_to_condition:
            return table_name, filename_to_condition[table_name]
        if strict:
            return None, None
        # Table name is unprefixed, mapping key is prefixed
        if table_name in aliases:
            return aliases[table_name]
        # Table name is prefixed, mapping key is not
        if "_" in table_name:
            stripped = table_name.split("_", 1)[1]
            if stripped in filename_to_condition:
                return stripped, filename_to_condition[stripped]
        return None, None

    return resolve


def main():
    args = parse_args()

    print("=" * 70)
    print("Preparing quantification table for ProtStatsWF::workflow_QC")
    print("=" * 70)

    ordered_conditions, file_to_cond, order_index = load_sample_conditions(
        args.sample_conditions
    )
    print(f"✓ Conditions file: {args.sample_conditions}")
    print(f"  {len(file_to_cond)} file(s) across {len(ordered_conditions)} condition(s): "
          f"{', '.join(ordered_conditions)}")

    resolve = build_name_resolver(file_to_cond, args.strict_names)

    with open(args.input_tsv, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"ERROR: {args.input_tsv} is empty")

        if header[0] != ID_COLUMN:
            print(f"  ! Expected first column '{ID_COLUMN}', found '{header[0]}' "
                  f"— using it as the ID column anyway")
        id_col_idx = 0

        # Locate the intensity block
        intensity_cols = [
            (idx, col[: -len(INTENSITY_SUFFIX)])
            for idx, col in enumerate(header)
            if col.endswith(INTENSITY_SUFFIX)
        ]
        if not intensity_cols:
            raise SystemExit(
                f"ERROR: no '*{INTENSITY_SUFFIX}' columns found in {args.input_tsv}"
            )
        print(f"✓ Input table: {args.input_tsv}")
        print(f"  {len(header)} column(s), {len(intensity_cols)} intensity column(s)")

        # Resolve every intensity column to a condition
        resolved = []      # (source_idx, table_name, mapping_key, condition)
        unresolved = []
        used_aliases = []
        for idx, table_name in intensity_cols:
            mapping_key, condition = resolve(table_name)
            if condition is None:
                unresolved.append(table_name)
                continue
            if mapping_key != table_name:
                used_aliases.append((table_name, mapping_key))
            resolved.append((idx, table_name, mapping_key, condition))

        if unresolved:
            raise SystemExit(
                "ERROR: these files are in the quantification table but not in "
                f"{args.sample_conditions}:\n  "
                + "\n  ".join(unresolved)
                + "\n\nAdd them to the conditions file (or re-run with matching names)."
            )

        if used_aliases:
            print(f"  ! Matched {len(used_aliases)} file(s) by stripping the condition "
                  f"prefix (table was generated before the RAW rename):")
            for table_name, mapping_key in used_aliases[:5]:
                print(f"      {table_name}  ->  {mapping_key}")
            if len(used_aliases) > 5:
                print(f"      ... and {len(used_aliases) - 5} more")

        mapped_keys = {mapping_key for _, _, mapping_key, _ in resolved}
        missing_from_table = [f for f in file_to_cond if f not in mapped_keys]
        if missing_from_table:
            print(f"  ! {len(missing_from_table)} file(s) in the conditions file are "
                  f"absent from this table (ignored): "
                  f"{', '.join(missing_from_table[:5])}"
                  + (" ..." if len(missing_from_table) > 5 else ""))

        # Order by condition (as listed in the conditions file), then by the file
        # order within that condition; assign replicate numbers accordingly.
        cond_rank = {c: i for i, c in enumerate(ordered_conditions)}
        resolved.sort(key=lambda r: (cond_rank[r[3]], order_index[r[2]]))

        per_condition_count = {}
        new_names = []
        column_map = {}
        for source_idx, table_name, mapping_key, condition in resolved:
            per_condition_count[condition] = per_condition_count.get(condition, 0) + 1
            new_name = f"{condition}_{per_condition_count[condition]:02d}"
            new_names.append((source_idx, new_name))
            column_map[new_name] = {
                "filename": mapping_key,
                "table_name": table_name,
                "condition": condition,
                "replicate": per_condition_count[condition],
                "source_column": header[source_idx],
            }

        out_header = [ID_COLUMN] + [name for _, name in new_names]
        source_indices = [idx for idx, _ in new_names]

        # Intensity columns are contiguous starting right after the ID column.
        intensity_start = 2
        intensity_end = 1 + len(new_names)

        os.makedirs(os.path.dirname(os.path.abspath(args.output_tsv)), exist_ok=True)
        n_rows = 0
        n_values = 0
        with open(args.output_tsv, "w", newline="") as out_fh:
            writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
            writer.writerow(out_header)
            for row in reader:
                if not row:
                    continue
                # Guard against ragged rows
                values = [
                    row[idx] if idx < len(row) else "" for idx in source_indices
                ]
                n_values += sum(1 for v in values if v not in ("", "None"))
                # Normalize Python's None placeholder to an empty field so R reads NA
                values = ["" if v == "None" else v for v in values]
                writer.writerow([row[id_col_idx] if row else ""] + values)
                n_rows += 1

    print(f"✓ Wrote: {args.output_tsv}")
    print(f"  {n_rows:,} row(s) x {len(out_header)} column(s)")
    print(f"  {n_values:,} non-empty intensity value(s)")
    print(f"  intensity_columns for workflow_QC: {intensity_start}:{intensity_end}")

    for condition in ordered_conditions:
        cols = [n for n, m in column_map.items() if m["condition"] == condition]
        if cols:
            print(f"    {condition}: {len(cols)} replicate(s) — {cols[0]} .. {cols[-1]}")

    column_map_path = args.output_column_map
    if column_map_path is None:
        base = os.path.splitext(args.output_tsv)[0]
        column_map_path = base + "_column_map.json"
    with open(column_map_path, "w") as fh:
        json.dump(
            {
                "id_column": ID_COLUMN,
                "intensity_column_start": intensity_start,
                "intensity_column_end": intensity_end,
                "n_intensity_columns": len(new_names),
                "conditions": ordered_conditions,
                "columns": column_map,
            },
            fh,
            indent=2,
        )
    print(f"✓ Wrote column map: {column_map_path}")


if __name__ == "__main__":
    main()
