#!/usr/bin/env python3
"""
Convert an OpenMS consensusXML into the components JSON that
analyze_pairwise_components.py consumes.

The original (OpenMS/FeatureLinker) workflow expresses its pairing result as a
consensusXML: one <consensusElement> per component, one <element> per member
feature. The clusterless workflow expresses the same thing as a components JSON.
Converting the former into the latter lets the pairwise within-component
distance analysis run on both, so the two workflows can be compared on identical
figures.

Field mapping (element attribute -> vertex field):
    mz  -> x_center      the member feature's m/z
    rt  -> y_center      the member feature's retention time, in seconds
    it  -> intensity
    map -> filename      resolved through <mapList>

The pairing cutoffs are read from the FeatureLinker parameters recorded in the
file's own dataProcessing block, so the figures are drawn against the cutoffs
that run actually used. --ppm/--rt override them.

RT CAVEAT: consensusXML members come from already RT-aligned featureXMLs, so the
retention times here are POST-alignment. analyze_pairwise_components.py notes
that clusterless y_center is exported PRE-alignment. The m/z panels compare
directly; the RT panels do not, and the difference is the alignment offset
(up to ~30 s per edge) rather than a property of either workflow.

OUTPUT LAYOUT is not free-form: the analyzer scans the JSON line-wise rather
than parsing it, so every field must sit on its own line (indent=2) and
"component_id" must precede its vertices - it is the record separator that
flushes the previous component.
"""

import argparse
import json
import os
import re
import sys
import xml.etree.ElementTree as ET


def _open_xml(path):
    """Open `path` positioned at the XML declaration.

    Tolerates junk before the declaration (a stray line typed into the file by
    an editor, for instance), which would otherwise abort the parse at line 1.
    """
    fh = open(path, "rb")
    head = fh.read(8192)
    idx = head.find(b"<?xml")
    if idx < 0:
        idx = head.find(b"<consensusXML")
    if idx < 0:
        fh.seek(0)
        return fh, 0
    fh.seek(idx)
    return fh, idx


def read_header(path):
    """Return (map_id -> filename, cutoffs) from the leading header.

    Stops at the consensus element list: everything needed sits before it, and
    the element list is the bulk of a multi-hundred-MB file.
    """
    maps = {}
    cutoffs = {}
    map_re = re.compile(r'<map id="(\d+)"[^>]*name="([^"]*)"')
    par_re = re.compile(r'name="parameter: algorithm:(distance_(?:MZ|RT):max_difference|'
                        r'distance_MZ:unit)" value="([^"]*)"')

    fh, _ = _open_xml(path)
    with fh:
        for raw in fh:
            line = raw.decode("latin-1")
            if "<consensusElementList" in line:
                break
            m = map_re.search(line)
            if m:
                name = os.path.basename(m.group(2))
                for ext in (".mzML", ".mzml", ".featureXML"):
                    if name.endswith(ext):
                        name = name[: -len(ext)]
                        break
                maps[m.group(1)] = name
            p = par_re.search(line)
            if p:
                cutoffs[p.group(1)] = p.group(2)
    return maps, cutoffs


def count_records(path):
    """Fast pre-pass for the metadata counts (byte scan, no XML parsing)."""
    n_comp = n_vert = 0
    fh, _ = _open_xml(path)
    with fh:
        for raw in fh:
            if b"<consensusElement id=" in raw:
                n_comp += 1
            elif b"<element map=" in raw:
                n_vert += 1
    return n_comp, n_vert


def iter_components(path):
    """Yield one component dict per <consensusElement>, streaming.

    Elements are cleared as they close, so peak memory stays flat regardless of
    file size.
    """
    fh, _ = _open_xml(path)
    with fh:
        for _, el in ET.iterparse(fh, events=("end",)):
            if el.tag != "consensusElement":
                continue

            vertices = []
            grouped = el.find("groupedElementList")
            if grouped is not None:
                for i, e in enumerate(grouped.findall("element")):
                    a = e.attrib
                    try:
                        mz = float(a["mz"])
                        rt = float(a["rt"])
                    except (KeyError, ValueError):
                        continue
                    try:
                        intensity = float(a.get("it", "nan"))
                    except ValueError:
                        intensity = None
                    try:
                        charge = int(a.get("charge", "0")) or None
                    except ValueError:
                        charge = None
                    vertices.append({
                        "filename": a.get("map", ""),
                        "openms_fid": a.get("id", ""),
                        "x_center": mz,
                        "y_center": rt,
                        "intensity": intensity,
                        "charge": charge,
                        "vertex_id": i,
                    })

            seqs = []
            for hit in el.iter("PeptideHit"):
                seq = hit.get("sequence")
                if seq and seq not in seqs:
                    seqs.append(seq)

            yield {
                "component_id": el.get("id", ""),
                "num_vertices": len(vertices),
                "pep_ident": seqs,
                "vertices": vertices,
            }
            el.clear()


def convert(xml_path, out_path, ppm_override=None, rt_override=None,
            progress_every=20000):
    maps, cut = read_header(xml_path)
    print(f"✓ {len(maps)} maps in <mapList>")

    ppm = ppm_override
    if ppm is None and cut.get("distance_MZ:max_difference"):
        # Only meaningful as a ppm cutoff when the unit says so
        if cut.get("distance_MZ:unit", "ppm").lower() == "ppm":
            ppm = float(cut["distance_MZ:max_difference"])
    rt_cutoff = rt_override
    if rt_cutoff is None and cut.get("distance_RT:max_difference"):
        rt_cutoff = float(cut["distance_RT:max_difference"])
    print(f"✓ cutoffs: ppm={ppm}  rt={rt_cutoff}")

    n_comp, n_vert = count_records(xml_path)
    print(f"✓ {n_comp:,} consensus elements / {n_vert:,} member features")

    metadata = {
        "source_consensusxml": os.path.basename(xml_path),
        "workflow": "original_openms_featurelinker",
        "total_components": n_comp,
        "total_vertices_in_components": n_vert,
        "cutoff_params": {
            "mz_cutoff": None,
            "ppm": ppm,
            "rt_cutoff": rt_cutoff,
            "edges_cutoff": None,
        },
        "rt_reference": "post_alignment",
    }

    written = paired = 0
    with open(out_path, "w") as out:
        out.write("{\n")
        out.write('  "metadata": ' + json.dumps(metadata, indent=2).replace("\n", "\n  "))
        out.write(",\n")
        out.write('  "components": {\n')

        first = True
        for comp in iter_components(xml_path):
            for v in comp["vertices"]:
                v["filename"] = maps.get(v["filename"], v["filename"])
            if not first:
                out.write(",\n")
            first = False
            body = json.dumps(comp, indent=2).replace("\n", "\n    ")
            out.write("    " + json.dumps(comp["component_id"]) + ": " + body)
            written += 1
            if comp["num_vertices"] >= 2:
                paired += 1
            if progress_every and written % progress_every == 0:
                print(f"  ... {written:,} components", flush=True)

        out.write("\n  }\n}\n")

    # The counting pass and the writing pass read the file separately, so a
    # disagreement means the consensusXML changed underneath us mid-conversion
    # (an editor saving it, an upstream job rewriting it). The output would be
    # a blend of two file states, so refuse it rather than hand back something
    # that merely looks complete.
    if written != n_comp:
        raise SystemExit(
            f"ERROR: counted {n_comp:,} consensus elements but wrote {written:,} "
            f"components. {xml_path} changed while it was being read - rerun "
            f"once nothing else is writing it. {out_path} is incomplete."
        )

    print(f"✓ wrote {written:,} components ({paired:,} with >=2 members) -> {out_path}")
    return written, paired


def main():
    p = argparse.ArgumentParser(
        description="Convert an OpenMS consensusXML into the components JSON "
                    "consumed by analyze_pairwise_components.py"
    )
    p.add_argument("--consensusxml", required=True, help="input consensusXML")
    p.add_argument("--out", required=True, help="output components JSON")
    p.add_argument("--ppm", type=float, default=None,
                   help="ppm cutoff to record (default: from the file's FeatureLinker params)")
    p.add_argument("--rt", type=float, default=None,
                   help="RT cutoff in seconds (default: from the file's FeatureLinker params)")
    args = p.parse_args()

    if not os.path.exists(args.consensusxml):
        sys.exit(f"ERROR: consensusXML not found: {args.consensusxml}")

    written, paired = convert(args.consensusxml, args.out, args.ppm, args.rt)
    if paired == 0:
        sys.exit("ERROR: no component has >=2 members - nothing for the pairwise analysis to measure")


if __name__ == "__main__":
    main()
