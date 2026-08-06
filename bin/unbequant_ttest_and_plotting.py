import argparse
import csv
import json
import math
import os
import sys

import tqdm
import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import false_discovery_control
import pandas as pd
from ast import literal_eval

csv.field_size_limit(sys.maxsize)

NONE_VALUES = {
    "nan",
    "inf",
    "#n/a",
    "n/a",
    "na",
    "",
}


def parse_ms2_scans(cells: str, origin_str: str):
    """ Parse the MS2 scans from a list of cell entries """
    ms2s = set()
    ms2_with_origin = []
    for cell, origin in zip(cells, origin_str):
        if cell:
            if cell.lower() in NONE_VALUES:
                ms2_with_origin.append(None)
                continue
            else:
                try:
                    ms2s = ",".join(str(x) for x in literal_eval(cell))
                    if ms2s:
                        ms2_with_origin.append(origin + ":" + ms2s)
                except (ValueError, SyntaxError):
                    # if the cell is not parsable it is a malformed list string element (we skip it)
                    ms2_with_origin.append(None)
    return ms2_with_origin


def parse_pep_prot_ident(cells: str):
    """ Parse the peptide/protein identification from a list of cell entries """
    peps_prots = set()
    for cell in cells:
        if cell:
            if cell.lower() in NONE_VALUES:
                continue
            else:
                try:
                    for c in literal_eval(cell):
                        peps_prots.add(c)
                except (ValueError, SyntaxError):
                    # if the cell is not parsable it is a malformed list string element (we skip it)
                    continue

    return peps_prots

def parse_charge(cell: str):
    """ Parse the charge value from a cell, returns None if not parsable"""
    if cell:
        try:
            return int(cell)
        except ValueError:
            # if the cell is not parsable (e.g. #N/A or similar, we return None)
            return None
    else:
        return None

def parse_intensity(cell: str):
    """ Parse the intensity value from a cell, returns None if the cell is empty or not parsable """
    if cell:
        if cell.lower() == "nan":
            return None
        if cell.lower() == "inf":
            return None
        try:
            return float(cell)
        except ValueError:
            # if the cell is not parsable (e.g. #N/A or similar, we return None)
            return None
    else:
        return None


def parse_args():
    """ Parse Arguments """
    parser = argparse.ArgumentParser(
        description="A small script to calculate the p-value in a CSV-file across two groups (same size), using column-headers."
    )

    # Input CSV
    parser.add_argument(
        "--input_unbequant_tsv", "-in", type=argparse.FileType("r"),
        help="The Output of Unbequant (intensities can be normalized before)."
    )

    parser.add_argument(
        "--groupA", "-gA", type=str, action="append",
        help="Group A. Specify here the column-headers, which belong to this group"
    )
    parser.add_argument(
        "--groupB", "-gB", type=str, action="append",
        help="Group B. Specify here the column-headers, which belong to this group"
    )

    parser.add_argument(
        "--percentage_of_missingness_per_group_allowed", "-p", type=float, default=0,
        help="A percentage (0 - 1), which is allowed to be missing in a group. Set it to 0.3 to allow 30 percent of missing values in a group"
    )

    parser.add_argument(
        "--log2_intensities_provided", "-lip", type=bool, default=True,
        help="Are there log2 intensities provided? If yes, set this to 1. "
    )

    # Output-CSV
    parser.add_argument(
        "--output", "-out", type=str, default=sys.stdout,
        help="Output-TSV containing the same columns and rows, with additional p_value, corrected p_value and fold_change"
        " Defaults to stdout."
    )

    # Output-CSV
    parser.add_argument(
        "--output_volcano", "-volcano", type=str, default=None,
        help="Output-HTML of the volcano plot. "
        " Defaults to None (--> No plots generated)."
    )

    return parser.parse_args()



def coloring_volcano(x):
    if x == 6:
        # is not identified and not significant
        return "Not Identified and not Significant"
    elif x == 7:
        # is identified and not significant
        return "Identified and not Significant"
    elif x == 11:
        # is not identified and significant
        return  "Not Identified and Significant"
    elif x == 12:
        # is identified and significant
        return "Identified and Significant"


VOLCANO_COLOR_MAP = {
    'Not Identified and not Significant': '#365575',
    'Identified and not Significant': '#4f7832',
    'Not Identified and Significant': '#076ddb',
    'Identified and Significant': '#59cc06',
}

# Every value a point in the volcano plot carries, in the order the expanded
# entries of the search result list show them.
VOLCANO_DETAIL_COLUMNS = [
    "openmsid",
    "fold_change_A_div_B",
    "ttest_ind_pvalue",
    "ttest_ind_corrected_pvalue",
    "-log10(pvalue)",
    "ttest_ind_statistic",
    "missing_values_in_a",
    "missing_values_in_b",
    "charge",
    "first_iso_global_min_mz",
    "first_iso_global_max_mz",
    "first_iso_global_min_rt",
    "first_iso_global_max_rt",
    "feature_global_min_mz",
    "feature_global_max_mz",
    "feature_global_min_rt",
    "feature_global_max_rt",
    "pep_ident",
    "raw_pep_ident",
    "prot_ident",
    "ms2s",
]

# The four values a collapsed entry of the search result list shows.
VOLCANO_SUMMARY_COLUMNS = [
    "fold_change_A_div_B",
    "ttest_ind_pvalue",
    "raw_pep_ident",
    "prot_ident",
]


VOLCANO_PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>__TITLE__</title>
<style>
  * { box-sizing: border-box; }
  html, body { height: 100%; margin: 0; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: 13px;
    color: #1f2933;
    background: #ffffff;
  }
  #app { display: flex; height: 100%; width: 100%; }
  #sidebar {
    width: 360px;
    min-width: 240px;
    max-width: 70%;
    flex: 0 0 auto;
    display: flex;
    flex-direction: column;
    border-right: 1px solid #d9e0e7;
    background: #f7f9fb;
    resize: horizontal;
    overflow: hidden;
  }
  #search-panel { padding: 12px 12px 8px 12px; border-bottom: 1px solid #e3e9ef; }
  #search-row { display: flex; align-items: center; gap: 10px; }
  #search-input {
    flex: 1 1 auto;
    min-width: 0;
    padding: 7px 9px;
    font-size: 13px;
    font-family: inherit;
    border: 1px solid #c3ced8;
    border-radius: 4px;
    background: #fff;
  }
  #search-input:focus { outline: none; border-color: #076ddb; box-shadow: 0 0 0 2px rgba(7,109,219,.15); }
  #unique-label {
    display: flex;
    align-items: center;
    gap: 5px;
    white-space: nowrap;
    cursor: pointer;
    user-select: none;
  }
  #unique-label input { cursor: pointer; margin: 0; }
  #status { margin-top: 8px; color: #5a6b7b; min-height: 16px; }
  #status .count { font-weight: 600; color: #1f2933; }
  /* min-height:0 lets the list shrink inside the flex column so that it scrolls */
  #hits { flex: 1 1 auto; min-height: 0; overflow-y: auto; overflow-x: hidden; padding: 6px; }
  .hit {
    border: 1px solid #e0e7ed;
    border-left-width: 4px;
    border-radius: 4px;
    background: #fff;
    margin-bottom: 5px;
  }
  .hit-head {
    display: flex;
    gap: 6px;
    padding: 6px 8px;
    cursor: pointer;
  }
  .hit-head:hover { background: #eef4fa; }
  .hit-caret { color: #8b9bab; flex: 0 0 auto; line-height: 1.5; display: inline-block; transition: transform .1s ease; }
  .hit.open > .hit-head > .hit-caret { transform: rotate(90deg); }
  .hit-main { flex: 1 1 auto; min-width: 0; }
  .hit-prot { font-weight: 600; word-break: break-word; }
  .hit-nums { color: #46586a; margin-top: 2px; }
  .hit-nums b { font-variant-numeric: tabular-nums; color: #1f2933; }
  .hit-nums span + span { margin-left: 10px; }
  .hit-pep { color: #6b7c8d; margin-top: 2px; word-break: break-all; font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 11.5px; }
  .hit-body { display: none; border-top: 1px solid #eef2f6; padding: 6px 8px 8px 8px; }
  .hit.open > .hit-body { display: block; }
  .hit-body table { width: 100%; border-collapse: collapse; }
  .hit-body th, .hit-body td { text-align: left; vertical-align: top; padding: 2px 4px; }
  .hit-body th { width: 42%; font-weight: 500; color: #5a6b7b; white-space: nowrap; }
  .hit-body td { word-break: break-all; font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 11.5px; }
  .hit-body tr:nth-child(even) { background: #f7f9fb; }
  #more {
    display: none;
    width: calc(100% - 4px);
    margin: 2px;
    padding: 6px;
    font: inherit;
    color: #076ddb;
    background: #fff;
    border: 1px dashed #c3ced8;
    border-radius: 4px;
    cursor: pointer;
  }
  #more:hover { background: #eef4fa; }
  #plot-area { flex: 1 1 auto; min-width: 0; height: 100%; }
  #plot-area > div { width: 100%; height: 100%; }
</style>
</head>
<body>
<div id="app">
  <aside id="sidebar">
    <div id="search-panel">
      <div id="search-row">
        <input id="search-input" type="search" placeholder="Search prot_ident, comma-separated&hellip;" autocomplete="off" spellcheck="false" title="One or more protein identifiers separated by commas; a point matches when any of them is found"/>
        <label id="unique-label" title="Only show points that map to exactly one protein accession, i.e. drop everything identified as several proteins or as an ambiguous protein group">
          <input id="unique-only" type="checkbox"/>Unique Only
        </label>
      </div>
      <div id="status"></div>
    </div>
    <div id="hits"></div>
    <button id="more" type="button"></button>
  </aside>
  <main id="plot-area">
__PLOT_DIV__
  </main>
</div>
<script id="volcano-data" type="application/json">__DATA_JSON__</script>
<script>
(function () {
  var DIM_OPACITY = 0.06;
  var HIT_OPACITY = 1.0;
  var PAGE_SIZE = 200;

  var data = JSON.parse(document.getElementById("volcano-data").textContent);
  var gd = document.getElementById("__PLOT_ID__");
  var searchInput = document.getElementById("search-input");
  var uniqueBox = document.getElementById("unique-only");
  var statusEl = document.getElementById("status");
  var hitsEl = document.getElementById("hits");
  var moreBtn = document.getElementById("more");

  var nPoints = data.n_points;
  var colIndex = {};
  data.columns.forEach(function (name, i) { colIndex[name] = i; });

  // Split the prot_ident cells into their accessions once, so that every
  // keystroke only walks pre-lowercased lists. ";" separates the identifications
  // of a point, "," the members of an ambiguous protein group within one of
  // them, and both make an identification non-unique.
  var protLists = new Array(nPoints);
  var protCells = data.values[colIndex[data.prot_column]];
  for (var i = 0; i < nPoints; i++) {
    var protCell = protCells[i];
    protLists[i] = protCell ? String(protCell).toLowerCase().split(/[;,]/).filter(function (p) { return p.length > 0; }) : [];
  }

  var traceIndices = data.trace_counts.map(function (_, i) { return i; });
  var baseOpacity = traceIndices.map(function (t) {
    var marker = (gd.data && gd.data[t] && gd.data[t].marker) || {};
    return marker.opacity === undefined ? null : marker.opacity;
  });

  var hits = [];
  var rendered = 0;
  var dimmed = false;

  function escapeHtml(value) {
    return String(value).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }

  function format(value) {
    if (value === null || value === undefined || value === "") return "\\u2013";
    if (typeof value !== "number") return String(value);
    if (!isFinite(value) || value === 0) return String(value);
    var magnitude = Math.abs(value);
    if (magnitude < 1e-3 || magnitude >= 1e5) return value.toExponential(3);
    return String(Math.round(value * 1e6) / 1e6);
  }

  function cell(row, column) {
    var value = data.values[colIndex[column]][row];
    return column === data.ms2s_column ? expandOrigins(value) : value;
  }

  // Puts the file names back in front of the ms2 scan lists
  function expandOrigins(value) {
    if (!value) return value;
    return String(value).split(";").map(function (part) {
      var colon = part.indexOf(":");
      var origin = data.ms2_origins[Number(part.slice(0, colon))];
      return (colon === -1 || origin === undefined) ? part : origin + ":" + part.slice(colon + 1);
    }).join(";");
  }

  // Several identifiers can be searched at once by separating them with a comma
  function parseQuery(raw) {
    return raw.toLowerCase().split(",").map(function (term) {
      return term.trim();
    }).filter(function (term) {
      return term.length > 0;
    });
  }

  function matches(row, terms, uniqueOnly) {
    var prots = protLists[row];
    if (uniqueOnly && prots.length !== 1) return false;
    for (var i = 0; i < prots.length; i++) {
      for (var j = 0; j < terms.length; j++) {
        if (prots[i].indexOf(terms[j]) !== -1) return true;
      }
    }
    return false;
  }

  function summaryHtml(row) {
    var nums = data.summary_columns.filter(function (c) { return typeof cell(row, c) === "number"; });
    var texts = data.summary_columns.filter(function (c) { return nums.indexOf(c) === -1; });
    var html = '<div class="hit-prot">' + escapeHtml(cell(row, data.prot_column) || "\\u2013") + "</div>";
    html += '<div class="hit-nums">' + nums.map(function (c) {
      return "<span>" + escapeHtml(data.summary_labels[c] || c) + " <b>" + escapeHtml(format(cell(row, c))) + "</b></span>";
    }).join("") + "</div>";
    texts.forEach(function (c) {
      if (c === data.prot_column) return;
      html += '<div class="hit-pep">' + escapeHtml(data.summary_labels[c] || c) + ": " + escapeHtml(cell(row, c) || "\\u2013") + "</div>";
    });
    return html;
  }

  function detailsHtml(row) {
    var rows = data.detail_columns.map(function (column) {
      return "<tr><th>" + escapeHtml(column) + "</th><td>" + escapeHtml(format(cell(row, column))) + "</td></tr>";
    }).join("");
    return "<table>" + rows + "</table>";
  }

  function renderPage() {
    var end = Math.min(hits.length, rendered + PAGE_SIZE);
    var html = "";
    for (var i = rendered; i < end; i++) {
      var row = hits[i];
      html += '<div class="hit" data-row="' + row + '" style="border-left-color:' +
        escapeHtml(data.trace_colors[data.point_trace[row]] || "#c3ced8") + '">' +
        '<div class="hit-head"><span class="hit-caret">\\u25b8</span><div class="hit-main">' +
        summaryHtml(row) + "</div></div><div class=\\"hit-body\\"></div></div>";
    }
    hitsEl.insertAdjacentHTML("beforeend", html);
    rendered = end;
    if (rendered < hits.length) {
      moreBtn.style.display = "block";
      moreBtn.textContent = "Show " + Math.min(PAGE_SIZE, hits.length - rendered) + " more (" + (hits.length - rendered) + " remaining)";
    } else {
      moreBtn.style.display = "none";
    }
  }

  function highlight(active) {
    if (!active) {
      // Redrawing a plot that is already undimmed only costs time
      if (dimmed) Plotly.restyle(gd, { "marker.opacity": baseOpacity }, traceIndices);
      dimmed = false;
      return;
    }
    dimmed = true;
    var opacities = data.trace_counts.map(function (count) {
      var arr = new Array(count);
      for (var i = 0; i < count; i++) arr[i] = DIM_OPACITY;
      return arr;
    });
    for (var i = 0; i < hits.length; i++) {
      opacities[data.point_trace[hits[i]]][data.point_index[hits[i]]] = HIT_OPACITY;
    }
    Plotly.restyle(gd, { "marker.opacity": opacities }, traceIndices);
  }

  function search() {
    var terms = parseQuery(searchInput.value);
    var uniqueOnly = uniqueBox.checked;

    hits = [];
    if (terms.length) {
      for (var row = 0; row < nPoints; row++) {
        if (matches(row, terms, uniqueOnly)) hits.push(row);
      }
    }

    hitsEl.innerHTML = "";
    rendered = 0;
    moreBtn.style.display = "none";

    if (!terms.length) {
      statusEl.innerHTML = "Type one or more protein identifiers, separated by commas, to highlight their points.";
    } else {
      statusEl.innerHTML = '<span class="count">' + hits.length + "</span> of " + nPoints + " points match" +
        (terms.length > 1 ? " any of " + terms.length + " identifiers" : "") +
        (uniqueOnly ? " (single-accession points only)" : "") + ".";
      renderPage();
    }
    highlight(terms.length > 0);
  }

  hitsEl.addEventListener("click", function (event) {
    var head = event.target.closest(".hit-head");
    if (!head) return;
    var hit = head.parentElement;
    var body = hit.querySelector(".hit-body");
    if (!body.innerHTML) body.innerHTML = detailsHtml(Number(hit.dataset.row));
    hit.classList.toggle("open");
  });

  var debounce = null;
  searchInput.addEventListener("input", function () {
    clearTimeout(debounce);
    debounce = setTimeout(search, 150);
  });
  uniqueBox.addEventListener("change", search);
  moreBtn.addEventListener("click", renderPage);
  // Keeps the plot fitted while the sidebar is dragged wider or narrower
  new ResizeObserver(function () { Plotly.Plots.resize(gd); }).observe(document.getElementById("plot-area"));

  search();
})();
</script>
</body>
</html>
"""


def json_safe(value):
    """ Convert a single pandas/numpy cell into something json.dumps emits as valid JSON """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        # NaN/Infinity are not valid JSON, so they are dropped to null
        return float(value) if math.isfinite(float(value)) else None
    if pd.isna(value):
        return None
    return str(value)


def map_points_to_traces(fig, df, color_column):
    """ Map every dataframe row to the (trace, point) it became in the plotly figure

    plotly express splits the points into one trace per color category and keeps
    the dataframe order within each trace, which is what is reconstructed here.
    Returns the per-row trace indices, the per-row point indices within that
    trace and the point count of each trace.
    """
    point_trace = [0] * len(df)
    point_index = [0] * len(df)
    trace_counts = []
    categories = df[color_column].to_numpy()
    for trace_idx, trace in enumerate(fig.data):
        positions = np.flatnonzero(categories == trace.name)
        if len(positions) != len(trace.x):
            raise RuntimeError(
                "Cannot map rows to trace '{}': {} rows but {} points".format(trace.name, len(positions), len(trace.x))
            )
        for point_idx, row in enumerate(positions):
            point_trace[row] = trace_idx
            point_index[row] = point_idx
        trace_counts.append(len(positions))
    return point_trace, point_index, trace_counts


def encode_ms2_origins(values):
    """ Replace the file name each ms2 scan list is prefixed with by a table index

    The file names repeat once per file and per point, which is the bulk of the
    ms2s column. Returns the encoded cells and the table to resolve them with.
    """
    origins = {}
    encoded = []
    for value in values:
        if not value:
            encoded.append(value)
            continue
        parts = []
        for part in value.split(";"):
            origin, colon, scans = part.partition(":")
            if not colon:
                # Not "<file>:<scans>", so there is no file name to fold away
                parts.append(part)
                continue
            if origin not in origins:
                origins[origin] = str(len(origins))
            parts.append(origins[origin] + ":" + scans)
        encoded.append(";".join(parts))
    return encoded, list(origins)


def write_volcano_html(fig, df, detail_columns, output_path, title=None):
    """ Write the volcano figure as an HTML page with a prot_ident search sidebar

    The sidebar filters the points on their prot_ident, dims everything that
    does not match and lists the matches with their full point information.
    """
    point_trace, point_index, trace_counts = map_points_to_traces(fig, df, "color")

    values = [[json_safe(v) for v in df[column].tolist()] for column in detail_columns]
    ms2s_column = detail_columns.index("ms2s")
    values[ms2s_column], ms2_origins = encode_ms2_origins(values[ms2s_column])

    payload = dict(
        n_points=len(df),
        columns=detail_columns,
        detail_columns=detail_columns,
        values=values,
        prot_column="prot_ident",
        ms2s_column="ms2s",
        ms2_origins=ms2_origins,
        # One color per trace, and every trace is exactly one color category, so
        # the category of a point does not have to be stored per point
        trace_colors=[trace.marker.color for trace in fig.data],
        summary_columns=VOLCANO_SUMMARY_COLUMNS,
        summary_labels={
            "fold_change_A_div_B": "log2(FC)",
            "ttest_ind_pvalue": "p",
            "raw_pep_ident": "raw_pep_ident",
        },
        point_trace=point_trace,
        point_index=point_index,
        trace_counts=trace_counts,
    )

    plot_id = "volcano-plot"
    plot_div = fig.to_html(full_html=False, include_plotlyjs=True, div_id=plot_id, default_width="100%", default_height="100%")

    page = VOLCANO_PAGE_TEMPLATE
    page = page.replace("__TITLE__", title or os.path.splitext(os.path.basename(output_path))[0])
    page = page.replace("__PLOT_ID__", plot_id)
    page = page.replace("__PLOT_DIV__", plot_div)
    # "</" inside a <script> block would end the block prematurely
    page = page.replace("__DATA_JSON__", json.dumps(payload, allow_nan=False).replace("</", "<\\/"))

    with open(output_path, "w") as out:
        out.write(page)



if __name__ == "__main__":
    # Parse arguments
    args = parse_args()

    # TODO should they be in the same size ?
    # if args.groupA is None or args.groupB is None or \
    #     len(args.groupA) != len(args.groupB) or len(args.groupA) == 0:
    #     raise Exception("Groups need to be specified and need to be of same length")

    # Initialize csv reader/writer
    tsv_in = csv.reader(args.input_unbequant_tsv, delimiter="\t")
    # Get indices of the columns
    header = next(tsv_in)
    gA_idx = [header.index(x) for x in args.groupA]
    gB_idx = [header.index(x) for x in args.groupB]
    openmsid_idx = header.index("openms_ceid")
    pep_ident_idx = [header.index(x) for x in header if x.endswith("_____l_pep_ident")]
    raw_pep_ident_idx = [header.index(x) for x in header if x.endswith("_____l_raw_pep_ident")]
    prot_ident_idx = [header.index(x) for x in header if x.endswith("_____l_prot_ident")]
    charge_idx = [header.index(x) for x in header if x.endswith("_____charge")]
    ms2_scans_idx = [header.index(x) for x in header if x.endswith("_____l_ms2_scans")]
    ms2_scans_origin = [x[:-len("_____l_ms2_scans")] for x in header if x.endswith("_____l_ms2_scans")]
    first_iso_global_min_mz_idx = header.index("first_iso_global_min_mz")
    first_iso_global_max_mz_idx = header.index("first_iso_global_max_mz")
    first_iso_global_min_rt_idx = header.index("first_iso_global_min_rt")
    first_iso_global_max_rt_idx = header.index("first_iso_global_max_rt")
    feature_global_min_mz_idx = header.index("feature_global_min_mz")
    feature_global_max_mz_idx = header.index("feature_global_max_mz")
    feature_global_min_rt_idx = header.index("feature_global_min_rt")
    feature_global_max_rt_idx = header.index("feature_global_max_rt")
    # Get additional headers, which were addded after the UnbeQuant analysis
    non_unbequant_headers = set([header.index(x) for x in header if "_____" not in x]) - set([openmsid_idx, first_iso_global_min_mz_idx, first_iso_global_max_mz_idx, first_iso_global_min_rt_idx, first_iso_global_max_rt_idx, feature_global_min_mz_idx, feature_global_max_mz_idx, feature_global_min_rt_idx, feature_global_max_rt_idx])

    # Iterate over each line and generate pd Frame
    dict_for_pd = {**dict(
        openmsid=[],
        ttest_ind_pvalue=[],
        ttest_ind_statistic=[],
        fold_change_A_div_B=[],
        missing_values_in_a=[],
        missing_values_in_b=[],
        charge=[],
        first_iso_global_min_mz=[],
        first_iso_global_max_mz=[],
        first_iso_global_min_rt=[],
        first_iso_global_max_rt=[],
        feature_global_min_mz=[],
        feature_global_max_mz=[],
        feature_global_min_rt=[],
        feature_global_max_rt=[],
        pep_ident=[],
        raw_pep_ident=[],
        prot_ident=[],
        ms2s=[],
    ), **{header[x]: [] for x in non_unbequant_headers}}
    additional_headers = [header[x] for x in non_unbequant_headers]
    for line in tqdm.tqdm(tsv_in):

        # Get Cell Contents
        ga = [parse_intensity(line[x]) for x in gA_idx]
        gb = [parse_intensity(line[x]) for x in gB_idx]

        # Check if we have enough values
        if ga.count(None)/len(ga) > args.percentage_of_missingness_per_group_allowed:
            # Too many missing values in group A
            continue
        if gb.count(None)/len(gb) > args.percentage_of_missingness_per_group_allowed:
            # Too many missing values in group B
            continue

        missing_ga = ga.count(None)/len(ga)
        missing_gb = gb.count(None)/len(gb)

        ga = [x for x in ga if x is not None]
        gb = [x for x in gb if x is not None]

        try:
            # Write pvalue and statistic
            res = ttest_ind(ga, gb)
            if pd.isna(res.pvalue):
                continue  # Skip if pvalue is nan. This case happens if a group  has 100% missing values
            # Add to the dictionary missing information and other static stuff
            dict_for_pd["openmsid"].append(line[openmsid_idx])
            dict_for_pd["pep_ident"].append(";".join(parse_pep_prot_ident([line[x] for x in pep_ident_idx])))
            dict_for_pd["raw_pep_ident"].append(";".join(parse_pep_prot_ident([line[x] for x in raw_pep_ident_idx])))
            dict_for_pd["prot_ident"].append(";".join(parse_pep_prot_ident([line[x] for x in prot_ident_idx])))
            dict_for_pd["charge"].append((set([parse_charge(line[x]) for x in charge_idx]) - {None}).pop())
            dict_for_pd["ms2s"].append(";".join(parse_ms2_scans([line[x] for x in ms2_scans_idx],  ms2_scans_origin)))
            dict_for_pd["ttest_ind_pvalue"].append(res.pvalue)
            dict_for_pd["ttest_ind_statistic"].append(res.statistic)
            if args.log2_intensities_provided:
                dict_for_pd["fold_change_A_div_B"].append(np.log2(np.mean([2**x for x in ga]) / np.mean([2**x for x in gb])))
            else:
                dict_for_pd["fold_change_A_div_B"].append(np.log2(np.mean(ga) / np.mean(gb)))
            dict_for_pd["missing_values_in_a"].append(missing_ga)
            dict_for_pd["missing_values_in_b"].append(missing_gb)
            dict_for_pd["first_iso_global_min_mz"].append(float(line[first_iso_global_min_mz_idx]))
            dict_for_pd["first_iso_global_max_mz"].append(float(line[first_iso_global_max_mz_idx]))
            dict_for_pd["first_iso_global_min_rt"].append(float(line[first_iso_global_min_rt_idx]))
            dict_for_pd["first_iso_global_max_rt"].append(float(line[first_iso_global_max_rt_idx]))
            dict_for_pd["feature_global_min_mz"].append(float(line[feature_global_min_mz_idx]))
            dict_for_pd["feature_global_max_mz"].append(float(line[feature_global_max_mz_idx]))
            dict_for_pd["feature_global_min_rt"].append(float(line[feature_global_min_rt_idx]))
            dict_for_pd["feature_global_max_rt"].append(float(line[feature_global_max_rt_idx]))
            for h in non_unbequant_headers:
                dict_for_pd[header[h]].append(line[h])
        except Exception as e:
            # Could not calculate write nones
            print("Could not calculate ttest_ind for rowkey: ", line[openmsid_idx], "Error: ", e)


    # Generate pandas dataframe
    df = pd.DataFrame(dict_for_pd)

    # Sort by ascending p value:
    df = df.sort_values(by=["ttest_ind_pvalue"], ascending=True)

    # Do p-value correction:
    df["ttest_ind_corrected_pvalue"] = false_discovery_control(df["ttest_ind_pvalue"], method="bh")

    # Calculate the p-value for volcano
    df["-log10(pvalue)"] = -np.log10(df["ttest_ind_pvalue"])

    # Save the dataframe
    df.to_csv(args.output, sep="\t", index=False)

    # Save the volcano plot
    if args.output_volcano:
        df["pep_ident_hovers"] = df["pep_ident"].apply(lambda x: x[:30] + (x[30:] and ".."))
        df["prot_ident_hovers"] = df["prot_ident"].apply(lambda x: x[:30] + (x[30:] and ".."))
        df["raw_pep_ident_hovers"] = df["raw_pep_ident"].apply(lambda x: x[:30] + (x[30:] and ".."))
        df["ms2s_hovers"] = df["ms2s"].apply(lambda x: x[:30] + (x[30:] and ".."))
        for h in additional_headers:
            df[h + "_hovers"] = df[h].apply(lambda x: str(x)[:30] + (str(x)[30:] and ".."))
        df["is_identified"] = df["pep_ident"].apply(lambda x: 2 if x != "" else 1)
        df["is_significat"] = (df["ttest_ind_corrected_pvalue"] < 0.05) * 5 + 5
        df["is_ident_significat"] = df["is_identified"] + df["is_significat"]
        df["color"] = df["is_ident_significat"].apply(coloring_volcano)

        import plotly.express as px
        fig = px.scatter(df, x="fold_change_A_div_B", y="-log10(pvalue)", color="color", color_discrete_map=VOLCANO_COLOR_MAP, hover_data={**{
            "color" : False,
            "openmsid" : True,
            "ttest_ind_pvalue" : True,
            "ttest_ind_statistic" : True,
            "fold_change_A_div_B" : True,
            "missing_values_in_a" : True,
            "missing_values_in_b" : True,
            "charge" : True,
            "first_iso_global_min_mz" : True,
            "first_iso_global_max_mz" : True,
            "first_iso_global_min_rt" : True,
            "first_iso_global_max_rt" : True,
            "feature_global_min_mz" : True,
            "feature_global_max_mz" : True,
            "feature_global_min_rt" : True,
            "feature_global_max_rt" : True,
            "pep_ident_hovers" : True,
            "raw_pep_ident_hovers" : True,
            "prot_ident_hovers" : True,
            "ms2s_hovers" : True,
        } , **{h + "_hovers": True for h in additional_headers}})
        fig.add_hline(y=-np.log10(0.05), line_width=1, line_dash='dash')
        fig.add_vline(x=np.log2(2), line_width=1, line_dash='dash')
        fig.add_vline(x=np.log2(0.5), line_width=1, line_dash='dash')

        # Write volcano, together with the prot_ident search sidebar
        write_volcano_html(fig, df, VOLCANO_DETAIL_COLUMNS + additional_headers, args.output_volcano)
