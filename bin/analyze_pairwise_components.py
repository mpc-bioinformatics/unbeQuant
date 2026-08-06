#!/usr/bin/env python3
"""
Pairwise within-component distance analysis for one finished clusterless run.

For every component, this measures EVERY pair of its members - a component with
members a,b,c,d contributes ab, ac, ad, bc, bd, cd - and reports how those
distances sit against the pairing cutoffs.

WHY EVERY PAIR AND NOT JUST THE EDGES
  The cutoffs are pairwise, but components are connected components, so
  membership is transitive: a-b and b-c can each pass while a-c is far outside.
  Nothing downstream re-checks. Plotting all pairs makes that drift visible and
  quantifies it, which the edge list alone cannot.

  Size-2 components are the control case: their single pair IS the tested edge,
  so a violation there means the cutoff was evaluated against different values
  than the ones stored here - not that the component drifted.

RT caveat: y_center is exported PRE-alignment (the loess corrections are applied
to copies inside pairing), so RT distances carry the alignment offset - measured
at up to ~30 s per edge. This is why size-2 components can appear to break the RT
cutoff by a few tens of seconds while never breaking the ppm one. m/z is exact.

OUTPUTS (into --outdir)
  <prefix>_mz_<theme>.png          |Delta m/z| histogram
  <prefix>_rt_<theme>.png          |Delta RT| histogram
  <prefix>_ppm_cutoff_<theme>.png  ppm pass/fail, split by component size
  <prefix>_rt_cutoff_<theme>.png   RT pass/fail, split by component size
  <prefix>_summary.json            the same numbers as data
"""
import argparse
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import LogLocator, NullFormatter  # noqa: E402


# Palette: one hue per metric, plus a shared within/above pair for the cutoff
# figures. Both themes are stepped for their own surface, not flipped.
THEMES = {
    "light": dict(
        surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", ink3="#8a8983",
        grid="#e3e2dd", mz="#2a78d6", rt="#1baf7a", rule="#e34948",
        within="#2a78d6", exceed="#eb6834",
    ),
    "dark": dict(
        surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", ink3="#8a8983",
        grid="#3a3a38", mz="#3987e5", rt="#199e70", rule="#e66767",
        within="#3987e5", exceed="#d95926",
    ),
}


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_metadata(path):
    """Parse the leading metadata block without loading the whole file.

    The components object that follows is hundreds of MB; only the header is
    needed to learn the cutoffs the run was configured with.
    """
    head = []
    with open(path) as fh:
        for line in fh:
            if '"components"' in line:
                break
            head.append(line)
    return json.loads("".join(head).rstrip().rstrip(",") + "}").get("metadata", {})


def stream_pairs(path, max_n=64, progress_every=20000):
    """Yield all within-component pairwise distances as flat arrays.

    Scanned line-wise rather than json.load()ed: the file runs to ~700 MB, and
    a full parse costs several GB of RAM for data that reduces to four arrays.
    This relies on the pretty-printed layout orchestrate_clusterless_cycles.py
    writes (json.dump(..., indent=2)), so each field sits on its own line.

    Returns (dmz, drt, ppm, size) - one entry per pair.
    """
    triu = {n: np.triu_indices(n, k=1) for n in range(2, max_n + 1)}
    dmz_parts, drt_parts, mzmean_parts, size_parts = [], [], [], []
    mz, rt = [], []
    n_components = 0

    def flush():
        n = len(mz)
        if n < 2:
            return
        if n not in triu:
            triu[n] = np.triu_indices(n, k=1)
        i, j = triu[n]
        a, b = np.asarray(mz), np.asarray(rt)
        dmz_parts.append(np.abs(a[i] - a[j]))
        drt_parts.append(np.abs(b[i] - b[j]))
        mzmean_parts.append(0.5 * (a[i] + a[j]))
        size_parts.append(np.full(i.size, n, dtype=np.int16))

    with open(path) as fh:
        for line in fh:
            if '"x_center"' in line:
                mz.append(float(line.split(":", 1)[1].strip().rstrip(",")))
            elif '"y_center"' in line:
                rt.append(float(line.split(":", 1)[1].strip().rstrip(",")))
            elif '"component_id"' in line:
                flush()
                mz, rt = [], []
                n_components += 1
                if progress_every and n_components % progress_every == 0:
                    print(f"  ... {n_components:,} components", flush=True)
    flush()

    if n_components == 0:
        raise SystemExit(
            f"ERROR: no components found in {path}. This parser expects the "
            "pretty-printed layout written by orchestrate_clusterless_cycles.py."
        )
    if not dmz_parts:
        raise SystemExit(
            f"ERROR: {n_components:,} components but none with >=2 members - "
            "nothing to pair."
        )

    dmz = np.concatenate(dmz_parts)
    drt = np.concatenate(drt_parts)
    ppm = dmz / np.concatenate(mzmean_parts) * 1e6
    return dmz, drt, ppm, np.concatenate(size_parts), n_components


# ---------------------------------------------------------------------------
# Shared chart furniture
# ---------------------------------------------------------------------------

def style(ax, t, xlabel, ylabel, title):
    ax.set_title(title, color=t["ink"], fontsize=12.5, fontweight="600",
                 loc="left", pad=12)
    ax.set_xlabel(xlabel, color=t["ink2"], fontsize=9.5, labelpad=8)
    ax.set_ylabel(ylabel, color=t["ink2"], fontsize=9.5, labelpad=8)
    ax.grid(axis="y", color=t["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["grid"])
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=t["ink2"], labelsize=9, length=3, width=0.8)


def footer(fig, t, main, label):
    fig.text(0.5, 0.028, main, color=t["ink3"], fontsize=8.5, ha="center")
    if label:
        fig.text(0.5, 0.008, label, color=t["ink3"], fontsize=8.5, ha="center")


# ---------------------------------------------------------------------------
# Figure 1/2: distance histograms
# ---------------------------------------------------------------------------

def plot_distribution(values, xmax, kind, t, out, n_components, label,
                      rt_cutoff=None, nbins=120):
    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    fig.patch.set_facecolor(t["surface"])
    ax.set_facecolor(t["surface"])

    inside = values[values <= xmax]
    n_over = values.size - inside.size
    counts, edges = np.histogram(inside, bins=nbins, range=(0.0, xmax))
    step = edges[1] - edges[0]
    color = t["mz"] if kind == "mz" else t["rt"]
    ax.bar(edges[:-1] + step / 2, counts, width=step * 0.88, color=color,
           linewidth=0, zorder=3)

    ax.set_yscale("log")
    ax.set_xlim(0, xmax)
    # Anchor the log floor just under the smallest bin: a decade-0 floor would
    # render every bar as one solid block and hide the shape.
    nz = counts[counts > 0]
    ax.set_ylim(bottom=10 ** np.floor(np.log10(nz.min())) if nz.size else 1,
                top=max(counts.max() * 2.2, 10))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_formatter(NullFormatter())

    if kind == "mz":
        style(ax, t,
              "|Δ m/z|  between two features in the same component  (m/z units, Th)",
              "pair count  (log scale)",
              "Pairwise m/z distance within components")
    else:
        style(ax, t,
              "|Δ RT|  between two features in the same component  (s)",
              "pair count  (log scale)",
              "Pairwise RT distance within components")
        if rt_cutoff is not None:
            ax.axvline(rt_cutoff, color=t["rule"], linewidth=1.2, zorder=5)
            ax.text(rt_cutoff + xmax * 0.02, 0.62, f"rt_cutoff {rt_cutoff:g} s",
                    transform=ax.get_xaxis_transform(), color=t["rule"],
                    fontsize=8.5, va="center", zorder=6,
                    bbox=dict(facecolor=t["surface"], edgecolor="none",
                              boxstyle="round,pad=0.35"))

    stats = (f"median {np.median(values):,.4g}    "
             f"p99 {np.percentile(values, 99):,.4g}    max {values.max():,.4g}\n"
             f"{n_over:,} pairs ({n_over / values.size * 100:.2f}%) beyond axis")
    ax.text(0.98, 0.95, stats, transform=ax.transAxes, ha="right", va="top",
            color=t["ink2"], fontsize=8.5, linespacing=1.6)

    fig.tight_layout(rect=[0, 0.075, 1, 1])
    footer(fig, t, f"All within-component feature pairs  ·  {values.size:,} "
                   f"pairs across {n_components:,} components", label)
    fig.savefig(out, dpi=150, facecolor=t["surface"])
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3/4: cutoff conformance by component size
# ---------------------------------------------------------------------------

def plot_cutoff(values, size, cutoff, unit, formula, t, out, n_components,
                label):
    sizes = np.arange(size.min(), size.max() + 1)
    bad = values > cutoff
    n_pairs = np.array([(size == s).sum() for s in sizes])
    n_exceed = np.array([bad[size == s].sum() for s in sizes])
    keep = n_pairs > 0
    sizes, n_pairs, n_exceed = sizes[keep], n_pairs[keep], n_exceed[keep]
    n_within = n_pairs - n_exceed
    clean = n_exceed.sum() == 0
    tag = f"{cutoff:g} {unit}"

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.8))
    fig.patch.set_facecolor(t["surface"])
    for ax in axes:
        ax.set_facecolor(t["surface"])
    x = np.arange(len(sizes))
    pct = n_exceed / n_pairs * 100

    # --- Panel A: violation rate per component size -----------------------
    ax = axes[0]
    ax.bar(x, pct, width=0.72, color=t["exceed"], linewidth=0, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(sizes, fontsize=8)
    style(ax, t, "component size  (vertices)",
          f"% of the component's pairs above {tag}",
          f"Share of pairs breaking the {tag} cutoff")
    ax.title.set_fontsize(12)
    if clean:
        # An all-zero series on an auto-scaled axis reads as a broken plot.
        ax.set_ylim(0, 1)
        ax.text(0.5, 0.5,
                f"No violations at any component size.\n\nEvery one of the "
                f"{values.size:,} pairs is within {tag} —\nincluding transitive "
                f"pairs, so the constraint\nholds across whole components, not "
                f"just edges.",
                transform=ax.transAxes, ha="center", va="center",
                color=t["ink2"], fontsize=10, linespacing=1.7)
    else:
        ax.set_ylim(0, pct.max() * 1.42)
        # Direct-label only the extremes; the axis carries the rest.
        for i in dict.fromkeys((0, int(np.argmax(pct)), len(sizes) - 1)):
            lab = "0" if pct[i] == 0 else f"{pct[i]:.2f}%"
            ax.text(x[i], pct[i] + pct.max() * 0.04, lab, ha="center",
                    va="bottom", color=t["ink2"], fontsize=8.5,
                    fontweight="600")
        # Size-2 is the control case - report what it actually does rather than
        # asserting it must be clean (for RT it often is not; see module docstring).
        n2 = int(n_exceed[0]) if sizes[0] == 2 else 0
        if sizes[0] == 2 and n2 == 0:
            note = ("size-2 components cannot violate:\n"
                    "their only pair is the edge itself")
        elif sizes[0] == 2:
            note = (f"{n2:,} size-2 pairs break the cutoff\n"
                    f"(max {values[size == 2].max():.1f} {unit}), yet their only\n"
                    f"pair IS the tested edge - so the test\n"
                    f"did not run on these values")
        else:
            note = ""
        if note:
            # Top-right: the left shoulder holds the small-size bars and labels.
            ax.text(0.98, 0.96, note, transform=ax.transAxes, ha="right",
                    va="top", color=t["ink3"], fontsize=8.5, linespacing=1.5)

    # --- Panel B: absolute pair counts ------------------------------------
    ax = axes[1]
    w = 0.4
    ax.bar(x - w / 2, n_within, width=w * 0.92, color=t["within"], linewidth=0,
           zorder=3, label=f"within {tag}")
    ax.bar(x + w / 2, n_exceed, width=w * 0.92, color=t["exceed"], linewidth=0,
           zorder=3, label=f"above {tag}")
    ax.set_yscale("log")
    ax.set_ylim(bottom=1, top=n_within.max() * 3)
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(x)
    ax.set_xticklabels(sizes, fontsize=8)
    style(ax, t, "component size  (vertices)", "pair count  (log scale)",
          "Pairs per component size")
    ax.title.set_fontsize(12)
    # Zeros cannot be drawn on a log axis - label them rather than show a gap.
    if clean:
        ax.text(0.035, 0.79, f"(no pairs above {tag} at any size)",
                transform=ax.transAxes, ha="left", va="top", color=t["ink2"],
                fontsize=9)
    else:
        for i in np.flatnonzero(n_exceed == 0):
            ax.text(x[i] + w / 2, 1.35, "0", ha="center", va="bottom",
                    color=t["ink2"], fontsize=8, fontweight="600")
    leg = ax.legend(frameon=False, fontsize=9, loc="upper left",
                    handlelength=1.1, handleheight=1.1)
    for txt in leg.get_texts():
        txt.set_color(t["ink2"])

    n_bad, n_all = int(n_exceed.sum()), int(n_pairs.sum())
    fig.suptitle(
        f"{n_bad:,} of {n_all:,} within-component pairs "
        f"({n_bad / n_all * 100:.2f}%) exceed the {tag} cutoff  ·  "
        f"{n_all - n_bad:,} ({(n_all - n_bad) / n_all * 100:.2f}%) stay within",
        color=t["ink"], fontsize=11.5, fontweight="600", y=0.985)
    # Only claim the edge/diameter gap when the data actually shows one.
    caveat = ("here it holds for whole components, not only edges" if clean
              else "the cutoff constrains graph edges, not component diameters")
    fig.tight_layout(rect=[0, 0.085, 1, 0.945])
    fig.text(0.5, 0.042,
             f"All pairs of features inside the same component, across "
             f"{n_components:,} components  ·  {formula}  ·  {caveat}",
             color=t["ink3"], fontsize=8.5, ha="center")
    if label:
        fig.text(0.5, 0.014, label, color=t["ink3"], fontsize=8.5, ha="center")
    fig.savefig(out, dpi=150, facecolor=t["surface"])
    plt.close(fig)


# ---------------------------------------------------------------------------

def summarize(dmz, drt, ppm, size, ppm_cutoff, rt_cutoff, n_components):
    ppm_bad, rt_bad = ppm > ppm_cutoff, drt > rt_cutoff
    m2 = size == 2
    q = lambda v, p: float(np.percentile(v, p))  # noqa: E731
    return {
        "components": int(n_components),
        "pairs": int(dmz.size),
        "ppm_cutoff": float(ppm_cutoff),
        "rt_cutoff": float(rt_cutoff),
        "ppm_exceed": int(ppm_bad.sum()),
        "ppm_exceed_pct": float(ppm_bad.mean() * 100),
        "rt_exceed": int(rt_bad.sum()),
        "rt_exceed_pct": float(rt_bad.mean() * 100),
        "mz_median": q(dmz, 50), "mz_p99": q(dmz, 99), "mz_max": float(dmz.max()),
        "rt_median": q(drt, 50), "rt_p99": q(drt, 99), "rt_max": float(drt.max()),
        "ppm_median": q(ppm, 50), "ppm_p99": q(ppm, 99), "ppm_max": float(ppm.max()),
        # Size-2 is the control case: its pair IS the tested edge.
        "size2_pairs": int(m2.sum()),
        "size2_ppm_exceed": int(ppm_bad[m2].sum()),
        "size2_rt_exceed": int(rt_bad[m2].sum()),
    }


TSV_COLUMNS = ["run", "components", "pairs", "ppm_cutoff", "rt_cutoff",
               "ppm_exceed", "ppm_exceed_pct", "rt_exceed", "rt_exceed_pct",
               "mz_median", "mz_p99", "mz_max", "rt_median", "rt_p99", "rt_max",
               "size2_ppm_exceed", "size2_rt_exceed"]


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--components_json", required=True,
                   help="final_components_all_cycles.json from the run")
    p.add_argument("--outdir", default=".", help="where the figures are written")
    p.add_argument("--prefix", default="pairwise", help="output filename prefix")
    p.add_argument("--ppm", type=float,
                   help="ppm cutoff to test against (default: from run metadata)")
    p.add_argument("--rt", type=float,
                   help="RT cutoff in seconds (default: from run metadata)")
    p.add_argument("--mz_xmax", type=float, default=0.012,
                   help="m/z histogram x limit, Th (default 0.012, ~p99.9)")
    p.add_argument("--rt_xmax", type=float, default=300.0,
                   help="RT histogram x limit, seconds (default 300)")
    p.add_argument("--theme", default="both", choices=["light", "dark", "both"])
    p.add_argument("--run", default="", help="run name, for the TSV row and captions")
    p.add_argument("--tsv_append",
                   help="append one row to this TSV (header written if new)")
    args = p.parse_args()

    meta = read_metadata(args.components_json)
    cut = meta.get("cutoff_params", {})
    # Cutoffs differ between runs, so they are read from the run's own metadata
    # unless explicitly overridden - never hardcoded.
    ppm_cutoff = args.ppm if args.ppm is not None else cut.get("ppm")
    rt_cutoff = args.rt if args.rt is not None else cut.get("rt_cutoff")
    if ppm_cutoff is None or rt_cutoff is None:
        raise SystemExit("ERROR: cutoffs not found in metadata; pass --ppm and --rt")

    print(f"Reading {args.components_json}")
    print(f"  cutoffs: {ppm_cutoff:g} ppm, {rt_cutoff:g} s")
    dmz, drt, ppm, size, n_components = stream_pairs(args.components_json)
    print(f"  {n_components:,} components -> {dmz.size:,} pairs")

    os.makedirs(args.outdir, exist_ok=True)
    themes = ["light", "dark"] if args.theme == "both" else [args.theme]
    label = args.run or os.path.basename(os.path.dirname(
        os.path.abspath(args.components_json)))

    for name in themes:
        t = THEMES[name]
        stem = os.path.join(args.outdir, args.prefix)
        plot_distribution(dmz, args.mz_xmax, "mz", t, f"{stem}_mz_{name}.png",
                          n_components, label)
        plot_distribution(drt, args.rt_xmax, "rt", t, f"{stem}_rt_{name}.png",
                          n_components, label, rt_cutoff=rt_cutoff)
        plot_cutoff(ppm, size, ppm_cutoff, "ppm",
                    "ppm = |Δ m/z| / mean(m/z) × 1e6", t,
                    f"{stem}_ppm_cutoff_{name}.png", n_components, label)
        plot_cutoff(drt, size, rt_cutoff, "s",
                    "RT distance = |Δ y_center| in seconds", t,
                    f"{stem}_rt_cutoff_{name}.png", n_components, label)
        print(f"  wrote 4 figures ({name})")

    m = summarize(dmz, drt, ppm, size, ppm_cutoff, rt_cutoff, n_components)
    with open(os.path.join(args.outdir, f"{args.prefix}_summary.json"), "w") as fh:
        json.dump(m, fh, indent=2)

    print(f"  ppm > {ppm_cutoff:g}  : {m['ppm_exceed']:,} "
          f"({m['ppm_exceed_pct']:.4f}%)   size-2: {m['size2_ppm_exceed']:,}")
    print(f"  rt  > {rt_cutoff:g}s : {m['rt_exceed']:,} "
          f"({m['rt_exceed_pct']:.4f}%)   size-2: {m['size2_rt_exceed']:,}")

    if args.tsv_append:
        new = (not os.path.exists(args.tsv_append)
               or os.path.getsize(args.tsv_append) == 0)
        row = dict(m, run=args.run)
        with open(args.tsv_append, "a") as fh:
            if new:
                fh.write("\t".join(TSV_COLUMNS) + "\n")
            fh.write("\t".join(
                f"{row[c]:.4f}" if isinstance(row[c], float) else str(row[c])
                for c in TSV_COLUMNS) + "\n")


if __name__ == "__main__":
    main()
