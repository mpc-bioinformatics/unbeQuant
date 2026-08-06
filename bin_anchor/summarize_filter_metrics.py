#!/usr/bin/env python3
"""
Component metrics for one finished clusterless run.

Reports the numbers the filter-mode grid exists to compare:

  LINKAGE      mean vertices per component. Wanted LARGE - it is how many files a
               feature was found in. A drop below ~7 on this dataset means the
               filter is discarding real matches.

  CONFLICTS    components holding two or more different peptides. Wanted LOW. The
               majority peptide (by vote over identified members) is taken as the
               component's true analyte; 'minority' counts the identified members
               that disagree with it.

  CONFORMANCE  components respecting the pairing cutoff. Connected components do
               NOT guarantee this: membership is transitive while the cutoff is
               pairwise, so A-B and B-C both inside the cutoff still permits A-C
               outside it, unboundedly. Two properties are separated:
                 PAIRWISE - every pair of members is inside the cutoff
                 STAR     - SOME member is inside the cutoff of all the others.
                            This is what an anchor filter in 'full' radius mode
                            guarantees, since the anchor is exactly such a member.
                            'full' therefore bounds drift at 2x the cutoff, not 1x.

RT caveat: y_center is exported PRE-alignment (the loess corrections are applied to
copies inside pairing), so RT spans carry the alignment offset - measured at up to
~30s per edge. m/z is exact.
"""
import argparse
import json
import os
import statistics
from collections import Counter


def measure(path, ppm, rt):
    with open(path) as f:
        data = json.load(f)
    comps = data['components']
    meta = data.get('metadata', {})

    sizes = []
    multi = []
    labelled = conflicting = labelled_verts = minority_verts = 0
    n_multi = pair_ok = star_ok = beyond_2x = 0
    worst_ppm = worst_rt = 0.0

    for comp in comps.values():
        verts = comp['vertices']
        sizes.append(len(verts))
        if len(verts) >= 2:
            multi.append(len(verts))
            n_multi += 1
            mz = [v['x_center'] for v in verts]
            rts = [v['y_center'] for v in verts]

            # Both axes are 1-D, so the worst pair is simply the extremes.
            span_ppm = (max(mz) - min(mz)) / min(mz) * 1e6
            span_rt = max(rts) - min(rts)
            worst_ppm = max(worst_ppm, span_ppm)
            worst_rt = max(worst_rt, span_rt)
            if span_ppm <= ppm and span_rt <= rt:
                pair_ok += 1
            if span_ppm > 2 * ppm or span_rt > 2 * rt:
                beyond_2x += 1

            # Star: tolerance taken from the candidate centre's own m/z, as pairing
            # does it. Short-circuits on the first centre that covers everything.
            for i in range(len(verts)):
                tol = mz[i] * (ppm / 1e6)
                if all(abs(mz[i] - mz[j]) <= tol and abs(rts[i] - rts[j]) <= rt
                       for j in range(len(verts))):
                    star_ok += 1
                    break

        idents = []
        for v in verts:
            pep = v.get('pep_ident') or []
            if isinstance(pep, str):
                pep = [pep]
            s = frozenset(x for x in pep if x)
            if s:
                idents.append(s)
        if not idents:
            continue
        labelled += 1
        labelled_verts += len(idents)
        votes = Counter()
        for s in idents:
            for pep in s:
                votes[pep] += 1
        if len(votes) > 1:
            conflicting += 1
            majority = votes.most_common(1)[0][0]
            minority_verts += sum(1 for s in idents if majority not in s)

    pct = lambda a, b: (100.0 * a / b) if b else 0.0  # noqa: E731
    return {
        'cycles': meta.get('total_cycles', ''),
        'components': len(comps),
        'vertices': sum(sizes),
        'avg_linkage': statistics.mean(sizes) if sizes else 0.0,
        'avg_linkage_ge2': statistics.mean(multi) if multi else 0.0,
        'median_size': statistics.median(sizes) if sizes else 0,
        'conflict_pct': pct(conflicting, labelled),
        'minority_pct': pct(minority_verts, labelled_verts),
        'pairwise_pct': pct(pair_ok, n_multi),
        'star_pct': pct(star_ok, n_multi),
        'worst_ppm_x': worst_ppm / ppm if ppm else 0.0,
        'worst_rt_x': worst_rt / rt if rt else 0.0,
        'beyond_2x': beyond_2x,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--components', required=True,
                   help='final_components_all_cycles.json from the run')
    p.add_argument('--ppm', type=float, required=True, help='m/z cutoff to test against')
    p.add_argument('--rt', type=float, required=True, help='RT cutoff to test against')
    p.add_argument('--run', default='', help='run name, for the TSV row')
    p.add_argument('--filter_method', default='', help='for the TSV row')
    p.add_argument('--radius', default='', help='for the TSV row')
    p.add_argument('--anchor_rt', default='', help='for the TSV row')
    p.add_argument('--tsv_append', help='append one row to this TSV (header written if new)')
    args = p.parse_args()

    m = measure(args.components, args.ppm, args.rt)

    print(f"  components      : {m['components']}  ({m['vertices']} vertices)")
    print(f"  avg linkage     : {m['avg_linkage']:.3f}   (>=2 members: {m['avg_linkage_ge2']:.3f})")
    print(f"  pep conflicts   : {m['conflict_pct']:.1f}%   minority vertices: {m['minority_pct']:.2f}%")
    print(f"  conformance     : pairwise {m['pairwise_pct']:.2f}%   star {m['star_pct']:.2f}%")
    print(f"  worst drift     : {m['worst_ppm_x']:.2f}x m/z, {m['worst_rt_x']:.2f}x RT   "
          f"(beyond 2x: {m['beyond_2x']})")

    if args.tsv_append:
        new = not os.path.exists(args.tsv_append) or os.path.getsize(args.tsv_append) == 0
        with open(args.tsv_append, 'a') as f:
            if new:
                f.write('run\tfilter_method\tradius\tanchor_rt\tcycles\tcomponents\t'
                        'avg_linkage\tavg_linkage_ge2\tconflict_pct\tminority_pct\t'
                        'pairwise_pct\tstar_pct\tworst_ppm_x\tworst_rt_x\tbeyond_2x\n')
            f.write(f"{args.run}\t{args.filter_method}\t{args.radius}\t{args.anchor_rt}\t"
                    f"{m['cycles']}\t{m['components']}\t{m['avg_linkage']:.3f}\t"
                    f"{m['avg_linkage_ge2']:.3f}\t{m['conflict_pct']:.2f}\t"
                    f"{m['minority_pct']:.2f}\t{m['pairwise_pct']:.2f}\t{m['star_pct']:.2f}\t"
                    f"{m['worst_ppm_x']:.2f}\t{m['worst_rt_x']:.2f}\t{m['beyond_2x']}\n")


if __name__ == '__main__':
    main()
