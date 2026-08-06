#!/usr/bin/env python3
"""
Pair matching features across multiple files (clusterless workflow).

Given per-file feature lists, this script finds nearest-neighbour matches between
features of the SAME charge across files using KD-trees, and emits a network graph
of "edges" (matched feature pairs) that survive the configured coordinate cutoffs.

Pipeline (see main()):
  1. Load features            - cycle 1 from *_feature_data.pkl, cycle 2+ from a
                                deleted_vertices JSON produced by the previous cycle.
  2. Optional RT alignment    - AlignTree (trafoXML) rewrites RTs up front; the fitted
                                JSON models are instead applied per file pair inside step 3.
  3. Pairwise matching        - feature_pairing_optimized(): charge-split KD-trees,
                                nearest-neighbour query, cutoff filtering at edge-creation.
  4. Save results             - edges pkl (consumed by build_network_graph_clusterless.py),
                                plus ident_lookup, unpaired vertices, and a metadata json.

Only edges are produced here; connected-component grouping happens downstream in
build_network_graph_clusterless.py.
"""

import argparse
import base64
import gc
import json
import os
import pickle
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

# The JSON RT-correction models may be pickled instances of these classes. Import the
# real definitions so unpickling works; fall back to minimal compatible stand-ins if the
# helper module is unavailable (e.g. when this script is run in isolation).
sys.path.insert(0, os.path.dirname(__file__))
try:
    from compare_retention_time_alignment import LOESSInterpolator, PiecewisePolynomial
except ImportError:
    class LOESSInterpolator:
        """LOESS interpolator that can be pickled (fallback)."""
        def __init__(self, x_data, y_data):
            self.x_data = np.array(x_data)
            self.y_data = np.array(y_data)

        def __call__(self, x_val):
            if isinstance(x_val, (list, np.ndarray)):
                return np.array([self._eval_scalar(x) for x in np.atleast_1d(x_val)])
            return self._eval_scalar(x_val)

        def _eval_scalar(self, x_val):
            if x_val <= self.x_data[0]:
                return self.y_data[0]
            if x_val >= self.x_data[-1]:
                return self.y_data[-1]
            idx = np.searchsorted(self.x_data, x_val)
            x0, x1 = self.x_data[idx - 1], self.x_data[idx]
            y0, y1 = self.y_data[idx - 1], self.y_data[idx]
            return y0 + (y1 - y0) * (x_val - x0) / (x1 - x0)

    class PiecewisePolynomial:
        """Piecewise polynomial model that can be pickled (fallback)."""
        pass


# ============================================================================
# RT alignment via AlignTree (OpenMS trafoXML files)
# ============================================================================

def build_trafoxmls_fid_dict(trafoxmls_dir: str, input_filenames: set = None) -> Dict[str, float]:
    """
    Build a mapping feature-ID -> transformed RT from OpenMS trafoXML files.

    trafoXML files contain <Pair note="<node>" to="<transformed_rt>"/> elements. The
    node attribute is converted to openms_fid format by prepending "f_".

    Args:
        trafoxmls_dir: Directory containing .trafoXML files (searched recursively).
        input_filenames: If given, only trafoXML files whose name starts with one of
                         these basenames are loaded, so unrelated files can't pollute
                         the mapping.

    Returns:
        Dict mapping "f_<node>" -> transformed_rt (float).
    """
    dir_path = Path(trafoxmls_dir).resolve()
    if not dir_path.exists():
        raise FileNotFoundError(f"TrafoXML directory does not exist: {trafoxmls_dir}")
    if not dir_path.is_dir():
        raise NotADirectoryError(f"TrafoXML path is not a directory: {trafoxmls_dir}")

    # Prefer a recursive search; fall back to non-recursive before giving up.
    trafo_files = sorted(dir_path.glob("**/*.trafoXML"))
    if not trafo_files:
        trafo_files = sorted(dir_path.glob("*.trafoXML"))
    if not trafo_files:
        raise FileNotFoundError(f"No .trafoXML files found in: {trafoxmls_dir}")

    # Keep only the trafoXML files relevant to this run's input files.
    if input_filenames:
        relevant = [f for f in trafo_files
                    if any(f.name.startswith(fname) for fname in input_filenames)]
        skipped = len(trafo_files) - len(relevant)
        if skipped:
            print(f"  Filtered trafoXML files: {len(relevant)} relevant, {skipped} skipped")
        trafo_files = relevant

    print(f"[PROGRESS] Building feature ID -> transformed RT mappings from {len(trafo_files)} trafoXML files")

    fid_to_rt_dict = {}
    for trafo_file in trafo_files:
        try:
            root = ET.parse(str(trafo_file)).getroot()
            transformation = root.find('Transformation')
            if transformation is None:
                continue
            pairs_elem = transformation.find('Pairs')
            if pairs_elem is None:
                continue

            pairs = pairs_elem.findall('Pair')
            for pair in pairs:
                node_str = pair.get('note')  # Node ID from trafoXML
                to_rt_str = pair.get('to')   # Transformed RT
                if node_str and to_rt_str:
                    fid_to_rt_dict[f"f_{node_str}"] = float(to_rt_str)
            print(f"  ✓ {trafo_file.name}: {len(pairs)} feature transformations loaded")
        except Exception as e:
            print(f"  ✗ Warning: Failed to parse {trafo_file.name}: {e}")

    print(f"  ✓ Total feature ID mappings: {len(fid_to_rt_dict)}")
    return fid_to_rt_dict


def apply_aligntree_rt_transformations(all_feature_data_lists: List[List[Dict]],
                                       fid_to_rt_dict: Dict[str, float]) -> List[List[Dict]]:
    """
    Apply AlignTree RT transformations by matching each feature's openms_fid.

    For every feature whose openms_fid is present in fid_to_rt_dict, y_center is replaced
    with the transformed RT value. Features are copied (inputs untouched).

    Returns the transformed feature lists.
    """
    transformed_lists = []
    total_matched = 0
    total_features = 0

    for file_idx, feature_list in enumerate(all_feature_data_lists):
        transformed_features = []
        matched = 0
        for feature in feature_list:
            total_features += 1
            feature_copy = feature.copy()
            fid = feature.get('openms_fid')
            if fid and fid in fid_to_rt_dict:
                feature_copy['y_center'] = fid_to_rt_dict[fid]
                matched += 1
                total_matched += 1
            transformed_features.append(feature_copy)
        transformed_lists.append(transformed_features)
        print(f"  File {file_idx + 1}: {matched}/{len(feature_list)} features transformed (y_center updated)")

    print(f"\n[PROGRESS] Total features transformed: {total_matched}/{total_features}")
    return transformed_lists


# ============================================================================
# RT alignment via fitted JSON models (polynomial / spline / rbf / piecewise / loess)
# ============================================================================

def build_rt_correction_model(correction_dict: Dict) -> Optional[Callable]:
    """
    Turn a fitted correction record into a callable model(rt) -> rt_diff.

    Supports polynomial (coefficients) and pickle-serialized callables
    (rbf / piecewise / loess / spline). Every supported model accepts a scalar or an
    array, so a whole file's RTs can be shifted in one call, and the deserialization
    happens here rather than per feature - one unpickle per file pair.

    Returns None when the record carries no usable fit, in which case the pair is left
    unaligned. Raises ValueError on an unsupported method or a deserialization failure.
    """
    method = correction_dict.get('method', 'polynomial')
    serialization = correction_dict.get('serialization', {})

    if method == 'polynomial':
        coefficients = serialization.get('coefficients', [])
        return np.poly1d(coefficients) if coefficients else None

    if method in ('rbf', 'piecewise', 'loess', 'spline'):
        pickle_str = serialization.get('model_pickle')
        if not pickle_str:
            return None
        try:
            return pickle.loads(base64.b64decode(pickle_str))
        except Exception as e:
            raise ValueError(f"Failed to deserialize correction model (method={method}): {e}")

    raise ValueError(f"Unsupported fitting method: {method}")


def load_rt_corrections_json(json_path: str) -> Dict[Tuple[str, str], Dict]:
    """
    Load RT corrections from a JSON file (supports all fitting methods).

    Returns a dict keyed by (file1_name, file2_name) tuples (with any .mzML extension
    stripped) mapping to the per-pair correction dict.

    Raises FileNotFoundError / ValueError on missing, empty, or malformed input.
    """
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"RT correction JSON file not found: {json_path}")
    if os.path.getsize(json_path) == 0:
        raise ValueError(f"RT correction JSON file is empty: {json_path}")

    with open(json_path, 'r') as f:
        content = f.read().strip()
    if not content:
        raise ValueError(f"RT correction JSON file contains no data: {json_path}")

    data = json.loads(content)
    if not isinstance(data, dict) or 'corrections' not in data:
        raise ValueError("JSON must have 'corrections' key")

    corrections = data['corrections']
    if not corrections:
        print("  ⚠ RT correction JSON has empty 'corrections' dict — no corrections will be applied")
        return {}

    result = {}
    for pair_key, correction_data in corrections.items():
        # pair_key format: "file1_vs_file2"
        parts = pair_key.split('_vs_')
        if len(parts) != 2:
            raise ValueError(f"Invalid pair key format: {pair_key}, expected 'file1_vs_file2'")
        file1 = parts[0].replace('.mzML', '')
        file2 = parts[1].replace('.mzML', '')
        result[(file1, file2)] = correction_data
    return result


def extract_and_save_unpaired_vertices(network_edges: Dict,
                                       all_feature_data_lists: List[List[Dict]],
                                       output_path: str = None) -> List[Dict]:
    """
    Find features that appear in no edge (unpaired) and optionally save them to JSON.

    A vertex is "paired" if its (file_idx, feature_idx) appears as an endpoint of some
    edge in network_edges; every other feature is returned as unpaired. If output_path is
    given and any unpaired vertices exist, they are written there.
    """
    paired_vertex_ids = set()
    for edge_list in network_edges.values():
        for edge in edge_list:
            paired_vertex_ids.add((edge['ref_file_idx'], edge['current_feature_idx']))
            paired_vertex_ids.add((edge['match_file_idx'], edge['matched_feature_idx']))

    unpaired_vertices = []
    for file_idx, features in enumerate(all_feature_data_lists):
        for feat_idx, feature in enumerate(features):
            if (file_idx, feat_idx) in paired_vertex_ids:
                continue
            unpaired_vertices.append({
                'vertex_id': f"{feature.get('filename', f'file_{file_idx}')}_{feat_idx}",
                'filename': feature.get('filename', f'file_{file_idx}'),
                'feature_idx': feat_idx,
                'file_idx': file_idx,
                'x_center': float(feature.get('x_center', 0)),
                'y_center': float(feature.get('y_center', 0)),
                'mz_start': float(feature.get('mz_start', 0)),
                'mz_end': float(feature.get('mz_end', 0)),
                'rt_start': float(feature.get('rt_start', 0)),
                'rt_end': float(feature.get('rt_end', 0)),
                'charge': int(feature.get('charge', 0)),
                'intensity': float(feature.get('intensity', 0)),
                'pep_ident': feature.get('pep_ident', []),
                'prot_ident': feature.get('prot_ident', []),
                'openms_fid': feature.get('openms_fid', ''),
                'ms2_scans': feature.get('ms2_scans', '')
            })

    if output_path and unpaired_vertices:
        try:
            with open(output_path, 'w') as f:
                json.dump({'total_unpaired': len(unpaired_vertices),
                           'unpaired_vertices': unpaired_vertices}, f, indent=2)
            print(f"\n[PROGRESS] Saved {len(unpaired_vertices)} unpaired vertices to: {os.path.basename(output_path)}")
        except Exception as e:
            print(f"✗ Warning: Could not save unpaired vertices to {output_path}: {e}")
    elif not unpaired_vertices:
        print("\n[PROGRESS] No unpaired vertices found (all features have at least one edge)")

    return unpaired_vertices


# ============================================================================
# Core pairing
# ============================================================================

def _resolve_rt_correction(query_file_name: str, ref_file_name: str,
                           rt_corrections: Dict[Tuple[str, str], Dict]
                           ) -> Tuple[Optional[Callable], float, str]:
    """
    Pick the model that maps `query_file_name`'s RTs onto `ref_file_name`'s RT scale.

    Curves are fitted by compare_retention_time_alignment.py as
        model_(A,B)(rt_A) = rt_B - rt_A
    i.e. the model's x-axis is FILE 1 (A) of the key pair. Two keys can therefore serve
    a given (ref, query) pair, and they are NOT equivalent:

      reverse key (query, ref):  model(rt_query) = rt_ref - rt_query
          -> rt_ref = rt_query + model(rt_query).  Evaluated on its OWN fitted domain
             (the query's RT), so this is EXACT.  Preferred.

      forward key (ref, query):  model(rt_ref) = rt_query - rt_ref
          -> rt_ref = rt_query - model(rt_ref), an implicit equation. Evaluating
             model(rt_query) instead is one fixed-point iteration, so it is APPROXIMATE
             (error ~ model'(rt) * drift; measured at ~0.1-4.4 s on real data).

    The reverse (exact) key is always tried first; the forward key is only a fallback for
    the rare case where a corrections file stores just one direction per pair.

    Returns (model, sign, status) where the RT shift is sign * model(rt) and status is
    one of 'exact' / 'approximate' / 'missing'. A pair with no key, or with a key whose
    record holds no usable fit, yields (None, 0.0, ...) and is matched unaligned.
    """
    for key, sign, status in (((query_file_name, ref_file_name), +1.0, 'exact'),
                              ((ref_file_name, query_file_name), -1.0, 'approximate')):
        if key in rt_corrections:
            model = build_rt_correction_model(rt_corrections[key])
            return (model, sign, status) if model is not None else (None, 0.0, 'missing')
    return None, 0.0, 'missing'


def _build_frame(features_by_charge: Dict[int, List[Dict]], scale_coords: Callable,
                 model: Optional[Callable] = None, sign: float = 0.0,
                 with_kdtree: bool = True) -> Dict[int, Dict]:
    """
    Express one file's features in one reference file's RT frame, split by charge.

    A "frame" is a reference file's RT coordinate system: that file keeps its native RTs
    and every other file is shifted onto them. The same file therefore has a different RT
    in each reference's frame, which is why the shifted values live here instead of being
    written back onto the feature dicts.

    Returns {charge: {'features', 'mz', 'rt', 'scaled', 'kdtree'}}. When `model` is given
    the RTs are shifted onto the reference frame by sign * model(rt) - one vectorized call
    per charge group. 'rt' always holds the values the cutoffs are checked against and
    'scaled' the normalized/weighted form the KD-tree metric uses. `with_kdtree` is False
    for the reference file itself, which is queried against but never queried into.
    """
    frame = {}
    for charge, charge_features in features_by_charge.items():
        mz = np.array([f['x_center'] for f in charge_features], dtype=float)
        rt = np.array([f['y_center'] for f in charge_features], dtype=float)
        if model is not None and len(rt) > 0:
            rt = rt + sign * np.asarray(model(rt), dtype=float).reshape(-1)
        coords = np.column_stack((mz, rt)) if len(mz) > 0 else np.empty((0, 2), dtype=float)
        scaled = scale_coords(coords)
        frame[charge] = {
            'features': charge_features,
            'mz': mz,
            'rt': rt,
            'scaled': scaled,
            'kdtree': cKDTree(scaled) if (with_kdtree and len(scaled) > 0) else None,
        }
    return frame


def feature_pairing_optimized(all_feature_data_lists: List[List[Dict]],
                              normalize_coordinates: bool = True,
                              distance_calc_before_scaling: bool = False,
                              rt_correction_json: str = None,
                              mz_rt_weight_ratio: float = 1.0,
                              mz_cutoff: float = None,
                              rt_cutoff: float = None,
                              edges_cutoff: float = None,
                              ppm: float = None,
                              max_match_candidates: int = 1) -> Tuple[dict, dict]:
    """
    Pair features across files with charge-split KD-trees and vectorized distances.

    For every feature, up to max_match_candidates nearest neighbours of the same charge in
    every other file are queried (ascending distance) and the first that passes the active
    cutoff becomes an edge. The cutoff is always evaluated in ORIGINAL coordinate space.
    These options are available for cutoffs:
        - ppm + rt_cutoff : per-feature m/z tolerance (mz * ppm / 1e6) plus an RT cutoff,
        - mz_cutoff + rt_cutoff : fixed m/z and RT cutoffs,
        - edges_cutoff : euclidean distance cutoff,
        - otherwise : no filtering.

    The stored edge 'distance' is the euclidean distance in the space selected by
    distance_calc_before_scaling (original if True, else the scaled/normalized space).

    Args:
        normalize_coordinates: Scale both axes into [0, 1] before building KD-trees.
        distance_calc_before_scaling: Compute the stored 'distance' in original space. #TODO obsolete, remove
        rt_correction_json: Path to fitted RT-correction models. When the file exists,
            every ordered file pair is aligned with the model fitted for that exact pair
            (see the per-reference-frame matching below).
        mz_rt_weight_ratio: RT weight relative to m/z in the scaled space (1.0 = equal).
        mz_cutoff / rt_cutoff / edges_cutoff / ppm: cutoff parameters (see above).
        max_match_candidates: how many nearest neighbours (in scaled space, ascending
            distance) to consider per (ref feature, query file) before giving up. The
            FIRST one that passes the cutoff becomes the edge. 1 = legacy behaviour
            (only the single nearest is tried, and the pair is dropped if it fails the
            cutoff). Larger values recover a valid in-cutoff match when the nearest
            neighbour falls outside the cutoff.

    Returns:
        (edges, file_index_map). edges has the form {ref_file_idx: [edge_dict, ...]}.
        Connected-component grouping is done downstream, not here.
    """
    if not all_feature_data_lists:
        return {}, {}

    # file_idx -> filename, used to restore filenames from indices downstream.
    file_index_map = {}
    for file_idx, features in enumerate(all_feature_data_lists):
        if features:
            file_index_map[file_idx] = features[0].get('filename', f'file_{file_idx}')

    # Single file: no cross-file edges possible.
    if len(all_feature_data_lists) == 1:
        return {0: []}, file_index_map

    print("\n[PROGRESS] Starting feature pairing (optimized, KD-tree nearest-neighbour)")
    # TODO: Push out, load after Aligntree loadin
    # --- Step 1: load JSON RT corrections (if available) -------------------
    rt_corrections = {}
    json_str = str(rt_correction_json).strip() if rt_correction_json else ""
    if json_str and os.path.exists(json_str):
        print(f"[PROGRESS] Loading RT corrections from: {json_str}")
        rt_corrections = load_rt_corrections_json(json_str)
        print(f"  ✓ Loaded {len(rt_corrections)} RT correction pairs")
    else:
        print("  ⚠ No RT correction JSON file provided - proceeding without corrections")

    # --- Step 2: coordinate normalization factors -------------------------
    # scale_coords maps (x, y) -> ((x-min_x)*scale_x, (y-min_y)*scale_y_weighted).
    # Note: Scaling occours before the rt-corrections are applied. Shouldnt make a large differance, but is still wrong in theory.
    # Best Practice would be to calculate the Normalization factors on corrected RT-Data
    # The one set of factors is shared by every reference frame, so the stored 'distance'
    # stays comparable across edges coming from different frames.
    scale_x = scale_y_weighted = 1.0
    min_x = min_y = 0.0
    if normalize_coordinates:
        all_x = [f['x_center'] for features in all_feature_data_lists for f in features]
        all_y = [f['y_center'] for features in all_feature_data_lists for f in features]
        if all_x and all_y:
            min_x, max_x = np.min(all_x), np.max(all_x)
            min_y, max_y = np.min(all_y), np.max(all_y)
            range_x = max_x - min_x if max_x > min_x else 1.0
            range_y = max_y - min_y if max_y > min_y else 1.0
            scale_x = 1.0 / range_x if range_x > 0 else 1.0
            scale_y = 1.0 / range_y if range_y > 0 else 1.0
            # mz_rt_weight_ratio > 1 emphasizes RT, < 1 emphasizes m/z.
            scale_y_weighted = scale_y * mz_rt_weight_ratio
            print(f"[PROGRESS] Coordinate normalization: X scale {scale_x:.6f}, "
                  f"Y scale {scale_y_weighted:.6f} (mz:RT ratio 1.0:{mz_rt_weight_ratio:.4f})")

    def scale_coords(coords):
        """Apply normalization + dimension weighting to an (N, 2) coordinate array."""
        if normalize_coordinates and len(coords) > 0:
            scaled = coords.astype(float).copy()
            scaled[:, 0] = (scaled[:, 0] - min_x) * scale_x
            scaled[:, 1] = (scaled[:, 1] - min_y) * scale_y_weighted
            return scaled
        return coords.astype(float)

    # --- Step 3: group features by charge, tagging original indices --------
    # file_features_by_charge[file_idx][charge] -> [feature_copy_with__file_idx/_feat_idx]
    file_features_by_charge = []
    for file_idx, features in enumerate(all_feature_data_lists):
        features_by_charge = defaultdict(list)
        for feat_idx, f in enumerate(features):
            f_copy = f.copy()
            f_copy['_feat_idx'] = feat_idx
            f_copy['_file_idx'] = file_idx
            features_by_charge[int(f['charge'])].append(f_copy)  # charge guaranteed present
        file_features_by_charge.append(dict(features_by_charge))

    print("  Charge distribution by file:")
    for file_idx, charge_dict in enumerate(file_features_by_charge):
        print(f"    File {file_idx}: {dict(sorted((c, len(f)) for c, f in charge_dict.items()))}")

    # --- Steps 4+5: per-pair RT alignment, KD-trees, nearest-neighbour matching ---
    # A "frame" is one file's features expressed in one reference file's RT coordinate
    # system: the reference keeps its native RTs and every other file is shifted onto them
    # by the model fitted for that exact pair. Each ordered pair is therefore matched in
    # the reference's own frame, using its own direct model rather than a correction
    # composed through some global reference file.
    #
    # A file appears in n-1 different frames with a different RT in each, so the alignment
    # must NOT be written back onto the feature dicts - every frame carries its own 'rt'
    # array, and the trees are rebuilt per reference file. That costs n*(n-1) tree builds
    # instead of n. Without corrections all frames are identical, so that case still builds
    # the trees once and reuses them.
    print("[PROGRESS] Performing pairwise feature matching...")
    if rt_corrections:
        print("  RT alignment: per file pair, in each reference file's own frame")
    edges = {}
    num_input_files = len(all_feature_data_lists)
    alignment_stats = defaultdict(int)

    static_frames = None
    if not rt_corrections:
        static_frames = [_build_frame(charge_dict, scale_coords)
                         for charge_dict in file_features_by_charge]

    for ref_file_idx in range(num_input_files):
        edges[ref_file_idx] = []
        ref_file_name = file_index_map.get(ref_file_idx)

        # Build this reference's frame: itself unshifted, every other file aligned onto it.
        if static_frames is not None:
            ref_frame = static_frames[ref_file_idx]
            query_frames = {j: static_frames[j]
                            for j in range(num_input_files) if j != ref_file_idx}
        else:
            ref_frame = _build_frame(file_features_by_charge[ref_file_idx], scale_coords,
                                     with_kdtree=False)
            query_frames = {}
            for query_file_idx in range(num_input_files):
                if query_file_idx == ref_file_idx:
                    continue
                query_file_name = file_index_map.get(query_file_idx)
                model, sign, status = _resolve_rt_correction(query_file_name, ref_file_name,
                                                             rt_corrections)
                alignment_stats[status] += 1
                if status == 'missing':
                    print(f"    No usable correction for {repr(query_file_name)} -> "
                          f"{repr(ref_file_name)}, matching unaligned")
                query_frames[query_file_idx] = _build_frame(
                    file_features_by_charge[query_file_idx], scale_coords, model, sign)

        for charge, ref_group in ref_frame.items():
            ref_features = ref_group['features']
            ref_mz, ref_rt, ref_scaled = ref_group['mz'], ref_group['rt'], ref_group['scaled']

            for ref_pos, ref_feature in enumerate(ref_features):
                current_feature_idx = ref_feature['_feat_idx']  # index within its file
                ref_point = ref_scaled[ref_pos]

                for query_file_idx, query_frame in query_frames.items():
                    # Only match against features of the same charge in the query file.
                    query_group = query_frame.get(charge)
                    if query_group is None or query_group['kdtree'] is None:
                        continue

                    query_features = query_group['features']
                    query_mz, query_rt = query_group['mz'], query_group['rt']

                    # Query the up-to-N nearest neighbours (ascending scaled distance) and
                    # keep the FIRST that passes the cutoff. N == 1 reproduces the legacy
                    # nearest-only behaviour; larger N recovers a valid in-cutoff match when
                    # the nearest neighbour(s) fall outside the cutoff.
                    k = min(max(1, max_match_candidates), len(query_features))
                    _, cand_idxs = query_group['kdtree'].query(ref_point, k=k)
                    cand_idxs = np.atleast_1d(np.asarray(cand_idxs).ravel())

                    for cand in cand_idxs:
                        idx = int(cand)
                        # scipy pads missing neighbours (k > n) with the tree size.
                        # TODO add warning if taht ahppens??
                        if idx < 0 or idx >= len(query_features):
                            continue

                        # Cutoffs are ALWAYS evaluated in unscaled coordinate space, on this
                        # frame's RTs (reference native vs query aligned onto it).
                        original_mz_distance = abs(ref_mz[ref_pos] - query_mz[idx])
                        original_rt_distance = abs(ref_rt[ref_pos] - query_rt[idx])

                        if ppm is not None and rt_cutoff is not None:
                            effective_mz_cutoff = ref_mz[ref_pos] * (ppm / 1_000_000)
                            if original_mz_distance > effective_mz_cutoff or original_rt_distance > rt_cutoff:
                                continue
                        elif mz_cutoff is not None and rt_cutoff is not None:
                            if original_mz_distance > mz_cutoff or original_rt_distance > rt_cutoff:
                                continue
                        elif edges_cutoff is not None:
                            original_exact_distance = float(
                                np.sqrt(original_mz_distance ** 2 + original_rt_distance ** 2))
                            if original_exact_distance > edges_cutoff:
                                continue

                        # Passed — store 'distance' in the selected space (original or
                        # scaled) and add the edge, then stop scanning worse candidates.
                        if distance_calc_before_scaling:
                            mz_sq = original_mz_distance ** 2
                            rt_sq = original_rt_distance ** 2
                        else:
                            matched_scaled = query_group['scaled'][idx]
                            mz_sq = (ref_point[0] - matched_scaled[0]) ** 2
                            rt_sq = (ref_point[1] - matched_scaled[1]) ** 2
                        exact_distance = float(np.sqrt(mz_sq + rt_sq))

                        edges[ref_file_idx].append({
                            'current_feature_idx': current_feature_idx,
                            'matched_feature_idx': query_features[idx]['_feat_idx'],
                            'distance': exact_distance,
                            'ref_file_idx': ref_file_idx,
                            'match_file_idx': query_file_idx,
                            'charge': charge,
                        })
                        break

        print(f"    File {ref_file_idx} complete: {len(edges[ref_file_idx])} edges")
        # This reference's frames are dead once its edges are built.
        del ref_frame, query_frames

    if rt_corrections:
        print(f"  File pairs aligned via exact (reverse) key: {alignment_stats['exact']}, "
              f"via approximate (forward) key: {alignment_stats['approximate']}, "
              f"unaligned: {alignment_stats['missing']}")

    total_edges = sum(len(edge_list) for edge_list in edges.values())
    print(f"[PROGRESS] Pairwise matching done: {total_edges} edges")

    # Match consolidation is intentionally not done here (handled downstream). Free the
    # large intermediate structures before returning.
    del file_features_by_charge, static_frames
    gc.collect()
    return edges, file_index_map


# ============================================================================
# Cycle 2+ input loading (from previous cycle's deleted vertices)
# ============================================================================

def load_vertices_from_json(json_path: str) -> List[Dict]:
    """
    Load vertices from a previous cycle's deleted_vertices.json for re-pairing.
    
    These vertices carry ORIGINAL (unaligned, unweighted) coordinates so this cycle can
    apply fresh RT alignment. Raises ValueError if required fields are missing.

    For the fitted-JSON RT paths (loess/polynomial/spline/rbf/piecewise) these vertices
    carry ORIGINAL (unaligned, unweighted) coordinates — those corrections are applied to
    copies during pairing, so each cycle re-aligns from scratch. For the aligntree path the
    carried y_center is whatever the previous cycle's trafoXML transform produced, but that
    is harmless: aligntree re-assigns y_center from the trafoXML table (keyed by the stable
    openms_fid) every cycle, so the transform is idempotent regardless of the incoming
    coordinate. Raises ValueError if required fields are missing.
    """
    with open(json_path, 'r') as f:
        data = json.load(f)
    vertices = data.get('deleted_vertices', [])
    print(f"  Loaded {len(vertices)} vertices from previous cycle for re-pairing")

    if vertices and isinstance(vertices[0], dict):
        required_fields = ['filename', 'feature_idx', 'x_center', 'y_center']
        missing = [f for f in required_fields if f not in vertices[0]]
        if missing:
            raise ValueError(f"Vertices missing required fields: {missing}")
    return vertices


def group_vertices_by_file(vertices: List[Dict]) -> List[List[Dict]]:
    """Group vertices into per-file feature lists, sorted by filename for consistency."""
    file_features = defaultdict(list)
    for vertex in vertices:
        file_features[vertex['filename']].append(vertex)

    result = [file_features[name] for name in sorted(file_features.keys())]
    print(f"  Grouped vertices into {len(result)} files")
    for idx, features in enumerate(result):
        print(f"    File {idx}: {len(features)} features")
    return result


# ============================================================================
# main
# ============================================================================

def _load_feature_data(args) -> Tuple[List[List[Dict]], bool]:
    """
    Load per-file feature lists for the current cycle.

    Cycle 1 reads *_feature_data.pkl (from --input_files, --input_dir, or the CWD).
    Cycle 2+ reads vertices from --input_vertices_json. Returns (feature_lists,
    has_charge_field) where has_charge_field is None if nothing could be loaded.
    """
    if args.input_vertices_json:
        print(f"\n[CYCLE {args.cycle_num}] Loading vertices from previous cycle: {args.input_vertices_json}")
        vertices = load_vertices_from_json(args.input_vertices_json)
        feature_lists = group_vertices_by_file(vertices)
        has_charge_field = 'charge' in vertices[0] if (vertices and isinstance(vertices[0], dict)) else None
        print(f"\n  Loaded {len(vertices)} vertices from {len(feature_lists)} files")
        return feature_lists, has_charge_field

    print(f"\n[CYCLE {args.cycle_num}] Loading features from pickle files...")
    if args.input_files:
        feature_files = sorted(Path(f) for f in args.input_files)
        print(f"  Using {len(feature_files)} explicitly provided feature file(s)")
    else:
        search_dir = args.input_dir if args.input_dir else '.'
        feature_files = sorted(Path(search_dir).glob("*_feature_data.pkl"))
        if not feature_files:
            print(f"✗ No feature pickle files found in {search_dir}!")
            print("  Use --input_dir <dir> or --input_files <file1> <file2> ...")
            return [], None
        print(f"  Found {len(feature_files)} feature file(s) in {search_dir}")

    feature_lists = []
    has_charge_field = None
    for pkl_file in feature_files:
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
        feature_lists.append(data)
        if has_charge_field is None and isinstance(data, list) and data and isinstance(data[0], dict):
            has_charge_field = 'charge' in data[0]
    return feature_lists, has_charge_field


def main():
    parser = argparse.ArgumentParser(
        description="Pair matching features across multiple files (clusterless version)")
    # Inputs
    parser.add_argument("--input_dir", help="Directory containing feature pickle files")
    parser.add_argument("--input_files", nargs='*', help="Explicit list of feature pickle files")
    parser.add_argument("--input_vertices_json", help="JSON file with vertices from previous cycle (cycle 2+)")
    parser.add_argument("--cycle_num", type=int, default=1, help="Current cycle number for tracking")
    # Outputs
    parser.add_argument("--output_pkl", required=True, help="Output pickle file path (edges + ident_lookup + metadata; consumed by build_network_graph_clusterless.py)")
    parser.add_argument("--output_json", required=True,
                        help="Path anchor for the derived ..._ident_lookup.json and ..._metadata.json outputs "
                             "(the full edges JSON itself is no longer written)")
    parser.add_argument("--output_unpaired_json", help="Output JSON file for unpaired vertices")
    # Cutoffs
    parser.add_argument("--edges_cutoff", type=float, default=None, help="Maximum euclidean distance cutoff for edges")
    parser.add_argument("--mz_cutoff", type=float, default=None, help="Maximum m/z difference for edge filtering")
    parser.add_argument("--ppm", type=float, default=None, help="Maximum m/z difference in ppm (mz*ppm/1e6). Mutually exclusive with --mz_cutoff")
    parser.add_argument("--rt_cutoff", type=float, default=None, help="Maximum RT difference for edge filtering")
    parser.add_argument("--match_cutoff", type=float, default=None, help="Minimum match score cutoff (recorded in metadata)")
    # Matching / normalization
    parser.add_argument("--max_match_candidates", type=int, default=1,
                        help="Nearest neighbours to try per (ref feature, query file), in ascending "
                             "distance, keeping the first that passes the cutoff. 1 = legacy "
                             "nearest-only. Larger values recover a valid in-cutoff match when the "
                             "nearest neighbour(s) fall outside the cutoff.")
    parser.add_argument("--distance_calc_before_scaling", default=False, action='store_true',
                        help="Store edge distance in original space (default: scaled space)")
    parser.add_argument("--normalize_coordinates", type=lambda x: x.lower() in ('true', '1', 'yes'),
                        nargs='?', const=True, default=None,
                        help="Enable/disable coordinate normalization before KD-tree matching. "
                             "Omitted => AUTO: OFF for coordinate-based cutoffs (mz/ppm + rt), "
                             "ON otherwise.")
    parser.add_argument("--mz_rt_weight_ratio", type=float, default=1.0,
                        help="RT weight relative to m/z (>1 emphasizes RT, <1 emphasizes m/z)")
    # RT correction
    parser.add_argument("--rt_correction_json", type=str, default=None,
                        help="Path to fitted RT correction models JSON. Each file pair is "
                             "aligned with the model fitted for that exact pair.")
    parser.add_argument("--trafoxmls_dir", type=str, default=None, help="Directory of trafoXML files for the AlignTree method")
    parser.add_argument("--rt_alignment_method", type=str, default='loess',
                        help="RT alignment method. 'loess' (and the other fitted methods: "
                             "polynomial, spline, rbf, piecewise) uses our own pairwise "
                             "alignment fitted on shared pep_idents, loaded from "
                             "--rt_correction_json and applied per file pair. "
                             "'aligntree' instead uses the global trafoXML transformations "
                             "produced by the OpenMS tree-guided aligner, read from "
                             "--trafoxmls_dir. 'noalign' applies no RT alignment at all.")

    args = parser.parse_args()

    if args.mz_cutoff is not None and args.ppm is not None:
        print("ERROR: --mz_cutoff and --ppm are mutually exclusive.")
        sys.exit(1)

    print(f"\n{'=' * 70}")
    print(f"Feature Pairing (Clusterless - Cycle {args.cycle_num})")
    print(f"{'=' * 70}")

    # --- Load feature data -------------------------------------------------
    all_feature_data_lists, has_charge_field = _load_feature_data(args)

    if has_charge_field is False:
        print("\n✗ ERROR: Feature pickle files are missing the required 'charge' field")
        sys.exit(1)
    if has_charge_field is None:
        print("\n✗ ERROR: Could not load any feature data")
        sys.exit(1)

    total_features_before_filtering = sum(len(fl) for fl in all_feature_data_lists)
    num_input_files = len(all_feature_data_lists)
    print(f"\nTotal features loaded: {total_features_before_filtering}")

    # aligntree requires a trafoXML directory.
    if args.rt_alignment_method == 'aligntree' and not args.trafoxmls_dir:
        print("\n✗ ERROR: RT alignment method 'aligntree' requires --trafoxmls_dir <path>")
        sys.exit(1)

    # --- Decide coordinate normalization -----------------------------------
    # Auto (when --normalize_coordinates omitted): OFF for coordinate-based cutoffs,
    # ON for euclidean-distance cutoff and standard matching.
    # TODO: Remove Euclidian Distance CUtoff and resuidual functions
    if args.normalize_coordinates is None:
        coord_cutoff = (args.mz_cutoff is not None or args.ppm is not None) and args.rt_cutoff is not None
        normalize_coordinates = not coord_cutoff
        print(f"  ℹ Coordinate normalization: AUTO ({'OFF' if coord_cutoff else 'ON'})")
    else:
        normalize_coordinates = args.normalize_coordinates
        print(f"  ℹ Coordinate normalization: MANUAL ({'ON' if normalize_coordinates else 'OFF'})")

    # --- Optional No RT transformation -------------------
    # 'noalign' disables every correction path: no trafoXML, no JSON models.
    if args.rt_alignment_method == 'noalign':
        print("\n[PROGRESS] RT alignment method 'noalign' - pairing without any RT corrections")
        args.rt_correction_json = None
        args.trafoxmls_dir = None
     # --- Optional AlignTree (trafoXML) RT transformation -------------------
    elif args.trafoxmls_dir and args.rt_alignment_method == 'aligntree':
        print(f"\n[PROGRESS] Loading AlignTree RT transformations from: {args.trafoxmls_dir}")
        try:
            input_filenames = {fl[0].get('filename', '') for fl in all_feature_data_lists
                               if fl and fl[0].get('filename', '')}
            fid_to_rt_dict = build_trafoxmls_fid_dict(args.trafoxmls_dir, input_filenames=input_filenames)
            all_feature_data_lists = apply_aligntree_rt_transformations(all_feature_data_lists, fid_to_rt_dict)
            print("  ✓ AlignTree transformations applied successfully")
            del fid_to_rt_dict
            # Corrections already applied via transformations; disable JSON correction.
            args.rt_correction_json = None
        except Exception as e:
            print(f"  ✗ ERROR applying AlignTree transformations: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    elif args.trafoxmls_dir:
        print(f"\n[WARNING] trafoXML directory provided but rt_alignment_method is "
              f"'{args.rt_alignment_method}' (not 'aligntree') - skipping trafoXML transformations")

    # --- Capture the coordinate-scaling constants --------------------------
    # feature_pairing_optimized derives these internally to build its scaled space,
    # but keeps them local. The anchor/window filters in build_network_graph need
    # the same space to rank candidates in, and cannot recompute it faithfully:
    # they only ever see PAIRED features, whereas these cover every input feature.
    # Computed here, immediately before the call and on the same data, so the values
    # are exactly the ones pairing goes on to use.
    #
    # Only the axis RANGES are exported. The origin cancels - the filters use these
    # solely to take differences between two scaled points - and the weight ratio is
    # passed to the graph step from the workflow parameter directly.
    #
    # Stored UNCONDITIONALLY, including when normalize_coordinates is off. That flag
    # governs whether PAIRING scales before building its KD-trees; the filters scale
    # regardless, because a raw-unit euclidean over (m/z in Da, RT in s) is dominated
    # by RT and would make the weight ratio meaningless.
    scaling_constants = {'range_x': None, 'range_y': None}
    _all_x = [f['x_center'] for fl in all_feature_data_lists for f in fl]
    _all_y = [f['y_center'] for fl in all_feature_data_lists for f in fl]
    if _all_x and _all_y:
        _min_x, _max_x = float(np.min(_all_x)), float(np.max(_all_x))
        _min_y, _max_y = float(np.min(_all_y)), float(np.max(_all_y))
        scaling_constants.update(
            range_x=(_max_x - _min_x) if _max_x > _min_x else 1.0,
            range_y=(_max_y - _min_y) if _max_y > _min_y else 1.0)
    del _all_x, _all_y
    print(f"[PROGRESS] Coordinate scaling constants: X range "
          f"{scaling_constants['range_x']}, Y range {scaling_constants['range_y']}")

    # --- Pair features -----------------------------------------------------
    network_edges, file_index_map = feature_pairing_optimized(
        all_feature_data_lists,
        normalize_coordinates=normalize_coordinates,
        distance_calc_before_scaling=args.distance_calc_before_scaling,
        rt_correction_json=args.rt_correction_json,
        mz_rt_weight_ratio=args.mz_rt_weight_ratio,
        mz_cutoff=args.mz_cutoff,
        rt_cutoff=args.rt_cutoff,
        edges_cutoff=args.edges_cutoff,
        ppm=args.ppm,
        max_match_candidates=args.max_match_candidates,
    )

    print(f"Built network graph with {sum(len(e) for e in network_edges.values())} edges")
    print(f"  File index mapping: {file_index_map}")

    # Count unique features (vertices) that survived cutoff filtering.
    unique_vertices = set()
    for edge_list in network_edges.values():
        for edge in edge_list:
            unique_vertices.add((edge['ref_file_idx'], edge['current_feature_idx']))
            unique_vertices.add((edge['match_file_idx'], edge['matched_feature_idx']))
    total_features_after_filter = len(unique_vertices)
    print(f"Total unique features in network edges after filtering: {total_features_after_filter}")

    # Unpaired vertices (features with no edges).
    extract_and_save_unpaired_vertices(network_edges, all_feature_data_lists, args.output_unpaired_json)

    # --- Build pep/prot identification lookup ------------------------------
    # Keyed "filename_featureidx"; built before feature data is freed.
    print("\nBuilding combined pep_ident/prot_ident lookup from all features...")
    ident_lookup = {}
    for file_idx, file_features in enumerate(all_feature_data_lists):
        for feat_idx, feature in enumerate(file_features):
            if not isinstance(feature, dict):
                continue
            filename = feature.get('filename', f'file_{file_idx}')
            ident_lookup[f"{filename}_{feat_idx}"] = {
                'pep_ident': feature.get('pep_ident', []) or [],
                'prot_ident': feature.get('prot_ident', []) or [],
                'x_center': feature.get('x_center', 0.0),
                'y_center': feature.get('y_center', 0.0),
                'intensity': feature.get('intensity', 0.0),
                'openms_fid': feature.get('openms_fid', ''),
                'ms2_scans': feature.get('ms2_scans', ''),
                'charge': feature.get('charge', 0),
            }
    print(f"Extracted identifications for {len(ident_lookup)} features")

    # Feature data is no longer needed; free it before serialization.
    del all_feature_data_lists
    gc.collect()

    # Flatten {file_idx: [edge, ...]} into a single list, popping as we go so the dict
    # and the full list never coexist in memory.
    flat_edges = []
    for key in list(network_edges.keys()):
        flat_edges.extend(network_edges.pop(key))
    del network_edges
    gc.collect()

    metadata_fields = {
        'total_features_before_filtering': total_features_before_filtering,
        'total_features_after_filtering': total_features_after_filter,
        'cutoff_params': {
            'mz_cutoff': args.mz_cutoff,
            'ppm': args.ppm,
            'rt_cutoff': args.rt_cutoff,
            'edges_cutoff': args.edges_cutoff,
        },
        'scaling_constants': scaling_constants,
        'pairing_parameters': {
            'pairing_method': 'kdtree',
            'best_match_only': True,
            'distance_calc_before_scaling': args.distance_calc_before_scaling,
            'normalize_coordinates': normalize_coordinates,
            'match_cutoff': args.match_cutoff,
            'max_match_candidates': args.max_match_candidates,
            'input_files': num_input_files,
        },
    }

    # --- Save primary output: paired_features_edges.pkl --------------------
    # This is the only edges serialization the workflow consumes: build_network_graph
    # reads it via --input_pkl (network_edges + file_index_map + cutoff_params). The
    # ident_lookup is also written separately below (that separate file is what
    # build_network_graph reads via --ident-lookup).
    output_data = {
        'network_edges': flat_edges,
        'ident_lookup': ident_lookup,
        'file_index_map': file_index_map,
        **metadata_fields,
    }

    print("Saving output (pickle)...")
    with open(args.output_pkl, 'wb') as f:
        pickle.dump(output_data, f)
    print(f"✓ Saved edges and ident_lookup (pickle): {os.path.basename(args.output_pkl)}")

    # Combined ident lookup, saved separately for fast loading downstream.
    ident_lookup_path = args.output_json.replace('.json', '_ident_lookup.json')
    with open(ident_lookup_path, 'w') as f:
        json.dump(ident_lookup, f)
    print(f"✓ Saved combined ident lookup: {os.path.basename(ident_lookup_path)}")

    del output_data, ident_lookup, flat_edges
    gc.collect()

    # --- Save metadata JSON (everything except edges + ident_lookup) -------
    metadata_json_path = args.output_json.replace('.json', '_metadata.json')
    try:
        with open(metadata_json_path, 'w') as f:
            json.dump(metadata_fields, f, indent=2)
        print(f"✓ Saved pairing metadata: {os.path.basename(metadata_json_path)}")
    except Exception as e:
        print(f"✗ Failed to save metadata JSON: {e}")


if __name__ == "__main__":
    main()
