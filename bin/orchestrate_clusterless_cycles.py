#!/usr/bin/env python3
"""
Orchestrate clusterless iterative component building.

This script manages cycles of:
1. Feature pairing
2. Component building and filtering
3. Collecting accepted components
4. Feeding deleted vertices to next cycle

Continues until convergence (no edges or no deleted vertices).
"""

import os
import sys
import json
import argparse
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Tuple


def load_json(path: str) -> Dict:
    """Load JSON file."""
    with open(path, 'r') as f:
        return json.load(f)


def save_json(data: Dict, path: str):
    """Save JSON file."""
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)


def count_vertices_in_pkls(pkl_paths: List[str]) -> int:
    """Count total features across all input PKL files (each PKL is a list of features)."""
    import pickle
    total = 0
    for path in pkl_paths:
        with open(path, 'rb') as f:
            data = pickle.load(f)
        total += len(data) if isinstance(data, list) else 0
    return total


def count_unpaired_in_json(unpaired_json: str) -> int:
    """Count unpaired vertices in a paired_features_unpaired.json file."""
    if not unpaired_json or not os.path.exists(unpaired_json):
        return 0
    data = load_json(unpaired_json)
    if isinstance(data, list):
        return len(data)
    return len(data.get('unpaired_vertices', []))


def vertex_key(v: Dict) -> tuple:
    """Return a stable hashable key for vertex deduplication across cycles.

    openms_fid is unique per feature within a file and is preserved unchanged
    across all cycle iterations, making it the correct stable identity key.
    """
    return (v.get('filename', ''), v.get('openms_fid', ''))


def count_edges_in_pkl(pkl_path: str) -> int:
    """Count total edges in edges pickle file."""
    import pickle
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)

    # Unwrap if stored as {network_edges: ..., ...}
    if isinstance(data, dict) and 'network_edges' in data:
        edges = data['network_edges']
    else:
        edges = data

    # Flat list format (new default)
    if isinstance(edges, list):
        return len(edges)

    # Legacy nested dict format {file_idx: [edge, ...]}
    if isinstance(edges, dict):
        return sum(len(v) for v in edges.values() if isinstance(v, list))

    return 0


def run_pairing_cycle(
    cycle_num: int,
    input_type: str,
    input_data: str,
    cycle_dir: str,
    config: Dict,
    script_dir: str
) -> Tuple[str, str, str, str]:
    """
    Run pairing for a single cycle.
    
    Args:
        cycle_num: Current cycle number
        input_type: 'tsv' or 'vertices_json'
        input_data: Path to input (TSV files glob or vertices JSON)
        cycle_dir: Output directory for this cycle
        config: Configuration dict with cutoffs, etc.
        script_dir: Directory containing scripts
    
    Returns:
        Tuple of (edges_pkl, edges_json, unpaired_json, ident_lookup_json)
    """
    print(f"\n  [Cycle {cycle_num}] Running feature pairing...")
    
    pair_script = os.path.join(script_dir, 'pair_features_1_clusterless.py')
    
    # Build command
    cmd = [
        pair_script,
        '--output_pkl', os.path.join(cycle_dir, 'paired_features_edges.pkl'),
        '--output_json', os.path.join(cycle_dir, 'paired_features_edges.json'),
        '--output_unpaired_json', os.path.join(cycle_dir, 'paired_features_unpaired.json'),
        '--cycle_num', str(cycle_num),
        '--skip-matchfinder',  # We only need edges for graph building
    ]
    
    # Input handling
    if input_type == 'tsv':
        # Cycle 1: Use input files
        cmd.extend(['--input_files'] + input_data.split())
    else:
        # Cycle 2+: Use vertices JSON
        cmd.extend(['--input_vertices_json', input_data])
    
    # Add cutoffs (ppm and mz_cutoff are mutually exclusive)
    if config.get('ppm') is not None:
        cmd.extend(['--ppm', str(config['ppm'])])
    elif config.get('mz_cutoff') is not None:
        cmd.extend(['--mz_cutoff', str(config['mz_cutoff'])])
    if config.get('rt_cutoff') is not None:
        cmd.extend(['--rt_cutoff', str(config['rt_cutoff'])])
    if config.get('edges_cutoff') is not None:
        cmd.extend(['--edges_cutoff', str(config['edges_cutoff'])])

    # RT alignment method (always pass so the pairing script knows which path to take)
    cmd.extend(['--rt_alignment_method', config.get('rt_alignment_method', 'loess')])

    # AlignTree: trafoXML-based corrections applied directly inside the pairing script
    if config.get('trafoxmls_dir'):
        cmd.extend(['--trafoxmls_dir', config['trafoxmls_dir']])

    # JSON-based RT correction (loess/polynomial/etc.)
    if config.get('rt_correction_json'):
        cmd.extend(['--rt_correction_json', config['rt_correction_json']])
        cmd.extend(['--rt_correction_mode', config.get('rt_correction_mode', 'none')])
        cmd.extend(['--rt_source', config.get('rt_source', 'y_center')])
    
    # Other parameters
    if config.get('normalize_coordinates') is not None:
        cmd.extend(['--normalize_coordinates', str(config['normalize_coordinates'])])
    if config.get('distance_calc_before_scaling'):
        cmd.append('--distance_calc_before_scaling')
    if config.get('mz_rt_weight_ratio', 1.0) != 1.0:
        cmd.extend(['--mz_rt_weight_ratio', str(config['mz_rt_weight_ratio'])])
    
    # Run pairing (stream output directly so it's visible in real-time)
    result = subprocess.run(cmd)

    if result.returncode != 0:
        raise RuntimeError(f"Pairing failed (exit code {result.returncode})")

    # Return output paths
    edges_pkl = os.path.join(cycle_dir, 'paired_features_edges.pkl')
    edges_json = os.path.join(cycle_dir, 'paired_features_edges.json')
    unpaired_json = os.path.join(cycle_dir, 'paired_features_unpaired.json')
    ident_lookup_json = os.path.join(cycle_dir, 'paired_features_edges_ident_lookup.json')
    
    return edges_pkl, edges_json, unpaired_json, ident_lookup_json


def run_graph_building_cycle(
    cycle_num: int,
    edges_pkl: str,
    unpaired_json: str,
    ident_lookup_json: str,
    cycle_dir: str,
    config: Dict,
    script_dir: str
) -> Tuple[str, str]:
    """
    Run graph building and filtering for a single cycle.
    
    Args:
        cycle_num: Current cycle number
        edges_pkl: Path to edges pickle file
        unpaired_json: Path to unpaired vertices JSON
        ident_lookup_json: Path to ident lookup JSON
        cycle_dir: Output directory for this cycle
        config: Configuration dict
        script_dir: Directory containing scripts
    
    Returns:
        Tuple of (accepted_components_json, deleted_vertices_json)
    """
    print(f"\n  [Cycle {cycle_num}] Building components and filtering...")
    
    graph_script = os.path.join(script_dir, 'build_network_graph_clusterless.py')
    
    # Build command
    cmd = [
        graph_script,
        '--input_pkl', edges_pkl,
        '--cycle_num', str(cycle_num),
        '--output_accepted_components', os.path.join(cycle_dir, 'accepted_components.json'),
        '--output_deleted_vertices', os.path.join(cycle_dir, 'deleted_vertices.json'),
        '--output_analysis', os.path.join(cycle_dir, 'graph_analysis.json'),
        '--output_composition', os.path.join(cycle_dir, 'graph_composition.json'),
    ]
    
    # Add cutoffs
    if config.get('mz_cutoff') is not None:
        cmd.extend(['--mz_cutoff', str(config['mz_cutoff'])])
    if config.get('rt_cutoff') is not None:
        cmd.extend(['--rt_cutoff', str(config['rt_cutoff'])])
    if config.get('edge_cutoff') is not None:
        cmd.extend(['--edge_cutoff', str(config['edge_cutoff'])])
    
    # Filtering
    if config.get('filter_duplicate_file_vertices', True):
        cmd.append('--filter-duplicate-file-vertices')
        cmd.extend(['--filter-method', config.get('filter_method', 'simplified-modularity')])
        cmd.extend(['--mega_complexity_threshold', str(config.get('mega_complexity_threshold', 5000))])
    
    # Ident lookup
    if ident_lookup_json and os.path.exists(ident_lookup_json):
        cmd.extend(['--ident-lookup', ident_lookup_json])
    
    # Unpaired vertices
    if unpaired_json and os.path.exists(unpaired_json):
        cmd.extend(['--unpaired-json', unpaired_json])
    
    # Multiprocessing
    if config.get('disable_multiprocessing'):
        cmd.append('--disable-multiprocessing')
    
    # Run graph building (stream output directly so it's visible in real-time)
    result = subprocess.run(cmd)

    if result.returncode != 0:
        raise RuntimeError(f"Graph building failed (exit code {result.returncode})")

    # Return output paths
    accepted_json = os.path.join(cycle_dir, 'accepted_components.json')
    deleted_json = os.path.join(cycle_dir, 'deleted_vertices.json')
    
    return accepted_json, deleted_json


def orchestrate_cycles(
    initial_input_files: List[str],
    output_dir: str,
    config: Dict,
    script_dir: str
) -> Dict:
    """
    Main orchestration loop for clusterless iterative component building.
    
    Args:
        initial_input_files: List of initial TSV feature files
        output_dir: Output directory for all results
        config: Configuration dictionary
        script_dir: Directory containing scripts
    
    Returns:
        Dictionary with final results metadata
    """
    print(f"\n{'='*70}")
    print("CLUSTERLESS ITERATIVE COMPONENT BUILDING")
    print(f"{'='*70}")
    print(f"Configuration:")
    if config.get('ppm') is not None:
        print(f"  m/z cutoff: {config['ppm']} ppm (per-feature: mz * ppm / 1000000)")
    else:
        print(f"  m/z cutoff: {config.get('mz_cutoff', 'None')}")
    print(f"  RT cutoff: {config.get('rt_cutoff', 'None')}")
    print(f"  Filter method: {config.get('filter_method', 'simplified-modularity')}")
    print(f"  Max cycles: {config.get('max_cycles', 10)}")
    print(f"{'='*70}\n")
    
    # Count total vertices before any pairing
    initial_vertices = count_vertices_in_pkls(initial_input_files)
    print(f"  Initial vertices (all input files): {initial_vertices:,}")

    # Initialize tracking
    cycle_num = 1
    current_input_type = 'tsv'
    current_input = ' '.join(initial_input_files)
    consecutive_empty_cycles = 0
    
    # Initialize final components structure
    final_components = {
        'metadata': {
            'total_cycles': 0,
            'total_components': 0,
            'total_vertices_in_components': 0,
            'total_unpaired_singletons': 0,
            'cutoff_params': {
                'mz_cutoff': config.get('mz_cutoff'),
                'ppm': config.get('ppm'),
                'rt_cutoff': config.get('rt_cutoff'),
                'edges_cutoff': config.get('edges_cutoff'),
            },
            'filter_params': {
                'filter_method': config.get('filter_method', 'simplified-modularity'),
                'mega_complexity_threshold': config.get('mega_complexity_threshold', 5000),
            }
        },
        'components': {},
        'unpaired_singletons': []
    }
    
    next_component_id = 0
    cycle_stats = []
    accumulated_singletons = {}  # vertex_key -> vertex dict; accumulates all unpaired+deleted across every cycle
    edges_count = 0  # safe default if cycle 1 fails before edges are counted

    # Main cycle loop
    while True:
        print(f"\n{'='*70}")
        print(f"CYCLE {cycle_num}")
        print(f"{'='*70}")
        
        cycle_start_time = time.time()
        
        # Create cycle output directory
        cycle_dir = os.path.join(output_dir, f'cycle_{cycle_num}')
        os.makedirs(cycle_dir, exist_ok=True)
        
        # Step 1: Pairing
        try:
            edges_pkl, edges_json, unpaired_json, ident_lookup_json = run_pairing_cycle(
                cycle_num=cycle_num,
                input_type=current_input_type,
                input_data=current_input,
                cycle_dir=cycle_dir,
                config=config,
                script_dir=script_dir
            )
        except Exception as e:
            print(f"\n✗ Pairing failed in cycle {cycle_num}: {e}")
            break
        
        # Check edges generated
        edges_count = count_edges_in_pkl(edges_pkl)
        print(f"\n  [Cycle {cycle_num}] Edges generated: {edges_count:,}")
        
        if edges_count == 0:
            print(f"\n  Convergence: No edges generated in cycle {cycle_num}")
            if os.path.exists(unpaired_json):
                unpaired_data = load_json(unpaired_json)
                unpaired_list = unpaired_data if isinstance(unpaired_data, list) else unpaired_data.get('unpaired_vertices', [])
                for v in unpaired_list:
                    accumulated_singletons[vertex_key(v)] = v
            break
        
        # Step 2: Graph building and filtering
        try:
            accepted_json, deleted_json = run_graph_building_cycle(
                cycle_num=cycle_num,
                edges_pkl=edges_pkl,
                unpaired_json=unpaired_json,
                ident_lookup_json=ident_lookup_json,
                cycle_dir=cycle_dir,
                config=config,
                script_dir=script_dir
            )
        except Exception as e:
            print(f"\n✗ Graph building failed in cycle {cycle_num}: {e}")
            break
        
        # Step 3: Merge accepted components into final JSON
        if os.path.exists(accepted_json):
            accepted_data = load_json(accepted_json)
            accepted_comps = accepted_data.get('components', [])
            
            vertices_in_accepted = 0
            for comp in accepted_comps:
                # Assign global component ID
                comp_id = f"component_{next_component_id}"
                comp['component_id'] = comp_id
                comp['cycle'] = cycle_num
                
                # Add to final components
                final_components['components'][comp_id] = comp
                next_component_id += 1
                vertices_in_accepted += comp.get('num_vertices', 0)
            
            print(f"\n  [Cycle {cycle_num}] Accepted components: {len(accepted_comps)}")
            print(f"  [Cycle {cycle_num}] Vertices in accepted components: {vertices_in_accepted:,}")
        else:
            accepted_comps = []
            vertices_in_accepted = 0
        
        # Step 4: Collect deleted vertices from graph building
        deleted_count = 0
        if os.path.exists(deleted_json):
            deleted_data = load_json(deleted_json)
            deleted_vertices = deleted_data.get('deleted_vertices', [])
            deleted_count = len(deleted_vertices)
            print(f"  [Cycle {cycle_num}] Deleted vertices: {deleted_count:,}")
        else:
            deleted_vertices = []

        # Collect unpaired vertices (no cross-file match within cutoff this cycle)
        unpaired_count = count_unpaired_in_json(unpaired_json)
        if unpaired_count > 0:
            unpaired_data = load_json(unpaired_json)
            unpaired_vertices_list = unpaired_data if isinstance(unpaired_data, list) else unpaired_data.get('unpaired_vertices', [])
        else:
            unpaired_vertices_list = []
        print(f"  [Cycle {cycle_num}] Unpaired vertices (no edge after pairing): {unpaired_count:,}")

        # Accumulate all unpaired and deleted vertices into the singleton tracking dict.
        # Using vertex_key as dict key deduplicates automatically if the same vertex
        # appears in multiple cycles (e.g. carried forward through combined input).
        for v in unpaired_vertices_list:
            accumulated_singletons[vertex_key(v)] = v
        for v in deleted_vertices:
            accumulated_singletons[vertex_key(v)] = v

        # Track cycle statistics
        cycle_elapsed = time.time() - cycle_start_time
        cycle_stats.append({
            'cycle': cycle_num,
            'edges_generated': edges_count,
            'components_accepted': len(accepted_comps),
            'vertices_in_components': vertices_in_accepted,
            'vertices_deleted': deleted_count,
            'unpaired_vertices': unpaired_count,
            'elapsed_seconds': round(cycle_elapsed, 2)
        })

        # Step 5: Build combined input for next cycle (deleted + unpaired), deduplicated
        # Deduplicate by stable (filename, x_center, y_center) key — deleted and unpaired
        # are disjoint within a cycle but guard here in case of unexpected overlap.
        combined_seen: Dict = {}
        for v in deleted_vertices + unpaired_vertices_list:
            combined_seen[vertex_key(v)] = v
        combined_vertices = list(combined_seen.values())
        combined_count = len(combined_vertices)
        raw_count = deleted_count + unpaired_count
        if combined_count != raw_count:
            print(f"  [Cycle {cycle_num}] Deduplicated combined input: {raw_count:,} → {combined_count:,} vertices")
        else:
            print(f"  [Cycle {cycle_num}] Combined input for next cycle: {combined_count:,} "
                  f"({deleted_count:,} deleted + {unpaired_count:,} unpaired)")

        # Step 6: Convergence check — stop after two consecutive cycles with no new components
        if len(accepted_comps) == 0:
            consecutive_empty_cycles += 1
            print(f"  [Cycle {cycle_num}] No new components accepted ({consecutive_empty_cycles}/2 consecutive)")
        else:
            consecutive_empty_cycles = 0

        if consecutive_empty_cycles >= 2:
            print(f"\n  Convergence: 2 consecutive cycles produced no new components")
            for v in combined_vertices:
                accumulated_singletons[vertex_key(v)] = v
            break

        if combined_count == 0:
            print(f"\n  Convergence: No vertices to re-pair after cycle {cycle_num}")
            break

        # Write combined vertices JSON (uses 'deleted_vertices' key — same format the pairing script expects)
        combined_json = os.path.join(cycle_dir, 'combined_vertices_for_next_cycle.json')
        save_json({'deleted_vertices': combined_vertices}, combined_json)

        # Step 7: Prepare next cycle
        current_input_type = 'vertices_json'
        current_input = combined_json
        cycle_num += 1

        # Safety check: max cycles
        if cycle_num > config.get('max_cycles', 10):
            print(f"\n  Max cycles ({config.get('max_cycles', 10)}) reached")
            for v in combined_vertices:
                accumulated_singletons[vertex_key(v)] = v
            break
    
    # Build singleton components from all accumulated unpaired/deleted vertices,
    # skipping any that were already accepted into a paired component.
    paired_keys = set()
    for comp in final_components['components'].values():
        for v in comp.get('vertices', []):
            paired_keys.add(vertex_key(v))

    singleton_count = 0
    for key, v in accumulated_singletons.items():
        if key in paired_keys:
            continue
        comp_id = f"component_{next_component_id}"
        next_component_id += 1
        singleton_count += 1
        final_components['components'][comp_id] = {
            'component_id': comp_id,
            'cycle': 'singleton',
            'num_vertices': 1,
            'pep_ident': v.get('pep_ident'),
            'prot_ident': v.get('prot_ident'),
            'vertices': [v]
        }

    print(f"\n  Singleton components created: {singleton_count:,}")
    final_components['unpaired_singletons'] = []  # now represented as components above

    # Update final metadata
    final_components['metadata']['total_cycles'] = cycle_num - 1 if edges_count > 0 else cycle_num
    final_components['metadata']['total_components'] = len(final_components['components'])
    final_components['metadata']['total_vertices_in_components'] = sum(
        c.get('num_vertices', 0) for c in final_components['components'].values()
    )
    final_components['metadata']['total_unpaired_singletons'] = len(final_components['unpaired_singletons'])
    
    # Save final results
    final_json_path = os.path.join(output_dir, 'final_components_all_cycles.json')
    save_json(final_components, final_json_path)
    print(f"\n✓ Saved final components: {final_json_path}")
    
    # Save cycle summary
    total_unpaired_accumulated = sum(c.get('unpaired_vertices', 0) for c in cycle_stats)
    summary_json_path = os.path.join(output_dir, 'cycle_history.json')
    save_json({
        'metadata': {
            'initial_vertices': initial_vertices,
            'total_unpaired_accumulated': total_unpaired_accumulated,
        },
        'cycles': cycle_stats,
    }, summary_json_path)
    print(f"✓ Saved cycle summary: {summary_json_path}")
    
    # Print final summary
    print(f"\n{'='*70}")
    print("CONVERGENCE SUMMARY")
    print(f"{'='*70}")
    print(f"  Total cycles: {final_components['metadata']['total_cycles']}")
    print(f"  Total components: {final_components['metadata']['total_components']:,}")
    print(f"  Vertices in components: {final_components['metadata']['total_vertices_in_components']:,}")
    print(f"  Unpaired singletons: {final_components['metadata']['total_unpaired_singletons']:,}")
    print(f"{'='*70}\n")
    
    return final_components['metadata']


def main():
    parser = argparse.ArgumentParser(
        description="Orchestrate clusterless iterative component building"
    )
    
    # Input
    parser.add_argument("--input_files", nargs='+', required=True,
                       help="Initial TSV feature files for cycle 1")
    parser.add_argument("--output_dir", default='.',
                       help="Output directory for all results (default: current directory)")
    
    # Cutoffs
    parser.add_argument("--mz_cutoff", type=float, default=None,
                       help="m/z (x-axis) coordinate cutoff for edge filtering (mutually exclusive with --ppm)")
    parser.add_argument("--ppm", type=float, default=None,
                       help="m/z cutoff as parts per million: cutoff = feature_mz * ppm / 1000000 (mutually exclusive with --mz_cutoff)")
    parser.add_argument("--rt_cutoff", type=float, default=None,
                       help="RT (y-axis) coordinate cutoff for edge filtering")
    parser.add_argument("--edges_cutoff", type=float, default=None,
                       help="Euclidean distance cutoff for edge filtering")
    
    # RT correction
    parser.add_argument("--rt_correction_json", type=str, default=None,
                       help="JSON file with RT correction models (loess/polynomial/spline methods)")
    parser.add_argument("--rt_correction_mode", type=str, default='none',
                       choices=['none', 'pre-kdtree', 'per-distance'],
                       help="RT correction mode")
    parser.add_argument("--rt_source", type=str, default='y_center',
                       help="RT field to use for correction")
    parser.add_argument("--rt_alignment_method", type=str, default='loess',
                       help="RT alignment method: 'loess', 'polynomial', 'spline', 'rbf', 'piecewise', 'aligntree'")
    parser.add_argument("--trafoxmls_dir", type=str, default=None,
                       help="Directory containing trafoXML files for aligntree method")
    
    # Pairing parameters
    parser.add_argument("--normalize_coordinates", type=lambda x: x.lower() in ('true', '1', 'yes'),
                       default=None, help="Normalize coordinates before pairing")
    parser.add_argument("--distance_calc_before_scaling", action='store_true',
                       help="Calculate distance before coordinate scaling")
    parser.add_argument("--mz_rt_weight_ratio", type=float, default=1.0,
                       help="Weight ratio for RT relative to m/z")
    
    # Filtering
    parser.add_argument("--filter_method", type=str, default='simplified-modularity',
                       choices=['simplified-modularity', 'modularity'],
                       help="Method for filtering duplicate vertices")
    parser.add_argument("--mega_complexity_threshold", type=float, default=5000,
                       help="Complexity threshold for heuristic filtering")
    parser.add_argument("--disable_multiprocessing", action='store_true',
                       help="Disable multiprocessing in filtering")
    
    # Cycle control
    parser.add_argument("--max_cycles", type=int, default=10,
                       help="Maximum number of cycles before stopping (safety limit)")
    
    args = parser.parse_args()

    if args.mz_cutoff is not None and args.ppm is not None:
        print("ERROR: --mz_cutoff and --ppm are mutually exclusive. Provide only one m/z cutoff parameter.")
        sys.exit(1)

    # Build configuration dictionary
    config = {
        'mz_cutoff': args.mz_cutoff,
        'ppm': args.ppm,
        'rt_cutoff': args.rt_cutoff,
        'edges_cutoff': args.edges_cutoff,
        'rt_correction_json': args.rt_correction_json,
        'rt_correction_mode': args.rt_correction_mode,
        'rt_source': args.rt_source,
        'rt_alignment_method': args.rt_alignment_method,
        'trafoxmls_dir': args.trafoxmls_dir,
        'normalize_coordinates': args.normalize_coordinates,
        'distance_calc_before_scaling': args.distance_calc_before_scaling,
        'mz_rt_weight_ratio': args.mz_rt_weight_ratio,
        'filter_duplicate_file_vertices': True,  # Always enabled
        'filter_method': args.filter_method,
        'mega_complexity_threshold': args.mega_complexity_threshold,
        'disable_multiprocessing': args.disable_multiprocessing,
        'max_cycles': args.max_cycles,
    }
    
    # Determine script directory (same as this script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run orchestration
    try:
        metadata = orchestrate_cycles(
            initial_input_files=args.input_files,
            output_dir=args.output_dir,
            config=config,
            script_dir=script_dir
        )
        print("\n✓ Orchestration completed successfully")
        return 0
    except Exception as e:
        print(f"\n✗ Orchestration failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
