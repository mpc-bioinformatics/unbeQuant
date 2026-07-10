#!/usr/bin/env python3

"""
Convert components from clusterless JSON format to ConsensusXML format.

This script takes components from final_components_all_cycles.json and converts them
to OpenMS ConsensusXML format compatible with downstream tools.

Each component becomes a consensusElement with:
- ID: component_id
- Centroid: average of all vertices in the component
- groupedElementList: references to all feature vertices in the component
"""

import json
import argparse
import sys
from datetime import datetime
import statistics
from typing import Dict, List


def bool_to_int(value: bool) -> int:
    """Convert boolean to integer for XML."""
    return 1 if value else 0


def load_json_components(json_file: str) -> Dict:
    """Load component data from JSON file."""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"✗ ERROR: Input JSON not found: {json_file}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"✗ ERROR: Failed to parse input JSON: {e}", file=sys.stderr)
        sys.exit(1)
    
    return data


def extract_unique_filenames(data: Dict) -> tuple[Dict[str, int], Dict[int, str]]:
    """
    Extract unique filenames from components and create bidirectional mapping.
    
    Returns:
        (filename_to_id, id_to_filename) - bidirectional mapping dictionaries
    """
    filenames = set()
    
    for component in data.get('components', {}).values():
        for vertex in component.get('vertices', []):
            filename = vertex.get('filename')
            if filename:
                filenames.add(filename)
    
    # Create sorted mapping
    filename_list = sorted(list(filenames))
    filename_to_id = {name: idx for idx, name in enumerate(filename_list)}
    id_to_filename = {idx: name for idx, name in enumerate(filename_list)}
    
    return filename_to_id, id_to_filename


def generate_unique_id(filename: str, map_id: int) -> str:
    """Generate a unique ID for each map entry."""
    combined = f"{filename}_{map_id}".encode()
    hash_val = hash(combined) & 0x7FFFFFFFFFFFFFFF  # Get positive 64-bit value
    return str(hash_val)


def create_consensus_xml(
    data: Dict,
    filename_to_id: Dict[str, int],
    id_to_filename: Dict[int, str],
    output_file: str,
    input_json_file: str = None,
    workflow_params: Dict = None
):
    """
    Create ConsensusXML content from component data.
    
    Args:
        data: Component data dictionary from final_components_all_cycles.json
        filename_to_id: Mapping of filenames to indices
        id_to_filename: Mapping of indices to filenames
        output_file: Output XML file path
        input_json_file: Input JSON file path
        workflow_params: Dictionary of workflow parameters to include in header
    """
    
    if workflow_params is None:
        workflow_params = {}
    
    # Collect statistics
    components = list(data.get('components', {}).values())
    total_components = len(components)
    total_vertices = sum(comp.get('num_vertices', 0) for comp in components)

    # Calculate component size statistics
    component_sizes = [comp.get('num_vertices', 0) for comp in components]
    avg_component_size = statistics.mean(component_sizes) if component_sizes else 0.0
    min_component_size = min(component_sizes) if component_sizes else 0
    max_component_size = max(component_sizes) if component_sizes else 0
    
    # Start building XML
    xml_lines = []
    
    # XML header
    xml_lines.append('<?xml version="1.0" encoding="ISO-8859-1"?>')
    xml_lines.append('<?xml-stylesheet type="text/xsl" href="https://www.openms.de/xml-stylesheet/ConsensusXML.xsl" ?>')
    
    # Root element with namespace
    xml_lines.append(
        '<consensusXML version="1.7" id="cm_clusterless_components" '
        'experiment_type="label-free" '
        'xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
    )
    
    # Data processing section with parameters
    timestamp = datetime.now().isoformat()
    xml_lines.append(f'    <dataProcessing completion_time="{timestamp}">')
    xml_lines.append('        <software name="clusters_json_to_consensusxml_clusterless" version="1.0" />')
    xml_lines.append('        <processingAction name="Component conversion from clusterless JSON" />')
    
    # Add input/output file parameters
    if input_json_file:
        xml_lines.append(f'        <UserParam type="string" name="parameter: in" value="{input_json_file}"/>')
    xml_lines.append(f'        <UserParam type="string" name="parameter: out" value="{output_file}"/>')
    
    # Add workflow parameters - Feature pairing
    xml_lines.append('        <!-- Feature pairing parameters -->')
    xml_lines.append(f'        <UserParam type="float" name="mmf_mz_cutoff" value="{workflow_params.get("mmf_mz_cutoff", 0.2)}"/>')
    xml_lines.append(f'        <UserParam type="float" name="mmf_rt_cutoff" value="{workflow_params.get("mmf_rt_cutoff", 100)}"/>')
    
    edges_cutoff = workflow_params.get('mmf_edges_cutoff')
    if edges_cutoff is not None:
        xml_lines.append(f'        <UserParam type="float" name="mmf_edges_cutoff" value="{edges_cutoff}"/>')
    
    # Add workflow parameters - Filtering
    xml_lines.append('        <!-- Filtering parameters -->')
    xml_lines.append(f'        <UserParam type="int" name="mmf_filter_duplicate_file_vertices" value="{bool_to_int(workflow_params.get("mmf_filter_duplicate_file_vertices", True))}"/>')
    xml_lines.append(f'        <UserParam type="string" name="mmf_filter_method" value="{workflow_params.get("mmf_filter_method", "simplified-modularity")}"/>')
    
    # Add component statistics
    xml_lines.append('        <!-- Component statistics -->')
    xml_lines.append(f'        <UserParam type="int" name="total_components" value="{total_components}"/>')
    xml_lines.append(f'        <UserParam type="int" name="total_vertices" value="{total_vertices}"/>')
    xml_lines.append(f'        <UserParam type="float" name="avg_component_size" value="{avg_component_size:.2f}"/>')
    xml_lines.append(f'        <UserParam type="int" name="min_component_size" value="{min_component_size}"/>')
    xml_lines.append(f'        <UserParam type="int" name="max_component_size" value="{max_component_size}"/>')
    
    xml_lines.append('    </dataProcessing>')
    
    # IdentificationRun section (placeholder - required for valid ConsensusXML)
    xml_lines.append('    <IdentificationRun id="PI_0" date="{}" search_engine="" search_engine_version="">'.format(timestamp))
    xml_lines.append('        <SearchParameters db="" db_version="" taxonomy="" mass_type="monoisotopic" charges="" enzyme="unknown_enzyme" missed_cleavages="0" precursor_peak_tolerance="0" precursor_peak_tolerance_ppm="false" peak_mass_tolerance="0" peak_mass_tolerance_ppm="false" >')
    xml_lines.append('        </SearchParameters>')
    xml_lines.append('        <ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0">')
    xml_lines.append('        </ProteinIdentification>')
    xml_lines.append('    </IdentificationRun>')
    
    # mapList section - list all unique files
    xml_lines.append('    <mapList count="{}">'.format(len(filename_to_id)))
    for filename, file_id in sorted(filename_to_id.items(), key=lambda x: x[1]):
        unique_id = generate_unique_id(filename, file_id)
        xml_lines.append(f'        <map id="{file_id}" name="{filename}" filename="{filename}" unique_id="{unique_id}" label="" size="0">')
        xml_lines.append('        </map>')
    xml_lines.append('    </mapList>')
    
    # consensusElementList section - one element per component
    xml_lines.append(f'    <consensusElementList>')
    
    consensus_id = 0
    for component in components:
        component_id = component.get('component_id', 'unknown')
        vertices = component.get('vertices', [])
        
        if not vertices:
            continue  # Skip empty components
        
        # Calculate centroid (average m/z and RT)
        mz_values = [v.get('x_center', v.get('mz', 0)) for v in vertices]
        rt_values = [v.get('y_center', v.get('rt', 0)) for v in vertices]

        centroid_mz = statistics.mean(mz_values) if mz_values else 0.0
        centroid_rt = statistics.mean(rt_values) if rt_values else 0.0
        centroid_intensity = sum(v.get('intensity', 0) for v in vertices)
        
        # Quality measure (average quality of vertices if available)
        qualities = [v.get('quality', 1.0) for v in vertices if 'quality' in v]
        avg_quality = statistics.mean(qualities) if qualities else 1.0
        
        # Start consensusElement
        xml_lines.append(
            f'        <consensusElement id="{consensus_id}" '
            f'quality="{avg_quality:.6f}" '
            f'charge="0">'
        )
        
        # Centroid coordinates
        xml_lines.append(f'            <centroid rt="{centroid_rt:.6f}" mz="{centroid_mz:.6f}" it="{centroid_intensity:.2f}"/>')
        
        # Grouped elements (vertices in this component)
        xml_lines.append(f'            <groupedElementList>')
        for vertex in vertices:
            filename = vertex.get('filename')
            if not filename:
                continue
            
            map_id = filename_to_id.get(filename)
            if map_id is None:
                continue
            
            # element id must be the numeric OpenMS unique feature ID so that
            # consensus_and_features_to_tsv.py can match it via "f_" + str(id) == openms_fid
            openms_fid = vertex.get('openms_fid', '')
            if openms_fid.startswith('f_'):
                feature_id = openms_fid[2:]
            elif openms_fid:
                feature_id = openms_fid
            else:
                feature_id = vertex.get('feature_idx', 0)
            mz = vertex.get('x_center', vertex.get('mz', 0))
            rt = vertex.get('y_center', vertex.get('rt', 0))
            intensity = vertex.get('intensity', 0)
            quality = vertex.get('quality', 1.0)
            charge = vertex.get('charge', 0)
            
            xml_lines.append(
                f'                <element map="{map_id}" '
                f'id="{feature_id}" '
                f'rt="{rt:.6f}" '
                f'mz="{mz:.6f}" '
                f'it="{intensity:.2f}" '
                f'charge="{charge}"/>'
            )
        
        xml_lines.append('            </groupedElementList>')
        
        # Add component metadata as UserParams
        xml_lines.append(f'            <UserParam type="string" name="component_id" value="{component_id}"/>')
        xml_lines.append(f'            <UserParam type="int" name="num_vertices" value="{len(vertices)}"/>')
        
        # Add cycle information if available
        cycle_num = component.get('cycle_num')
        if cycle_num is not None:
            xml_lines.append(f'            <UserParam type="int" name="cycle_num" value="{cycle_num}"/>')
        
        xml_lines.append('        </consensusElement>')
        
        consensus_id += 1
    
    xml_lines.append('    </consensusElementList>')
    
    # Close root element
    xml_lines.append('</consensusXML>')
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(xml_lines))
    
    print(f"✓ Created ConsensusXML with {consensus_id} consensus elements")
    print(f"  Components: {total_components}")
    print(f"  Total vertices: {total_vertices}")
    print(f"  Average component size: {avg_component_size:.2f}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert clusterless components JSON to ConsensusXML format"
    )
    parser.add_argument(
        "--input_json",
        required=True,
        help="Input JSON file with components (final_components_all_cycles.json)"
    )
    parser.add_argument(
        "--output_xml",
        required=True,
        help="Output ConsensusXML file"
    )
    
    # Workflow parameters (optional, for metadata)
    parser.add_argument("--mmf_mz_cutoff", type=float, default=0.2)
    parser.add_argument("--mmf_rt_cutoff", type=float, default=100)
    parser.add_argument("--mmf_edges_cutoff", type=float, default=None)
    parser.add_argument("--mmf_filter_duplicate_file_vertices", type=bool, default=True)
    parser.add_argument("--mmf_filter_method", type=str, default="simplified-modularity")
    
    return parser.parse_args()


def main():
    args = parse_args()
    
    print("\n" + "="*70)
    print("Converting Components to ConsensusXML (Clusterless)")
    print("="*70 + "\n")
    
    # Load component data
    print(f"Loading components from: {args.input_json}")
    data = load_json_components(args.input_json)
    
    # Extract unique filenames
    print("Extracting unique filenames...")
    filename_to_id, id_to_filename = extract_unique_filenames(data)
    print(f"  Found {len(filename_to_id)} unique files")
    
    # Collect workflow parameters
    workflow_params = {
        'mmf_mz_cutoff': args.mmf_mz_cutoff,
        'mmf_rt_cutoff': args.mmf_rt_cutoff,
        'mmf_edges_cutoff': args.mmf_edges_cutoff,
        'mmf_filter_duplicate_file_vertices': args.mmf_filter_duplicate_file_vertices,
        'mmf_filter_method': args.mmf_filter_method,
    }
    
    # Generate ConsensusXML
    print(f"\nGenerating ConsensusXML: {args.output_xml}")
    create_consensus_xml(
        data,
        filename_to_id,
        id_to_filename,
        args.output_xml,
        args.input_json,
        workflow_params
    )
    
    print("\n" + "="*70)
    print("✅ Conversion complete!")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
