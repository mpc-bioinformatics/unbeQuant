#!/usr/bin/env python3
"""
Generate vertex coordinates from paired_features.json for coordinate-based Graphviz layout
Uses m/z (normalized) as X and RT (normalized) as Y coordinates
"""

import json
import sys
from pathlib import Path

def generate_coordinates_from_features(features_json_path, output_coords_json):
    """
    Extract m/z and RT from paired features and generate coordinates
    """
    try:
        with open(features_json_path, 'r') as f:
            data = json.load(f)
        
        # Handle different data structures
        features = data
        if isinstance(data, dict) and 'features' in data:
            features = data['features']
        
        if not isinstance(features, list):
            print(f"❌ Error: Expected list of features, got {type(features)}")
            return False
        
        print(f"✓ Loaded {len(features)} features")
        
        # Extract m/z and RT values
        mz_values = []
        rt_values = []
        coordinates = {}
        
        for i, feature in enumerate(features):
            if isinstance(feature, dict):
                # Try to find m/z and RT fields
                mz = feature.get('mz') or feature.get('m/z') or feature.get('mass_to_charge')
                rt = feature.get('rt') or feature.get('RT') or feature.get('retention_time')
                
                if mz is not None and rt is not None:
                    mz_values.append(float(mz))
                    rt_values.append(float(rt))
                    coordinates[str(i)] = [float(mz), float(rt)]
        
        if not coordinates:
            print(f"❌ Error: No m/z and RT values found in features")
            return False
        
        # Normalize coordinates to range [0, 10]
        min_mz = min(mz_values)
        max_mz = max(mz_values)
        min_rt = min(rt_values)
        max_rt = max(rt_values)
        
        mz_range = max_mz - min_mz if max_mz > min_mz else 1
        rt_range = max_rt - min_rt if max_rt > min_rt else 1
        
        for vid, (mz, rt) in coordinates.items():
            norm_mz = ((mz - min_mz) / mz_range * 10) if mz_range > 0 else 5
            norm_rt = ((rt - min_rt) / rt_range * 10) if rt_range > 0 else 5
            coordinates[vid] = [norm_mz, norm_rt]
        
        # Save coordinates
        output = {
            "description": f"Vertex coordinates for {len(coordinates)} features",
            "source": features_json_path,
            "mz_range": [min_mz, max_mz],
            "rt_range": [min_rt, max_rt],
            "normalized_range": [0, 10],
            "coordinates": coordinates
        }
        
        with open(output_coords_json, 'w') as f:
            json.dump(output, f, indent=2)
        
        print(f"✓ Generated coordinates for {len(coordinates)} vertices")
        print(f"  m/z range: {min_mz:.2f} - {max_mz:.2f}")
        print(f"  RT range: {min_rt:.2f} - {max_rt:.2f}")
        print(f"✓ Saved to: {output_coords_json}")
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <features.json> <output_coordinates.json>")
        print(f"Example: {sys.argv[0]} work/ae/55ff3b4ce9f787e88c38951a6241c0/paired_features.json results/clustering_samples/coordinates.json")
        sys.exit(1)
    
    features_file = sys.argv[1]
    output_file = sys.argv[2]
    
    success = generate_coordinates_from_features(features_file, output_file)
    sys.exit(0 if success else 1)
