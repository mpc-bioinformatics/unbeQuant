#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.map_alignment_heatmaps_dir = "$PWD/results"  // Directory containing the heatmap PNG files
params.map_alignment_features_tsv_dir = "/workspaces/unbeQuant/results/quantifications/features_with_annotated_identifications"  // Directory containing features TSV files
params.map_alignment_outdir = "$PWD/results/heatmaps_with_features"  // Output-Directory where annotated heatmaps will be stored
params.map_alignment_export_data = "true"  // Boolean, if true, will export data into map_alignment_outdir

// Optional Parameters
params.map_alignment_box_color = "red"  // Color of the drawn boxes
params.map_alignment_box_width = 2  // Width of the box lines in pixels
params.map_alignment_box_alpha = 0.7  // Transparency of the boxes (0-1)


// Standalone Workflow
workflow {
    // Get all TSV feature files
    features_tsvs = Channel.fromPath(params.map_alignment_features_tsv_dir + "/*.tsv")
    
    // For each TSV file, find its corresponding heatmap and dictionaries
    heatmap_images = features_tsvs.map { tsv_file ->
        // Extract the base name without extension
        base_name = tsv_file.baseName
        // Look for the corresponding heatmap
        heatmap_file = file(params.map_alignment_heatmaps_dir + "/${base_name}.png")
        if (heatmap_file.exists()) {
            return heatmap_file
        } else {
            error("No corresponding heatmap found for features file: ${tsv_file.name}")
        }
    }
    
    mz_dicts = Channel.fromPath(params.map_alignment_heatmaps_dir + "/*_mz_dict.pkl")
    rt_dicts = Channel.fromPath(params.map_alignment_heatmaps_dir + "/*_rt_dict.pkl")
    
    map_alignment(heatmap_images, features_tsvs, mz_dicts, rt_dicts)
}

// Importable Workflow
workflow map_alignment {
    take:
        // Takes heatmap PNG files from map_mzml
        heatmap_images
        // Takes features TSV files
        features_tsvs
        // Takes mz_dict pickle files
        mz_dicts
        // Takes rt_dict pickle files
        rt_dicts
    main:
        // Create a mapping of heatmap basenames to their files
        heatmap_channel = heatmap_images.map { file -> tuple(file.baseName, file) }
        features_channel = features_tsvs.map { file -> tuple(file.baseName, file) }
        mz_dict_channel = mz_dicts.map { file -> tuple(file.baseName.replaceFirst(/_mz_dict$/, ''), file) }
        rt_dict_channel = rt_dicts.map { file -> tuple(file.baseName.replaceFirst(/_rt_dict$/, ''), file) }

        features_channel = features_channel
            .join(mz_dict_channel, by: 0)
            .join(rt_dict_channel, by: 0)
            .map { tuple(it[0], it[1], it[2], it[3]) }

        // Combine heatmaps with their corresponding features TSV files and dictionaries
        combined = heatmap_channel
            .join(features_channel, by: 0)
        
        // Draw features on heatmaps
        draw_features_on_heatmaps(
            combined.map { tuple(it[1], it[2], it[3], it[4]) },
            params.map_alignment_box_color,
            params.map_alignment_box_width,
            params.map_alignment_box_alpha
        )
    emit:
        // Return annotated heatmap images
        annotated_heatmaps = draw_features_on_heatmaps.out.annotated_images
}


process draw_features_on_heatmaps {
    label 'high_memory'
    container 'python:3.10'
    memory '23 GB'
    cpus 10
    
    publishDir "${params.map_alignment_outdir}/", mode:'copy', enabled:"${params.map_alignment_export_data}"

    input:
    tuple path(heatmap_image), path(features_tsv), path(mz_dict_file), path(rt_dict_file)
    val box_color
    val box_width
    val box_alpha
    
    output:
    path "*_with_features_3.png", emit: annotated_images
    path "*_interactive.html", emit: interactive_html
    
    script:
    """
    pip install -q pillow numpy pandas matplotlib base64 2>/dev/null
    
    python3 << 'EOF'
    import pandas as pd
    import numpy as np
    from PIL import Image, ImageDraw
    import ast
    import os
    import gc
    import pickle
    import base64
    
    # Disable PIL's decompression bomb check for large images
    Image.MAX_IMAGE_PIXELS = None
    
    heatmap_file = "${heatmap_image}"
    features_file = "${features_tsv}"
    mz_dict_file = "${mz_dict_file}"
    rt_dict_file = "${rt_dict_file}"
    
    # Resolve any symlinks to absolute paths
    import os
    heatmap_file = os.path.realpath(heatmap_file)
    features_file = os.path.realpath(features_file)
    mz_dict_file = os.path.realpath(mz_dict_file)
    rt_dict_file = os.path.realpath(rt_dict_file)
    
    print(f"Processing heatmap: {heatmap_file}")
    print(f"Loading features from: {features_file}")
    print(f"Loading dictionaries from: {mz_dict_file}, {rt_dict_file}")
    
    # Load the coordinate dictionaries from map_mzml
    print("Loading coordinate dictionaries...")
    with open(mz_dict_file, 'rb') as f:
        mz_dict = pickle.load(f)
    with open(rt_dict_file, 'rb') as f:
        rt_dict = pickle.load(f)
    print(f"Loaded mz_dict with {len(mz_dict)} entries")
    print(f"Loaded rt_dict with {len(rt_dict)} entries")

    rt_total_arr = list(rt_dict)

    # Round RT values to 8 decimal places to ensure consistency
    rt_dict = {round(rt, 8): y for rt, y in rt_dict.items()}
    
    # Create reverse dictionaries for pixel-to-coordinate mapping
    reverse_mz_dict = {v: k for k, v in mz_dict.items()}
    reverse_rt_dict = {v: k for k, v in rt_dict.items()}

    # Load the features TSV data
    features_data = pd.read_csv(features_file, sep='\t')
    print(f"Loaded {len(features_data)} features from TSV data")
    
    # Function to parse list strings from TSV
    def parse_list(s):
        if isinstance(s, str):
            try:
                parsed = ast.literal_eval(s)
                # Convert strings to floats
                return [float(x) for x in parsed]
            except:
                return []
        return []
    
    # Load the heatmap image
    print("Loading heatmap image...")
    heatmap = Image.open(heatmap_file)
    img_width, img_height = heatmap.size
    print(f"Heatmap size: {img_width} x {img_height}")
    
    # Convert to RGB directly to avoid keeping multiple color spaces in memory
    heatmap_rgb = heatmap.convert('RGB')
    heatmap = None  # Release original image
    gc.collect()
    
    # Create overlay only for drawing
    overlay = Image.new('RGBA', (img_width, img_height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    
    # Parse box color
    box_alpha = ${box_alpha}
    box_color = "${box_color}"
    box_width = ${box_width}
    
    color_map = {
        'red': (255, 0, 0, int(255 * box_alpha)),
        'blue': (0, 0, 255, int(255 * box_alpha)),
        'green': (0, 255, 0, int(255 * box_alpha)),
        'yellow': (255, 255, 0, int(255 * box_alpha)),
        'cyan': (0, 255, 255, int(255 * box_alpha)),
        'magenta': (255, 0, 255, int(255 * box_alpha)),
        'white': (255, 255, 255, int(255 * box_alpha)),
        'black': (0, 0, 0, int(255 * box_alpha))
    }
    box_color_rgba = color_map.get(box_color.lower(), (255, 0, 0, int(255 * box_alpha)))
    
    # Store feature boxes for interactive HTML
    feature_boxes = []
    
    # Draw boxes for each feature
    print("Drawing feature boxes...")
    box_count = 0
    skipped_count = 0
    for idx, row in features_data.iterrows():
        mz_starts = parse_list(row['l_mz_start'])
        mz_ends = parse_list(row['l_mz_end'])
        rt_starts = parse_list(row['l_rt_start'])
        rt_ends = parse_list(row['l_rt_end'])
        
        # Get the min and max values for this feature
        if mz_starts and mz_ends and rt_starts and rt_ends:
            mz_start = min(mz_starts + mz_ends)
            mz_end = max(mz_starts + mz_ends)
            rt_start = min(rt_starts + rt_ends)
            rt_end = max(rt_starts + rt_ends)
            
            # Round to matching precision (same as heatmap generation)
            # Determine precision by checking dictionary key format
            mz_precision = len(str(list(mz_dict.keys())[0]).split('.')[-1]) if mz_dict else 2
            
            mz_start_rounded = round(mz_start, mz_precision)
            mz_end_rounded = round(mz_end, mz_precision)
            
            # Look up pixel coordinates using dictionaries from map_mzml
            try:
                x_start = mz_dict[mz_start_rounded]
                x_end = mz_dict[mz_end_rounded]
                y_start = rt_dict[round(rt_start, 8)]
                y_end = rt_dict[round(rt_end, 8)]
                
                # Draw rectangle (box)
                draw.rectangle(
                    [(x_start, y_start), (x_end, y_end)],
                    outline=box_color_rgba,
                    width=box_width
                )
                
                # Store feature box info for interactive HTML
                feature_boxes.append({
                    'x_start': int(min(x_start, x_end)),
                    'x_end': int(max(x_start, x_end)),
                    'y_start': int(min(y_start, y_end)),
                    'y_end': int(max(y_start, y_end)),
                    'mz_start': mz_start,
                    'mz_end': mz_end,
                    'rt_start': rt_start,
                    'rt_end': rt_end
                })
                box_count += 1
                
                if box_count % 100 == 0:
                    print(f"  Drew {box_count} boxes...", end='\\r')
            except KeyError as e:
                if skipped_count < 5:
                    print(f"Warning: Feature {idx} coordinate out of bounds: {e}")
                skipped_count += 1
                continue
        else:
            skipped_count += 1
    
    print(f"\\nFinished drawing {box_count} feature boxes (skipped {skipped_count})")
    
    # Composite the overlay onto the RGB image
    print("Compositing overlay onto heatmap...")
    heatmap_rgba = heatmap_rgb.convert('RGBA')
    heatmap_rgb = None  # Release RGB version
    gc.collect()
    
    result = Image.alpha_composite(heatmap_rgba, overlay)
    overlay = None  # Release overlay
    heatmap_rgba = None  # Release RGBA version
    gc.collect()
    
    result_rgb = result.convert('RGB')
    result = None  # Release RGBA result
    gc.collect()
    
    # Save the annotated image
    print("Saving annotated heatmap...")
    output_filename = heatmap_file.replace('.png', '_with_features_3.png')
    result_rgb.save(output_filename, format='PNG', optimize=False)
    print(f"Saved annotated heatmap to {output_filename}")
    ################################################
    
    """