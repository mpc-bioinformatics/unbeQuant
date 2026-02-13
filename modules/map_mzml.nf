#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.map_mzml_input_dir = "$PWD/mzmls"  // Directory containing .mzML files
params.map_mzml_outdir = "$PWD/results"  // Output-Directory where the heatmap images will be stored
params.map_mzml_export_data = "true"  // Boolean, if true, will export data into map_mzml_outdir

// Optional Parameters
params.map_mzml_round_up_to = 2  // Number of decimal places to round m/z values //Carefull when adjusting this parameter, every increase of the used decimals increases the used memory tenfold. Check if you have enough RAM first. 
params.map_mzml_log_scale = true  // If True, apply logarithmic scaling to intensities
params.map_mzml_scale_colors = true  // If True, scale colors based on min/max intensities
params.map_mzml_invert_colors = true  // If True, invert color scale (high intensity = dark)
params.map_mzml_add_features = false  // If True, add feature markers to heatmap (not yet implemented)

// Standalone Workflow
workflow {
    // Map mzML files to heatmaps
    mzml_files = Channel.fromPath(params.map_mzml_input_dir + "/*.mzML")
    map_mzml(mzml_files)
}

// Importable Workflow
workflow map_mzml {
    take:
        // Takes .mzML files
        mzml_files
    main:
        // Create RGB intensity dictionary if needed (only runs once)
        if (params.map_mzml_scale_colors) {
            create_rgb_intensity_dict()
            rgb_dict = create_rgb_intensity_dict.out.dict_file
        } else {
            rgb_dict = Channel.empty()
        }
        
        // Map mzML files to heatmap images
        generate_heatmaps(
            mzml_files,
            params.map_mzml_round_up_to,
            params.map_mzml_log_scale,
            params.map_mzml_scale_colors,
            params.map_mzml_invert_colors,
            rgb_dict
        )
    emit:
        // Return each heatmap image and optional npy array
        heatmap_images = generate_heatmaps.out.heatmap_images
        npy_array = generate_heatmaps.out.npy_array
        mz_dicts = generate_heatmaps.out.mz_dicts
        rt_dicts = generate_heatmaps.out.rt_dicts
}


process create_rgb_intensity_dict {
    label 'high_memory'
    container 'python:3.10'
    
    output:
    path "rgb_intensity_dict.pkl", emit: dict_file
    
    script:
    """
    pip install -q numpy 2>/dev/null
    
    python3 << 'EOF'
    import numpy as np
    import pickle
    
    print("Precomputing RGB intensity dictionary...")
    
    def intensity_to_rgb(i):
        r = 0
        g = 0
        b = 0
        
        # Distribute 24 bits round-robin: bits 0,3,6,9... → B, bits 1,4,7,10... → G, bits 2,5,8,11... → R
        for bit_pos in range(24):
            if i & (1 << bit_pos):
                channel = bit_pos % 3
                bit_index = bit_pos // 3
                
                if channel == 0:
                    b |= (1 << bit_index)
                elif channel == 1:
                    g |= (1 << bit_index)
                else:
                    r |= (1 << bit_index)
        
        return (r, g, b)
    
    rgb_intensity_dict = {i: intensity_to_rgb(i) for i in range(16777216)}
    
    with open('rgb_intensity_dict.pkl', 'wb') as f:
        pickle.dump(rgb_intensity_dict, f)
    
    print("RGB intensity dictionary created and saved.")
    EOF
    """
}


process generate_heatmaps {
    label 'high_memory'
    container 'python:3.10'
    
    publishDir "${params.map_mzml_outdir}/", mode:'copy', enabled:"${params.map_mzml_export_data}"

    input:
    path mzml_file
    val round_up_to
    val log_scale
    val scale_colors
    val invert_colors
    path rgb_dict_file, stageAs: 'rgb_dict.pkl'
    
    output:
    path "*.png", emit: heatmap_images
    path "img_frame.npy", emit: npy_array, optional: true
    path "*_mz_dict.pkl", emit: mz_dicts, optional: true
    path "*_rt_dict.pkl", emit: rt_dicts, optional: true
    
    script:
    """
    pip install -q pyopenms numpy pillow psutil matplotlib 2>/dev/null
    
    python3 << 'EOF'
    import pyopenms as oms
    from urllib.request import urlretrieve
    import itertools
    import numpy as np
    from PIL import Image
    from collections import Counter
    import pickle
    import time
    import psutil
    import gc
    import os
    import matplotlib.pyplot as plt
    
    # Configuration
    round_up_to = ${round_up_to}
    testing = False
    log_scale = ${log_scale ? 'True' : 'False'}
    scale_colors = ${scale_colors ? 'True' : 'False'}
    invert_colors = ${invert_colors ? 'True' : 'False'}
    mzml_file = "${mzml_file}"
    
    # Start timer
    start_time = time.time()
    
    # RGB intensity mapping
    if scale_colors == True:
        print("Loading precomputed RGB intensity dictionary...")
        with open('rgb_dict.pkl', 'rb') as f:
            rgb_intensity_dict = pickle.load(f)
        print("RGB intensity dictionary loaded.")
        
        def intesity_to_RGB(rgbint):
            return (rgb_intensity_dict[rgbint])
    else:
        def rgb_int2tuple(rgbint):
            return (rgbint // 256 // 256 % 256, rgbint // 256 % 256, rgbint % 256)
    
    # Load mzML file
    exp = oms.MSExperiment()
    oms.MzMLFile().load(mzml_file, exp)
    
    spec_total = []
    RTINSECONDS_arr = []
    mz_arr = []
    counter = 0
    
    for spec in exp:
        if spec.getMSLevel() == 1:
            spec_curr = []
            RT = spec.getRT()
            mz, intensity = spec.get_peaks()
            counter += 1
            RTINSECONDS_arr.append(RT)
            print(f"Spectrum {counter}", end='\\r')
            if len(mz) == len(intensity) and len(mz) > 0:
                for i in range(0, len(mz)):
                    mz_i = mz[i]
                    spec_curr = [RT, mz_i, intensity[i]]
                    spec_total.append(spec_curr)
                    mz_arr.append(mz_i)
            else:
                print(f"Error: Mismatch in lengths of m/z and intensity arrays or empty spectrum at RT = {RT}")
                exit()
    
    spec_total = np.array(spec_total)
    mz_arr = np.array(mz_arr)
    mz_arr = np.round(mz_arr, round_up_to)
    mz_total_num = (max(mz_arr) - min(mz_arr)) * (10 ** round_up_to)
    mz_total_arr = np.arange(min(mz_arr), max(mz_arr) + (1 / (10 ** round_up_to)), (1 / (10 ** round_up_to)))
    mz_total_arr = np.round(mz_total_arr, round_up_to)
    mz_total_arr = list(mz_total_arr)
    
    # Calculate predicted image size in bytes
    img_height = int((len(RTINSECONDS_arr) + 10))
    img_width = int(mz_total_num) + 10
    img_size_bytes = img_height * img_width * 2
    img_size_gb = img_size_bytes / (1024 ** 3)
    
    # Get available RAM in bytes
    available_ram_bytes = psutil.virtual_memory().available
    available_ram_gb = available_ram_bytes / (1024 ** 3)
    
    # Use memmap only if predicted size exceeds 80% of available RAM
    use_memmap = img_size_bytes > (available_ram_bytes * 0.8)
    
    print(f"Predicted image size: {img_size_gb:.2f} GB")
    print(f"Available RAM: {available_ram_gb:.2f} GB")
    print(f"Using memmap (disk-based): {use_memmap}")
    
    # Create array (in RAM or on disk depending on available memory)
    if use_memmap:
        img_filepath = "img_frame.npy"
        img_frame = np.memmap(img_filepath, dtype=np.uint8, mode='w+',
                                shape=(img_height, img_width, 3))
        img_frame[:] = 255
    else:
        img_frame = np.zeros((img_height, img_width, 3), dtype=np.uint8)
        img_frame[:] = 255
    
    counter = 0
    min_spec = spec_total[:, 2].min()
    max_spec = spec_total[:, 2].max()
    spec_diff = max_spec - min_spec
    
    # Apply logarithmic normalization to handle extreme values
    if log_scale == True:
        log_min = np.log1p(min_spec)
        log_max = np.log1p(max_spec)
        log_diff = log_max - log_min
    else:
        log_min = min_spec
        log_max = max_spec
        log_diff = spec_diff
    
    # Determine base filename once for downstream artifacts
    base_filename = os.path.splitext(os.path.basename(mzml_file))[0]

    # Create lookup dictionaries for O(1) access
    mz_dict = {v: i for i, v in enumerate(mz_total_arr)}
    rt_dict = {v: i for i, v in enumerate(RTINSECONDS_arr)}
    
    # Save dictionaries as pickle files for downstream reuse
    print("Saving coordinate dictionaries...")
    with open(f'{base_filename}_mz_dict.pkl', 'wb') as f:
        pickle.dump(mz_dict, f)
    with open(f'{base_filename}_rt_dict.pkl', 'wb') as f:
        pickle.dump(rt_dict, f)
    print(f"Dictionaries saved: {base_filename}_mz_dict.pkl, {base_filename}_rt_dict.pkl")


    for i in range(0, len(spec_total) - 1):
        counter += 1
        index_mz = mz_dict[round(spec_total[i][1], round_up_to)]
        index_sec = rt_dict[spec_total[i][0]]
        
        if log_scale == True:
            log_intensity = np.log1p(spec_total[i][2])
        else:
            log_intensity = spec_total[i][2]
        
        if invert_colors == True:
            norm_log_intensity = 1 - ((log_intensity - log_min) / log_diff)
        else:
            norm_log_intensity = ((log_intensity - log_min) / log_diff)
        
        log_intesity_as_RGB = intesity_to_RGB(int(round((norm_log_intensity * 16777215))))
        img_frame[index_sec][index_mz] = log_intesity_as_RGB
        
        if counter % 10000 == 0:
            print(f"Mapping spectrum point {counter} out of {len(spec_total)}", end='\\r')
    
    # Flush memmap to disk if using disk-based storage
    if use_memmap:
        img_frame.flush()
    
    # Save as RGB PNG image
    base_filename = os.path.splitext(os.path.basename(mzml_file))[0]
    image_filepath = f"{base_filename}.png"
    print(f"\\nSaving image to {image_filepath}...")
    
    if use_memmap:
        # For memmap, save chunk-by-chunk to avoid loading entire array
        img_pil = Image.new('RGB', (img_width, img_height))
        pixels = img_pil.load()
        
        for row_idx in range(img_height):
            if row_idx % 100 == 0:
                print(f"  Saving row {row_idx}/{img_height}", end='\\r')
            for col_idx in range(img_width):
                pixel = tuple(img_frame[row_idx, col_idx])
                pixels[col_idx, row_idx] = pixel
        
        img_pil.save(image_filepath)
    else:
        # For in-memory array, direct save is fine
        img_pil = Image.fromarray(img_frame, mode='RGB')
        img_pil.save(image_filepath, format='PNG')
    
    print(f"\\nimg_frame saved to {image_filepath}")
    
    # Cleanup
    del spec_total
    del mz_arr
    del RTINSECONDS_arr
    del mz_total_arr
    del exp
    del img_frame
    plt.close('all')
    
    # Manually trigger garbage collection
    collected = gc.collect()
    print(f"Garbage collector collected {collected} objects.")
    
    # End timer
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"\\nScript execution time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
    EOF
    """
}
