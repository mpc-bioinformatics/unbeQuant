#!/usr/bin/env python3
"""
Create heatmap image from processed mzML spectral data.
Refactored from map_mzml_batch_tsv.py
"""

import argparse
import pickle
import numpy as np
import os
import gc
import h5py
from PIL import Image
from typing import Dict, List

#TODO: Create an 3D-Heatmap
def initialize_rgb_dict(scale_colors: bool = True):
    """Precompute RGB dictionary once for all files."""
    if scale_colors:
        print("Precomputing RGB intensity dictionary...")
        
        def intensity_to_rgb(i):
            r = g = b = 0
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
        
        rgb_dict = {i: intensity_to_rgb(i) for i in range(16777216)}
        print("RGB dictionary precomputed and ready")
        return rgb_dict
    return {}


def intensity_to_rgb_lookup(intensity_int: int, rgb_dict: Dict, scale_colors: bool = True):
    """Look up RGB values using precomputed dictionary."""
    if scale_colors and rgb_dict:
        return rgb_dict[intensity_int]
    return rgb_int2tuple(intensity_int)


def rgb_int2tuple(rgbint: int):
    """Convert RGB integer to tuple for non-scaled colors."""
    return (rgbint // 256 // 256 % 256, rgbint // 256 % 256, rgbint % 256)


def save_heatmap_to_hdf5(img_frame: np.ndarray, img_raw: np.ndarray, output_path: str):
    """Save heatmap image data to HDF5 format for efficient storage and retrieval."""
    print(f"  Saving heatmap to HDF5: {os.path.basename(output_path)}")
    with h5py.File(output_path, 'w') as hf:
        # Store RGB image data with compression
        hf.create_dataset('img_frame', data=img_frame, compression='gzip', compression_opts=4, dtype=np.uint8)
        # Store raw intensity data with compression
        hf.create_dataset('img_raw', data=img_raw, compression='gzip', compression_opts=4, dtype=np.float32)
        # Store metadata
        hf.attrs['format'] = 'heatmap_image_v1'
        hf.attrs['img_height'] = img_frame.shape[0]
        hf.attrs['img_width'] = img_frame.shape[1]


def load_heatmap_from_hdf5(hdf5_path: str):
    """Load heatmap image data from HDF5 format."""
    with h5py.File(hdf5_path, 'r') as hf:
        img_frame = hf['img_frame'][:]
        img_raw = hf['img_raw'][:]
    return img_frame, img_raw


def export_heatmap_to_png(hdf5_path: str, png_output_path: str):
    """Export heatmap from HDF5 to PNG format."""
    img_frame, _ = load_heatmap_from_hdf5(hdf5_path)
    print(f"  Exporting PNG from HDF5: {os.path.basename(png_output_path)}")
    img_pil = Image.fromarray(img_frame, mode='RGB')
    img_pil.save(png_output_path, format='PNG')
    img_pil = None
    gc.collect()


def create_heatmap_image(spec_total: np.ndarray, RTINSECONDS_arr: List[float], 
                         mz_total_arr: List[float], mz_dict: Dict, rt_dict: Dict,
                         output_path: str, log_scale: bool = True, scale_colors: bool = True,
                         invert_colors: bool = False, save_png: bool = False):
    """
    Create and save heatmap image directly to HDF5 datasets for minimal memory usage.
    Writes data incrementally to HDF5 instead of creating large arrays in memory.
    Args:
        save_png: If True, also exports a PNG copy from the HDF5 file
    Returns: HDF5 file path
    """
    print(f"  Creating heatmap image...")
    
    rgb_dict = initialize_rgb_dict(scale_colors)
    
    img_height = int(len(RTINSECONDS_arr) + 10)
    mz_total_num = (max(mz_total_arr) - min(mz_total_arr)) * (10 ** 2)  # Assuming round_up_to=2
    img_width = int(mz_total_num) + 10
    
    # Prepare HDF5 output path
    hdf5_output = output_path if output_path.endswith('.h5') else output_path.replace('.png', '.h5')
    
    # Create HDF5 file and datasets directly (no large arrays in memory)
    print(f"  Initializing HDF5 file: {os.path.basename(hdf5_output)}")
    with h5py.File(hdf5_output, 'w') as hf:
        # Create datasets with compression and initialize with default values
        img_frame_dset = hf.create_dataset(
            'img_frame', 
            shape=(img_height, img_width, 3), 
            dtype=np.uint8, 
            compression='gzip', 
            compression_opts=4
        )
        img_raw_dset = hf.create_dataset(
            'img_raw', 
            shape=(img_height, img_width), 
            dtype=np.float32, 
            compression='gzip', 
            compression_opts=4
        )
        
        # Initialize img_frame with white background
        img_frame_dset[:] = 255
        
        # Store metadata
        hf.attrs['format'] = 'heatmap_image_v1'
        hf.attrs['img_height'] = img_height
        hf.attrs['img_width'] = img_width
        
        # Normalize intensity values
        min_spec = spec_total[:, 2].min()
        max_spec = spec_total[:, 2].max()
        
        if log_scale:
            log_min = np.log1p(min_spec)
            log_max = np.log1p(max_spec)
            log_diff = log_max - log_min
        else:
            log_min = min_spec
            log_max = max_spec
            log_diff = max_spec - min_spec
        
        # Map spectrum points to HDF5 datasets incrementally
        print(f"  Mapping {len(spec_total)} spectrum points to HDF5 datasets...")
        for i in range(len(spec_total)):
            if i % 50000 == 0:
                print(f"    Mapped {i}/{len(spec_total)}", end='\r')
            
            try:
                index_mz = mz_dict[round(spec_total[i][1], 2)]
                index_sec = rt_dict[round(spec_total[i][0], 8)]
                intensity_cur = spec_total[i][2]
                
                if log_scale:
                    log_intensity = np.log1p(intensity_cur)
                else:
                    log_intensity = intensity_cur
                
                if invert_colors:
                    norm_log_intensity = 1 - ((log_intensity - log_min) / log_diff)
                else:
                    norm_log_intensity = ((log_intensity - log_min) / log_diff)
                
                rgb_tuple = intensity_to_rgb_lookup(int(round(norm_log_intensity * 16777215)), rgb_dict, scale_colors)
                # Write directly to HDF5 dataset (minimal memory footprint)
                img_frame_dset[index_sec, index_mz] = rgb_tuple
                img_raw_dset[index_sec, index_mz] = intensity_cur
            except (KeyError, IndexError):
                pass
        
        print(f"\n  Successfully saved heatmap to HDF5: {os.path.basename(hdf5_output)}")
    
    # Optionally export PNG from HDF5
    if save_png:
        png_output = output_path if output_path.endswith('.png') else output_path.replace('.h5', '.png')
        export_heatmap_to_png(hdf5_output, png_output)
    
    gc.collect()
    
    return hdf5_output


def main():
    parser = argparse.ArgumentParser(
        description="Create heatmap image from processed mzML spectral data"
    )
    parser.add_argument("--spectrum_pkl", help="Path to spectrum pickle file")
    parser.add_argument("--spectrum_hdf5", help="Path to spectrum HDF5 file")
    parser.add_argument("--output_h5", help="Output HDF5 heatmap file path (recommended)")
    parser.add_argument("--output_png", help="Output PNG heatmap file path (optional)")
    parser.add_argument("--export_png", action="store_true", help="Also export PNG from existing HDF5")
    parser.add_argument("--hdf5_input", help="Input HDF5 to export as PNG")
    parser.add_argument("--log_scale", type=bool, default=True, help="Apply logarithmic scaling")
    parser.add_argument("--scale_colors", type=bool, default=True, help="Scale colors based on intensities")
    parser.add_argument("--invert_colors", type=bool, default=False, help="Invert color scale")
    parser.add_argument("--round_up_to", type=int, default=2, help="Decimal places for m/z rounding")
    
    args = parser.parse_args()
    
    print(f"\n{'='*70}")
    print(f"Creating Heatmap Image")
    print(f"{'='*70}")
    
    # Mode 1: Export existing HDF5 to PNG
    if args.export_png and args.hdf5_input and args.output_png:
        export_heatmap_to_png(args.hdf5_input, args.output_png)
        print(f"Heatmap successfully exported to {args.output_png}")
        return
    
    # Mode 2: Load spectrum data and create heatmap
    # Load spectrum data from pickle or HDF5
    if args.spectrum_pkl:
        print("Loading spectrum data from pickle file...")
        with open(args.spectrum_pkl, 'rb') as f:
            data = pickle.load(f)
    elif args.spectrum_hdf5:
        print("Loading spectrum data from HDF5 file...")
        with h5py.File(args.spectrum_hdf5, 'r') as hf:
            data = {
                'spec_total': hf['spec_total'][:],
                'RTINSECONDS_arr': hf['RTINSECONDS_arr'][:],
                'mz_total_arr': hf['mz_total_arr'][:],
                'mz_dict': dict(hf['mz_dict'].attrs),
                'rt_dict': dict(hf['rt_dict'].attrs)
            }
    else:
        print("Error: Please provide either --spectrum_pkl or --spectrum_hdf5")
        return
    
    spec_total = data['spec_total']
    RTINSECONDS_arr = data['RTINSECONDS_arr']
    mz_total_arr = data['mz_total_arr']
    mz_dict = data['mz_dict']
    rt_dict = data['rt_dict']
    
    # Determine output paths
    if not args.output_h5 and not args.output_png:
        print("Error: Please provide either --output_h5 or --output_png")
        return
    
    # Default to HDF5 if no format specified
    output_path = args.output_h5 if args.output_h5 else args.output_png
    save_png = bool(args.output_png)
    
    # Create and save heatmap
    create_heatmap_image(
        spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict,
        output_path,
        log_scale=args.log_scale,
        scale_colors=args.scale_colors,
        invert_colors=args.invert_colors,
        save_png=save_png
    )
    
    print(f"Heatmap successfully created and saved to {output_path}")


if __name__ == "__main__":
    main()
