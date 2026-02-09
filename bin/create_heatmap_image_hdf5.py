#!/usr/bin/env python3
"""
Create heatmap image from processed mzML spectral data.
Refactored from map_mzml_batch_tsv.py
Optimized for fast processing with vectorized operations and batched HDF5 writes.
"""

import argparse
import pickle
import numpy as np
import os
import gc
import h5py
from PIL import Image
from typing import Dict, List, Tuple

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
        
        # Vectorized RGB conversion
        intensity_array = np.arange(16777216, dtype=np.uint32)
        rgb_array = np.zeros((16777216, 3), dtype=np.uint8)
        
        for bit_pos in range(24):
            mask = (intensity_array >> bit_pos) & 1
            channel = bit_pos % 3
            bit_index = bit_pos // 3
            rgb_array[:, channel] |= (mask << bit_index).astype(np.uint8)
        
        rgb_dict = {i: tuple(rgb_array[i]) for i in range(0, 16777216, 1)}
        print("RGB dictionary precomputed and ready")
        return rgb_array  # Return array for faster lookup
    return None


def intensity_to_rgb_lookup(intensity_int: int, rgb_array: np.ndarray, scale_colors: bool = True):
    """Look up RGB values using precomputed array."""
    if scale_colors and rgb_array is not None:
        return tuple(rgb_array[intensity_int])
    return rgb_int2tuple(intensity_int)


def intensity_to_rgb_vectorized(intensities: np.ndarray, rgb_array: np.ndarray, scale_colors: bool = True):
    """Vectorized RGB lookup for arrays of intensities."""
    if scale_colors and rgb_array is not None:
        return rgb_array[intensities.astype(np.uint32)]
    # Fallback for non-scaled colors
    return np.stack([
        (intensities // 256 // 256) % 256,
        (intensities // 256) % 256,
        intensities % 256
    ], axis=-1).astype(np.uint8)


def rgb_int2tuple(rgbint: int):
    """Convert RGB integer to tuple for non-scaled colors."""
    return (rgbint // 256 // 256 % 256, rgbint // 256 % 256, rgbint % 256)


def save_heatmap_to_hdf5(img_frame: np.ndarray, output_path: str):
    """Save heatmap image data to HDF5 format for efficient storage and retrieval."""
    print(f"  Saving heatmap to HDF5: {os.path.basename(output_path)}")
    with h5py.File(output_path, 'w') as hf:
        # Store RGB image data with compression
        hf.create_dataset('img_frame', data=img_frame, compression='gzip', compression_opts=4, dtype=np.uint8)
        # Store metadata
        hf.attrs['format'] = 'heatmap_image_v1'
        hf.attrs['img_height'] = img_frame.shape[0]
        hf.attrs['img_width'] = img_frame.shape[1]


def load_heatmap_from_hdf5(hdf5_path: str):
    """Load heatmap image data from HDF5 format."""
    with h5py.File(hdf5_path, 'r') as hf:
        img_frame = hf['img_frame'][:]
    return img_frame


def export_heatmap_to_png(hdf5_path: str, png_output_path: str):
    """Export heatmap from HDF5 to PNG format."""
    img_frame = load_heatmap_from_hdf5(hdf5_path)
    print(f"  Exporting PNG from HDF5: {os.path.basename(png_output_path)}")
    img_pil = Image.fromarray(img_frame, mode='RGB')
    img_pil.save(png_output_path, format='PNG')
    img_pil = None
    gc.collect()


def create_heatmap_image(spec_total: np.ndarray, RTINSECONDS_arr: List[float], 
                         mz_total_arr: List[float], mz_dict: Dict, rt_dict: Dict,
                         output_path: str, log_scale: bool = True, scale_colors: bool = True,
                         invert_colors: bool = False, save_png: bool = True, batch_size: int = 100000,
                         row_batch_size: int = 10, compression_level: int = None):
    """
    Create and save heatmap image directly to HDF5 datasets for minimal memory usage.
    Writes data in batches to HDF5 for optimal performance.
    Args:
        batch_size: Deprecated (kept for backwards compatibility)
        row_batch_size: Number of rows to prepare and write together (default: 10)
                       Higher values = faster writes but more RAM (each row ~3.2MB)
                       Lower values = slower writes but less RAM
        compression_level: GZIP compression level (None=no compression, 1=light, 4=default, 9=max)
                         None=fastest (1min, 36MB), 1=balanced (7min, 0.36MB), 4=original (13min, 0.2MB)
        save_png: If True, also exports a PNG copy from the HDF5 file
    Returns: HDF5 file path
    """
    print(f"  Creating heatmap image...")
    
    rgb_array = initialize_rgb_dict(scale_colors)
    
    img_height = int(len(RTINSECONDS_arr) + 10)
    mz_total_num = (max(mz_total_arr) - min(mz_total_arr)) * (10 ** 2)  # Assuming round_up_to=2
    img_width = int(mz_total_num) + 10
    
    # Prepare HDF5 output path
    hdf5_output = output_path if output_path.endswith('.h5') else output_path.replace('.png', '.h5')
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(hdf5_output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"  Created output directory: {output_dir}")
    
    # Remove existing file if it exists to avoid lock issues
    if os.path.exists(hdf5_output):
        try:
            os.remove(hdf5_output)
            print(f"  Removed existing file: {hdf5_output}")
        except Exception as e:
            print(f"  Warning: Could not remove existing file: {e}")
    
    # Create HDF5 file and datasets directly (no large arrays in memory)
    print(f"  Initializing HDF5 file: {os.path.basename(hdf5_output)}")
    try:
        with h5py.File(hdf5_output, 'w', libver='latest', swmr=False) as hf:
            # Create datasets with optional compression
            print(f"  Creating HDF5 datasets with shape ({img_height}, {img_width})...")
            if compression_level is not None:
                print(f"  Using GZIP compression level {compression_level}")
                img_frame_dset = hf.create_dataset(
                    'img_frame', 
                    shape=(img_height, img_width, 3), 
                    dtype=np.uint8, 
                    compression='gzip',
                    compression_opts=compression_level,
                    chunks=True,
                    fillvalue=255
                )
            else:
                print(f"  No compression (fastest write speed)")
                img_frame_dset = hf.create_dataset(
                    'img_frame', 
                    shape=(img_height, img_width, 3), 
                    dtype=np.uint8, 
                    chunks=True,
                    fillvalue=255
                )
            print(f"  Created 'img_frame' dataset")
            # Store metadata
            hf.attrs['format'] = 'heatmap_image_v1'
            hf.attrs['img_height'] = img_height
            hf.attrs['img_width'] = img_width
            print(f"  HDF5 file initialized and ready for data writing")
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
            print(f"  Normalized intensity range: min={log_min}, max={log_max}")
            # Vectorized normalization: map spectrum points
            print(f"  Mapping {len(spec_total)} spectrum points to HDF5 datasets...")
            
            # Vectorize dictionary lookups
            try:
                index_mz_arr = np.array([mz_dict.get(round(mz, 2), None) for mz in spec_total[:, 1]])
                index_rt_arr = np.array([rt_dict.get(round(rt, 8), None) for rt in spec_total[:, 0]])
            except:
                # Fallback if dictionaries have issues
                index_mz_arr = np.full(len(spec_total), None, dtype=object)
                index_rt_arr = np.full(len(spec_total), None, dtype=object)
                for i in range(len(spec_total)):
                    index_mz_arr[i] = mz_dict.get(round(spec_total[i][1], 2), None)
                    index_rt_arr[i] = rt_dict.get(round(spec_total[i][0], 8), None)
            print(f"  Completed mapping spectrum points to indices, filtering valid entries...")
            # Filter valid entries
            valid_mask = (index_mz_arr != None) & (index_rt_arr != None)
            valid_indices = np.where(valid_mask)[0]
            print(f"  Found {len(valid_indices)} valid spectrum points to process")
            if len(valid_indices) == 0:
                print(f"  Warning: No valid spectrum points found")
                return hdf5_output
            
            # Pre-extract indices for chunked processing
            batch_index_mz = index_mz_arr[valid_indices].astype(int)
            batch_index_rt = index_rt_arr[valid_indices].astype(int)
            
            print(f"  Ready for chunked processing with vectorized operations")
            
            # All-at-once processing with row-wise HDF5 writes
            # One row at a time: extract, fill, write, discard (minimal memory)
            print(f"  Processing {len(valid_indices)} points with row-by-row strategy...")
            
            # Extract all valid data at once for efficient filtering
            all_rt = index_rt_arr[valid_indices].astype(int)
            all_mz = index_mz_arr[valid_indices].astype(int)
            all_intensity = spec_total[valid_indices, 2]
            
            # Compute RGB for all points (vectorized - very fast)
            all_log_intensity = np.log1p(all_intensity) if log_scale else all_intensity
            
            if invert_colors:
                all_norm = 1 - ((all_log_intensity - log_min) / log_diff)
            else:
                all_norm = (all_log_intensity - log_min) / log_diff
            
            all_norm = np.clip(all_norm, 0, 1)
            all_rgb_indices = (all_norm * 16777215).astype(np.uint32)
            all_rgb = intensity_to_rgb_vectorized(all_rgb_indices, rgb_array, scale_colors)
            
            # Get unique RT indices
            unique_rt = np.unique(all_rt)
            print(f"  Found {len(unique_rt)} RT indices to write")
            print(f"  Processing {len(unique_rt)} rows in batches of {row_batch_size}...")
            
            # Write rows in batches: prepare multiple rows, write together, discard
            total_written = 0
            num_rows = len(unique_rt)
            
            for batch_start_idx in range(0, num_rows, row_batch_size):
                batch_end_idx = min(batch_start_idx + row_batch_size, num_rows)
                batch_rt_indices = unique_rt[batch_start_idx:batch_end_idx]
                
                # Prepare all rows in this batch
                batch_rows_rgb = {}
                batch_rows_raw = {}
                
                for rt_idx in batch_rt_indices:
                    # Find indices of points for this RT
                    rt_mask = (all_rt == rt_idx)
                    rt_point_indices = np.where(rt_mask)[0]
                    
                    # Create row for this RT
                    row_rgb = np.full((img_width, 3), 255, dtype=np.uint8)
                    
                    # Fill m/z positions for this RT
                    for point_idx in rt_point_indices:
                        mz_idx = all_mz[point_idx]
                        row_rgb[mz_idx] = all_rgb[point_idx]
                    
                    batch_rows_rgb[rt_idx] = row_rgb
                
                # Write all rows in batch at once
                for rt_idx in batch_rt_indices:
                    img_frame_dset[rt_idx, :] = batch_rows_rgb[rt_idx]
                    total_written += np.sum(all_rt == rt_idx)
                
                # Progress counter
                progress_pct = 100 * batch_end_idx / num_rows
                print(f"  Writing rows: {batch_end_idx}/{num_rows} ({progress_pct:.1f}%) - {total_written} points written", end='\r')
                
                # Explicitly delete batch to free memory immediately
                del batch_rows_rgb
            
            print(f"\n  Successfully written all {total_written} points to HDF5")
            
            print(f"\n  Successfully saved heatmap to HDF5: {os.path.basename(hdf5_output)}")
    except Exception as e:
        print(f"  Error creating HDF5 file: {e}")
        print(f"  Attempting to clean up...")
        if os.path.exists(hdf5_output):
            try:
                os.remove(hdf5_output)
            except:
                pass
        raise
    
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
    parser.add_argument("--batch_size", type=int, default=100000, help="Deprecated parameter (ignored)")
    parser.add_argument("--row_batch_size", type=int, default=10, help="Number of rows to prepare and write together (higher=faster but more RAM, each row ~3.2MB)")
    parser.add_argument("--compression_level", type=int, default=None, help="GZIP compression level: None=no compression (fastest, ~1min), 1=light (7min, 0.36MB), 4=balanced (13min, 0.2MB), 9=max (29min, 0.2MB)")
    parser.add_argument("--disable_batching", action="store_true", help="Deprecated parameter (ignored)")
    
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
    output_h5 = create_heatmap_image(
        spec_total, RTINSECONDS_arr, mz_total_arr, mz_dict, rt_dict,
        output_path,
        log_scale=args.log_scale,
        scale_colors=args.scale_colors,
        invert_colors=args.invert_colors,
        save_png=save_png,
        row_batch_size=args.row_batch_size,
        compression_level=args.compression_level
    )
    
    print(f"Heatmap successfully created and saved to {output_h5}")


if __name__ == "__main__":
    main()
