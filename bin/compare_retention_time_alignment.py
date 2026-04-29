#!/usr/bin/env python3
"""
Compare retention time alignment across multiple feature files.

For each pair of files, matches features by pep_ident and rounded m/z center,
then plots how retention time differences vary across the m/z range.

Useful for analyzing retention time drifts between LC runs or alignment quality.
"""

import argparse
import json
import glob
import os
from typing import Dict, List, Tuple
from collections import defaultdict
from datetime import datetime
import pickle
import base64
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from scipy.interpolate import UnivariateSpline, Rbf
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
try:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF as SKLearnRBF, ConstantKernel as C
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


def _format_coefficient(coeff: float, decimal_points: int = 35) -> str:
    """Format coefficient with specified decimal precision, stripping trailing zeros.
    
    Note: Maximum precision capped at 40 decimals to keep output readable.
    Coefficients smaller than 1e-40 will still be formatted but with reduced precision.
    """
    # Cap decimal points at 40 maximum for readability
    max_decimals = min(decimal_points, 40)
    format_str = f"{{:+.{max_decimals}f}}"
    return format_str.format(coeff).rstrip('0').rstrip('.')


def calculate_percentile_threshold(values: List[float], percentile: float) -> Tuple[float, float]:
    """
    Calculate outlier threshold using percentile-based approach (IQR method).
    
    Uses the specified percentile to define the acceptable data range.
    For example, percentile=0.95 keeps the central 95% of data, excluding the outer 5%.
    
    Args:
        values: List of absolute RT differences (should be absolute values)
        percentile: Percentile to keep (0.0-1.0), e.g., 0.95 for 95th percentile
                   If 0 or None, returns no filtering (threshold of infinity)
    
    Returns:
        Tuple of (threshold_value, actual_percentile_used)
        threshold_value: The cutoff for abs(diff) to include in fitting
        actual_percentile_used: The percentile that was calculated (for reporting)
    """
    if not values or percentile <= 0:
        # No filtering - return infinity
        return float('inf'), 0.0
    
    if percentile >= 1.0:
        # Use full data
        return float('inf'), 100.0
    
    # Calculate the percentile value using IQR-based approach
    # For percentile=0.95, calculate the 97.5th percentile to exclude top 5% outliers
    # This is equivalent to keeping the central (percentile*100)% of data
    lower_percentile = (1.0 - percentile) / 2.0 * 100  # e.g., 2.5 for percentile=0.95
    upper_percentile = (100 - lower_percentile)          # e.g., 97.5 for percentile=0.95
    
    lower_bound = np.percentile(values, lower_percentile)
    upper_bound = np.percentile(values, upper_percentile)
    
    # Use the maximum of the two bounds (should be symmetric around 0 for RT diffs)
    threshold = max(abs(lower_bound), abs(upper_bound))
    
    return threshold, percentile * 100.0


def fit_polynomial_equation(x_values: List[float], y_values: List[float], degree: int = 8, decimal_points: int = 35) -> Tuple[np.poly1d, str, int]:
    """
    Fit a polynomial to data and return the polynomial object and equation string.
    
    Args:
        x_values: X coordinates
        y_values: Y coordinates
        degree: Polynomial degree (default 8)
        decimal_points: Decimal places for formatting coefficients (default 35)
    
    Returns:
        Tuple of (polynomial object, equation string, actual_degree_used)
    """
    if len(x_values) < degree + 1:
        # Not enough points, reduce degree
        degree = max(1, len(x_values) - 1)
    
    # Fit polynomial
    coeffs = np.polyfit(x_values, y_values, degree)
    poly = np.poly1d(coeffs)
    
    # Generate equation string
    terms = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        coeff_str = _format_coefficient(coeff, decimal_points)
        if power == 0:
            terms.append(coeff_str)
        elif power == 1:
            terms.append(f"{coeff_str}x")
        else:
            terms.append(f"{coeff_str}x<sup>{power}</sup>")
    
    equation = "y = " + " ".join(terms).replace("+ ", "+ ").replace("- ", "- ")
    equation = equation.replace("+-", "-")
    
    return poly, equation, degree


def fit_polynomial_equation_plain(x_values: List[float], y_values: List[float], degree: int = 8, decimal_points: int = 35) -> Tuple[np.poly1d, str, int]:
    """
    Fit a polynomial to data and return plain-text equation (no HTML tags).
    
    Args:
        x_values: X coordinates
        y_values: Y coordinates
        degree: Polynomial degree (default 8)
        decimal_points: Decimal places for formatting coefficients (default 35)
    
    Returns:
        Tuple of (polynomial object, plain-text equation string, actual_degree_used)
    """
    if len(x_values) < degree + 1:
        degree = max(1, len(x_values) - 1)
    
    coeffs = np.polyfit(x_values, y_values, degree)
    poly = np.poly1d(coeffs)
    
    # Generate plain-text equation string
    terms = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        coeff_str = _format_coefficient(coeff, decimal_points)
        if power == 0:
            terms.append(coeff_str)
        elif power == 1:
            terms.append(f"{coeff_str}*x")
        else:
            terms.append(f"{coeff_str}*x^{power}")
    
    equation = "y = " + " ".join(terms).replace("+ ", "+ ").replace("- ", "- ")
    equation = equation.replace("+-", "-")
    
    return poly, equation, degree


def fit_spline_equation(x_values: List[float], y_values: List[float], smoothing_factor: float = 100, degree: int = 3, decimal_points: int = 35) -> Tuple[object, str, str]:
    """
    Fit a cubic spline to data using UnivariateSpline.
    
    Args:
        x_values: X coordinates (will be sorted internally)
        y_values: Y coordinates (corresponding to x_values)
        smoothing_factor: Smoothing factor for spline (higher = smoother)
        degree: Polynomial degree of spline (default 3 for cubic)
        decimal_points: Decimal places for knot display
    
    Returns:
        Tuple of (spline_object, html_equation_string, plain_text_equation_string)
    """
    # UnivariateSpline requires sorted and unique x values
    # Sort both x and y by x values
    x_sorted = np.array(x_values)
    y_sorted = np.array(y_values)
    sort_idx = np.argsort(x_sorted)
    x_sorted = x_sorted[sort_idx]
    y_sorted = y_sorted[sort_idx]
    
    # Remove duplicate x values (keep first occurrence)
    unique_idx = np.concatenate(([True], x_sorted[1:] != x_sorted[:-1]))
    x_unique = x_sorted[unique_idx]
    y_unique = y_sorted[unique_idx]
    
    if len(x_unique) < 4:
        raise ValueError(f"Not enough unique x values for spline (need ≥4, got {len(x_unique)})")
    
    # Fit spline with adjusted degree based on available points
    adjusted_degree = min(degree, len(x_unique) - 1)
    spline = UnivariateSpline(x_unique, y_unique, s=smoothing_factor, k=adjusted_degree)
    
    # Generate knot information for display
    knots = spline.get_knots()
    
    html_eq = f"Cubic Spline (k={adjusted_degree}, s={smoothing_factor:.2f}), {len(knots)} knots"
    plain_eq = f"Cubic Spline (k={adjusted_degree}, s={smoothing_factor:.2f}), {len(knots)} knots"
    
    return spline, html_eq, plain_eq


def fit_loess_equation(x_values: List[float], y_values: List[float], smoothing_fraction: float = 0.3, decimal_points: int = 35) -> Tuple[object, str, str]:
    """
    Fit LOESS (Locally Weighted Regression) to data using statsmodels.
    
    Args:
        x_values: X coordinates  
        y_values: Y coordinates
        smoothing_fraction: Fraction of data to use for local regression (0-1)
        decimal_points: Decimal places for coefficient display
    
    Returns:
        Tuple of (loess_object, html_equation_string, plain_text_equation_string)
    """
    if not HAS_STATSMODELS:
        raise ImportError("statsmodels required for LOESS fitting. Install: pip install statsmodels")
    
    # Sort data for LOESS
    sorted_indices = np.argsort(x_values)
    x_sorted = np.array(x_values)[sorted_indices]
    y_sorted = np.array(y_values)[sorted_indices]
    
    # Fit LOESS model
    loess_model = lowess(y_sorted, x_sorted, frac=smoothing_fraction, it=3)
    
    # loess_model is an ndarray: [[x1, y1_fitted], [x2, y2_fitted], ...]
    # Use module-level LOESSInterpolator class (supports pickling)
    loess_interpolator = LOESSInterpolator(loess_model[:, 0], loess_model[:, 1])
    
    html_eq = f"LOESS (frac={smoothing_fraction:.3f}, it=3), {len(loess_model)} fitted points"
    plain_eq = f"LOESS (frac={smoothing_fraction:.3f}, it=3), {len(loess_model)} fitted points"
    
    return loess_interpolator, html_eq, plain_eq


def fit_rbf_equation(x_values: List[float], y_values: List[float], rbf_function: str = 'multiquadric', epsilon: float = 1.0, decimal_points: int = 35) -> Tuple[object, str, str]:
    """
    Fit Radial Basis Function to data using scipy.interpolate.Rbf.
    
    Args:
        x_values: X coordinates (will be sorted internally)
        y_values: Y coordinates (corresponding to x_values)
        rbf_function: RBF function type (multiquadric, inverse_multiquadric, gaussian, etc.)
        epsilon: RBF shape parameter
        decimal_points: Decimal places for parameter display
    
    Returns:
        Tuple of (rbf_object, html_equation_string, plain_text_equation_string)
    """
    # Sort by x values
    x_array = np.array(x_values)
    y_array = np.array(y_values)
    sort_idx = np.argsort(x_array)
    x_sorted = x_array[sort_idx]
    y_sorted = y_array[sort_idx]
    
    # RBF can handle duplicate x values, but remove them for better fit
    unique_idx = np.concatenate(([True], x_sorted[1:] != x_sorted[:-1]))
    x_unique = x_sorted[unique_idx]
    y_unique = y_sorted[unique_idx]
    
    rbf_model = Rbf(x_unique, y_unique, function=rbf_function, epsilon=epsilon)
    
    html_eq = f"RBF ({rbf_function}, ε={epsilon:.4f}), {len(x_unique)} basis centers"
    plain_eq = f"RBF ({rbf_function}, epsilon={epsilon:.4f}), {len(x_unique)} basis centers"
    
    return rbf_model, html_eq, plain_eq


# Define LOESSInterpolator class at module level so it can be pickled
class LOESSInterpolator:
    """LOESS interpolator that can be pickled."""
    def __init__(self, x_data, y_data):
        self.x_data = np.array(x_data)
        self.y_data = np.array(y_data)
    
    def __call__(self, x_val):
        """Evaluate LOESS interpolator at x_val (scalar or array)."""
        # Handle both scalar and array inputs
        if isinstance(x_val, (list, np.ndarray)):
            return np.array([self._eval_scalar(x) for x in np.atleast_1d(x_val)])
        else:
            return self._eval_scalar(x_val)
    
    def _eval_scalar(self, x_val):
        """Evaluate at a single scalar value."""
        # Linear interpolation between LOESS fitted points
        if x_val <= self.x_data[0]:
            return self.y_data[0]
        if x_val >= self.x_data[-1]:
            return self.y_data[-1]
        idx = np.searchsorted(self.x_data, x_val)
        # Linear interpolation
        x0, x1 = self.x_data[idx-1], self.x_data[idx]
        y0, y1 = self.y_data[idx-1], self.y_data[idx]
        return y0 + (y1 - y0) * (x_val - x0) / (x1 - x0)


# Define PiecewisePolynomial class at module level so it can be pickled
class PiecewisePolynomial:
    """Piecewise polynomial model that can be pickled."""
    def __init__(self, segments):
        self.segments = segments
    
    def __call__(self, x_val):
        """Evaluate piecewise polynomial at x_val (scalar or array)."""
        # Handle both scalar and array inputs
        if isinstance(x_val, (list, np.ndarray)):
            return np.array([self._eval_scalar(x) for x in np.atleast_1d(x_val)])
        else:
            return self._eval_scalar(x_val)
    
    def _eval_scalar(self, x_val):
        """Evaluate at a single scalar value."""
        for seg in self.segments:
            x_min, x_max = seg['x_range']
            if x_min <= x_val <= x_max:
                return seg['poly'](x_val)
        # Extrapolate using nearest segment
        if x_val < self.segments[0]['x_range'][0]:
            return self.segments[0]['poly'](x_val)
        return self.segments[-1]['poly'](x_val)


def fit_piecewise_polynomial(x_values: List[float], y_values: List[float], num_segments: int = 5, poly_degree: int = 2, decimal_points: int = 35) -> Tuple[object, str, str]:
    """
    Fit piecewise polynomials (split x range into segments, fit polynomial to each).
    
    Args:
        x_values: X coordinates
        y_values: Y coordinates
        num_segments: Number of segments to split into
        poly_degree: Polynomial degree for each segment
        decimal_points: Decimal places for coefficient display
    
    Returns:
        Tuple of (piecewise_object, html_equation_string, plain_text_equation_string)
    """
    x_array = np.array(x_values)
    y_array = np.array(y_values)
    
    # Sort by x
    sorted_indices = np.argsort(x_array)
    x_sorted = x_array[sorted_indices]
    y_sorted = y_array[sorted_indices]
    
    # Split into segments
    segments = []
    x_min, x_max = x_sorted[0], x_sorted[-1]
    
    for i in range(num_segments):
        seg_start = x_min + (i / num_segments) * (x_max - x_min)
        seg_end = x_min + ((i + 1) / num_segments) * (x_max - x_min)
        
        # Get points in this segment (with small overlap)
        overlap_ext = 0.05 * (x_max - x_min)  # 5% overlap
        mask = (x_sorted >= seg_start - overlap_ext) & (x_sorted <= seg_end + overlap_ext)
        
        if np.sum(mask) >= poly_degree + 1:
            x_seg = x_sorted[mask]
            y_seg = y_sorted[mask]
            
            # Fit polynomial to this segment
            coeffs = np.polyfit(x_seg, y_seg, poly_degree)
            poly = np.poly1d(coeffs)
            
            segments.append({
                'x_range': (seg_start, seg_end),
                'coefficients': coeffs.tolist(),
                'poly': poly
            })
    
    # Use module-level PiecewisePolynomial class (supports pickling)
    pw_model = PiecewisePolynomial(segments)
    
    html_eq = f"Piecewise Polynomial ({num_segments} segments, degree={poly_degree})"
    plain_eq = f"Piecewise Polynomial ({num_segments} segments, degree={poly_degree})"
    
    return pw_model, html_eq, plain_eq


def fit_retention_time_model(x_values: List[float], y_values: List[float], 
                            fitting_method: str = 'polynomial', 
                            polynomial_degree: int = 8,
                            spline_degree: int = 3,
                            spline_smoothing: float = 100,
                            loess_fraction: float = 0.3,
                            loess_iterations: int = 3,
                            decimal_points: int = 35) -> Tuple[object, str, str, Dict]:
    """
    Unified function to fit RT correction model using specified method.
    
    All methods fit RT vs RT-Difference (RT on x-axis, RT difference on y-axis).
    This represents how RT correction varies across the retention time spectrum.
    
    Args:
        x_values: Retention time values from reference file (x-axis input for correction model)
        y_values: RT differences between query and reference files (y-axis, correction amount)
        fitting_method: Method to use ('polynomial', 'spline', 'loess', 'rbf', 'piecewise')
        polynomial_degree: Degree for polynomial fitting (only used if method='polynomial')
        spline_degree: Degree for spline fitting (only used if method='spline', default 3)
        spline_smoothing: Smoothing factor for splines (only used if method='spline', default 100)
        loess_fraction: LOESS smoothing fraction - proportion of data for local regression (default 0.3, range 0.01-1.0)
        loess_iterations: LOESS iterations for robustness (default 3, used if method='loess')
        decimal_points: Decimal places for display (only used if method='polynomial')
    
    Returns:
        Tuple of (model_object, html_equation, plain_equation, serialization_dict)
    """
    x_array = np.array(x_values)
    y_array = np.array(y_values)
    
    fit_params = {}
    serialization_data = {}
    
    if fitting_method == 'polynomial':
        model, html_eq, degree = fit_polynomial_equation(x_values, y_values, degree=polynomial_degree, decimal_points=decimal_points)
        plain_eq = html_eq  # Same for polynomial
        fit_params['degree'] = degree
        fit_params['decimal_points'] = decimal_points
        # Round coefficients to 40 decimals max for storage and calculations
        def round_coefficient(c: float) -> float:
            return round(float(c), 40)
        serialization_data['coefficients'] = [round_coefficient(c) for c in model.coefficients]
        serialization_data['degree'] = int(degree)
    
    elif fitting_method == 'spline':
        model, html_eq, plain_eq = fit_spline_equation(x_values, y_values, smoothing_factor=spline_smoothing, degree=spline_degree, decimal_points=decimal_points)
        fit_params['smoothing_factor'] = spline_smoothing
        fit_params['degree'] = spline_degree
        # Spline serialized as pickle (like RBF, LOESS, piecewise)
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')
    
    elif fitting_method == 'loess':
        if not HAS_STATSMODELS:
            raise ImportError("statsmodels required for LOESS. Install: pip install statsmodels")
        model, html_eq, plain_eq = fit_loess_equation(x_values, y_values, smoothing_fraction=loess_fraction, decimal_points=decimal_points)
        fit_params['smoothing_fraction'] = loess_fraction
        fit_params['iterations'] = loess_iterations
        # For LOESS, we'll serialize as pickle since it has a complex structure
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')
    
    elif fitting_method == 'rbf':
        model, html_eq, plain_eq = fit_rbf_equation(x_values, y_values, rbf_function='multiquadric', epsilon=1.0, decimal_points=decimal_points)
        fit_params['rbf_function'] = 'multiquadric'
        fit_params['epsilon'] = 1.0
        # RBF also serialized as pickle
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')
    
    elif fitting_method == 'piecewise':
        model, html_eq, plain_eq = fit_piecewise_polynomial(x_values, y_values, num_segments=5, poly_degree=2, decimal_points=decimal_points)
        fit_params['num_segments'] = 5
        fit_params['poly_degree'] = 2
        # Piecewise serialized as pickle
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')
    
    else:
        raise ValueError(f"Unknown fitting method: {fitting_method}")
    
    return model, html_eq, plain_eq, serialization_data


def load_feature_json(filepath: str) -> List[Dict]:
    """Load feature data from JSON file."""
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"  ✗ Error loading {os.path.basename(filepath)}: {e}")
        return []


def round_with_precision(value: float, decimals: int) -> float:
    """Round a value to specified decimal places."""
    return round(value, decimals)


def find_decimal_precision(value: float, max_decimals: int = 4) -> Tuple[int, float]:
    """
    Find the minimum decimal precision needed for a value.
    Tries from 0 decimals up to max_decimals.
    """
    for decimals in range(max_decimals + 1):
        rounded = round_with_precision(value, decimals)
        # Check if this precision distinguishes the value
        if abs(rounded - value) < 10 ** (-(decimals + 1)):
            return decimals, rounded
    return max_decimals, round_with_precision(value, max_decimals)


def get_coordinate(feature: Dict, coord_type: str, source: str) -> float:
    """
    Extract a coordinate value from a feature based on type and source.
    
    Args:
        feature: Feature dictionary
        coord_type: 'mz' or 'rt'
        source: 'geo_center', 'center', or 'start'
    
    Returns:
        Coordinate value
    """
    if source == 'geo_center':
        key = 'x_center_geo' if coord_type == 'mz' else 'y_center_geo'
    elif source == 'center':
        key = 'x_center' if coord_type == 'mz' else 'y_center'
    elif source == 'start':
        key = 'mz_start' if coord_type == 'mz' else 'rt_start'
    else:
        raise ValueError(f"Unknown coordinate source: {source}")
    
    return feature.get(key, 0.0)


def resolve_mz_collisions(features_by_pep: Dict[str, Dict], mz_source: str = 'geo_center') -> Dict[str, Tuple[int, float]]:
    """
    Resolve m/z collisions by reducing decimal precision when needed.
    
    Args:
        mz_source: Which m/z value to use ('geo_center', 'center', or 'start')
    
    Returns: Dict mapping pep_ident -> (decimal_places_used, rounded_mz)
    """
    pep_to_mz = {}
    assigned_peps = set()
    
    # Try decimal places from 0 to 4
    for target_decimals in range(5):
        occupied_mz = {}  # maps rounded_mz -> list of pep_idents
        
        for pep_ident, pep_data in features_by_pep.items():
            if pep_ident in assigned_peps:
                continue  # Already assigned
            
            mz = get_coordinate(pep_data, 'mz', mz_source)
            rounded_mz = round_with_precision(mz, target_decimals)
            
            if rounded_mz not in occupied_mz:
                occupied_mz[rounded_mz] = []
            occupied_mz[rounded_mz].append(pep_ident)
        
        # Assign non-conflicting peps and mark conflicting ones
        for rounded_mz, peps in occupied_mz.items():
            if len(peps) == 1:
                pep_ident = peps[0]
                pep_to_mz[pep_ident] = (target_decimals, rounded_mz)
                assigned_peps.add(pep_ident)
        
        # If all assigned, we're done
        if len(assigned_peps) == len(features_by_pep):
            break
    
    # Assign any remaining peps with max decimals
    for pep_ident in features_by_pep.keys():
        if pep_ident not in assigned_peps:
            mz = get_coordinate(features_by_pep[pep_ident], 'mz', mz_source)
            rounded_mz = round_with_precision(mz, 4)
            pep_to_mz[pep_ident] = (4, rounded_mz)
    
    return pep_to_mz


def match_features_by_closest_mz_and_rt(features1: List[Dict], features2: List[Dict], 
                                        mz_weight: float = 1.0, rt_weight: float = 1.0,
                                        mz_source: str = 'geo_center', rt_source: str = 'geo_center') -> List[Tuple[Dict, Dict]]:
    """
    Match features from file1 to file2 by finding closest neighbors in (m/z, RT) space.
    Each feature in file2 is used only once.
    
    Uses weighted Euclidean distance to avoid false matches:
    - Features with same m/z but very different RT are penalized
    - Features with similar RT and m/z are preferred
    
    Args:
        features1: List of feature dicts from file1, all with same pep_ident
        features2: List of feature dicts from file2, all with same pep_ident
        mz_weight: Weight for m/z difference in distance calculation
        rt_weight: Weight for RT difference in distance calculation
        mz_source: Which m/z value to use ('geo_center', 'center', or 'start')
        rt_source: Which RT value to use ('geo_center', 'center', or 'start')
    
    Returns:
        List of (feature1, feature2) tuples representing matched pairs
    """
    if not features1 or not features2:
        return []
    
    # Sort features by m/z for initial organization
    f1_sorted = sorted(features1, key=lambda f: get_coordinate(f, 'mz', mz_source))
    f2_sorted = sorted(features2, key=lambda f: get_coordinate(f, 'mz', mz_source))
    
    matches = []
    used_f2_indices = set()
    
    # For each feature in file1, find closest neighbor in 2D space (m/z, RT)
    for f1 in f1_sorted:
        min_distance = float('inf')
        best_f2_idx = None
        
        for f2_idx, f2 in enumerate(f2_sorted):
            if f2_idx in used_f2_indices:
                continue
            
            # Calculate 2D distance in (m/z, RT) space
            f1_mz = get_coordinate(f1, 'mz', mz_source)
            f1_rt = get_coordinate(f1, 'rt', rt_source)
            f2_mz = get_coordinate(f2, 'mz', mz_source)
            f2_rt = get_coordinate(f2, 'rt', rt_source)
            
            mz_diff = abs(f1_mz - f2_mz)
            rt_diff = abs(f1_rt - f2_rt)
            
            # Weighted Euclidean distance
            distance = ((mz_weight * mz_diff) ** 2 + (rt_weight * rt_diff) ** 2) ** 0.5
            
            if distance < min_distance:
                min_distance = distance
                best_f2_idx = f2_idx
        
        # Match if we found an unused feature in file2
        if best_f2_idx is not None:
            matches.append((f1, f2_sorted[best_f2_idx]))
            used_f2_indices.add(best_f2_idx)
    
    return matches



def compare_files(file1_data: List[Dict], file2_data: List[Dict], 
                  file1_name: str, file2_name: str, one_to_one_only: bool = False,
                  mz_source: str = 'geo_center', rt_source: str = 'geo_center') -> Tuple[Dict[float, float], Dict]:
    """
    Compare two feature files and return RT differences keyed by m/z, plus statistics.
    Matches multiple features with the same pep_ident by closest m/z.
    
    Args:
        one_to_one_only: If True, only include pep_idents that appear exactly once in both files
        mz_source: Which m/z value to use ('geo_center', 'center', or 'start')
        rt_source: Which RT value to use ('geo_center', 'center', or 'start')
    
    Returns: 
        Tuple of:
        - [{x: rounded_mz, y: rt_diff, pep_ident: pep_id}, ...]
        - {num_matches, mean_rt_diff, std_rt_diff}
    """
    # Group features by pep_ident for each file
    file1_by_pep = defaultdict(list)
    file2_by_pep = defaultdict(list)
    
    for feature in file1_data:
        pep_idents = feature.get('l_raw_pep_ident', [])
        if pep_idents:
            # l_raw_pep_ident can be a list of strings
            if isinstance(pep_idents, list):
                for pep in pep_idents:
                    file1_by_pep[pep].append(feature)
            else:
                file1_by_pep[pep_idents].append(feature)
    
    for feature in file2_data:
        pep_idents = feature.get('l_raw_pep_ident', [])
        if pep_idents:
            if isinstance(pep_idents, list):
                for pep in pep_idents:
                    file2_by_pep[pep].append(feature)
            else:
                file2_by_pep[pep_idents].append(feature)
    
    # Find common pep_idents
    common_peps = set(file1_by_pep.keys()) & set(file2_by_pep.keys())
    
    # Filter to 1:1 matches if requested
    if one_to_one_only:
        common_peps = {pep for pep in common_peps 
                      if len(file1_by_pep[pep]) == 1 and len(file2_by_pep[pep]) == 1}
    
    if not common_peps:
        print(f"    ⚠ No common pep_idents between {file1_name} and {file2_name}")
        return [], {}
    
    # Match features by closest m/z for each pep_ident
    # First, create a dict of single features for m/z collision resolution
    file1_first_per_pep = {pep: file1_by_pep[pep][0] for pep in common_peps}
    file2_first_per_pep = {pep: file2_by_pep[pep][0] for pep in common_peps}
    
    # Resolve m/z collisions
    pep_to_mz = resolve_mz_collisions(file1_first_per_pep, mz_source=mz_source)
    
    # Collect all matched pairs across all pep_idents
    all_matches = []
    
    for pep_ident in common_peps:
        decimals = pep_to_mz[pep_ident][0]
        
        # Get features for this pep_ident from both files
        f1_features = file1_by_pep[pep_ident]
        f2_features = file2_by_pep[pep_ident]
        
        # Match by closest point in (m/z, RT) space
        # This prevents false matches where m/z is close but RT is very different
        matches = match_features_by_closest_mz_and_rt(f1_features, f2_features,
                                                       mz_source=mz_source, rt_source=rt_source)
        
        for f1, f2 in matches:
            # Round m/z with determined precision
            f1_mz = get_coordinate(f1, 'mz', mz_source)
            f2_mz = get_coordinate(f2, 'mz', mz_source)
            rounded_mz = round_with_precision(f1_mz, decimals)
            
            # In 1:1 mode, accept all matches since each pep_ident appears exactly once
            # In multi-match mode, verify f2's rounded m/z also matches to avoid false pairings
            if one_to_one_only:
                all_matches.append((f1, f2, rounded_mz, pep_ident))
            else:
                f2_rounded_mz = round_with_precision(f2_mz, decimals)
                if f2_rounded_mz == rounded_mz:
                    all_matches.append((f1, f2, rounded_mz, pep_ident))
    
    # Build result list and calculate statistics
    result = []
    rt_differences = []
    
    for f1, f2, rounded_mz, pep_ident in all_matches:
        f1_rt = get_coordinate(f1, 'rt', rt_source)
        f2_rt = get_coordinate(f2, 'rt', rt_source)
        rt_diff = f2_rt - f1_rt
        result.append({
            'x': rounded_mz,
            'y': rt_diff,
            'rt_file1': f1_rt,  # RT from file 1
            'pep_ident': pep_ident
        })
        rt_differences.append(rt_diff)
    
    # Calculate statistics
    stats = {}
    if rt_differences:
        stats['num_matches'] = len(rt_differences)
        stats['mean_rt_diff'] = float(np.mean(rt_differences))
        stats['std_rt_diff'] = float(np.std(rt_differences))
    else:
        stats['num_matches'] = 0
        stats['mean_rt_diff'] = 0.0
        stats['std_rt_diff'] = 0.0
    
    return result, stats



def generate_curve_approximation(model: object, x_values: List[float], method: str = 'polynomial') -> Tuple[List[float], List[float]]:
    """
    Generate curve approximation points from a fitted model.
    Uses n_points = int(max(x_values)) for approximation count.
    
    Args:
        model: Fitted model object (poly1d, UnivariateSpline, Rbf, LOESSInterpolator, etc.)
        x_values: Input x values (used to determine range and point count)
        method: Fitting method name for debugging
    
    Returns:
        Tuple of (x_approximation, y_approximation) lists
    """
    try:
        x_array = np.array(x_values)
        min_x = float(np.min(x_array))
        max_x = float(np.max(x_array))
        n_points = int(max_x)
        
        # Ensure at least some points
        n_points = max(50, n_points)
        
        # Generate approximation points
        x_fit = np.linspace(min_x, max_x, n_points)
        y_fit = model(x_fit)
        
        return x_fit.tolist(), y_fit.tolist()
    except Exception as e:
        raise ValueError(f"Failed to generate curve approximation for {method}: {str(e)}")


def create_comparison_plot(file1_name: str, comparisons: Dict[str, list], outlier_threshold: float = 0.95, polynomial_degree: int = 8, decimal_points: int = 35, fitting_method: str = 'polynomial', corrections_dict: Dict = None) -> Tuple[go.Figure, Dict[str, str]]:
    """
    Create a Plotly figure showing RT differences for one file vs all others.
    Includes fitted curves for each comparison.
    
    Args:
        file1_name: Base file name
        comparisons: {file2_name: [{'x': mz, 'y': rt_diff, 'pep_ident': pep_id}, ...]}
        outlier_threshold: Percentile threshold for outlier exclusion (0.0-1.0).
                          0.95 = keep central 95% of data, exclude outer 5% (default 0.95).
                          0 or None = no filtering (keep all data).
        polynomial_degree: Polynomial degree for fitting (default 8)
        decimal_points: Decimal places for coefficient formatting (default 35)
    
    Returns:
        Tuple of (Plotly Figure object, {file2_name: equation_string})
    """
    fig = go.Figure()
    equations = {}
    
    # Add a trace for each file comparison
    for file2_name, data_points in sorted(comparisons.items()):
        if not data_points:
            continue
        
        # Sort by m/z for line plot
        sorted_data = sorted(data_points, key=lambda d: d['x'])
        mz_values = [d['x'] for d in sorted_data]
        rt_diffs = [d['y'] for d in sorted_data]
        
        # Create custom hover text with pep_ident
        hover_text = [
            f"<b>{file2_name}</b><br>M/Z: {d['x']:.4f}<br>RT Diff: {d['y']:.2f}s<br>Pep ID: {d['pep_ident']}"
            for d in sorted_data
        ]
        
        # Add data trace
        fig.add_trace(go.Scatter(
            x=mz_values,
            y=rt_diffs,
            mode='lines+markers',
            name=f'vs {file2_name}',
            visible=True,
            hovertext=hover_text,
            hoverinfo='text'
        ))
        
        # Fit polynomial and add fitted curve
        if len(mz_values) >= 3:  # At least 3 points needed for fitting
            # Calculate percentile-based outlier threshold using IQR method
            abs_rt_diffs = [abs(d['y']) for d in sorted_data]
            threshold_value, percentile_used = calculate_percentile_threshold(abs_rt_diffs, outlier_threshold)
            
            # Filter out extreme outliers for fitting (but keep all points displayed)
            filtered_mz = [d['x'] for d in sorted_data if abs(d['y']) <= threshold_value]
            filtered_rt = [d['y'] for d in sorted_data if abs(d['y']) <= threshold_value]
            
            if len(filtered_mz) >= 3:
                poly, equation, deg = fit_polynomial_equation(filtered_mz, filtered_rt, degree=polynomial_degree, decimal_points=decimal_points)
                percentile_note = f" ({percentile_used:.1f}% data, outliers excluded)" if percentile_used > 0 else ""
                equations[file2_name] = f"{equation} <span style='font-size: 0.85em; color: #666;'>(degree {deg}){percentile_note}</span>"
            else:
                equations[file2_name] = "<span style='color: #999;'>(Not enough points after filtering)</span>"
                poly = None
            
            # Generate fitted curve points (only if we have a valid polynomial)
            if poly is not None:
                x_fit = np.linspace(min(mz_values), max(mz_values), 200)
                y_fit = poly(x_fit)
                
                # Convert numpy arrays to lists for Plotly serialization
                x_fit_list = x_fit.tolist()
                y_fit_list = y_fit.tolist()
                
                # Add fitted curve trace (more visible: thicker, solid line, bright color)
                colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8']
                color_idx = list(sorted(comparisons.keys())).index(file2_name) % len(colors)
                
                fig.add_trace(go.Scatter(
                    x=x_fit_list,
                    y=y_fit_list,
                    mode='lines',
                    name=f'Fit: {file2_name}',
                    line=dict(width=4, color=colors[color_idx]),
                    visible=False,
                    hoverinfo='skip',
                    customdata=[f"fit_{file2_name}"] * len(x_fit_list)
                ))
    
    # Update layout
    fig.update_layout(
        title={
            'text': f'Retention Time Alignment – {file1_name}',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18}
        },
        xaxis_title='M/Z (m/z)',
        yaxis_title='RT Difference (seconds)',
        hovermode='x unified',
        height=600,
        template='plotly_white',
        font=dict(size=12),
        legend=dict(
            title='File Comparisons',
            yanchor='top',
            y=0.99,
            xanchor='right',
            x=0.99
        )
    )
    
    return fig, equations


def create_rt_vs_rtdiff_plot(file1_name: str, comparisons: Dict[str, list], outlier_threshold: float = 0.95, polynomial_degree: int = 8, decimal_points: int = 35, fitting_method: str = 'polynomial', spline_degree: int = 3, spline_smoothing: float = 100, loess_fraction: float = 0.3, loess_iterations: int = 3, corrections_dict: Dict = None, equations_plain_dict: Dict = None) -> Tuple[go.Figure, Dict[str, str], Dict[str, str]]:
    """
    Create a Plotly figure showing RT difference vs RT from file1.
    This reveals if the RT drift is uniform across the LC run or varies by retention time.
    Includes fitted curves for each comparison using the specified fitting method.
    
    Args:
        file1_name: Base file name
        comparisons: {file2_name: [{'x': mz, 'y': rt_diff, 'rt_file1': rt1, 'pep_ident': pep_id}, ...]}
        outlier_threshold: Percentile threshold for outlier exclusion (0.0-1.0).
                          0.95 = keep central 95% of data, exclude outer 5% (default 0.95).
                          0 or None = no filtering (keep all data).
        polynomial_degree: Polynomial degree for fitting (only used if method='polynomial')
        decimal_points: Decimal places for coefficient formatting (default 35)
        fitting_method: Method used for fitting ('polynomial', 'spline', 'rbf', 'loess', 'piecewise')
        spline_degree: Degree for spline fitting (only used if method='spline', default 3)
        spline_smoothing: Smoothing factor for splines (only used if method='spline', default 100)
        loess_fraction: LOESS smoothing fraction (only used if method='loess', default 0.3)
        loess_iterations: LOESS iterations (only used if method='loess', default 3)
        corrections_dict: Pre-fitted corrections {'file1': {'file2': {...model_data...}}}
        equations_plain_dict: Plain-text equations {'file1': {'file2': equation_text}}
    
    Returns:
        Tuple of (Plotly Figure object, {file2_name: equation_string_html}, {file2_name: equation_string_plain})
    """
    fig = go.Figure()
    equations = {}
    equations_plain = {}
    
    if corrections_dict is None:
        corrections_dict = {}
    if equations_plain_dict is None:
        equations_plain_dict = {}
    
    # Add a trace for each file comparison
    for file2_name, data_points in sorted(comparisons.items()):
        if not data_points:
            continue
        
        # Sort by RT from file 1 for line plot
        sorted_data = sorted(data_points, key=lambda d: d['rt_file1'])
        rt_values = [d['rt_file1'] for d in sorted_data]
        rt_diffs = [d['y'] for d in sorted_data]
        
        # Create custom hover text with pep_ident
        hover_text = [
            f"<b>{file2_name}</b><br>RT (File1): {d['rt_file1']:.2f}s<br>RT Diff: {d['y']:.2f}s<br>Pep ID: {d['pep_ident']}"
            for d in sorted_data
        ]
        
        # Add data trace
        fig.add_trace(go.Scatter(
            x=rt_values,
            y=rt_diffs,
            mode='lines+markers',
            name=f'vs {file2_name}',
            visible=True,
            hovertext=hover_text,
            hoverinfo='text'
        ))
        
        # Prepare equation text for button
        equation_text = None
        model_obj = None
        
        # Try to get pre-fitted model from corrections_dict
        if file1_name in corrections_dict and file2_name in corrections_dict[file1_name]:
            correction_data = corrections_dict[file1_name][file2_name]
            
            # Get plain-text equation if available
            if equations_plain_dict and file1_name in equations_plain_dict and file2_name in equations_plain_dict[file1_name]:
                equation_text = equations_plain_dict[file1_name][file2_name]
            
            # Try to deserialize the model
            try:
                if fitting_method == 'polynomial' and 'coefficients' in correction_data:
                    coefficients = correction_data['coefficients']
                    model_obj = np.poly1d(coefficients)
                    if not equation_text:
                        equation_text = f"Polynomial (degree {len(coefficients)-1})"
                else:
                    # For other methods, try to deserialize pickled model
                    if 'model' in correction_data:
                        model_str = correction_data['model']
                        if model_str.startswith('pikle:'):
                            pickling_data = model_str.replace('pikle:', '')
                            model_obj = pickle.loads(base64.b64decode(pickling_data))
                        if not equation_text:
                            equation_text = f"fitted curve for {file2_name}"
            except Exception as e:
                pass
        
        # If no pre-fitted model, fit one now for backward compatibility
        if model_obj is None:
            # Calculate percentile-based outlier threshold using IQR method
            abs_rt_diffs = [abs(d['y']) for d in sorted_data]
            threshold_value, percentile_used = calculate_percentile_threshold(abs_rt_diffs, outlier_threshold)
            
            # Filter data for fitting
            filtered_rt = [d['rt_file1'] for d in sorted_data if abs(d['y']) <= threshold_value]
            filtered_diff = [d['y'] for d in sorted_data if abs(d['y']) <= threshold_value]
            
            if len(filtered_rt) >= 3:
                try:
                    if fitting_method == 'polynomial':
                        poly, equation, deg = fit_polynomial_equation(filtered_rt, filtered_diff, degree=polynomial_degree, decimal_points=decimal_points)
                        _, equation_plain, _ = fit_polynomial_equation_plain(filtered_rt, filtered_diff, degree=polynomial_degree, decimal_points=decimal_points)
                        model_obj = poly
                        equation_text = f"{equation} <span style='font-size: 0.85em; color: #666;'>(degree {deg})</span>"
                    elif fitting_method == 'spline':
                        model_obj, _, equation_plain = fit_spline_equation(filtered_rt, filtered_diff, smoothing_factor=spline_smoothing, degree=spline_degree, decimal_points=decimal_points)
                        equation_text = f"fitted curve for {file2_name}"
                    elif fitting_method == 'rbf':
                        model_obj, _, equation_plain = fit_rbf_equation(filtered_rt, filtered_diff, rbf_function='multiquadric', epsilon=1.0, decimal_points=decimal_points)
                        equation_text = f"fitted curve for {file2_name}"
                    elif fitting_method == 'loess':
                        model_obj, _, equation_plain = fit_loess_equation(filtered_rt, filtered_diff, smoothing_fraction=loess_fraction, decimal_points=decimal_points)
                        equation_text = f"fitted curve for {file2_name}"
                    elif fitting_method == 'piecewise':
                        model_obj, _, equation_plain = fit_piecewise_polynomial(filtered_rt, filtered_diff, num_segments=5, poly_degree=2, decimal_points=decimal_points)
                        equation_text = f"fitted curve for {file2_name}"
                    else:
                        equation_text = f"(Unknown method: {fitting_method})"
                        equation_plain = "(Unknown method)"
                    
                    equations_plain[file2_name] = equation_plain
                except Exception as e:
                    equation_text = f"(Fitting failed: {str(e)})"
                    equation_plain = f"(Fitting failed: {str(e)})"
                    equations_plain[file2_name] = equation_plain
                    model_obj = None
            else:
                equation_text = "<span style='color: #999;'>(Not enough points after filtering)</span>"
                equation_plain = "(Not enough points after filtering)"
                equations_plain[file2_name] = equation_plain
                model_obj = None
        
        # Add equation text
        if equation_text:
            equations[file2_name] = equation_text if isinstance(equation_text, str) else equation_text.decode() if isinstance(equation_text, bytes) else str(equation_text)
        
        # Generate fitted curve if model is available
        if model_obj is not None:
            try:
                x_fit, y_fit = generate_curve_approximation(model_obj, rt_values, fitting_method)
                
                # Add fitted curve trace (more visible: thicker, solid line, bright color)
                colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8']
                color_idx = list(sorted(comparisons.keys())).index(file2_name) % len(colors)
                
                fig.add_trace(go.Scatter(
                    x=x_fit,
                    y=y_fit,
                    mode='lines',
                    name=f'Fit: {file2_name}',
                    line=dict(width=4, color=colors[color_idx]),
                    visible=False,
                    hoverinfo='skip',
                    customdata=[f"fit_{file2_name}"] * len(x_fit)
                ))
            except Exception as e:
                print(f"    Warning: Could not generate curve for {file2_name}: {e}")
    
    # Update layout
    fig.update_layout(
        title={
            'text': f'RT Drift Pattern – {file1_name}',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18}
        },
        xaxis_title='Retention Time (seconds)',
        yaxis_title='RT Difference (seconds)',
        hovermode='x unified',
        height=600,
        template='plotly_white',
        font=dict(size=12),
        legend=dict(
            title='File Comparisons',
            yanchor='top',
            y=0.99,
            xanchor='right',
            x=0.99
        )
    )
    
    return fig, equations, equations_plain


def create_rt_transformation_plots(file1_name: str, all_files_data: Dict[str, List[Dict]], 
                                    corrections_dict: Dict[Tuple[str, str], Dict]) -> Dict[str, go.Figure]:
    """
    Create transformation plots showing Original RT vs Transformed RT for each file against a reference.
    
    For each non-reference file, creates a plot with Original RT (x-axis) vs Transformed RT (y-axis),
    using the fitted correction model from that file to the reference file.
    
    Args:
        file1_name: Reference file name (usually the first file)
        all_files_data: {filename: [list of features]}
        corrections_dict: {(file1, file2): {forward/reverse: fitted_model_data}}
    
    Returns:
        {filename: Plotly Figure} - one figure per non-reference file
    """
    transformation_figs = {}
    
    # Get all files except reference
    other_files = [f for f in all_files_data.keys() if f != file1_name]
    
    for query_file in other_files:
        fig = go.Figure()
        
        # Get the correction model from query_file to reference (file1_name)
        correction_key = (query_file, file1_name)
        
        if correction_key not in corrections_dict:
            # Try reverse direction
            correction_key = (file1_name, query_file)
            if correction_key not in corrections_dict:
                print(f"  ⚠ No correction model found for {query_file} → {file1_name}")
                continue
            # Use reverse direction (need to negate the correction)
            use_reverse = True
        else:
            use_reverse = False
        
        correction_data = corrections_dict[correction_key]
        
        # Determine the direction to use
        if use_reverse:
            model_data = correction_data.get('reverse', {})
        else:
            model_data = correction_data.get('forward', {})
        
        if not model_data or 'model_pickle' not in model_data:
            print(f"  ⚠ No fitted model in correction data for {query_file} → {file1_name}")
            continue
        
        # Unpickle the model
        try:
            model_b64 = model_data.get('model_pickle', '')
            model_bytes = base64.b64decode(model_b64.encode('utf-8'))
            model = pickle.loads(model_bytes)
        except Exception as e:
            print(f"  ✗ Failed to load model for {query_file}: {e}")
            continue
        
        # Get features from query file
        features = all_files_data.get(query_file, [])
        if not features:
            continue
        
        # Extract original RT values
        original_rts = []
        for feature in features:
            if 'y_center' in feature:
                original_rts.append(feature['y_center'])
        
        if not original_rts:
            continue
        
        original_rts = np.array(original_rts)
        
        # Apply correction to get transformed RTs
        try:
            if use_reverse:
                # Negate for reverse direction
                transformed_rts = original_rts - model(original_rts)
            else:
                transformed_rts = original_rts + model(original_rts)
        except Exception as e:
            print(f"  ✗ Failed to apply correction model for {query_file}: {e}")
            continue
        
        # Add scatter plot
        fig.add_trace(go.Scatter(
            x=original_rts,
            y=transformed_rts,
            mode='markers',
            name=f'{query_file} → {file1_name}',
            marker=dict(
                size=6,
                opacity=0.6,
                line=dict(width=0.5, color='white')
            ),
            hovertemplate='<b>Original RT:</b> %{x:.2f}<br>' +
                         '<b>Transformed RT:</b> %{y:.2f}<br>' +
                         '<b>Δ RT:</b> %{customdata:.2f}<extra></extra>',
            customdata=transformed_rts - original_rts
        ))
        
        # Add perfect alignment line (y=x)
        rt_min, rt_max = original_rts.min(), original_rts.max()
        fig.add_trace(go.Scatter(
            x=[rt_min, rt_max],
            y=[rt_min, rt_max],
            mode='lines',
            name='Perfect Alignment (y=x)',
            line=dict(dash='dash', color='gray', width=2),
            hoverinfo='skip',
            showlegend=True,
            visible=True
        ))
        
        # Update layout
        fig.update_layout(
            title=f'RT Transformation: {query_file} → Reference ({file1_name})',
            xaxis_title='Original RT (query file)',
            yaxis_title='Transformed RT (after correction)',
            hovermode='closest',
            width=900,
            height=600,
            template='plotly_white',
            font=dict(size=12),
            legend=dict(
                yanchor='top',
                y=0.99,
                xanchor='left',
                x=0.01
            )
        )
        
        # Make axes equal for better visualization
        fig.update_xaxes(scaleanchor='y', scaleratio=1)
        fig.update_yaxes(scaleanchor='x', scaleratio=1)
        
        transformation_figs[query_file] = fig
    
    return transformation_figs


def serialize_model_to_json(method: str, model: object, x_vals: np.ndarray, y_vals: np.ndarray, 
                           fit_params: Dict[str, any]) -> Dict[str, any]:
    """
    Serialize a fitted model to JSON-compatible format.
    
    Args:
        method: Fitting method name ('polynomial', 'spline', 'loess', 'rbf', 'piecewise')
        model: The fitted model object
        x_vals: Original X values used for fitting
        y_vals: Original Y values used for fitting
        fit_params: Additional parameters specific to the method
    
    Returns:
        Dictionary with 'method', 'serialization', and 'metadata' keys
    """
    serialization_data = {
        'method': method,
        'metadata': {
            'n_points': len(x_vals),
            'x_min': float(np.min(x_vals)),
            'x_max': float(np.max(x_vals)),
            'y_min': float(np.min(y_vals)),
            'y_max': float(np.max(y_vals)),
            'fit_params': fit_params
        }
    }
    
    if method == 'polynomial':
        # Polynomial coefficients are simple floats
        serialization_data['coefficients'] = [float(c) for c in model.coefficients]
        serialization_data['degree'] = int(model.order)
    
    elif method == 'spline':
        # UnivariateSpline: serialize knots and coefficients
        knots = model.get_knots().tolist()
        coeffs = model.get_coeffs().tolist()
        serialization_data['knots'] = knots
        serialization_data['coefficients'] = coeffs
        serialization_data['degree'] = int(model.k)
    
    elif method in ['loess', 'rbf', 'piecewise']:
        # For complex methods, serialize the model as base64-encoded pickle
        import pickle
        import base64
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')
    
    return serialization_data


def export_fitted_corrections_json(corrections_dict: Dict[str, Dict[str, any]], 
                                  output_path: str, fitting_method: str):
    """
    Export all fitted corrections to JSON file.
    
    Args:
        corrections_dict: Dict[file_pair_name][direction] = serialized_model_dict
        output_path: Path to output JSON file
        fitting_method: Name of fitting method used
    """
    output_data = {
        'fitting_method': fitting_method,
        'corrections': corrections_dict,
        'timestamp': datetime.now().isoformat()
    }
    
    output_dir = os.path.dirname(output_path)
    if output_dir:  # Only create directory if it's not empty
        os.makedirs(output_dir, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"✓ Exported corrections to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Compare retention time alignment across multiple feature files'
    )
    parser.add_argument('--input_dir', default='results/feature_analysis/feature_data_lists',
                       help='Directory containing feature_data.json files')
    parser.add_argument('--input_pattern', default='*_feature_data.json',
                       help='Glob pattern for feature data files')
    parser.add_argument('--output_html', default='results/retention_time_alignment_report.html',
                       help='Output HTML report file')
    parser.add_argument('--output_json', default='results/feature_analysis/rt_alignment_comparison/fitted_corrections.json',
                       help='Output JSON file with fitted corrections (default: results/feature_analysis/rt_alignment_comparison/fitted_corrections.json)')
    parser.add_argument('--fitting_method', default='polynomial', 
                       choices=['polynomial', 'spline', 'loess', 'rbf', 'piecewise'],
                       help='Method for fitting RT correction curves (default: polynomial)')
    parser.add_argument('--outlier_threshold', type=float, default=0.95,
                       help='Percentile threshold for outlier exclusion (0.0-1.0). 0.95 = keep central 95%% of data, exclude outer 5%% (default: 0.95). Set to 0 to disable filtering.')
    parser.add_argument('--one_to_one_only', action='store_true',
                       help='Only include pep_idents that appear exactly once in both files')
    parser.add_argument('--mz_source', default='geo_center', choices=['geo_center', 'center', 'start'],
                       help='Which m/z value to use: geo_center (default), center, or start')
    parser.add_argument('--rt_source', default='geo_center', choices=['geo_center', 'center', 'start'],
                       help='Which RT value to use: geo_center (default), center, or start')
    parser.add_argument('--polynomial_degree', type=int, default=8,
                       help='Polynomial degree for fitting RT correction curves (default: 8)')
    parser.add_argument('--spline_degree', type=int, default=8,
                       help='Spline degree for curve fitting (default: 3, used with spline method)')
    parser.add_argument('--spline_smoothing', type=float, default=1000,
                       help='Smoothing factor for splines (default: 100, used with spline method)')
    parser.add_argument('--loess_fraction', type=float, default=0.3,
                       help='LOESS smoothing fraction - proportion of data for local regression (default: 0.3, range 0.01-1.0)')
    parser.add_argument('--loess_iterations', type=int, default=3,
                       help='LOESS iterations for robustness (default: 3)')
    parser.add_argument('--decimal_points', type=int, default=35,
                       help='Decimal places for polynomial coefficients (default: 35)')
    parser.add_argument('--show-plots', action='store_true',
                       help='Display plots interactively in browser during execution (requires browser access)')
    
    args = parser.parse_args()
    
    # Find all feature data files
    if not os.path.isdir(args.input_dir):
        print(f"✗ Input directory not found: {args.input_dir}")
        return
    
    pattern = os.path.join(args.input_dir, args.input_pattern)
    json_files = sorted(glob.glob(pattern))
    
    if not json_files:
        print(f"✗ No files matching pattern: {pattern}")
        return
    
    print(f"\n{'='*70}")
    print(f"Retention Time Alignment Analysis")
    print(f"{'='*70}")
    print(f"Found {len(json_files)} feature data files\n")
    
    # Load all files - use REAL FILENAMES as keys, NOT basenames
    file_data = {}
    file_mapping = {}  # Maps basename to real filename
    for filepath in json_files:
        basename = os.path.basename(filepath).replace('_feature_data.json', '')
        print(f"Loading {basename}...", end='', flush=True)
        
        data = load_feature_json(filepath)
        if data:
            # Extract the real filename from the first feature's filename field
            if data and isinstance(data, list) and len(data) > 0 and 'filename' in data[0]:
                real_filename = data[0]['filename'].replace('_feature_data.json', '')
            else:
                # Fallback: use basename if no filename field
                real_filename = basename
            
            # Use REAL FILENAME as the dictionary key (not basename)
            file_data[real_filename] = data
            file_mapping[basename] = real_filename
            print(f" ✓ ({len(data)} features, real name: '{real_filename}')", flush=True)
        else:
            print(f" ✗", flush=True)
    
    print(f"\nFile mapping (JSON basename → real filename): {file_mapping}\n", flush=True)
    
    if len(file_data) < 2:
        print(f"\n✗ Need at least 2 files for comparison")
        return
    
    print(f"\n{'='*70}")
    print(f"Comparing files ({len(file_data)} files, {len(file_data)*(len(file_data)-1)//2} pairs)")
    print(f"{'='*70}\n")
    
    # Create comparisons: base_file -> {comparison_file: {mz: rt_diff}}
    all_comparisons = {}
    all_stats = {}  # Store statistics for each comparison
    
    for i, (file1_name, file1_data) in enumerate(sorted(file_data.items())):
        comparisons = {}
        stats_dict = {}
        
        for file2_name, file2_data in sorted(file_data.items()):
            if file1_name == file2_name:
                continue
            
            print(f"Comparing {file1_name} vs {file2_name}...", end='')
            result, stats = compare_files(file1_data, file2_data, file1_name, file2_name, 
                                         one_to_one_only=args.one_to_one_only,
                                         mz_source=args.mz_source, rt_source=args.rt_source)
            
            if result:
                comparisons[file2_name] = result
                stats_dict[file2_name] = stats
                print(f" ✓ ({len(result)} matches)")
            else:
                print(f" ⚠ (no matches)")
        
        all_comparisons[file1_name] = comparisons
        all_stats[file1_name] = stats_dict
    
    # Fit RT corrections using specified method
    print(f"\n{'='*70}")
    print(f"Fitting RT corrections ({args.fitting_method} method)...")
    print(f"{'='*70}\n")
    
    # First, generate the RT vs RT-diff equations (these will be reused in the corrections)
    print(f"Generating RT drift pattern equations...")
    all_equations_plain = {}  # {file1_name: {file2_name: equation_string}}
    
    for file1_name in sorted(all_comparisons.keys()):
        comparisons = all_comparisons[file1_name]
        equations_plain_for_file = {}
        
        for file2_name, data_points in sorted(comparisons.items()):
            if not data_points or len(data_points) < 3:
                continue
            
            # Extract RT and RT-diff for this comparison
            rt_values = [d['rt_file1'] for d in data_points]
            rt_diffs = [d['y'] for d in data_points]
            
            # Filter by percentile-based outlier threshold
            abs_rt_diffs = [abs(d['y']) for d in data_points]
            threshold_value, percentile_used = calculate_percentile_threshold(abs_rt_diffs, args.outlier_threshold)
            filtered_rt = [d['rt_file1'] for d in data_points if abs(d['y']) <= threshold_value]
            filtered_diff = [d['y'] for d in data_points if abs(d['y']) <= threshold_value]
            
            if len(filtered_rt) >= 3:
                try:
                    # Generate equation for RT vs RT-diff using the same method as the corrections
                    if args.fitting_method == 'polynomial':
                        _, equation_plain, _ = fit_polynomial_equation_plain(filtered_rt, filtered_diff, degree=args.polynomial_degree, decimal_points=args.decimal_points)
                    elif args.fitting_method == 'spline':
                        _, _, equation_plain = fit_spline_equation(filtered_rt, filtered_diff, smoothing_factor=args.spline_smoothing, degree=args.spline_degree, decimal_points=args.decimal_points)
                    elif args.fitting_method == 'loess':
                        _, _, equation_plain = fit_loess_equation(filtered_rt, filtered_diff, smoothing_fraction=args.loess_fraction, decimal_points=args.decimal_points)
                    elif args.fitting_method == 'rbf':
                        _, _, equation_plain = fit_rbf_equation(filtered_rt, filtered_diff, rbf_function='multiquadric', epsilon=1.0, decimal_points=args.decimal_points)
                    elif args.fitting_method == 'piecewise':
                        _, _, equation_plain = fit_piecewise_polynomial(filtered_rt, filtered_diff, num_segments=5, poly_degree=2, decimal_points=args.decimal_points)
                    else:
                        equation_plain = f"(Unknown method: {args.fitting_method})"
                    equations_plain_for_file[file2_name] = equation_plain
                except Exception as e:
                    equations_plain_for_file[file2_name] = f"(Fitting failed: {str(e)})"
        
        all_equations_plain[file1_name] = equations_plain_for_file
    
    print(f"  ✓ Generated {sum(len(v) for v in all_equations_plain.values())} equations\n")
    
    # Now fit corrections using the pre-generated equations
    corrections_dict = {}  # Store corrections for JSON export
    
    for file1_name in sorted(all_comparisons.keys()):
        comparisons = all_comparisons[file1_name]
        
        for file2_name, data_points in sorted(comparisons.items()):
            if not data_points:
                continue
            
            # Extract RT (x-axis for drift pattern) and RT diff (y-axis) from data points
            x_values = [d['rt_file1'] for d in data_points]  # RT from file 1 (for drift pattern equation)
            y_values = [d['y'] for d in data_points]  # RT difference
            
            # Filter by percentile-based outlier threshold for fitting
            abs_y_values = [abs(y) for y in y_values]
            threshold_value, percentile_used = calculate_percentile_threshold(abs_y_values, args.outlier_threshold)
            filtered_pairs = [(d['rt_file1'], d['y']) for d in data_points if abs(d['y']) <= threshold_value]
            
            if len(filtered_pairs) < 3:
                print(f"  Skipping {file1_name} → {file2_name}: insufficient data after filtering")
                continue
            
            filtered_x = [p[0] for p in filtered_pairs]
            filtered_y = [p[1] for p in filtered_pairs]
            
            try:
                print(f"  Fitting {file1_name} → {file2_name}...", end='', flush=True)
                
                # Fit model using specified method
                model, html_eq, plain_eq, serialization_data = fit_retention_time_model(
                    filtered_x, filtered_y,
                    fitting_method=args.fitting_method,
                    polynomial_degree=args.polynomial_degree,
                    spline_degree=args.spline_degree,
                    spline_smoothing=args.spline_smoothing,
                    loess_fraction=args.loess_fraction,
                    loess_iterations=args.loess_iterations,
                    decimal_points=args.decimal_points
                )
                
                # Get the pre-generated equation (RT vs RT-diff, not M/Z vs RT-diff)
                equation_to_save = all_equations_plain.get(file1_name, {}).get(file2_name, plain_eq)
                
                # Store correction in dict
                pair_key = f"{file1_name}_vs_{file2_name}"
                corrections_dict[pair_key] = {
                    'method': args.fitting_method,
                    'direction': 'forward',
                    'file1': file_mapping.get(file1_name, file1_name),
                    'file2': file_mapping.get(file2_name, file2_name),
                    'n_points': len(data_points),
                    'n_points_fitted': len(filtered_pairs),
                    'n_points_filtered': len(data_points) - len(filtered_pairs),
                    'serialization': serialization_data,
                    'equation': equation_to_save
                }
                
                print(f" ✓ ({len(filtered_pairs)}/{len(data_points)} points)")
                
            except Exception as e:
                print(f" ✗ (Error: {str(e)})")
    
    # Export corrections to JSON if output_json parameter is specified
    # Always write JSON (even if empty) to ensure file exists and is valid
    if args.output_json:
        try:
            if corrections_dict:
                export_fitted_corrections_json(corrections_dict, args.output_json, args.fitting_method)
            else:
                # Write skeleton JSON with empty corrections dict
                output_data = {
                    'fitting_method': args.fitting_method,
                    'corrections': {},
                    'timestamp': datetime.now().isoformat(),
                    'note': 'No corrections were found - all file pairs either had insufficient data or fitting failed'
                }
                output_dir = os.path.dirname(args.output_json)
                if output_dir:  # Only create directory if it's not empty
                    os.makedirs(output_dir, exist_ok=True)
                with open(args.output_json, 'w') as f:
                    json.dump(output_data, f, indent=2)
                print(f"⚠ No corrections found - wrote skeleton JSON to: {args.output_json}")
        except Exception as e:
            print(f"⚠ Failed to export corrections JSON: {e}")
    
    # Create HTML report with plots
    print(f"\n{'='*70}")
    print(f"Generating HTML report...")
    print(f"{'='*70}\n")
    
    figures = []
    all_equations = {}  # Store equations for each plot
    
    # Create RT vs RT-diff plots (drift patterns) with fitted curves
    # Note: corrections_dict and all_equations_plain are already populated above
    for file1_name in sorted(all_comparisons.keys()):
        comparisons = all_comparisons[file1_name]
        if comparisons:
            fig, equations, equations_plain = create_rt_vs_rtdiff_plot(
                file1_name, comparisons, 
                outlier_threshold=args.outlier_threshold, 
                polynomial_degree=args.polynomial_degree, 
                decimal_points=args.decimal_points,
                fitting_method=args.fitting_method,
                spline_degree=args.spline_degree,
                spline_smoothing=args.spline_smoothing,
                loess_fraction=args.loess_fraction,
                loess_iterations=args.loess_iterations,
                corrections_dict=corrections_dict,
                equations_plain_dict=all_equations_plain
            )
            mapped_file1_name = file_mapping.get(file1_name, file1_name)
            
            # Remap equations dict keys to use real filenames instead of basenames
            remapped_equations = {}
            for file2_name, equation in equations.items():
                mapped_file2_name = file_mapping.get(file2_name, file2_name)
                remapped_equations[mapped_file2_name] = equation
            
            figures.append((fig, f'RT Drift Pattern – {mapped_file1_name}', remapped_equations))
            all_equations[f'RT Drift Pattern – {mapped_file1_name}'] = remapped_equations
    
    if not figures:
        print("✗ No valid comparisons to plot")
        return
    
    # Display plots interactively if requested
    if args.show_plots:
        print(f"\n{'='*70}")
        print(f"Displaying plots interactively ({len(figures)} plots)...")
        print(f"{'='*70}")
        print(f"Each plot will open in your default browser. Close the browser tab to continue to the next plot.\n")
        
        for i, (fig, title, equations) in enumerate(figures, 1):
            print(f"Opening plot {i}/{len(figures)}: {title}")
            try:
                fig.show()
                print(f"  ✓ Plot {i} displayed\n")
            except Exception as e:
                print(f"  ✗ Failed to display plot {i}: {e}\n")
    
    # Create HTML report with all figures
    outlier_threshold = args.outlier_threshold
    outlier_threshold_percent = outlier_threshold * 100.0 if outlier_threshold > 0 else 0.0
    outlier_threshold_display = f"{outlier_threshold_percent:.1f}th percentile" if outlier_threshold > 0 else "No filtering"
    mz_source = args.mz_source
    rt_source = args.rt_source
    one_to_one_only = args.one_to_one_only
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Retention Time Alignment Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background-color: #f5f5f5;
                color: #333;
                line-height: 1.6;
            }}
            .header {{
                background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%);
                color: #333;
                padding: 40px 20px;
                text-align: center;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
            .header h1 {{
                font-size: 2.5em;
                margin-bottom: 10px;
            }}
            .header p {{
                font-size: 1.1em;
                opacity: 0.9;
            }}
            .container {{
                max-width: 1400px;
                margin: 0 auto;
                padding: 30px 20px;
            }}
            .plot-section {{
                background: white;
                border-radius: 8px;
                padding: 20px;
                margin-bottom: 30px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            }}
            .plot-container {{
                width: 100%;
                height: 650px;
            }}
            .stats {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
                margin-bottom: 20px;
            }}
            .stat-box {{
                background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%);
                color: #333;
                padding: 15px;
                border-radius: 6px;
                text-align: center;
            }}
            .stat-box .value {{
                font-size: 2em;
                font-weight: bold;
            }}
            .stat-box .label {{
                font-size: 0.9em;
                opacity: 0.9;
                margin-top: 5px;
            }}
            .info-box {{
                background: #FFFACD;
                border-left: 4px solid #FFD700;
                padding: 15px;
                margin-bottom: 20px;
                border-radius: 4px;
            }}
            .info-box p {{
                margin: 5px 0;
                font-size: 0.95em;
            }}
            .stats-table {{
                width: 100%;
                border-collapse: collapse;
                background: white;
                border-radius: 6px;
                overflow: hidden;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin-bottom: 30px;
            }}
            .stats-table th {{
                background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%);
                color: #333;
                padding: 15px;
                text-align: left;
                font-weight: 600;
            }}
            .stats-table td {{
                padding: 12px 15px;
                border-bottom: 1px solid #eee;
            }}
            .stats-table tr:hover {{
                background-color: #f9f9f9;
            }}
            .stats-table tr:last-child td {{
                border-bottom: none;
            }}
            .footer {{
                text-align: center;
                padding: 20px;
                color: #666;
                font-size: 0.9em;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Retention Time Alignment Analysis</h1>
            <p>Comparing feature retention times across multiple LC runs</p>
        </div>
        
        <div class="container">
            <div class="info-box">
                <p><strong>Analysis Description:</strong></p>
                <p>These plots show how retention time (RT) differences between files vary across the m/z range. 
                Each line represents a comparison against another file. Lines are interactive—hover over points for details, 
                and click legend entries to toggle visibility of individual comparisons.</p>
            </div>
            
            <div class="info-box" style="background: #F0F8FF; border-left-color: #4169E1;">
                <p><strong>Analysis Configuration:</strong></p>
                <p>Outlier threshold: <code>{outlier_threshold_display}</code></p>
                <p style="font-size: 0.9em; color: #666; margin-top: 8px;">Using percentile-based outlier detection (IQR method). Points outside the {outlier_threshold_percent:.1f}% central range are excluded from fitting but remain visible in plots.</p>
                <p style="margin-top: 12px; padding-top: 12px; border-top: 1px solid #ddd;"><strong>Coordinate Sources:</strong> M/Z = <code>{mz_source}</code>, RT = <code>{rt_source}</code></p>
                {f'<p style="margin-top: 12px;"><strong>Pep ID Filtering:</strong> Only 1:1 matches used (pep_idents appearing exactly once in both files)</p>' if one_to_one_only else ''}
            </div>
    """
    
    # Add stats
    total_comparisons = sum(len(v) for v in all_comparisons.values())
    total_matches = sum(sum(len(comp) for comp in v.values()) for v in all_comparisons.values())
    
    html_content += f"""
            <div class="stats">
                <div class="stat-box">
                    <div class="value">{len(file_data)}</div>
                    <div class="label">Input Files</div>
                </div>
                <div class="stat-box">
                    <div class="value">{total_comparisons}</div>
                    <div class="label">File Comparisons</div>
                </div>
                <div class="stat-box">
                    <div class="value">{total_matches}</div>
                    <div class="label">Total Matches</div>
                </div>
            </div>
    """
    
    # Add statistics table
    html_content += """
            <div style="background: white; border-radius: 8px; padding: 20px; margin-bottom: 30px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);">
                <h2 style="margin-bottom: 15px; color: #333;">Comparison Statistics</h2>
                <table class="stats-table">
                    <thead>
                        <tr>
                            <th>File 1</th>
                            <th>File 2</th>
                            <th>Matches</th>
                            <th>Mean RT Diff (sec)</th>
                            <th>Std Dev (sec)</th>
                        </tr>
                    </thead>
                    <tbody>
    """
    
    # Add statistics rows
    for file1_name in sorted(all_stats.keys()):
        stats_dict = all_stats[file1_name]
        for file2_name in sorted(stats_dict.keys()):
            stats = stats_dict[file2_name]
            num_matches = stats['num_matches']
            mean_rt = stats['mean_rt_diff']
            std_rt = stats['std_rt_diff']
            
            # Use mapped filenames for display
            mapped_file1_name = file_mapping.get(file1_name, file1_name)
            mapped_file2_name = file_mapping.get(file2_name, file2_name)
            
            html_content += f"""
                        <tr>
                            <td><strong>{mapped_file1_name}</strong></td>
                            <td><strong>{mapped_file2_name}</strong></td>
                            <td style="text-align: center;">{num_matches}</td>
                            <td style="text-align: right;">{mean_rt:.2f}</td>
                            <td style="text-align: right;">{std_rt:.2f}</td>
                        </tr>
            """
    
    html_content += """
                    </tbody>
                </table>
            </div>
    """
    
    # Add plots with equations
    for i, (fig, title, equations) in enumerate(figures):
        plot_html = fig.to_html(include_plotlyjs=False, div_id=f"plot_{i}")
        
        # Build equations HTML with clickable styling
        equations_html = ""
        if equations:
            equations_html = "<div style='background: #FFFACD; padding: 15px; border-radius: 6px; margin-top: 15px;'>"
            equations_html += "<p style='margin-bottom: 10px; font-weight: bold; color: #333;'>Fitted Equations (click to toggle):</p>"
            for file2_name, equation in sorted(equations.items()):
                equations_html += f"""<p style='margin: 8px 0; font-family: monospace; color: #555; cursor: pointer; padding: 8px; border-radius: 4px; transition: background 0.2s;' 
                   class='equation-toggle' 
                   data-plot-id='plot_{i}' 
                   data-fit-name='Fit: {file2_name}'
                   onmouseover="this.style.background='#FFE680'"
                   onmouseout="this.style.background='transparent'">
                   <strong>{file2_name}:</strong> {equation}</p>"""
            equations_html += "</div>"
        
        html_content += f"""
            <div class="plot-section">
                <div class="plot-container">
                    {plot_html}
                </div>
                {equations_html}
            </div>
        """
    
    html_content += """
        </div>
        
        <div class="footer">
            <p>Generated by Retention Time Alignment Analysis Tool</p>
        </div>
        
        <script>
        // Wait for Plotly to be ready, then bind toggle event handlers
        function initializeToggles() {
            document.querySelectorAll('.equation-toggle').forEach(element => {
                element.addEventListener('click', function(e) {
                    e.preventDefault();
                    e.stopPropagation();
                    
                    const plotId = this.getAttribute('data-plot-id');
                    const fitName = this.getAttribute('data-fit-name');
                    
                    // Find the Plotly plot div
                    const plotDiv = document.getElementById(plotId);
                    if (!plotDiv) {
                        console.warn(`Plot ${plotId} not found`);
                        return;
                    }
                    
                    // Get the plot data from Plotly (stored in properties when rendering)
                    // The data is stored as a property on the element
                    let plotData = plotDiv.data;
                    if (!plotData) {
                        // Try to access via Plotly's internal structure
                        // When using to_html(), data is stored differently
                        console.warn(`No data property on plot ${plotId}, trying Plotly.getPlotData()`);
                        try {
                            // Force Plotly to expose the internal plot
                            // This gets the current state of traces
                            if (window.Plotly && plotDiv._fullLayout) {
                                plotData = plotDiv._fullData || [];
                            } else {
                                return;
                            }
                        } catch (err) {
                            console.warn(`Could not access plot data: ${err}`);
                            return;
                        }
                    }
                    
                    // Find the fitted curve trace with this exact name
                    let traceIndex = -1;
                    for (let i = 0; i < plotData.length; i++) {
                        if (plotData[i].name === fitName) {
                            traceIndex = i;
                            break;
                        }
                    }
                    
                    if (traceIndex === -1) {
                        console.warn(`Trace "${fitName}" not found in plot ${plotId}. Available traces: ${plotData.map(t => t.name).join(', ')}`);
                        return;
                    }
                    
                    // Toggle visibility
                    const newVisible = !plotData[traceIndex].visible;
                    Plotly.restyle(plotDiv, {visible: newVisible}, [traceIndex]);
                    
                    // Highlight the equation when toggled
                    this.style.backgroundColor = newVisible ? '#FFE680' : 'transparent';
                    this.style.opacity = newVisible ? '0.6' : '1';
                });
                
                // Set initial styling (curves off by default)
                element.style.opacity = '0.6';
                element.style.backgroundColor = 'transparent';
            });
        }
        
        // Initialize after Plotly is loaded and DOM is fully ready
        function waitForPlotly() {
            if (window.Plotly && document.readyState !== 'loading') {
                // Give Plotly a moment to initialize all plots
                setTimeout(initializeToggles, 100);
            } else {
                setTimeout(waitForPlotly, 100);
            }
        }
        
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', waitForPlotly);
        } else {
            waitForPlotly();
        }
        </script>
    </body>
    </html>
    """
    
    # Write HTML file
    output_dir = os.path.dirname(args.output_html)
    if output_dir:  # Only create directory if a path is specified
        os.makedirs(output_dir, exist_ok=True)
    with open(args.output_html, 'w') as f:
        f.write(html_content)
    
    print(f"✓ Report generated: {args.output_html}")
    
    print(f"\nSummary:")
    print(f"  Files analyzed: {len(file_data)}")
    print(f"  Coordinate sources: M/Z = {args.mz_source}, RT = {args.rt_source}")
    print(f"  File comparisons: {total_comparisons}")
    print(f"  Total matches: {total_matches}")
    print(f"  Plots generated: {len(figures)}")
    if args.show_plots:
        print(f"  Interactive display: ENABLED (plots were shown)")
    print(f"  Outlier threshold: {args.outlier_threshold} seconds")
    if args.one_to_one_only:
        print(f"  1:1 pep_ident filter: ENABLED (only pep_idents appearing once in both files)")


if __name__ == "__main__":
    main()
