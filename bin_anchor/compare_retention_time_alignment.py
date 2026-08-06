#!/usr/bin/env python3
"""
Compare retention time alignment across multiple feature files.

For every ordered pair of files, features are matched by shared pep_ident (and
closest m/z / RT), the RT difference `f2_rt - f1_rt` is computed, and a correction
model `model(rt_file1) = rt_diff` is fitted per pair. The fitted models are written
to `fitted_corrections.json` (consumed downstream by the pairing step) alongside a
human-facing interactive HTML report.

Useful for analyzing retention time drifts between LC runs or alignment quality.

Sections (in file order):
  1. Formatting & outlier helpers
  2. Curve-fitting models      - per-method fitters, picklable wrappers, dispatcher
  3. Feature loading & coordinate access
  4. Feature matching & pairwise comparison
  5. HTML report (plots)
  6. Corrections export (JSON)
  7. Entry point / orchestration (main)
"""

import argparse
import json
import glob
import os
from typing import Any, Dict, List, Tuple
from collections import defaultdict
from datetime import datetime
import pickle
import base64
import plotly.graph_objects as go
import numpy as np
from scipy.interpolate import UnivariateSpline, Rbf
try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


# =============================================================================
# 1. FORMATTING & OUTLIER HELPERS
# Pure helpers: format polynomial coefficients for display, and compute a
# symmetric percentile-trim threshold that excludes RT-diff outliers before fitting.
# =============================================================================

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
    Calculate an outlier threshold by symmetric percentile trimming.

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
    
    # Trim symmetrically: split the excluded fraction between both tails.
    # For percentile=0.95, calculate the 97.5th percentile to exclude top 5% outliers
    # This is equivalent to keeping the central (percentile*100)% of data
    lower_percentile = (1.0 - percentile) / 2.0 * 100  # e.g., 2.5 for percentile=0.95
    upper_percentile = (100 - lower_percentile)          # e.g., 97.5 for percentile=0.95
    
    lower_bound = np.percentile(values, lower_percentile)
    upper_bound = np.percentile(values, upper_percentile)
    
    # Use the maximum of the two bounds (should be symmetric around 0 for RT diffs)
    threshold = max(abs(lower_bound), abs(upper_bound))
    
    return threshold, percentile * 100.0


# =============================================================================
# 2. CURVE-FITTING MODELS
# One fitter per method (polynomial / spline / loess / rbf / piecewise), two small
# picklable model wrappers (LOESS + piecewise have no native callable form), and
# fit_retention_time_model(), the dispatcher that also serializes the fit for JSON.
# Each fitter takes (x=rt_file1, y=rt_diff); all but the polynomial one return
# (model, html_eq, plain_eq), where those two strings are identical except for RBF.
#
# Fitting is the expensive part of this script, so every fit happens exactly once:
# Phase 5 fits each ordered file pair and the resulting model is reused for the
# JSON equation, the pickled payload, and the report curve. Rendering is kept
# separate from fitting - format_polynomial_equation() turns already-fitted
# coefficients into either markup form without refitting.
# =============================================================================

def format_polynomial_equation(coefficients, degree: int, decimal_points: int = 35, html: bool = True) -> str:
    """
    Render fitted polynomial coefficients as an equation string.

    Two renderings of the same numbers: `html=True` gives `x<sup>2</sup>` for the
    report, `html=False` gives `*x^2` for the JSON and console. Taking the
    coefficients as an argument means a caller that needs both forms fits once.

    Args:
        coefficients: Polynomial coefficients, highest power first (np.poly1d order)
        degree: Degree the fit actually used
        decimal_points: Decimal places for formatting coefficients
        html: Superscript markup if True, plain ASCII if False

    Returns:
        Equation string, e.g. "y = 1.5*x^2 - 0.2*x + 3.0"
    """
    terms = []
    for i, coeff in enumerate(coefficients):
        power = degree - i
        coeff_str = _format_coefficient(coeff, decimal_points)
        if power == 0:
            terms.append(coeff_str)
        elif power == 1:
            terms.append(f"{coeff_str}x" if html else f"{coeff_str}*x")
        else:
            terms.append(f"{coeff_str}x<sup>{power}</sup>" if html else f"{coeff_str}*x^{power}")

    return ("y = " + " ".join(terms)).replace("+-", "-")


def fit_polynomial_equation(x_values: List[float], y_values: List[float], degree: int = 8, decimal_points: int = 35) -> Tuple[np.poly1d, str, int]:
    """
    Fit a polynomial to data and return the polynomial object and equation string.

    Args:
        x_values: X coordinates
        y_values: Y coordinates
        degree: Polynomial degree (default 8)
        decimal_points: Decimal places for formatting coefficients (default 35)

    Returns:
        Tuple of (polynomial object, HTML equation string, actual_degree_used)
    """
    if len(x_values) < degree + 1:
        # Not enough points, reduce degree
        degree = max(1, len(x_values) - 1)

    coeffs = np.polyfit(x_values, y_values, degree)
    poly = np.poly1d(coeffs)

    return poly, format_polynomial_equation(coeffs, degree, decimal_points, html=True), degree


def fit_spline_equation(x_values: List[float], y_values: List[float], smoothing_factor: float = 100, degree: int = 3) -> Tuple[object, str, str]:
    """
    Fit a spline of varying degree to data using UnivariateSpline.
    
    Args:
        x_values: X coordinates (will be sorted internally)
        y_values: Y coordinates (corresponding to x_values)
        smoothing_factor: Smoothing factor for spline (higher = smoother)
        degree: Polynomial degree of spline (default 3 for cubic)
    
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

    # Cap the degree at what the data allows and at UnivariateSpline's hard limit of 5.
    adjusted_degree = min(degree, 5, len(x_unique) - 1)
    spline = UnivariateSpline(x_unique, y_unique, s=smoothing_factor, k=adjusted_degree)

    # Generate knot information for display
    knots = spline.get_knots()

    html_eq = f"Spline (k={adjusted_degree}, s={smoothing_factor:.2f}), {len(knots)} knots"
    plain_eq = f"Spline (k={adjusted_degree}, s={smoothing_factor:.2f}), {len(knots)} knots"
    
    return spline, html_eq, plain_eq


def fit_loess_equation(x_values: List[float], y_values: List[float], smoothing_fraction: float = 0.3, iterations: int = 3) -> Tuple[object, str, str]:
    """
    Fit LOESS (Locally Weighted Regression) to data using statsmodels.

    Args:
        x_values: X coordinates
        y_values: Y coordinates
        smoothing_fraction: Fraction of data used for each local regression (0-1); smaller = more local
        iterations: Number of robustifying iterations passed to lowess (default 3)

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
    loess_model = lowess(y_sorted, x_sorted, frac=smoothing_fraction, it=iterations)

    # loess_model is an ndarray: [[x1, y1_fitted], [x2, y2_fitted], ...]
    # Use module-level LOESSInterpolator class (supports pickling)
    loess_interpolator = LOESSInterpolator(loess_model[:, 0], loess_model[:, 1])

    html_eq = f"LOESS (frac={smoothing_fraction:.3f}, it={iterations}), {len(loess_model)} fitted points"
    plain_eq = f"LOESS (frac={smoothing_fraction:.3f}, it={iterations}), {len(loess_model)} fitted points"

    return loess_interpolator, html_eq, plain_eq


def fit_rbf_equation(x_values: List[float], y_values: List[float], rbf_function: str = 'multiquadric', epsilon: float = 1.0) -> Tuple[object, str, str]:
    """
    Fit Radial Basis Function to data using scipy.interpolate.Rbf.
    
    Args:
        x_values: X coordinates (will be sorted internally)
        y_values: Y coordinates (corresponding to x_values)
        rbf_function: RBF function type (multiquadric, inverse_multiquadric, gaussian, etc.)
        epsilon: RBF shape parameter
    
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
# --- Picklable model wrappers -------------------------------------------------
# LOESS and piecewise fits have no single evaluable/serializable form of their own,
# so we wrap them in module-level classes. They MUST stay at module scope with
# stable names/interfaces: the pairing script imports them to unpickle the models
# stored in fitted_corrections.json.

class LOESSInterpolator:
    """LOESS interpolator that can be pickled."""
    def __init__(self, x_data, y_data):
        # lowess() returns only sampled (x, fitted_y) points, not a callable model.
        # Store them (already sorted ascending by x) so we can evaluate anywhere.
        self.x_data = np.array(x_data)
        self.y_data = np.array(y_data)

    def __call__(self, x_val):
        """Evaluate LOESS interpolator at x_val (scalar or array)."""
        # Dispatch: array/list in -> array out (HTML curve); scalar in -> scalar out (pairing).
        if isinstance(x_val, (list, np.ndarray)):
            return np.array([self._eval_scalar(x) for x in np.atleast_1d(x_val)])
        else:
            return self._eval_scalar(x_val)

    def _eval_scalar(self, x_val):
        """Evaluate at a single scalar value."""
        # Outside the fitted range: clamp to the endpoint (flat extrapolation, no runaway slope).
        if x_val <= self.x_data[0]:
            return self.y_data[0]
        if x_val >= self.x_data[-1]:
            return self.y_data[-1]
        # Inside: find the bracketing pair x_data[idx-1] < x_val <= x_data[idx]...
        idx = np.searchsorted(self.x_data, x_val)
        # ...and linearly interpolate between those two fitted points.
        x0, x1 = self.x_data[idx-1], self.x_data[idx]
        y0, y1 = self.y_data[idx-1], self.y_data[idx]
        # Duplicate fitted x's (lowess does not dedupe) would divide by zero here.
        if x1 == x0:
            return y0
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


def fit_piecewise_polynomial(x_values: List[float], y_values: List[float], num_segments: int = 5, poly_degree: int = 2) -> Tuple[object, str, str]:
    """
    Fit piecewise polynomials (split x range into segments, fit polynomial to each).
    
    Args:
        x_values: X coordinates
        y_values: Y coordinates
        num_segments: Number of segments to split into
        poly_degree: Polynomial degree for each segment
    
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


# --- Unified fitting dispatcher -----------------------------------------------
# Picks the fitter for `fitting_method`, and builds the JSON-serializable payload:
# polynomial -> rounded coefficient list; everything else -> base64 pickle of the model.

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
    serialization_data = {}

    if fitting_method == 'polynomial':
        model, html_eq, degree = fit_polynomial_equation(x_values, y_values, degree=polynomial_degree, decimal_points=decimal_points)
        # Same coefficients, ASCII rendering - no refit needed.
        plain_eq = format_polynomial_equation(model.coefficients, degree, decimal_points, html=False)
        # Polynomial serializes as its plain coefficient list (rounded to 40 decimals max).
        def round_coefficient(c: float) -> float:
            return round(float(c), 40)
        serialization_data['coefficients'] = [round_coefficient(c) for c in model.coefficients]
        serialization_data['degree'] = int(degree)

    elif fitting_method == 'spline':
        model, html_eq, plain_eq = fit_spline_equation(x_values, y_values, smoothing_factor=spline_smoothing, degree=spline_degree)
        # Spline/RBF/LOESS/piecewise serialize as a base64-encoded pickle of the model.
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')

    elif fitting_method == 'loess':
        model, html_eq, plain_eq = fit_loess_equation(x_values, y_values, smoothing_fraction=loess_fraction, iterations=loess_iterations)
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')

    elif fitting_method == 'rbf':
        model, html_eq, plain_eq = fit_rbf_equation(x_values, y_values, rbf_function='multiquadric', epsilon=1.0)
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')

    elif fitting_method == 'piecewise':
        model, html_eq, plain_eq = fit_piecewise_polynomial(x_values, y_values, num_segments=5, poly_degree=2)
        model_bytes = pickle.dumps(model)
        serialization_data['model_pickle'] = base64.b64encode(model_bytes).decode('utf-8')

    else:
        raise ValueError(f"Unknown fitting method: {fitting_method}")

    return model, html_eq, plain_eq, serialization_data


# =============================================================================
# 3. FEATURE LOADING & COORDINATE ACCESS
# Read the per-file feature JSON, and resolve which numeric field a coordinate
# comes from. get_coordinate() maps the 'start'/'center'/'geo_center' source (set
# by the .nf) to the actual feature key; resolve_mz_collisions() finds, per pep,
# the minimum m/z decimal precision that keeps peptides distinguishable.
# =============================================================================

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


def get_coordinate(feature: Dict, coord_type: str, source: str) -> float:
    """
    Extract a coordinate value from a feature based on type and source.

    Args:
        feature: Feature dictionary
        coord_type: 'mz' or 'rt'
        source: 'geo_center', 'center', or 'start'

    Returns:
        Coordinate value

    Raises:
        ValueError: if the requested key is missing or the value is 0.0. A real m/z or RT
        is always > 0, so a missing/zero coordinate means malformed input, and silently
        defaulting it would corrupt the fitted corrections; fail loudly instead.
    """
    if source == 'geo_center':
        key = 'x_center_geo' if coord_type == 'mz' else 'y_center_geo'
    elif source == 'center':
        key = 'x_center' if coord_type == 'mz' else 'y_center'
    elif source == 'start':
        key = 'mz_start' if coord_type == 'mz' else 'rt_start'
    else:
        raise ValueError(f"Unknown coordinate source: {source}")

    value = feature.get(key)
    if value is None or value == 0.0:
        raise ValueError(
            f"Missing or zero {coord_type} coordinate (key '{key}', source '{source}') "
            f"for feature {feature.get('openms_fid', feature.get('filename', '<unknown>'))}"
        )
    return value


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


# =============================================================================
# 4. FEATURE MATCHING & PAIRWISE COMPARISON
# match_features_by_closest_mz_and_rt() greedily pairs same-pep features across two
# files by weighted (m/z, RT) distance; compare_files() runs that over all shared
# peptides and returns the matched RT differences (+ stats) that feed the fit.
# =============================================================================

def match_features_by_closest_mz_and_rt(features1: List[Dict], features2: List[Dict],
                                        mz_weight: float = 1.0, rt_weight: float = 0.001,
                                        mz_source: str = 'geo_center', rt_source: str = 'geo_center') -> List[Tuple[Dict, Dict]]:
    """
    Match features from file1 to file2 by finding closest neighbors in (m/z, RT) space.
    Each feature in file2 is used only once.

    Uses a weighted Euclidean distance. RT is down-weighted by default (rt_weight << mz_weight,
    mirroring params.mmf_mz_rt_weight_ratio) because the RT axis spans a much larger range and
    drifts far more than m/z between runs, so matching is driven mainly by m/z.

    Args:
        features1: List of feature dicts from file1, all with same pep_ident
        features2: List of feature dicts from file2, all with same pep_ident
        mz_weight: Weight for the m/z difference in the distance
        rt_weight: Weight for the RT difference in the distance (small by default; see above)
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
                  file1_name: str, file2_name: str, one_to_one_only: bool = True,
                  mz_source: str = 'geo_center', rt_source: str = 'geo_center') -> Tuple[List[Dict], Dict]:
    """
    Compare two feature files and return a list of matched RT-difference records, plus statistics.
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
    
    #filling filex_by_pep with pep_idents
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
    
    # Filter to 1:1 matches if requested (default)
    if one_to_one_only:
        common_peps = {pep for pep in common_peps 
                      if len(file1_by_pep[pep]) == 1 and len(file2_by_pep[pep]) == 1}
    
    if not common_peps:
        print(f"    ⚠ No common pep_idents between {file1_name} and {file2_name}")
        return [], {}
    
    # Pick one representative feature per pep from file1, then find the minimum m/z decimal
    # precision that keeps those peptides distinguishable (used to round the display m/z).
    file1_first_per_pep = {pep: file1_by_pep[pep][0] for pep in common_peps}
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



# =============================================================================
# 5. HTML REPORT (PLOTS)
# Build the interactive Plotly report. create_rt_vs_rtdiff_plot() draws, per file,
# the raw RT-vs-RT-diff scatter for each other file plus a (re-fitted) correction
# curve overlay; generate_curve_approximation() samples a model into that curve.
# This section is display-only and never touches fitted_corrections.json.
# =============================================================================

# Number of points used to render a fitted curve as a smooth line in the HTML report.
# This is a display-only resolution — plenty for a smooth line, and independent of the
# RT magnitude (previously n_points = int(max_x) produced thousands of points and bloated
# the report by ~50x for no visible benefit).
CURVE_RENDER_POINTS = 200


def generate_curve_approximation(model: object, x_values: List[float], method: str = 'loess') -> Tuple[List[float], List[float]]:
    """
    Sample a fitted model over the x range to draw its curve in the HTML report.

    Args:
        model: Fitted model object (poly1d, UnivariateSpline, Rbf, LOESSInterpolator, etc.)
        x_values: Input x values (used only to determine the min/max range)
        method: Fitting method name for debugging

    Returns:
        Tuple of (x_approximation, y_approximation) lists
    """
    try:
        x_array = np.array(x_values)
        min_x = float(np.min(x_array))
        max_x = float(np.max(x_array))

        # Fixed render resolution (see CURVE_RENDER_POINTS).
        x_fit = np.linspace(min_x, max_x, CURVE_RENDER_POINTS)
        y_fit = model(x_fit)

        return x_fit.tolist(), y_fit.tolist()
    except Exception as e:
        raise ValueError(f"Failed to generate curve approximation for {method}: {str(e)}")


def create_rt_vs_rtdiff_plot(file1_name: str, comparisons: Dict[str, list], fitted_models: Dict[str, Tuple] = None, fitting_method: str = 'polynomial') -> Tuple[go.Figure, Dict[str, str]]:
    """
    Create a Plotly figure showing RT difference vs RT from file1.
    This reveals if the RT drift is uniform across the LC run or varies by retention time.
    Includes fitted curves for each comparison using the specified fitting method.
    
    Draws the curves from models fitted in Phase 5; it does not fit anything itself.

    Args:
        file1_name: Base file name
        comparisons: {file2_name: [{'x': mz, 'y': rt_diff, 'rt_file1': rt1, 'pep_ident': pep_id}, ...]}
        fitted_models: {file2_name: (model, html_equation, degree)} from Phase 5.
            A pair is absent when it had too few points to fit; model is None when
            the fit raised, and html_equation then carries the error text.
        fitting_method: Method the models were fitted with - selects how the curve
            is evaluated in generate_curve_approximation()

    Returns:
        Tuple of (Plotly Figure object, {file2_name: equation_string_html})
    """
    fig = go.Figure()
    equations = {}

    # Add a trace for each file comparison
    for plot_idx, (file2_name, data_points) in enumerate(sorted(comparisons.items())):
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
        
        # Reuse the model Phase 5 already fitted for this pair. A pair missing from
        # fitted_models had too few points after outlier filtering to fit.
        entry = fitted_models.get(file2_name) if fitted_models else None

        if entry is None:
            model_obj = None
            equation_text = "<span style='color: #999;'>(Not enough points after filtering)</span>"
        else:
            model_obj, html_equation, fitted_degree = entry
            if model_obj is None:
                equation_text = html_equation          # carries the fitting error
            elif fitting_method == 'polynomial':
                equation_text = f"{html_equation} <span style='font-size: 0.85em; color: #666;'>(degree {fitted_degree})</span>"
            else:
                equation_text = f"fitted curve for {file2_name}"

        # Add equation text
        if equation_text:
            equations[file2_name] = equation_text
        
        # Generate fitted curve if model is available
        if model_obj is not None:
            try:
                x_fit, y_fit = generate_curve_approximation(model_obj, rt_values, fitting_method)
                
                # Add fitted curve trace (more visible: thicker, solid line, bright color)
                colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8']
                color_idx = plot_idx % len(colors)
                
                # Curve starts hidden; the report's equation-toggle JS shows/hides it by
                # matching this trace's `name`. (No customdata needed — it was 1 repeated
                # string per point and nothing read it.)
                fig.add_trace(go.Scatter(
                    x=x_fit,
                    y=y_fit,
                    mode='lines',
                    name=f'Fit: {file2_name}',
                    line=dict(width=4, color=colors[color_idx]),
                    visible=False,
                    hoverinfo='skip'
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
    
    return fig, equations


# =============================================================================
# 6. CORRECTIONS EXPORT (JSON)
# Write the fitted per-pair correction models to disk. This is the pipeline-consumed
# artifact: the pairing step reads it back to RT-align features before pairing.
# =============================================================================

def export_fitted_corrections_json(corrections_dict: Dict[str, Dict[str, Any]],
                                  output_path: str, fitting_method: str):
    """
    Export all fitted corrections to JSON file.

    Args:
        corrections_dict: {file_pair_name: correction_record} (record holds method + serialization)
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


# =============================================================================
# 7. ENTRY POINT / ORCHESTRATION
# Parse args, load all feature files, compare every ordered file pair, fit the
# correction models, write fitted_corrections.json, and render the HTML report.
# =============================================================================

def main():
    # ----- Phase 1: parse CLI arguments -----
    # Defaults here are the standalone defaults; the .nf overrides most of them
    # explicitly (e.g. --fitting_method loess, --mz_source/--rt_source start).
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
    parser.add_argument('--fitting_method', default='loess', 
                       choices=['polynomial', 'spline', 'loess', 'rbf', 'piecewise'],
                       help='Method for fitting RT correction curves (default: polynomial)')
    parser.add_argument('--outlier_threshold', type=float, default=0.95,
                       help='Percentile threshold for outlier exclusion (0.0-1.0). 0.95 = keep central 95%% of data, exclude outer 5%% (default: 0.95). Set to 0 to disable filtering. Best set between 0.95-0.99 depending on data quality (lower RT shift -> higher value).')
    parser.add_argument('--one_to_one_only', action='store_true',
                       help='Only include pep_idents that appear exactly once in both files. True by default. Strongly advised to set True - only set to False if measurement has very few (>200) pep_idents')
    parser.add_argument('--mz_source', default='geo_center', choices=['geo_center', 'center', 'start'],
                       help='Which m/z value to use: geo_center (default), center, or start')
    parser.add_argument('--rt_source', default='geo_center', choices=['geo_center', 'center', 'start'],
                       help='Which RT value to use: geo_center (default), center, or start')
    parser.add_argument('--polynomial_degree', type=int, default=8,
                       help='Polynomial degree for fitting RT correction curves (default: 8)')
    parser.add_argument('--spline_degree', type=int, default=3,
                       help='Spline degree for curve fitting (default: 3, used with spline method; UnivariateSpline requires 1-5)')
    parser.add_argument('--spline_smoothing', type=float, default=100,
                       help='Smoothing factor for splines (default: 100, used with spline method)')
    parser.add_argument('--loess_fraction', type=float, default=0.3,
                       help='LOESS smoothing fraction - proportion of data for local regression (default: 0.3, range 0.01-1.0)')
    parser.add_argument('--loess_iterations', type=int, default=3,
                       help='LOESS iterations for robustness (default: 3)')
    parser.add_argument('--decimal_points', type=int, default=35,
                       help='Decimal places for polynomial coefficients (default: 35)')

    args = parser.parse_args()

    # ----- Phase 2: discover feature-data files -----
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
    
    # ----- Phase 3: load features -----
    # Key file_data by the REAL filename (read from each feature's 'filename' field),
    # not the on-disk basename, so downstream keys/labels match the rest of the pipeline.
    file_data = {}
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
            print(f" ✓ ({len(data)} features, real name: '{real_filename}')", flush=True)
        else:
            print(f" ✗", flush=True)

    if len(file_data) < 2:
        print(f"\n✗ Need at least 2 files for comparison")
        return
    
    print(f"\n{'='*70}")
    print(f"Comparing files ({len(file_data)} files, {len(file_data)*(len(file_data)-1)//2} pairs)")
    print(f"{'='*70}\n")
    
    # ----- Phase 4: compare every ordered file pair -----
    # For each file1, match features against every other file2 and collect the RT
    # differences. Ordered pairs (both directions), so a file is never compared to itself.
    # all_comparisons: file1 -> {file2: [ {x, y=rt_diff, rt_file1, pep_ident}, ... ]}
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
    
    # ----- Phase 5: fit correction models -----
    # ONE fit per ordered file pair. The fitted model is then reused for all three
    # things that need it: the plain-text equation stored in the JSON, the pickled
    # payload the pairing step consumes, and the curve drawn in the HTML report
    # (handed to create_rt_vs_rtdiff_plot below via fitted_models).
    #
    # All models fit RT-vs-RT-difference on the outlier-filtered points.
    print(f"\n{'='*70}")
    print(f"Fitting RT corrections ({args.fitting_method} method)...")
    print(f"{'='*70}\n")

    corrections_dict = {}   # "{file1}_vs_{file2}" -> record written to the JSON
    fitted_models = {}      # {file1_name: {file2_name: (model, html_equation, degree)}}

    for file1_name in sorted(all_comparisons.keys()):
        comparisons = all_comparisons[file1_name]
        models_for_file = {}

        for file2_name, data_points in sorted(comparisons.items()):
            if not data_points:
                continue

            # x = RT in the reference file, y = RT difference to the query file
            y_values = [d['y'] for d in data_points]

            # Drop outliers by percentile before fitting, so a few bad matches
            # cannot drag the correction curve.
            abs_y_values = [abs(y) for y in y_values]
            threshold_value, _ = calculate_percentile_threshold(abs_y_values, args.outlier_threshold)
            filtered_pairs = [(d['rt_file1'], d['y']) for d in data_points if abs(d['y']) <= threshold_value]

            if len(filtered_pairs) < 3:
                print(f"  Skipping {file1_name} → {file2_name}: insufficient data after filtering")
                continue

            filtered_x = [p[0] for p in filtered_pairs]
            filtered_y = [p[1] for p in filtered_pairs]

            try:
                print(f"  Fitting {file1_name} → {file2_name}...", end='', flush=True)

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

                # (model, html label, degree) - degree is only set for polynomial
                models_for_file[file2_name] = (model, html_eq, serialization_data.get('degree'))

                # Only 'method' + 'serialization' are consumed downstream (by pairing);
                # the other fields are inspection-only metadata.
                pair_key = f"{file1_name}_vs_{file2_name}"
                corrections_dict[pair_key] = {
                    'method': args.fitting_method,
                    'file1': file1_name,
                    'file2': file2_name,
                    'n_points': len(data_points),
                    'n_points_fitted': len(filtered_pairs),
                    'n_points_filtered': len(data_points) - len(filtered_pairs),
                    'serialization': serialization_data,
                    'equation': plain_eq
                }

                print(f" ✓ ({len(filtered_pairs)}/{len(data_points)} points)")

            except Exception as e:
                print(f" ✗ (Error: {str(e)})")
                models_for_file[file2_name] = (None, f"(Fitting failed: {str(e)})", None)

        fitted_models[file1_name] = models_for_file

    # ----- Phase 6: export corrections JSON (the pipeline-consumed artifact) -----
    # Always write the file (even with zero corrections) so the downstream path() input
    # in the .nf stays valid.
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
    
    # ----- Phase 7: render the HTML report (display-only) -----
    # Build one drift-pattern figure per file, then assemble them (with a stats table
    # and the fitted-equation toggles) into a single self-contained HTML page.
    print(f"\n{'='*70}")
    print(f"Generating HTML report...")
    print(f"{'='*70}\n")

    figures = []

    # Build one drift-pattern figure (RT vs RT-diff, with fitted-curve overlays) per file
    for file1_name in sorted(all_comparisons.keys()):
        comparisons = all_comparisons[file1_name]
        if comparisons:
            fig, equations = create_rt_vs_rtdiff_plot(
                file1_name, comparisons,
                fitted_models=fitted_models.get(file1_name),
                fitting_method=args.fitting_method
            )
            figures.append((fig, f'RT Drift Pattern – {file1_name}', equations))

    if not figures:
        print("✗ No valid comparisons to plot")
        return
    
    # Assemble the final HTML page. Below: the page shell + header/config info-box
    # (this f-string), then the stats table, then each figure with its equation
    # toggles, then the closing JS that wires the toggles to Plotly. Local vars here
    # are just short aliases used inside the f-string.
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
                <p style="font-size: 0.9em; color: #666; margin-top: 8px;">Using percentile-based outlier detection. Points outside the {outlier_threshold_percent:.1f}% central range are excluded from fitting but remain visible in plots.</p>
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

            html_content += f"""
                        <tr>
                            <td><strong>{file1_name}</strong></td>
                            <td><strong>{file2_name}</strong></td>
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
    print(f"  Outlier threshold: {args.outlier_threshold} (fraction of central data kept)")
    if args.one_to_one_only:
        print(f"  1:1 pep_ident filter: ENABLED (only pep_idents appearing once in both files)")


if __name__ == "__main__":
    main()
