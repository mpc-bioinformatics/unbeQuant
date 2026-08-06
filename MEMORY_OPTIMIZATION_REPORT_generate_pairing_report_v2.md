# Memory Optimization Report: generate_pairing_report_v2.py

**Date:** May 30, 2026  
**Script Size:** 3,701 lines  
**Purpose:** Generate statistical summary and HTML report of feature pairing results  
**Current Issue:** Memory consumption scales linearly with input file count (particularly edges_json size)  

---

## 1. SCRIPT OVERVIEW

### Purpose & Workflow
The script processes feature pairing outputs to generate:
- Summary text report (statistics)
- CSV statistics file
- HTML visualization report (optional)

### Input Files (from Nextflow workflow)
1. **edges_json** - Network edges with distances (10.5M edges = ~1-2GB) ✓ REQUIRED for histogram
2. **graph_composition_json** - Aggregated component statistics (~10-100MB)
3. **metadata_json** - Pairing parameters and feature counts (~1-5MB)
4. **ident_lookup_json** - Feature ID → identifications mapping (~100MB-1GB)
5. **clusters_json** - Cluster data (optional, ~200MB-2GB)
6. **filtered_clusters_json** - Filtered cluster data (optional, ~200MB-2GB)
7. **resolution_iterations_json** - Optimization iterations (optional, ~10-50MB)
8. **deleted_vertices_json** - Deletion report (optional, ~50-200MB)

---

## 2. CRITICAL MEMORY BOTTLENECKS

### 🔴 **CRITICAL: Duplicate edges_json Loading**
**Location:** Lines 357-371 AND lines 377-393  
**Impact:** 100-200% memory overhead for edges_json  
**Details:**
- Line 357-371: `distances = []` → Loads full edges_json, extracts distances
- Line 377-393: `intra_group_distances = []` → **Loads the SAME edges_json AGAIN**, extracts distances

Both variables store **identical data** (all edge distances). This means:
- For a 2GB edges_json: **4GB memory consumed** (2GB + 2GB) instead of 1GB
- For 10.5M edges: **Wasteful reload of 10.5M edge objects**

**Root Cause:** Copy-paste code without deduplication.

**Optimization:** Reuse single loading - only one `distances` list needed (both histogram and stats use same data).

---

### 🔴 **CRITICAL: Complete edges_json in Memory**
**Location:** Line 362 (`json.load(f)`)  
**Impact:** 1-2GB for large datasets  
**Details:**
- Entire edges_json file loaded into memory at once
- Structure: `{charge_key: [{edge_dict}, {edge_dict}, ...]}`
- For 10.5M edges across 15 charges: Could be 2GB+

**Problem:** JSON loading doesn't stream; holds entire file in memory.

**Optimization Opportunity:** 
- Stream parsing instead of full load
- OR only keep distances, not full edge objects
- OR extract distance statistics without storing all values (online statistics)

---

### 🟡 **HIGH: ident_lookup_json Complete Load**
**Location:** Lines 422-433  
**Impact:** 100MB-1GB depending on dataset  
**Details:**
- Entire ident_lookup dict loaded with `json.load()`
- Full structure held in memory even after counting identifications
- Example structure: `{feature_id: {pep_ident: [...], prot_ident: [...]}}`

**Problem:** Full dict kept in memory for entire function execution.

**Optimization:** Stream counting or generator pattern - only need count, not full dict.

---

### 🟡 **HIGH: graph_composition Complete Load**
**Location:** Lines 139-152  
**Impact:** 10-100MB  
**Details:**
- Full graph composition dict held in memory throughout function
- Passed to `generate_html_report()` (line 547) - kept in memory during HTML generation
- Multiple accesses to same structure (lines 321, 327, 437, etc.)

**Problem:** Large dict held throughout entire function + passed to another function.

**Optimization:** Extract only needed statistics, not full structure. Pass extracted summary instead.

---

### 🟡 **HIGH: clusters_data & filtered_clusters_data Complete Load**
**Location:** Lines 201-217  
**Impact:** 200MB-2GB each  
**Details:**
- Both loaded completely with `json.load()`
- Passed to `generate_html_report()` (line 547)
- Kept in memory during HTML generation
- Never modified after loading

**Problem:** Large structures held in memory longer than necessary.

**Optimization:** Only load if HTML output is required. Stream processing if possible.

---

### 🟡 **HIGH: All Distances Stored in List**
**Location:** Line 367 (distances.append), Line 385 (intra_group_distances.append)  
**Impact:** 40-200MB depending on edge count  
**Details:**
- Entire distances list kept in memory
- For 10.5M edges with float64: ~80MB just for distances
- Entire edges_data dict ALSO kept in memory: ~1-2GB
- Both exist simultaneously during this block

**Problem:** Complete edge objects kept alongside extracted distances.

**Optimization:** 
- Stream distances directly to histogram function
- OR use online statistics calculation (no need to store all distances)
- OR use generator expressions instead of lists

---

### 🟠 **MEDIUM: HTML Generation Memory Accumulation**
**Location:** Lines 545-557 (`generate_html_report()` call)  
**Impact:** All loaded data kept in memory during HTML generation  
**Details:**
```python
generate_html_report(
    output_html,
    ... (scalars),
    graph_composition,           # ← 10-100MB dict
    clusters_data,               # ← 200MB-2GB dict
    filtered_clusters_data,      # ← 200MB-2GB dict
    intra_group_distances,       # ← 40-200MB list
    ...
)
```

**Problem:** All large structures kept in memory simultaneously while HTML is being generated (~30+ lines of complex HTML construction).

**Optimization:** 
- Pass only summarized data, not full dicts
- Stream HTML generation instead of building full strings
- Delete structures after extract summary

---

### 🟠 **MEDIUM: resolution_iterations Complete Load**
**Location:** Lines 223-243  
**Impact:** 10-50MB  
**Details:**
- Full iteration history loaded with `json.load()`
- List of 5-10 iterations with detailed metrics
- Passed to HTML generation (line 553)

**Optimization:** Extract only best iteration data, not full history (unless needed for visualization).

---

### 🟠 **MEDIUM: matplotlib Figure Objects Not Cleaned**
**Location:** `create_distance_histogram_base64()` (lines 728-760)  
**Impact:** Memory spikes during histogram generation  
**Details:**
- Creates matplotlib figure with 10.5M data points
- Figure objects can consume significant memory
- `plt.close(fig)` at line 757 helps but memory cleanup may be delayed

**Optimization:** Explicit resource cleanup, smaller datasets for visualization (sampling).

---

### 🟠 **MEDIUM: Multiple Passes Over Large Data**
**Location:** Various locations  
**Details:**
- `graph_composition` iterated multiple times:
  - Line 327: Get total_features
  - Line 437: Extract file_match_distribution
  - Line 442: Get file_match_uniqueness
  
**Problem:** Same dict accessed repeatedly, could be consolidated.

---

## 3. UNUSED OR UNNECESSARY VARIABLES

### Unused Variables:
1. **`paired_features`** (Line 315) - Initialized as `[]` but never populated
   - Comment says "component building was removed"
   - Never used, but declared anyway
   - Wasted assignment

2. **`pep_ident_stats`** (Line 316) - Initialized as `{}` but never updated
   - Never used throughout function
   - Only passed to HTML generation as empty dict

3. **`pairing_parameters`** (Line 318) - Initialized but never populated from metadata
   - Comment says "Extract pairing parameters (not in graph_composition, use defaults)"
   - Could be removed

4. **`cutoff_params`** (Line 319) - Initialized but largely unused
   - Only populated with values, never actually used in calculations

5. **`files_matched`** (Line 374) - Initialized as `[]` for comment purposes only
   - Never used

6. **`total_deleted`** (Line 348) - Set to 0, never updated
   - Passed to HTML generation (line 544)
   - Calculated differently from optimization results

### Partially Used Variables:
1. **`edges_data`** (Line 362) - Kept in scope for two separate iterations
   - Could be used in single pass
   - Garbage collected only after line 393

---

## 4. RUNTIME BOTTLENECKS

### Slow Operations:
1. **JSON parsing** (8-10 files × json.load)
   - JSON parsing is CPU-intensive for large files
   - Some files parsed multiple times

2. **String building in HTML generation**
   - Large string concatenations in `generate_html_report()`
   - Could use StringIO for efficiency

3. **ident_lookup iteration** (Line 425-431)
   - Full dict iteration over millions of features
   - Only counting, could use streaming

4. **matplotlib histogram creation**
   - Plotting 10.5M data points
   - Very slow for visualization
   - Consider sampling

---

## 5. OPTIMIZATION OPPORTUNITIES (Ranked by Impact)

### Priority 1: CRITICAL (40-50% memory savings)
1. **Eliminate duplicate edges_json loading**
   - Merge lines 357-371 and 377-393 into single pass
   - Single `distances` list used for both histogram and stats
   - **Expected Savings:** 1-2GB

2. **Stream ident_lookup instead of full load**
   - Count identifications without storing full dict
   - **Expected Savings:** 100MB-1GB
   
3. **Only extract needed data from graph_composition**
   - Create summary dict with only needed fields
   - Don't pass full 100MB dict to HTML generation
   - **Expected Savings:** 10-100MB

### Priority 2: HIGH (20-30% memory savings)
4. **Use online statistics instead of storing all distances**
   - Calculate mean, min, max, percentiles without storing array
   - Only store sample for histogram (e.g., 100k samples)
   - **Expected Savings:** 40-200MB

5. **Lazy load clusters_data/filtered_clusters_data**
   - Only load if HTML output requested
   - Stream processing instead of full load
   - **Expected Savings:** 200MB-2GB

6. **Delete/unref large structures after use**
   - Explicitly delete edges_data after line 393
   - Delete ident_lookup after counting
   - Delete clusters_data/filtered_clusters_data after summary extraction
   - **Expected Savings:** 1-3GB (cumulative memory reduction during execution)

### Priority 3: MEDIUM (10-20% memory savings)
7. **Sample distances for matplotlib histogram**
   - Instead of plotting 10.5M points, plot random sample of 100k
   - Histogram shape identical, much faster rendering
   - **Expected Savings:** Memory spike reduction + runtime

8. **Use StringIO for HTML generation**
   - Replace string concatenation with StringIO
   - **Expected Savings:** 10-50MB

9. **Consolidate multiple graph_composition accesses**
   - Single pass to extract total_features, file_match_distribution, file_match_uniqueness
   - **Expected Savings:** 5-10% runtime improvement

---

## 6. SPECIFIC CODE LOCATIONS & FIXES

### Fix #1: Eliminate Duplicate edges_json Load
```python
# CURRENT (LINES 357-393):
distances = []
if edges_json and os.path.exists(edges_json):
    with open(edges_json, 'r') as f:
        edges_data = json.load(f)
    for charge_key, edges_list in edges_data.items():
        for edge in edges_list:
            if 'distance' in edge:
                distances.append(edge['distance'])  # ← Load #1

intra_group_distances = []
if edges_json and os.path.exists(edges_json):
    with open(edges_json, 'r') as f:
        edges_data = json.load(f)  # ← Load #2 - DUPLICATE!
    for charge_key, edges_list in edges_data.items():
        for edge in edges_list:
            if 'distance' in edge:
                intra_group_distances.append(edge['distance'])  # Same data

# OPTIMIZED:
distances = intra_group_distances = []
total_paired_groups = 0
if edges_json and os.path.exists(edges_json):
    with open(edges_json, 'r') as f:
        edges_data = json.load(f)
    for charge_key, edges_list in edges_data.items():
        for edge in edges_list:
            if 'distance' in edge:
                distances.append(edge['distance'])  # Single pass
            total_paired_groups += 1
    # Delete edges_data immediately after extraction
    del edges_data
    import gc
    gc.collect()
```

**Expected Memory Savings:** 1-2GB

---

### Fix #2: Stream ident_lookup Instead of Full Load
```python
# CURRENT (LINES 422-433):
if ident_lookup_json and os.path.exists(ident_lookup_json):
    with open(ident_lookup_json, 'r') as f:
        ident_lookup = json.load(f)  # ← Full dict in memory
    for feature_id, feature_data in ident_lookup.items():
        if isinstance(feature_data, dict):
            pep_ident = feature_data.get('pep_ident', [])
            prot_ident = feature_data.get('prot_ident', [])
            if (pep_ident and len(pep_ident) > 0) or (prot_ident and len(prot_ident) > 0):
                features_with_ident += 1

# OPTIMIZED (STREAMING):
features_with_ident = 0
if ident_lookup_json and os.path.exists(ident_lookup_json):
    try:
        with open(ident_lookup_json, 'r') as f:
            for line in f:
                feature_data = json.loads(line)
                # ... process line ...
                features_with_ident += 1
    except: pass  # Handle streaming errors

# Or use ijson for streaming JSON parsing
```

**Expected Memory Savings:** 100MB-1GB

---

### Fix #3: Extract Summary from graph_composition
```python
# CURRENT (LINE 321, 437, 442):
composition_summary = graph_composition.get('composition_summary', {})
total_features = composition_summary.get('total_features', 0)

# Later: full graph_composition passed to HTML generation

# OPTIMIZED:
# Extract needed fields once
graph_summary = {
    'total_features': graph_composition.get('composition_summary', {}).get('total_features', 0),
    'file_match_distribution': extract_from_composition(...),
    'file_match_uniqueness': graph_composition.get('file_match_uniqueness', {}),
}
# Delete full graph_composition after extraction
del graph_composition
gc.collect()

# Pass summary instead to generate_html_report
```

**Expected Memory Savings:** 10-100MB

---

## 7. WORKFLOW SCALABILITY ANALYSIS

### Current Memory Profile (Estimated)
| Input Size | Edges JSON | ident_lookup | clusters | Total Memory |
|-----------|-----------|-----------|---------|------------|
| Small (1-5 files) | 100MB | 50MB | 50MB | ~300MB |
| Medium (10-20 files) | 500MB | 300MB | 300MB | ~1.5GB |
| **Large (50+ files)** | **2GB** | **1GB** | **2GB** | **~7-8GB** |

### With Optimizations (Estimated)
| Input Size | Edges JSON | ident_lookup | clusters | Total Memory |
|-----------|-----------|-----------|---------|------------|
| Small (1-5 files) | 50MB | 5MB | 0MB* | ~100MB |
| Medium (10-20 files) | 250MB | 50MB | 0MB* | ~400MB |
| **Large (50+ files)** | **1GB** | **100MB** | **0MB*** | **~1.2GB** |

*Only loaded if HTML output requested

### Scalability Impact:
- **Current:** Memory usage grows exponentially with file count (7-8GB for 50+ files)
- **Optimized:** Near-linear growth (1-1.2GB for 50+ files)
- **Improvement:** 6-7x memory reduction for large datasets

---

## 8. SUMMARY TABLE: Optimization Recommendations

| Issue | Type | Memory Impact | Runtime Impact | Effort | Priority |
|-------|------|---------------|-----------------|--------|----------|
| Duplicate edges_json load | CRITICAL | 1-2GB | -5-10s | LOW | P1 |
| Stream ident_lookup | HIGH | 100MB-1GB | +5-15s* | MEDIUM | P1 |
| Graph composition summary | HIGH | 10-100MB | -2-5s | LOW | P1 |
| Online statistics (no storage) | MEDIUM | 40-200MB | +10-20s** | MEDIUM | P2 |
| Lazy load clusters data | HIGH | 200MB-2GB | -5-10s | MEDIUM | P2 |
| Explicit memory cleanup | MEDIUM | +500MB-1GB peak | negligible | LOW | P2 |
| Histogram sampling | MEDIUM | 200MB spike reduction | -10-30s | LOW | P3 |
| StringIO for HTML | LOW | 10-50MB | +2-5s | LOW | P3 |

*Streaming may be slightly slower than direct iteration  
**Online stats calculation may require more CPU (marginal)

---

## 9. RECOMMENDATIONS

### Immediate Actions (No Breaking Changes):
1. ✅ **Fix duplicate edges_json loading** → 1-2GB saved
2. ✅ **Add explicit cleanup** → 500MB-1GB peak reduction
3. ✅ **Extract graph_composition summary** → 10-100MB saved

### Short-term (Medium Effort):
4. ✅ **Implement streaming for ident_lookup** → 100MB-1GB saved
5. ✅ **Lazy-load clusters data** → 200MB-2GB conditional savings
6. ✅ **Sample distances for histogram** → Runtime improvement

### Long-term (Architecture Change):
7. ⏰ **Stream JSON parsing with ijson** → Per-file memory bounded
8. ⏰ **Generator-based report generation** → Constant memory footprint
9. ⏰ **Split into separate processes** → Parallel processing of independent analyses

---

## 10. EXPECTED OUTCOMES

### After Priority 1 Optimizations:
- **Memory Usage:** 7-8GB → 2-3GB (65-75% reduction)
- **Peak Memory:** During edge loading (reduced from 3-4GB to 1GB)
- **Runtime:** Slight improvement (eliminated duplicate load)
- **Scalability:** Handles 100+ input files comfortably

### After All Optimizations:
- **Memory Usage:** 7-8GB → 500MB-1GB (85-95% reduction)
- **Peak Memory:** Bounded to 1-2GB regardless of input size
- **Runtime:** 10-20% faster overall
- **Scalability:** Linear scaling with input count

---

## Conclusion

**The script's memory problems are primarily structural, not algorithmic.** Key issues are:
1. Duplicate file loading (easy fix)
2. Keeping entire data structures in memory when only summaries needed (medium effort)
3. Missing garbage collection between processing stages (easy fix)

**Estimated fixes required:** 2-3 hours for Priority 1-2 optimizations  
**Expected improvement:** 65-90% memory reduction
