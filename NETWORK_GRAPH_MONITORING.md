# Network Graph Build Monitoring

## Command
```bash
python3 bin/build_network_graph.py \
  --input_pkl results/feature_analysis/feature_data_lists/edges.pkl \
  --enable-clustering \
  --clustering-method louvain \
  --output-clusters-json results/clustering_samples/clusters_louvain.json \
  --output_image results/clustering_samples/network_louvain.svg \
  --use-graphviz-clusters \
  --layout_engine sfdp \
  --layout_use_weights
```

## Monitoring Schedule
- Check interval: Every 4 hours
- Token budget: 100k (max 50% of 200k)
- Retry strategy: Halve `--fraction` flag on crash (0.5, 0.25, 0.125, etc.)

## Status Log

### Check 1 - Initial (2026-02-24 18:06)
- **Status**: ✓ RUNNING
- **PID**: 64817
- **CPU**: 27.5%
- **Memory**: 22.7% (5.7GB)
- **Runtime**: ~1m 37s
- **Notes**: Script actively working on graph rendering. Phase 2/3 (rendering SVG)

---

### Next Check: 2026-02-24 22:06 (in ~4 hours)

---

## Retry History
(None yet - monitoring only on first run)

---

## Output Files
- Target: `results/clustering_samples/clusters_louvain.json`
- Target: `results/clustering_samples/network_louvain.svg`
- Status: Building...
