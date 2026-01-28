# MolMason - Chemical Gate Benchmarks

**Phase-scoped gates** protect chemistry quality without blocking early development.

## Phase-Scoped Gate Summary

| Phase | Benchmark | Metric | Gate Type | Threshold |
|------:|-----------|--------|-----------|-----------|
| 1 | MOSES | Validity | **Hard fail** | ≥99.0% |
| 1 | MOSES | Uniqueness | **Hard fail** | ≥95.0% |
| 1 | MOSES | Novelty | Warn | ≥0.70 |
| 1 | GuacaMol | Goal-directed | Warn | ≥0.50 |
| 2+ | MOSES | Novelty | **Hard fail** | ≥0.80 |
| 2+ | GuacaMol | Goal-directed | **Hard fail** | ≥0.70 |
| 2+ | Property | LogP/TPSA RMSE | **Hard fail** | ≤0.20 |

## Running Benchmarks

```bash
./tools/benchmarks/run_gates.sh generated.sdf --phase 1
```

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | All gates pass |
| 1 | MOSES hard fail |
| 2 | GuacaMol hard fail |
| 3 | Property hard fail |
| 10 | Warning |
