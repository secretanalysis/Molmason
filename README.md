# MolMason

**Local-first molecule generation + Pareto decision engine with uncertainty, audit logs, and schema-first route plans.**

> MolMason is the new name for the project formerly known as ChemOracle. The rename reflects the product's focus: **crafting molecules worth making** and producing **auditable, uncertainty-aware decisions** from target → candidate set → route plan.

---

## What is MolMason?

MolMason is a **constrained multi-objective molecular design desktop tool** for daily lab use. It selects *what to make next* and produces **schema-first route plans** under **Pareto optimization** with **uncertainty bands** and **closed-loop learning**.

### Key Features

- **Multi-objective optimization** — Pareto-ranked candidates via pymoo (NSGA-II)
- **Uncertainty quantification** — P10/P50/P90 bands carried end-to-end
- **Schema-first protocols** — RoutePlan/StepRecord as canonical data
- **Immutable audit trail** — Hash-chained JSONL for reproducibility
- **Air-gapped safe** — stdio transport by default, no network required
- **Lab-ready exports** — SDF, CSV, ORD JSON with full provenance

---

## Architecture

```
┌─────────────────┐
│   Qt/QML UI     │  Presentation only
└────────┬────────┘
         │ FFI
┌────────▼────────┐
│  Rust Core      │  Protocol/supervision/router/audit
└────────┬────────┘
         │ JSON-RPC (stdio)
┌────────▼────────┐
│  Python Engine  │  RDKit + generators + scorers + pymoo
└─────────────────┘
```

### Design Principles

| Principle | Implementation |
|-----------|----------------|
| **Stdio primary** | Zero ports by default (pharma/air-gapped safe) |
| **Schema-first** | RoutePlan + StepRecord are canonical |
| **Bandit routing** | One generator per action (no SMILES blending) |
| **Ranking in Python** | pymoo owns Pareto logic; Rust orchestrates only |
| **Immutable audit** | Hash-chained JSONL (no SQLite for audit) |

---

## Quick Start

### Prerequisites

- Python 3.10+
- RDKit (via conda: `conda install -c conda-forge rdkit`)
- pymoo, numpy

### Installation

```bash
cd python-engine
pip install -e ".[full]"
```

### Test RPC

```bash
# Ping
echo '{"jsonrpc":"2.0","id":1,"method":"ping","params":{}}' | python -m molmason_engine.server

# Generate molecules
echo '{"jsonrpc":"2.0","id":2,"method":"generate","params":{"seeds":["CCO"]}}' | python -m molmason_engine.server
```

---

## Directory Structure

```
MolMason/
├── CLAUDE.md                 # Authoritative spec
├── BENCHMARKS.md             # Chemical gate thresholds
├── python-engine/            # Python engine (primary)
│   └── molmason_engine/
├── crates/                   # Rust workspace
├── apps/desktop-qt/          # Qt 6 UI
├── schemas/v2/               # JSON schemas
└── benchmarks/               # Chemical gate data
```

---

## RPC Methods

| Method | Description |
|--------|-------------|
| `ping` | Health check |
| `generate` | Create molecules from seeds |
| `rank` | Pareto sort candidates |
| `cancel` | Stop running request |

---

## Chemical Gates

| Benchmark | Metric | Threshold |
|-----------|--------|-----------|
| MOSES | Validity | ≥99.0% |
| MOSES | Uniqueness | ≥95.0% |
| MOSES | Novelty | ≥0.80 |
| GuacaMol | Goal-directed | ≥0.70 |
| Property | LogP/TPSA RMSE | ≤0.20 |

---

## License

MIT
