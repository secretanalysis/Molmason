# MolMason Architecture

## Overview

MolMason is a **local-first molecule generation + decision engine**.

## Stack

- **UI**: Qt 6 / QML (presentation only)
- **Core**: Rust (orchestration, audit, schema validation)
- **Engine**: Python (RDKit, pymoo, generators, scorers)

## Communication

Rust ↔ Python via JSON-RPC over stdio (air-gapped safe).

## Key Principles

1. **Ranking lives in Python** (pymoo) - Rust must not implement ranking
2. **Audit log is authoritative** - hash-chained JSONL
3. **Schema-first** - RoutePlan/StepRecord are canonical
4. **One generator per action** - no SMILES blending
