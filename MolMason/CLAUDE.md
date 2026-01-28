# CLAUDE.md — MolMason

> Authoritative specification for Claude when working on MolMason.
> Behavioral guidelines + project architecture. Follow strictly.

---

## Part 1: Behavioral Guidelines

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

### 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them—don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

**MolMason-specific:**
- If a change might affect audit integrity, stop and confirm.
- If a generator modification could produce invalid SMILES, flag it.
- If unclear whether something belongs in Rust or Python, ask.

### 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

**MolMason-specific:**
- One generator per bandit action. No SMILES blending.
- pymoo for ranking. Don't reimplement NSGA-II.
- RDKit for validation. Don't write custom sanitization.

### 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it—don't delete it.

**MolMason-specific:**
- Don't modify audit log format without explicit approval.
- Don't change RPC contract without updating both Rust and Python.
- Don't touch chemical gate thresholds (BENCHMARKS.md is frozen).

### 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

**MolMason-specific success criteria:**
- Generator change → MOSES validity ≥99% still passes
- RPC change → `echo '{"jsonrpc":"2.0","id":1,"method":"ping","params":{}}' | python -m molmason_engine.server` returns ok
- Scoring change → Property RMSE ≤0.20 still holds

---

## Part 2: Project Architecture

### Overview

**MolMason** is a **local-first molecule generation + decision engine** that selects *what to make next* and produces **schema-first route plans** under **multi-objective (Pareto-first) optimization** with **uncertainty** and **closed-loop learning**.

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

| Principle | Implementation | Rationale |
|-----------|----------------|-----------|
| **Stdio primary** | Zero ports by default | Pharma/air-gapped safe |
| **Schema-first** | RoutePlan + StepRecord canonical | Markdown is renderer, not source |
| **Bandit routing** | One generator per action | No SMILES blending across generators |
| **Pareto + ranking in Python** | Python engine is source of truth (pymoo/NSGA-style). Rust must not implement ranking. | Avoid duplicate implementations |
| **Immutable audit** | Hash-chained JSONL | No SQLite dual-write bugs |
| **Uncertainty bands** | P10/P50/P90 carried end-to-end | Honest confidence intervals |

### Directory Layout

```
MolMason/
├── CLAUDE.md                 # THIS FILE (authoritative)
├── BENCHMARKS.md             # Chemical gate thresholds (frozen)
├── python-engine/            # PRIMARY Python engine
│   └── molmason_engine/
├── crates/                   # Rust workspace
├── apps/desktop-qt/          # Qt 6 UI
├── schemas/v2/               # JSON schemas (source of truth)
└── benchmarks/               # Chemical gate data
```

---

## Part 3: Priorities

```
P0  config/observability     Logging, env vars, user config
P1  RPC hello-molecule       ping + generate working E2E
P2  immutable audit log      Hash-chained JSONL, crash recovery
P3  chemical gate            MOSES/GuacaMol/property, CI hard-fail
P4  exports                  SDF, CSV, ORD JSON with provenance
P5  UI safety                Error handling, validation feedback
P6  CLI                      Command-line for scripting
P7  packaging + wizards      Installers, first-run setup
```

**Gate features by priority.** If asked about P5 while P2 is incomplete, flag it.

---

## Part 4: RPC Contract

### Transport

| Mode | Binding | Use Case |
|------|---------|----------|
| **stdio** (default) | N/A | Production, air-gapped |
| **HTTP** (dev only) | 127.0.0.1 ONLY | Debugging |

HTTP is **disabled** if `MOLMASON_RELEASE=1` unless `MOLMASON_HTTP_ENABLE=1`.

### Methods

| Method | Progress | Cancellable | Description |
|--------|----------|-------------|-------------|
| `ping` | No | No | Health check |
| `generate` | Yes | Yes | Create molecules from seeds |
| `rank` | Yes | No | pymoo Pareto sorting |
| `cancel` | No | N/A | Stop running request |

---

## Part 5: Chemical Gates

### Phase-Scoped Gates

| Phase | Benchmark | Metric | Gate Type | Threshold |
|------:|-----------|--------|----------|-----------|
| 1 | MOSES | Validity | **Hard fail** | ≥99.0% |
| 1 | MOSES | Uniqueness | **Hard fail** | ≥95.0% |
| 1 | MOSES | Novelty | **Warn** | ≥0.70 |
| 2+ | MOSES | Novelty | **Hard fail** | ≥0.80 |
| 2+ | GuacaMol | Goal-directed | **Hard fail** | ≥0.70 |

---

## Part 6: Audit System

**Only durable store:** per-session hash-chained JSONL (append-only).

### SQLite Policy (Allowed Only as Derived Cache)

SQLite may be used only for **derived, reconstructible** data (UI indexes, caches). **Prohibited** for authoritative state.

### Rules

- **Never** modify existing log entries.
- **Always** include prev_hash for chain integrity.
- **Invariant:** every `generate` action must record exactly one `generator_id`.

---

## Part 7: Quick Reference

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `MOLMASON_RELEASE` | `0` | Release mode |
| `MOLMASON_HTTP_ENABLE` | `0` | Force HTTP in release |
| `MOLMASON_LOG_DIR` | `~/.molmason/logs` | Audit logs |
| `MOLMASON_DEBUG` | `0` | Verbose logging |

### Key Commands

```bash
# Start server (stdio)
python -m molmason_engine.server

# Test ping
echo '{"jsonrpc":"2.0","id":1,"method":"ping","params":{}}' | python -m molmason_engine.server

# Run tests
cd python-engine && pytest tests/ -v

# Run chemical gates
./tools/benchmarks/run_gates.sh molecules.sdf
```

---

**End of CLAUDE.md**
