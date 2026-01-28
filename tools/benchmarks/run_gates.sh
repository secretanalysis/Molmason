#!/bin/bash
# MolMason - Run Chemical Gates
set -e
INPUT="${1:-generated.sdf}"
PHASE="${MOLMASON_GATE_PHASE:-1}"

while [[ $# -gt 0 ]]; do
    case $1 in
        --phase) PHASE="$2"; shift 2 ;;
        *) INPUT="$1"; shift ;;
    esac
done

echo "MolMason Chemical Gates"
echo "Input: $INPUT, Phase: $PHASE"

# Run benchmarks
python -m molmason_engine.benchmarks.moses --input "$INPUT" --phase "$PHASE" || exit 1
echo "Gates complete"
