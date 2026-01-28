#!/bin/bash
# MolMason - Audit Log Invariant Checker
set -e
LOG_FILE="${1:-$HOME/.molmason/logs/latest.jsonl}"

if [ ! -f "$LOG_FILE" ]; then
    echo "ERROR: Log file not found: $LOG_FILE"
    exit 10
fi

echo "Checking audit invariants: $LOG_FILE"
python3 << EOFPY
import json
import sys

log_file = "$LOG_FILE"
errors = []
open_requests = {}

with open(log_file, 'r') as f:
    for line_num, line in enumerate(f, 1):
        if not line.strip():
            continue
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            errors.append(f"Line {line_num}: Invalid JSON")
            continue
        
        entry_type = entry.get('type')
        entry_id = entry.get('id')
        
        if entry_type == 'START':
            open_requests[entry_id] = entry
            if entry.get('action') == 'generate' and not entry.get('generator_id'):
                errors.append(f"Line {line_num}: generate missing generator_id")
        elif entry_type in ('END', 'ABORTED'):
            open_requests.pop(entry_id, None)

for req_id in open_requests:
    errors.append(f"Unresolved START: {req_id}")

if errors:
    print("ERRORS:")
    for e in errors:
        print(f"  ✗ {e}")
    sys.exit(1)
else:
    print("✓ All invariants pass")
EOFPY
