"""
MolMason - Progress Event Emitter
==================================
"""

import sys
import json
import threading
from typing import TextIO


class ProgressEmitter:
    """Emit JSON-RPC progress notifications."""
    
    def __init__(self, output: TextIO = None):
        self.output = output or sys.stdout
        self._lock = threading.Lock()
    
    def emit(self, request_id: str, percent: float, message: str):
        event = {
            "jsonrpc": "2.0",
            "method": "progress",
            "params": {
                "request_id": request_id,
                "percent": min(100.0, max(0.0, percent)),
                "message": message
            }
        }
        with self._lock:
            self.output.write(json.dumps(event) + "\n")
            self.output.flush()
