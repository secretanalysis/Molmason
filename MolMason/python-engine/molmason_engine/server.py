"""
MolMason Engine - JSON-RPC Server
==================================
Primary transport: stdio (default)
Secondary: HTTP on localhost (dev only, disabled in release)
"""

import sys
import json
import time
import os
import logging
import argparse
import threading
from typing import Dict, Any, List
from dataclasses import dataclass, field

LOG_DIR = os.environ.get("MOLMASON_LOG_DIR", os.path.expanduser("~/.molmason/logs"))
DEBUG = os.environ.get("MOLMASON_DEBUG", "0") == "1"

logging.basicConfig(
    level=logging.DEBUG if DEBUG else logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
)
logger = logging.getLogger("molmason.server")


@dataclass
class ActiveRequest:
    id: str
    method: str
    cancelled: bool = False
    checkpoint: Dict = field(default_factory=dict)


class ProgressEmitter:
    def __init__(self, output_stream=None):
        self.output = output_stream or sys.stdout
        self._lock = threading.Lock()
    
    def emit(self, request_id: str, percent: float, message: str):
        event = {
            "jsonrpc": "2.0",
            "method": "progress",
            "params": {
                "request_id": request_id,
                "percent": percent,
                "message": message
            }
        }
        with self._lock:
            self.output.write(json.dumps(event) + "\n")
            self.output.flush()


class MolMasonServer:
    def __init__(self):
        self.progress = ProgressEmitter()
        self.active_requests: Dict[str, ActiveRequest] = {}
        self._lock = threading.Lock()
        
        self.rdkit_available = False
        self.pymoo_available = False
        
        try:
            from rdkit import Chem
            self.rdkit_available = True
        except ImportError:
            pass
        
        try:
            import pymoo
            self.pymoo_available = True
        except ImportError:
            pass
    
    def handle_request(self, request: Dict) -> Dict:
        method = request.get("method")
        params = request.get("params", {})
        request_id = request.get("id")
        
        handlers = {
            "ping": self._handle_ping,
            "generate": self._handle_generate,
            "rank": self._handle_rank,
            "cancel": self._handle_cancel,
        }
        
        handler = handlers.get(method)
        if not handler:
            return self._error_response(request_id, -32601, f"Method not found: {method}")
        
        try:
            result = handler(request_id, params)
            return {"jsonrpc": "2.0", "id": request_id, "result": result}
        except Exception as e:
            logger.exception(f"Error handling {method}")
            return self._error_response(request_id, -32000, str(e))
    
    def _handle_ping(self, request_id: str, params: Dict) -> Dict:
        return {
            "status": "ok",
            "product_name": "MolMason",
            "version": "2.0.0",
            "rdkit_available": self.rdkit_available,
            "pymoo_available": self.pymoo_available,
            "timestamp": time.time()
        }
    
    def _handle_generate(self, request_id: str, params: Dict) -> Dict:
        seeds = params.get("seeds", [])
        if not seeds:
            raise ValueError("No seeds provided")
        
        with self._lock:
            self.active_requests[request_id] = ActiveRequest(id=request_id, method="generate")
        
        try:
            candidates = []
            for i, seed in enumerate(seeds):
                with self._lock:
                    req = self.active_requests.get(request_id)
                    if req and req.cancelled:
                        return {"cancelled": True, "checkpoint": {"completed_seeds": i}, "candidates": candidates}
                
                percent = (i / len(seeds)) * 100
                self.progress.emit(request_id, percent, f"Processing seed {i+1}/{len(seeds)}")
                
                seed_candidates = self._generate_from_seed(seed, params)
                candidates.extend(seed_candidates)
            
            self.progress.emit(request_id, 100, "Generation complete")
            return {"candidates": candidates, "total": len(candidates), "seeds_processed": len(seeds)}
        finally:
            with self._lock:
                self.active_requests.pop(request_id, None)
    
    def _generate_from_seed(self, seed: str, params: Dict) -> List[Dict]:
        candidates = []
        
        if self.rdkit_available:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(seed)
            if mol is None:
                return []
            seed = Chem.MolToSmiles(mol)
        
        # Placeholder generation
        modifications = [("", seed), ("C", seed + "C"), ("O", seed + "O"), ("N", seed + "N")]
        
        for mod_name, smiles in modifications:
            if self.rdkit_available:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                smiles = Chem.MolToSmiles(mol)
            
            candidates.append({
                "smiles": smiles,
                "seed": seed,
                "modification": mod_name,
                "source": "placeholder_generator"
            })
        
        return candidates
    
    def _handle_rank(self, request_id: str, params: Dict) -> Dict:
        candidates = params.get("candidates", [])
        if not candidates:
            raise ValueError("No candidates provided")
        
        if not self.pymoo_available:
            return {"ranked": candidates, "fronts": [list(range(len(candidates)))], "method": "simple"}
        
        from .pareto import rank_with_pymoo
        self.progress.emit(request_id, 10, "Computing Pareto fronts")
        ranked, fronts = rank_with_pymoo(candidates, params.get("objectives", []))
        self.progress.emit(request_id, 100, "Ranking complete")
        return {"ranked": ranked, "fronts": fronts, "method": "pymoo"}
    
    def _handle_cancel(self, request_id: str, params: Dict) -> Dict:
        target_id = params.get("request_id")
        if not target_id:
            raise ValueError("No request_id provided")
        
        with self._lock:
            req = self.active_requests.get(target_id)
            if req:
                req.cancelled = True
                return {"cancelled": True, "checkpoint": req.checkpoint}
            return {"cancelled": False, "reason": "Request not found"}
    
    def _error_response(self, request_id: Any, code: int, message: str) -> Dict:
        return {"jsonrpc": "2.0", "id": request_id, "error": {"code": code, "message": message}}


def run_stdio(server: MolMasonServer):
    logger.info("MolMason Engine starting (stdio mode)")
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            request = json.loads(line)
            response = server.handle_request(request)
            print(json.dumps(response), flush=True)
        except json.JSONDecodeError as e:
            print(json.dumps({"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": f"Parse error: {e}"}}), flush=True)


def run_http(server: MolMasonServer, port: int):
    from http.server import HTTPServer, BaseHTTPRequestHandler
    
    class Handler(BaseHTTPRequestHandler):
        def do_POST(self):
            content_length = int(self.headers.get("Content-Length", 0))
            body = self.rfile.read(content_length).decode("utf-8")
            try:
                request = json.loads(body)
                response = server.handle_request(request)
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps(response).encode("utf-8"))
            except Exception as e:
                self.send_response(500)
                self.end_headers()
        
        def log_message(self, format, *args):
            pass
    
    is_release = os.environ.get("MOLMASON_RELEASE", "0") == "1"
    http_forced = os.environ.get("MOLMASON_HTTP_ENABLE", "0") == "1"
    
    if is_release and not http_forced:
        logger.error("HTTP mode disabled in release. Set MOLMASON_HTTP_ENABLE=1 to override.")
        sys.exit(1)
    
    httpd = HTTPServer(("127.0.0.1", port), Handler)
    logger.info(f"MolMason Engine starting (HTTP mode on 127.0.0.1:{port})")
    httpd.serve_forever()


def main():
    parser = argparse.ArgumentParser(description="MolMason Engine")
    parser.add_argument("--port", type=int, help="Run HTTP server on port (dev only)")
    args = parser.parse_args()
    
    server = MolMasonServer()
    if args.port:
        run_http(server, args.port)
    else:
        run_stdio(server)


if __name__ == "__main__":
    main()
