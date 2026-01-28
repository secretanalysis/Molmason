"""Tests for MolMason server."""

import json
import pytest
from molmason_engine.server import MolMasonServer


@pytest.fixture
def server():
    return MolMasonServer()


def test_ping(server):
    request = {"jsonrpc": "2.0", "id": 1, "method": "ping", "params": {}}
    response = server.handle_request(request)
    
    assert response["id"] == 1
    assert "result" in response
    assert response["result"]["status"] == "ok"
    assert response["result"]["product_name"] == "MolMason"


def test_generate_simple(server):
    request = {
        "jsonrpc": "2.0",
        "id": 2,
        "method": "generate",
        "params": {"seeds": ["C"]}
    }
    response = server.handle_request(request)
    
    assert response["id"] == 2
    assert "result" in response
    assert response["result"]["total"] > 0


def test_generate_no_seeds(server):
    request = {"jsonrpc": "2.0", "id": 3, "method": "generate", "params": {}}
    response = server.handle_request(request)
    
    assert "error" in response


def test_unknown_method(server):
    request = {"jsonrpc": "2.0", "id": 4, "method": "unknown", "params": {}}
    response = server.handle_request(request)
    
    assert "error" in response
    assert response["error"]["code"] == -32601
