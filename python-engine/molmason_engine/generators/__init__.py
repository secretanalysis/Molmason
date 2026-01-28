"""
MolMason - Molecule Generators
===============================
One generator per bandit action. No SMILES blending.
"""

from typing import List, Dict, Protocol


class Generator(Protocol):
    """Generator interface."""
    
    def generate(self, seed: str, params: Dict) -> List[Dict]:
        """Generate candidates from seed."""
        ...


class PlaceholderGenerator:
    """Placeholder generator for P1 testing."""
    
    def generate(self, seed: str, params: Dict) -> List[Dict]:
        return [{"smiles": seed, "source": "placeholder"}]
