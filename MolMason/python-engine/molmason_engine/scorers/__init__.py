"""
MolMason - Property Scorers
============================
"""

from typing import Dict


def score_candidate(smiles: str) -> Dict[str, float]:
    """Score a candidate molecule."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
        }
    except ImportError:
        return {}
