"""
MolMason - Pareto Ranking via pymoo
====================================
Python-side ranking (Rust must NOT implement ranking).
"""

from typing import List, Dict, Tuple
import logging

logger = logging.getLogger("molmason.pareto")


def rank_with_pymoo(
    candidates: List[Dict],
    objectives: List[str],
    directions: List[str] = None
) -> Tuple[List[Dict], List[List[int]]]:
    """
    Rank candidates using pymoo NSGA-II style non-dominated sorting.
    
    Returns:
        ranked: Candidates sorted by front, then crowding distance
        fronts: List of fronts, each front is list of candidate indices
    """
    if not candidates:
        return [], []
    
    try:
        import numpy as np
        from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
    except ImportError:
        logger.warning("pymoo not available - returning unsorted")
        return candidates, [list(range(len(candidates)))]
    
    if not objectives:
        objectives = ["score"]
    
    if directions is None:
        directions = ["max"] * len(objectives)
    
    # Build objective matrix
    n_candidates = len(candidates)
    n_objectives = len(objectives)
    
    F = np.zeros((n_candidates, n_objectives))
    
    for i, cand in enumerate(candidates):
        scores = cand.get("scores", {})
        for j, obj in enumerate(objectives):
            val = scores.get(obj, 0.0)
            # pymoo minimizes, so negate for maximization
            if directions[j] == "max":
                F[i, j] = -val
            else:
                F[i, j] = val
    
    # Non-dominated sorting
    nds = NonDominatedSorting()
    fronts = nds.do(F)
    
    # Convert to list of lists
    fronts_list = [list(front) for front in fronts]
    
    # Add front info to candidates
    for front_idx, front in enumerate(fronts_list):
        for cand_idx in front:
            candidates[cand_idx]["front"] = front_idx
    
    # Sort by front
    ranked = sorted(candidates, key=lambda c: c.get("front", 999))
    
    return ranked, fronts_list
