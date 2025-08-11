"""
Simulation models for CRISPR-Cas9 editing and delivery.
"""
from typing import Dict, Any, Optional, Tuple, List
import random

# Import at function level to avoid circular imports
DELIVERY_EFFICIENCY = None

class ChromatinContext:
    """
    Models the chromatin context of a target site, which can affect editing efficiency.
    """
    
    def __init__(self, accessibility_score: float = 1.0, histone_marks: Optional[Dict[str, float]] = None) -> None:
        """
        Initialize the chromatin context.
        
        Args:
            accessibility_score: Chromatin accessibility (0-1, where 1 is most accessible)
            histone_marks: Dictionary of histone marks and their levels
        """
        self.accessibility = max(0.0, min(1.0, accessibility_score))
        self.active_marks = histone_marks or {}
    
    def modify_cutting_efficiency(self, base_efficiency: float) -> float:
        """
        Modify the base cutting efficiency based on chromatin state.
        
        Args:
            base_efficiency: Base cutting efficiency (0-1)
            
        Returns:
            Modified efficiency (0-1)
        """
        mod = 1.0
        
        # Apply histone mark effects
        if self.active_marks.get("H3K27ac", 0) > 0:
            # Active mark - increases efficiency
            mod *= 1.2
        if self.active_marks.get("H3K9me3", 0) > 0:
            # Repressive mark - decreases efficiency
            mod *= 0.7
        
        # Apply accessibility effect
        return base_efficiency * self.accessibility * mod

def cas9_activity(temp_c: float, protein_conc: float) -> float:
    """
    Calculate Cas9 activity based on temperature and protein concentration.
    
    Args:
        temp_c: Temperature in Celsius
        protein_conc: Protein concentration in nM
        
    Returns:
        Activity multiplier (0-1)
    """
    # Activity optimal at 37C, drops at extremes
    temp_factor = max(0.0, 1.0 - abs(temp_c - 37) / 20)
    
    # Protein concentration effect (saturating)
    conc_factor = min(1.0, protein_conc / 100.0)
    
    return temp_factor * conc_factor

def repair_pathway_prob(cell_cycle: str) -> float:
    """
    Get the probability of HDR vs NHEJ based on cell cycle phase.
    
    Args:
        cell_cycle: Cell cycle phase ("G1", "S/G2", or "Other")
        
    Returns:
        Probability of HDR (0-1)
    """
    if cell_cycle == "S/G2":
        return 0.7  # High HDR in S/G2
    elif cell_cycle == "G1":
        return 0.2  # Low HDR in G1
    else:
        return 0.4  # Intermediate in other phases

def simulate_delivery(organ: str, method: str) -> bool:
    """
    Simulate delivery of CRISPR components to target cells.
    
    Args:
        organ: Target organ/tissue
        method: Delivery method
        
    Returns:
        True if delivery was successful, False otherwise
    """
    # Import here to avoid circular imports and support local runs
    try:
        from CRISPRcas9_simV3.config import DELIVERY_EFFICIENCY
    except ImportError:  # Fallback when running without package install
        from config import DELIVERY_EFFICIENCY
    
    if organ not in DELIVERY_EFFICIENCY or method not in DELIVERY_EFFICIENCY[organ]:
        return False
    
    # Simple probability check based on delivery efficiency
    success_rate = DELIVERY_EFFICIENCY[organ][method]
    return random.random() < success_rate

def delivery_efficiency(base: float, dose: float, route: str) -> float:
    """
    Calculate delivery efficiency based on dose and route.
    
    Args:
        base: Base efficiency (0-1)
        dose: Dose amount
        route: Administration route ("IV", "IM", or "IP")
        
    Returns:
        Effective delivery efficiency (0-1)
    """
    # Route efficiency factors
    route_factor = {"IV": 1.0, "IM": 0.7, "IP": 0.5}.get(route, 0.5)
    
    # Dose-response: log-scale, saturating
    dose_factor = 1 - 0.5 * (1 / (1 + dose))
    
    return base * dose_factor * route_factor
