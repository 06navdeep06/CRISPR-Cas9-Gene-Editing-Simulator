"""
CRISPR-Cas9 Core Model

Contains the Cas9 class and related functions for target finding and editing simulation.
"""
from typing import List, Dict, Any, Optional, Tuple
from Bio.Seq import Seq
import random

# Import at function level to avoid circular imports
# These will be imported when needed inside the methods
IUPAC_CODES = None
PAM_PATTERNS = None
is_valid_pam = None

class Cas9:
    """
    A class representing the Cas9 protein for CRISPR gene editing simulation.
    """
    
    def __init__(self, pam_pattern: str = "NGG", max_mismatches: int = 2) -> None:
        """
        Initialize the Cas9 model.
        
        Args:
            pam_pattern: The PAM pattern this Cas9 variant recognizes (IUPAC format)
            max_mismatches: Maximum number of mismatches allowed in gRNA binding
        """
        self.pam_pattern: str = pam_pattern.upper()
        self.max_mismatches: int = max_mismatches
        self.pam_sites_log: List[Dict[str, Any]] = []

    def find_pam_sites(self, seq: str) -> List[int]:
        """
        Find all PAM sites in the given sequence.
        
        Args:
            seq: DNA sequence to search for PAM sites
            
        Returns:
            List of start positions of PAM sites
        """
        # Import here to avoid circular imports
        try:
            from CRISPRcas9_simV3.utils.sequence_utils import is_valid_pam
        except ImportError:  # Fallback when running without package install
            from utils.sequence_utils import is_valid_pam
        
        seq = seq.upper()
        pam_len = len(self.pam_pattern)
        sites = []
        for i in range(len(seq) - pam_len + 1):
            pam_seq = seq[i:i+pam_len]
            if is_valid_pam(pam_seq, self.pam_pattern):
                sites.append(i)
        return sites

    def find_targets(self, seq: str, guide_rna: str) -> List[Dict[str, Any]]:
        """
        Find all potential target sites for the given guide RNA in the sequence.
        
        Args:
            seq: DNA sequence to search in
            guide_rna: Guide RNA sequence (without PAM)
            
        Returns:
            List of target site dictionaries with position, sequence, and strand info
        """
        # Import here to avoid circular imports
        try:
            from CRISPRcas9_simV3.utils.sequence_utils import is_valid_pam, reverse_complement
        except ImportError:  # Fallback when running without package install
            from utils.sequence_utils import is_valid_pam, reverse_complement
        
        if not guide_rna:
            return []
            
        guide_rna = guide_rna.upper()
        seq = seq.upper()
        targets = []
        
        # Search on forward strand
        for i in range(len(seq) - len(guide_rna) + 1):
            target_seq = seq[i:i+len(guide_rna)]
            if target_seq == guide_rna:
                # Check for PAM at the 3' end (NGG for standard SpCas9)
                pam_pos = i + len(guide_rna)
                if pam_pos <= len(seq) - len(self.pam_pattern):
                    pam_seq = seq[pam_pos:pam_pos+len(self.pam_pattern)]
                    if is_valid_pam(pam_seq, self.pam_pattern):
                        targets.append({
                            'position': i,
                            'sequence': target_seq,
                            'pam_sequence': pam_seq,
                            'strand': '+',
                            'mismatches': 0
                        })
        
        # Search on reverse complement strand
        rev_guide = reverse_complement(guide_rna)
        for i in range(len(seq) - len(rev_guide) + 1):
            target_seq = seq[i:i+len(rev_guide)]
            if target_seq == rev_guide:
                # Check for PAM at the 5' end (CCN for standard SpCas9 on reverse strand)
                pam_pos = i - len(self.pam_pattern)
                if pam_pos >= 0:
                    pam_seq = seq[pam_pos:pam_pos+len(self.pam_pattern)]
                    # For reverse strand, we need to check the reverse complement of the PAM
                    rev_pam = reverse_complement(pam_seq)
                    if is_valid_pam(rev_pam, self.pam_pattern):
                        targets.append({
                            'position': i,
                            'sequence': target_seq,
                            'pam_sequence': pam_seq,
                            'strand': '-',
                            'mismatches': 0
                        })
        
        return targets

    def simulate_edit(
        self,
        seq: str,
        target: Dict[str, Any],
        edit_type: str,
        donor_template: Optional[str] = None
    ) -> str:
        """
        Simulate a gene edit at the target site.
        
        Args:
            seq: Original DNA sequence
            target: Target site information
            edit_type: Type of edit ("Knockout (NHEJ)" or "Knock-in/Edit (HDR)")
            donor_template: Donor template for HDR (required for knock-in)
            
        Returns:
            Edited DNA sequence
        """
        seq = str(seq)
        if edit_type == "Knockout (NHEJ)":
            # Simulate NHEJ by introducing small deletions
            cut_site = target['position'] + len(target['sequence']) // 2
            del_len = random.randint(1, 3)
            edited = seq[:cut_site] + seq[cut_site+del_len:]
            return edited
        elif edit_type == "Knock-in/Edit (HDR)":
            if donor_template:
                # Simulate HDR by inserting the donor template at the cut site
                cut_site = target['position'] + len(target['sequence']) // 2
                edited = seq[:cut_site] + donor_template + seq[cut_site:]
                return edited
            else:
                return seq
        else:
            return seq
