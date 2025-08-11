"""
CRISPR-Cas9 Gene Editing Simulator

A Streamlit application for simulating CRISPR-Cas9 gene editing experiments.
"""
from models import calculate_on_target_score
import streamlit as st
from Bio import SeqIO
import io
import random
import difflib
from typing import Optional, Dict, Any, List

# Standard library imports
import os
import sys
import re
import random
import difflib
from typing import Optional, Dict, Any, List, Tuple

# Third-party imports
import streamlit as st
import numpy as np
import pandas as pd
from Bio.Seq import Seq

# Add the project root to the Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Local application imports
try:
    from CRISPRcas9_simV3.constants import IUPAC_CODES
    from CRISPRcas9_simV3.config import PAM_PATTERNS, ORGANS, DELIVERY_EFFICIENCY, DELIVERY_METHODS
    from CRISPRcas9_simV3.models import (
        Cas9,
        calculate_on_target_score,
        off_target_score,
        pam_score,
        secondary_structure_penalty,
        ChromatinContext,
        cas9_activity,
        repair_pathway_prob,
        simulate_delivery,
        delivery_efficiency
    )
    from CRISPRcas9_simV3.utils.sequence_utils import is_valid_pam, reverse_complement, find_pam_sites
    from CRISPRcas9_simV3.utils.visualization import plot_editing_results, plot_delivery_efficiency, plot_sequence_alignment
except ImportError:
    # Fallback to relative imports if the package isn't installed
    from constants import IUPAC_CODES
    from config import PAM_PATTERNS, ORGANS, DELIVERY_EFFICIENCY, DELIVERY_METHODS
    from models import (
        Cas9,
        calculate_on_target_score,
        off_target_score,
        pam_score,
        secondary_structure_penalty,
        ChromatinContext,
        cas9_activity,
        repair_pathway_prob,
        simulate_delivery,
        delivery_efficiency
    )
    from utils.sequence_utils import is_valid_pam, reverse_complement, find_pam_sites
    from utils.visualization import plot_editing_results, plot_delivery_efficiency, plot_sequence_alignment

# --- Page Configuration ---
st.set_page_config(
    page_title="CRISPR-Cas9 Simulator",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Session State Initialization ---
if 'seq_record' not in st.session_state:
    st.session_state.seq_record = None
if 'guide_rna' not in st.session_state:
    st.session_state.guide_rna = ""
if 'edited_seq' not in st.session_state:
    st.session_state.edited_seq = None
if 'targets' not in st.session_state:
    st.session_state.targets = []

# --- Helper Functions ---
def reset_simulation():
    """Reset the simulation state."""
    st.session_state.edited_seq = None
    st.session_state.targets = []

def load_example_sequence():
    """Load an example DNA sequence."""
    example_seq = """
    >Example_Sequence
    ATGCCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
    """
    st.session_state.seq_record = "".join(example_seq.split())
    st.session_state.guide_rna = "TAGCTAGCTAGCTAGCTAGC"  # Example gRNA
    reset_simulation()

# --- Sidebar ---
st.sidebar.title("CRISPR-Cas9 Simulator")
st.sidebar.markdown("---")

# 1. Sequence Input
st.sidebar.header("1. Input DNA Sequence")
input_mode = st.sidebar.radio("Input mode", ["Paste sequence", "Upload FASTA", "Example"], 
                            on_change=reset_simulation)

if input_mode == "Paste sequence":
    raw_seq = st.sidebar.text_area("Paste DNA sequence here", height=8, 
                                 on_change=reset_simulation)
    if raw_seq:
        st.session_state.seq_record = raw_seq.replace("\n", "").replace(" ", "").upper()
        
elif input_mode == "Upload FASTA":
    fasta_file = st.sidebar.file_uploader("Upload FASTA file", type=["fasta", "fa"], 
                                        on_change=reset_simulation)
    if fasta_file:
        fasta_obj = next(SeqIO.parse(io.StringIO(fasta_file.getvalue().decode()), "fasta"))
        st.session_state.seq_record = str(fasta_obj.seq).upper()
        st.session_state.fasta_header = fasta_obj.description
        
elif input_mode == "Example":
    if st.sidebar.button("Load Example"):
        load_example_sequence()

# 2. Guide RNA
st.sidebar.header("2. Guide RNA")
guide_rna = st.sidebar.text_input("Guide RNA sequence (20 nt)", 
                                value=st.session_state.guide_rna,
                                max_chars=25,
                                on_change=reset_simulation).upper()
st.session_state.guide_rna = guide_rna

# 3. Editing Options
st.sidebar.header("3. Editing Options")
edit_type = st.sidebar.selectbox("Edit type", 
                                ["Knockout (NHEJ)", "Knock-in/Edit (HDR)"],
                                on_change=reset_simulation)

donor_template = None
if edit_type == "Knock-in/Edit (HDR)":
    donor_template = st.sidebar.text_input("Donor template (for HDR)", 
                                         max_chars=30,
                                         on_change=reset_simulation).upper()

# 4. Off-target & PAM
st.sidebar.header("4. Off-target & PAM")
max_mismatches = st.sidebar.slider("Max mismatches allowed", 0, 2, 2,
                                  on_change=reset_simulation)

# 5. Delivery Parameters
st.sidebar.header("5. Delivery Parameters")
organ = st.sidebar.selectbox("Target organ/tissue", ORGANS)
delivery_method = st.sidebar.selectbox("Delivery method", DELIVERY_METHODS)
dose = st.sidebar.slider("Dose (arbitrary units)", 1, 100, 10)
route = st.sidebar.selectbox("Administration route", ["IV", "IM", "IP"])

# 6. Experimental Parameters
st.sidebar.header("6. Experimental Parameters")
temp_c = st.sidebar.slider("Temperature (°C)", 20, 42, 37)
protein_conc = st.sidebar.slider("Cas9 protein conc (nM)", 10, 200, 100)
cell_cycle = st.sidebar.selectbox("Cell cycle phase", ["G1", "S/G2", "Other"])

# 7. Chromatin Context
st.sidebar.header("7. Chromatin Context")
chromatin_accessibility = st.sidebar.slider(
    "Chromatin accessibility (0=closed, 1=open)", 
    0.0, 1.0, 1.0, 0.05
)
histone_selected = st.sidebar.multiselect(
    "Histone marks", 
    ["H3K27ac", "H3K9me3"]
)
histone_marks = {mark: 1 for mark in histone_selected}
chromatin_ctx = ChromatinContext(
    accessibility_score=chromatin_accessibility, 
    histone_marks=histone_marks
)

# 8. PAM Input and Cas Variant
st.sidebar.header("8. PAM Sequence / Cas Variant")
cas_options = list(PAM_PATTERNS.keys()) + ["Custom"]
cas_type = st.sidebar.selectbox(
    "Select Cas protein (Cas variant)",
    cas_options,
    key="pam_variant_selectbox"
)

if cas_type == "Custom":
    pam_pattern = st.sidebar.text_input(
        "Enter custom PAM pattern (IUPAC, e.g. NGG, TTTV)", 
        value="NGG"
    ).upper()
else:
    pam_pattern = PAM_PATTERNS[cas_type]

# --- Main Content ---
st.title("CRISPR-Cas9 Gene Editing Simulator")

# Sequence Info
if st.session_state.seq_record:
    st.subheader("Sequence Information")
    col1, col2 = st.columns(2)
    
    with col1:
        st.metric("Sequence Length", f"{len(st.session_state.seq_record)} bp")
        
        # Check for ambiguous bases
        ambiguous = set(st.session_state.seq_record) - set("ACGT")
        if ambiguous:
            st.warning(f"Sequence contains ambiguous bases: {', '.join(ambiguous)}. "
                     f"This may affect PAM finding.")
    
    with col2:
        # Check for PAM sites
        test_cas9 = Cas9(pam_pattern=pam_pattern)
        pam_sites = test_cas9.find_pam_sites(st.session_state.seq_record)
        st.metric(f"PAM sites found ({pam_pattern})", len(pam_sites))
        
        if len(pam_sites) == 0:
            st.info("No PAM sites found. Try a different PAM or check your sequence.")
    
    # gRNA troubleshooting
    if guide_rna and len(guide_rna) >= 15:
        st.subheader("Guide RNA Analysis")
        
        # Find targets
        cas9 = Cas9(pam_pattern=pam_pattern, max_mismatches=max_mismatches)
        st.session_state.targets = cas9.find_targets(st.session_state.seq_record, guide_rna)
        
        # Calculate scores
        on_targets = [t for t in st.session_state.targets if t['mismatches'] == 0]
        off_targets = [t for t in st.session_state.targets if t['mismatches'] > 0]
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Total Target Sites", len(st.session_state.targets))
        with col2:
            st.metric("On-target Sites", len(on_targets))
        with col3:
            st.metric("Off-target Sites", len(off_targets))
        
        if not st.session_state.targets:
            st.warning("No target sites found for this gRNA. "
                      "Try changing the gRNA or PAM, or check the sequence orientation.")
            
            st.markdown("""
            **Tips for successful editing:**
            - Ensure your gRNA matches a region in the DNA sequence adjacent to a PAM site
            - Try using a gRNA sequence present in your DNA (copy a 20 nt region from your sequence)
            - Verify your DNA sequence is in the correct orientation (5'→3')
            - Use the 'Custom' PAM option if your Cas variant is different
            """)
        
        # Simulation button
        if st.button("Simulate CRISPR Editing"):
            with st.spinner("Running simulation..."):
                # Initialize edited sequence
                st.session_state.edited_seq = st.session_state.seq_record
                
                # Simulate delivery
                base_eff = DELIVERY_EFFICIENCY[organ][delivery_method]
                deliv_eff = delivery_efficiency(base_eff, dose, route)
                
                # Temperature/protein/cell cycle effects
                cas9_eff = cas9_activity(temp_c, protein_conc)
                hdr_prob = repair_pathway_prob(cell_cycle)
                
                # Modify by chromatin accessibility and histone marks
                effective_eff = chromatin_ctx.modify_cutting_efficiency(deliv_eff * cas9_eff)
                
                # Simulate editing for each target
                for target in on_targets:
                    if random.random() < effective_eff:  # Successful delivery and cutting
                        if edit_type == "Knock-in/Edit (HDR)" and donor_template and random.random() < hdr_prob:
                            # HDR
                            st.session_state.edited_seq = cas9.simulate_edit(
                                st.session_state.edited_seq, 
                                target, 
                                "Knock-in/Edit (HDR)", 
                                donor_template
                            )
                        else:
                            # NHEJ
                            st.session_state.edited_seq = cas9.simulate_edit(
                                st.session_state.edited_seq, 
                                target, 
                                "Knockout (NHEJ)"
                            )
                
                st.success("Simulation completed!")
    
        # Show results if available
    if st.session_state.edited_seq is not None:
        st.subheader("Simulation Results")
        
        # Delivery efficiency
        st.markdown("### Delivery Efficiency")
        base_eff = DELIVERY_EFFICIENCY[organ][delivery_method]
        deliv_eff = delivery_efficiency(base_eff, dose, route)
        effective_eff = chromatin_ctx.modify_cutting_efficiency(deliv_eff * cas9_activity(temp_c, protein_conc))
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Base Delivery Efficiency", f"{base_eff*100:.1f}%")
        with col2:
            st.metric("With Delivery Parameters", f"{deliv_eff*100:.1f}%")
        with col3:
            st.metric("Effective Editing Rate", f"{effective_eff*100:.1f}%")
        
        # Visualize delivery efficiency
        plot_delivery_efficiency(effective_eff)
        
        # Show sequence changes
        st.markdown("### Sequence Changes")
        plot_sequence_alignment(st.session_state.seq_record, st.session_state.edited_seq)
        
        # Show editing results
        if st.session_state.targets:
            st.markdown("### Editing Results")
            plot_editing_results(
                st.session_state.seq_record,
                st.session_state.edited_seq,
                st.session_state.targets
            )
        
        # Download button
        st.download_button(
            label="Download Edited Sequence",
            data=st.session_state.edited_seq,
            file_name="edited_sequence.fasta",
            mime="text/plain"
        )
else:
    st.info("Please input a DNA sequence to begin.")

# --- Footer ---
st.markdown("---")
st.markdown(
    """
    **How it works:**  
    - Finds all PAM sites (both strands) based on selected Cas variant  
    - Locates gRNA binding regions (configurable mismatches)  
    - Simulates NHEJ (knockout) or HDR (knock-in) editing  
    - Models delivery efficiency by organ/method  
    - Reports off-targets and shows edited sequence  
    
    *This is a simulation tool for educational and research purposes only.*
    """
)

        
       
        
       
