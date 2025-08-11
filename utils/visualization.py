"""
Visualization utilities for CRISPR-Cas9 simulation results.
"""
from typing import Dict, Any, List, Optional
import matplotlib.pyplot as plt
import streamlit as st

def plot_editing_results(seq_record: str, edited_seq: str, targets: List[Dict[str, Any]]) -> None:
    """
    Plot the editing results showing the original and edited sequences.
    
    Args:
        seq_record: Original DNA sequence
        edited_seq: Edited DNA sequence
        targets: List of target sites
    """
    try:
        from dna_features_viewer import GraphicFeature, GraphicRecord
        
        # Create features for visualization
        features = []
        for i, t in enumerate(targets):
            # Color on-targets red, off-targets blue
            color = "#ff6b6b" if t['mismatches'] == 0 else "#4dabf7"
            # Our targets come from Cas9.find_targets() with keys:
            # 'position' (int), 'sequence' (str), 'pam_sequence' (str), 'strand' ('+'|'-'), 'mismatches' (int)
            start = t.get('position', 0)
            end = start + len(t.get('sequence', ''))

            features.append(GraphicFeature(
                start=start,
                end=end,
                strand=1 if t['strand'] == "+" else -1,
                color=color,
                label=f"{'On' if t['mismatches'] == 0 else 'Off'}-target {i+1}"
            ))
        
        # Create the plot
        record = GraphicRecord(sequence_length=len(seq_record), features=features)
        ax, _ = record.plot(figure_width=10)
        
        # Customize the plot
        ax.set_title("CRISPR-Cas9 Editing Results")
        ax.set_xlabel("Position (bp)")
        
        # Show the plot in Streamlit
        st.pyplot(ax.figure)
        
    except ImportError:
        st.warning("Install dna_features_viewer for graphical visualization.")
        
        # Fallback to text-based visualization
        st.subheader("Editing Results")
        st.write("Original sequence length:", len(seq_record))
        st.write("Edited sequence length:", len(edited_seq))
        
        # Show differences
        import difflib
        diff = difflib.ndiff(seq_record, edited_seq)
        diff_text = ''.join(diff)
        st.text_area("Sequence changes (diff)", diff_text, height=200)

def plot_delivery_efficiency(effective_eff: float) -> None:
    """
    Plot the delivery efficiency as a bar chart.
    
    Args:
        effective_eff: Effective delivery efficiency (0-1)
    """
    fig, ax = plt.subplots(figsize=(6, 2))
    
    # Create a simple bar chart
    ax.bar(["Success", "Failure"], 
           [effective_eff, 1-effective_eff], 
           color=["#4caf50", "#f44336"])
    
    ax.set_ylim(0, 1)
    ax.set_ylabel("Probability")
    ax.set_title("Delivery Success Probability")
    
    # Add value labels
    for i, v in enumerate([effective_eff, 1-effective_eff]):
        ax.text(i, v + 0.02, f"{v:.1%}", ha='center')
    
    # Show the plot in Streamlit
    st.pyplot(fig)

def plot_sequence_alignment(original: str, edited: str) -> None:
    """
    Generate a side-by-side comparison of original and edited sequences.
    
    Args:
        original: Original DNA sequence
        edited: Edited DNA sequence
    """
    st.subheader("Sequence Comparison")
    
    # Truncate long sequences for display
    max_display = 500
    if len(original) > max_display or len(edited) > max_display:
        st.warning(f"Sequences truncated to first {max_display} bases for display.")
        original = original[:max_display]
        edited = edited[:max_display]
    
    # Create a two-column layout
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Original Sequence**")
        st.code(original, language="text")
    
    with col2:
        st.markdown("**Edited Sequence**")
        st.code(edited, language="text")
    
    # Highlight differences
    st.markdown("**Changes (red=deleted, green=added, yellow=changed):**")
    
    import difflib
    diff = difflib.ndiff(original, edited)
    diff_text = ''.join(diff)
    
    # Convert diff to HTML with colors
    html_diff = []
    for line in diff_text.split('\n'):
        line_html = []
        for char in line:
            if char.startswith('-'):
                line_html.append(f'<span style="color: red;">{char[-1]}</span>')
            elif char.startswith('+'):
                line_html.append(f'<span style="color: green;">{char[-1]}</span>')
            elif char.startswith('?'):
                line_html.append(f'<span style="background-color: yellow;">^</span>')
            else:
                line_html.append(char)
        html_diff.append(''.join(line_html))
    
    st.markdown(''.join(html_diff), unsafe_allow_html=True)
