

# CRISPR-Cas9 Gene Editing Simulator

> **Note**: This is a refactored version of the original project by [Nirajan2005](https://github.com/Nirajan2005). This version by 06navdeep06 includes alot of code improvements and bug fixes while maintaining the original functionality.

A comprehensive simulation tool for CRISPR-Cas9 gene editing experiments, now with a modular and maintainable codebase.

## Features

- **Interactive Web Interface**: Built with Streamlit for an intuitive user experience
- **Multiple Cas Variants**: Supports SpCas9, SaCas9, Cpf1, and custom PAM sequences
- **Realistic Simulation**: Models NHEJ and HDR repair pathways
- **Delivery Modeling**: Simulates different delivery methods and their efficiencies
- **Chromatin Context**: Accounts for chromatin accessibility and histone modifications
- **Visualization**: Interactive visualization of editing results

## Project Structure

```
CRISPRcas9_simV3/
├── app.py                  # Main application entry point
├── config.py               # Configuration and constants
├── models/                 # Data models and core logic
│   ├── __init__.py
│   ├── cas9.py            # Cas9 class and target finding
│   ├── scoring.py         # On/off-target scoring
│   └── simulation.py      # Delivery and editing simulation
├── utils/                  # Utility functions
│   ├── __init__.py
│   ├── sequence_utils.py  # Sequence manipulation
│   └── visualization.py   # Plotting and visualization
└── requirements.txt        # Dependencies
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/CRISPRcas9_simV3.git
   cd CRISPRcas9_simV3
   ```

2. Create and activate a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Run the Streamlit application:
   ```bash
   streamlit run app.py
   ```

2. Open your web browser and navigate to the URL shown in the terminal (usually http://localhost:8501)

3. Use the sidebar to:
   - Input or upload a DNA sequence
   - Configure guide RNA and PAM settings
   - Set experimental parameters
   - Run the simulation

## Dependencies

- Python 3.8+
- Streamlit
- Biopython
- dna-features-viewer
- Matplotlib
- NumPy

## License

MIT License - see [LICENSE](LICENSE) for details.
