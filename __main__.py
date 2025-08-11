"""
Package entry point for the CRISPR-Cas9 Simulator.

Run with:
    python -m CRISPRcas9_simV3

This will launch the Streamlit app.
"""

import os
import sys
import subprocess


def _run_streamlit_app() -> int:
    project_root = os.path.dirname(os.path.abspath(__file__))
    app_path = os.path.join(project_root, "app.py")
    cmd = [sys.executable, "-m", "streamlit", "run", app_path]
    # Inherit the current environment and stdio
    return subprocess.call(cmd)


if __name__ == "__main__":
    sys.exit(_run_streamlit_app())
