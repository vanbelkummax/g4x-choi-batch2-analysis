#!/usr/bin/env python3
"""
Centralized configuration for Spatial Hackathon 2026 analysis.
Single source of truth for paths, metadata, and settings.

Import this module in all analysis scripts to ensure consistency.
"""
from pathlib import Path

# =============================================================================
# Paths
# =============================================================================

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "outputs" / "adata" / "annotated"
FIGURE_DIR = PROJECT_ROOT / "outputs" / "figures" / "showcase"
TABLE_DIR = PROJECT_ROOT / "outputs" / "tables"
POLYMATHIC_DIR = PROJECT_ROOT / "outputs" / "adata" / "polymathic"

# Ensure directories exist
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
POLYMATHIC_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# Statistical Settings
# =============================================================================

P_VALUE_CUTOFF = 0.05
FDR_METHOD = 'fdr_bh'  # Benjamini-Hochberg
RANDOM_SEED = 42

# Dual statistical testing settings (added 2026-01-20)
DUAL_STATS = True           # Report both MWU and Welch's t-test
SHOW_FDR = True             # Show FDR (q-values) alongside raw p-values
DEEP_DIVE_THRESHOLD = 0.10  # Trigger deep-dive analysis for p < 0.1
BOOTSTRAP_N = 1000          # Number of bootstrap iterations for CI

# Legacy setting (kept for backward compatibility)
STATISTICAL_TEST = 'welch'  # 'welch' or 'mwu' - Welch has better power for small n

# =============================================================================
# Sample Metadata (SINGLE SOURCE OF TRUTH)
# =============================================================================

# YP04A excluded: 134 cells after QC - too sparse for reliable analysis
PDAC_METADATA = {
    "YP03A": {"patient": "YP03", "timepoint": "Pre", "response": "NR"},
    "YP03C": {"patient": "YP03", "timepoint": "Post", "response": "NR"},
    # "YP04A" excluded - 134 cells after QC (too sparse)
    "YP04C": {"patient": "YP04", "timepoint": "Post", "response": "NR"},
    "YP12A": {"patient": "YP12", "timepoint": "Pre", "response": "R"},
    "YP12C": {"patient": "YP12", "timepoint": "Post", "response": "R"},
    "YP15A": {"patient": "YP15", "timepoint": "Pre", "response": "R"},
    "YP15C": {"patient": "YP15", "timepoint": "Post", "response": "R"},
}

# Excluded samples (documented for reference)
EXCLUDED_SAMPLES = {
    "YP04A": "134 cells after QC - too sparse for reliable analysis"
}

# Sample groupings for comparative analysis
SAMPLE_GROUPS = {
    'R_Pre': ['YP12A', 'YP15A'],
    'R_Post': ['YP12C', 'YP15C'],
    'NR_Pre': ['YP03A'],  # YP04A excluded
    'NR_Post': ['YP03C', 'YP04C'],
}

# =============================================================================
# Visualization
# =============================================================================

RESPONSE_COLORS = {'R': '#2ecc71', 'NR': '#e74c3c'}
TIMEPOINT_COLORS = {'Pre': '#3498db', 'Post': '#9b59b6'}

# Figure settings
FIGURE_DPI = 300
FIGURE_FORMAT = ['png', 'pdf']

# =============================================================================
# G4X Gating Thresholds
# =============================================================================

# Use percentile-based thresholds instead of hardcoded values
G4X_GATING_PERCENTILE = 90  # Top 10% are considered positive

# =============================================================================
# Statistical Power Notes
# =============================================================================

POWER_NOTES = """
STATISTICAL POWER LIMITATIONS:
- Sample size: n=4 Responders, n=3 Non-Responders (7 total)
- With Welch's t-test, p-value floor is ~0.03 (not 0.057 like MWU)
- Power for detecting large effects (d=1.5): ~50%
- Power for detecting moderate effects (d=0.8): ~20%
- All findings with p < 0.1 are HYPOTHESIS-GENERATING
- Validation in larger cohorts is required
"""

# =============================================================================
# Helper Functions
# =============================================================================

def get_response(sample_name: str) -> str:
    """Get response status for a sample."""
    return PDAC_METADATA.get(sample_name, {}).get('response', 'Unknown')

def get_timepoint(sample_name: str) -> str:
    """Get timepoint for a sample."""
    return PDAC_METADATA.get(sample_name, {}).get('timepoint', 'Unknown')

def is_excluded(sample_name: str) -> bool:
    """Check if a sample is excluded from analysis."""
    return sample_name in EXCLUDED_SAMPLES

def get_samples_by_response(response: str) -> list:
    """Get all samples with a given response status."""
    return [s for s, m in PDAC_METADATA.items() if m['response'] == response]

def get_samples_by_timepoint(timepoint: str) -> list:
    """Get all samples with a given timepoint."""
    return [s for s, m in PDAC_METADATA.items() if m['timepoint'] == timepoint]
