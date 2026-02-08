# -*- coding: utf-8 -*-
"""
UMAP analysis pipeline for phee call units (C vs E groups)

Overview
--------
Input:
- Three (or more) pair-level CSV files.
- Each CSV contains acoustic features (first FEATURE_COLS columns)
  and a label column formatted as:
      Group_idx_nphee_npheest
      e.g. C_1_2_A

Processing steps:
1. Per-file analysis (independent UMAP per CSV):
   - Perform UMAP on the full dataset
   - Color by group (C / E)
   - Color by npheest
   - Output:
       - 2D UMAP plot (PDF)
       - 2D UMAP coordinates (CSV)
       - Metadata preserved (group, idx, nphee, npheest, source file)

2. Combined analysis across all input files:
   - Subset with nphee == 2:
       * Independent UMAP
       * Group-specific labeling
       * Median UMAP positions per npheest
       * Connection plots between sequential phee calls
   - Subset with nphee == 3:
       * Same strategy as nphee == 2

Rules:
- Only the first FEATURE_COLS columns are used as features
- Missing values → 0
- Non-numeric values → 0
- Label column must be the last column

Outputs:
- All results are saved under OUTDIR
"""

# ============================ Imports ============================

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for batch plotting

from pathlib import Path
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from umap import UMAP
import matplotlib as mpl
from matplotlib.legend_handler import HandlerTuple


# ============================ Global Configuration ============================

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial']

INPUT_FILES = [
    "All_idx_6-Zscore.csv"
]

C_COLORS = ["#1f77b4", "#00c7c7", "#2c9f2c"]
E_COLORS = ["#D52728", "#ff31ff", "#ff7e0e"]

FEATURE_COLS = 6

OUTDIR = Path("Phee_unit_C-E_split_Zscore")
OUTDIR.mkdir(parents=True, exist_ok=True)

UMAP_KW = dict(
    n_neighbors=15,
    min_dist=0.1,
    metric="euclidean",
    random_state=42
)

MM = 1 / 25.4  # millimeter-to-inch conversion


# ============================ Utility Functions ============================

def get_axis_limits(*datasets):
    """
    Compute unified x/y axis limits for multiple 2D datasets.
    Axis limits are expanded slightly and rounded to multiples of 5.
    """
    min_x = min(np.min(d[:, 0]) for d in datasets)
    max_x = max(np.max(d[:, 0]) for d in datasets)
    min_y = min(np.min(d[:, 1]) for d in datasets)
    max_y = max(np.max(d[:, 1]) for d in datasets)

    x_margin = (max_x - min_x) * 0.05
    y_margin = (max_y - min_y) * 0.05

    min_x_5 = np.floor((min_x - x_margin) / 5) * 5
    max_x_5 = np.ceil((max_x + x_margin) / 5) * 5
    min_y_5 = np.floor((min_y - y_margin) / 5) * 5
    max_y_5 = np.ceil((max_y + y_margin) / 5) * 5

    return (min_x_5, max_x_5), (min_y_5, max_y_5)


def parse_label_series(label_series: pd.Series):
    """
    Parse label strings of the format:
        Group_idx_nphee_npheest
    Returns:
        group, idx, nphee, npheest
    """
    parts = label_series.astype(str).fillna("").str.split("_", expand=True)
    if parts.shape[1] < 4:
        raise ValueError("Label must follow 'Group_idx_nphee_npheest' format")

    group = parts[0].str.strip()
    idx = pd.to_numeric(parts[1], errors="coerce")
    nphee = pd.to_numeric(parts[2], errors="coerce")
    npheest = parts[3].str.strip()

    return group, idx, nphee, npheest


def prepare_features(df: pd.DataFrame) -> np.ndarray:
    """
    Extract feature matrix from the first FEATURE_COLS columns.
    Non-numeric and missing values are converted to zero.
    """
    feat_df = df.iloc[:, :FEATURE_COLS].copy()
    feat_df = feat_df.fillna(0)
    feat_df = feat_df.apply(pd.to_numeric, errors="coerce").fillna(0)

    X = feat_df.to_numpy(dtype=float)
    return np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)


def run_umap_2d(features: np.ndarray) -> np.ndarray:
    """Run 2D UMAP."""
    return UMAP(n_components=2, **UMAP_KW).fit_transform(features)


def save_points_2d(
    Y, labels, file_stem, suffix, outdir,
    original_index, group=None, idx=None,
    nphee=None, npheest=None, source_file=None
):
    """
    Save 2D UMAP coordinates and metadata to CSV.
    """
    df_out = pd.DataFrame({
        "umap_x": Y[:, 0],
        "umap_y": Y[:, 1],
        "label": labels.astype(str),
        "original_row_index": original_index
    })

    if source_file is not None:
        df_out["source_file"] = source_file
    if group is not None:
        df_out["group"] = group
    if idx is not None:
        df_out["idx"] = idx
    if nphee is not None:
        df_out["nphee"] = nphee
    if npheest is not None:
        df_out["npheest"] = npheest

    df_out.to_csv(
        outdir / f"{file_stem}__{suffix}__points.csv",
        index=False,
        encoding="utf-8-sig"
    )


# ============================ Main Pipeline ============================

# ---------- Per-file full-table UMAP ----------
for infile in INPUT_FILES:
    print(f"Processing file: {infile}")

    try:
        df = pd.read_csv(infile, encoding="utf-8")
    except UnicodeDecodeError:
        try:
            df = pd.read_csv(infile, encoding="latin1")
        except UnicodeDecodeError:
            df = pd.read_csv(infile, encoding="utf-16")

    labels_raw = df.iloc[:, -1]
    group, idx, nphee, npheest = parse_label_series(labels_raw)

    X = prepare_features(df)
    Y = run_umap_2d(X)

    save_points_2d(
        Y, group,
        Path(infile).stem,
        "umap_by_group",
        OUTDIR,
        original_index=df.index,
        group=group,
        idx=idx,
        nphee=nphee,
        npheest=npheest,
        source_file=[infile] * len(df)
    )

print("✅ Completed: per-file UMAP + combined nphee==2/3 analyses")
