"""Exploratory quality-control summaries for a gene-by-sample count matrix.

This module is not a differential-expression analysis pipeline. It is designed
for transparent input validation and exploratory QC using synthetic examples.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_count_matrix(path: str | Path) -> pd.DataFrame:
    """Load a CSV count matrix with gene IDs in its first column."""
    counts = pd.read_csv(path, index_col=0)
    if counts.empty or counts.shape[0] < 2 or counts.shape[1] < 2:
        raise ValueError("Count matrix must contain at least two genes and two samples.")
    if counts.index.isna().any() or (counts.index.astype(str).str.strip() == "").any():
        raise ValueError("Gene identifiers must be present and non-empty.")
    if not counts.index.is_unique:
        raise ValueError("Gene identifiers must be unique.")
    if counts.columns.isna().any() or (counts.columns.astype(str).str.strip() == "").any():
        raise ValueError("Sample names must be present and non-empty.")
    if not counts.columns.is_unique:
        raise ValueError("Sample names must be unique.")
    numeric = counts.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        raise ValueError("All count values must be numeric and non-missing.")
    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Count values must be finite.")
    if (values < 0).any():
        raise ValueError("Count values cannot be negative.")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError("Raw count values must be whole numbers.")
    return numeric.astype("int64")


def sample_qc(counts: pd.DataFrame) -> pd.DataFrame:
    """Return basic library-size and gene-detection metrics per sample."""
    library_size = counts.sum(axis=0)
    if (library_size <= 0).any():
        bad = library_size[library_size <= 0].index.tolist()
        raise ValueError(f"Samples with zero library size are not supported: {bad}")
    return pd.DataFrame({
        "sample": counts.columns.astype(str),
        "library_size": library_size.to_numpy(dtype="int64"),
        "genes_detected": (counts > 0).sum(axis=0).to_numpy(dtype="int64"),
        "zero_fraction": (counts == 0).mean(axis=0).to_numpy(dtype=float),
        "median_raw_count": counts.median(axis=0).to_numpy(dtype=float),
    })


def log_cpm(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    """Calculate log2 counts per million for exploratory visualization only."""
    if pseudocount <= 0:
        raise ValueError("pseudocount must be greater than zero")
    library_size = counts.sum(axis=0)
    if (library_size <= 0).any():
        raise ValueError("All samples must have positive library sizes.")
    cpm = counts.div(library_size, axis=1) * 1_000_000
    return np.log2(cpm + pseudocount)


def pca_coordinates(counts: pd.DataFrame) -> pd.DataFrame:
    """Return two sample PCA coordinates using centered log-CPM and SVD."""
    transformed = log_cpm(counts).T.to_numpy(dtype=float)
    centered = transformed - transformed.mean(axis=0, keepdims=True)
    u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    n_components = min(2, u.shape[1], len(singular_values))
    scores = u[:, :n_components] * singular_values[:n_components]
    if n_components == 1:
        scores = np.column_stack([scores[:, 0], np.zeros(scores.shape[0])])
    return pd.DataFrame({
        "sample": counts.columns.astype(str),
        "PC1": scores[:, 0],
        "PC2": scores[:, 1],
    })


def run(counts_path: str | Path, outdir: str | Path) -> dict[str, str | int]:
    """Run the exploratory QC and write tables plus a small JSON summary."""
    counts = load_count_matrix(counts_path)
    qc = sample_qc(counts)
    pca = pca_coordinates(counts)
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    qc_path = output / "sample_qc.csv"
    pca_path = output / "pca_coordinates.csv"
    summary_path = output / "run_summary.json"
    qc.to_csv(qc_path, index=False)
    pca.to_csv(pca_path, index=False)
    summary = {
        "input_file": Path(counts_path).name,
        "genes": int(counts.shape[0]),
        "samples": int(counts.shape[1]),
        "analysis_scope": "exploratory QC only; not differential-expression analysis",
        "outputs": [qc_path.name, pca_path.name],
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return {"genes": counts.shape[0], "samples": counts.shape[1], "outdir": str(output)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True, help="Gene-by-sample count matrix CSV")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()
    result = run(args.counts, args.outdir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
