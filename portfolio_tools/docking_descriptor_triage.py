"""Transparent sorting and heuristic descriptor flags for docking result tables.

This utility does not perform docking or predict ADMET. Thresholds are
illustrative and must be reviewed for the target, compound class, and study.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED = ["compound_id", "docking_score", "molecular_weight", "logp", "hbd", "hba", "tpsa"]


def triage_table(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    missing = [column for column in REQUIRED if column not in table.columns]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    if table.empty:
        raise ValueError("Input table has no compounds.")
    if table["compound_id"].isna().any() or (table["compound_id"].astype(str).str.strip() == "").any():
        raise ValueError("compound_id values must be present and non-empty.")
    if table["compound_id"].duplicated().any():
        raise ValueError("compound_id values must be unique.")
    numeric_cols = REQUIRED[1:]
    for column in numeric_cols:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    if table[numeric_cols].isna().any().any():
        raise ValueError("Docking score and descriptor columns must contain numeric values.")
    if not np.isfinite(table[numeric_cols].to_numpy(dtype=float)).all():
        raise ValueError("Numeric values must be finite.")
    if (table[["molecular_weight", "hbd", "hba", "tpsa"]] < 0).any().any():
        raise ValueError("Molecular weight, H-bond counts, and TPSA cannot be negative.")

    # These flags are transparent heuristics, not ADMET predictions.
    table["mw_rule"] = table["molecular_weight"] <= 500
    table["logp_rule"] = table["logp"] <= 5
    table["hbd_rule"] = table["hbd"] <= 5
    table["hba_rule"] = table["hba"] <= 10
    table["tpsa_rule"] = table["tpsa"] <= 140
    rule_cols = ["mw_rule", "logp_rule", "hbd_rule", "hba_rule", "tpsa_rule"]
    table["passes_illustrative_descriptor_rules"] = table[rule_cols].all(axis=1)
    table = table.sort_values(["docking_score", "compound_id"], ascending=[True, True]).reset_index(drop=True)
    table.insert(0, "rank_by_input_docking_score", range(1, len(table) + 1))
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="CSV with docking scores and descriptors")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()
    result = triage_table(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / "docking_descriptor_triage.csv"
    result.to_csv(output_path, index=False)
    print(f"Wrote {len(result)} rows to {output_path}")
    print("Reminder: docking rank and descriptor flags are not proof of binding or ADMET properties.")


if __name__ == "__main__":
    main()
