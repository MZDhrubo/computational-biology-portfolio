"""Summarize a protein PDB file without claiming full structural validation."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def inspect_pdb(path: str | Path) -> dict:
    """Return basic chain and residue inventory for ATOM records in a PDB file."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    chains: dict[str, dict] = defaultdict(lambda: {"residues": {}, "atoms": 0, "bfactors": [], "occupancies": []})
    total_atoms = 0
    for line_number, line in enumerate(source.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
        if not line.startswith("ATOM  "):
            continue
        if len(line) < 54:
            raise ValueError(f"Malformed ATOM record at line {line_number}: expected fixed-width PDB fields.")
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip().upper()
        chain_id = line[21].strip() or "_"
        residue_number = line[22:26].strip()
        insertion_code = line[26].strip()
        altloc = line[16].strip()
        try:
            occupancy = float(line[54:60].strip() or 0)
            bfactor = float(line[60:66].strip() or 0)
        except ValueError as exc:
            raise ValueError(f"Invalid occupancy or B-factor at line {line_number}.") from exc
        residue_key = (residue_number, insertion_code, residue_name)
        chain = chains[chain_id]
        residue = chain["residues"].setdefault(residue_key, {"atoms": set(), "has_ca": False, "resname": residue_name})
        # Atom name set avoids double-counting alternate conformations for coverage.
        residue["atoms"].add(atom_name)
        if atom_name == "CA" and altloc in ("", "A", "1"):
            residue["has_ca"] = True
        chain["atoms"] += 1
        chain["bfactors"].append(bfactor)
        chain["occupancies"].append(occupancy)
        total_atoms += 1
    if total_atoms == 0:
        raise ValueError("No ATOM records found in the PDB file.")

    chain_summaries = []
    for chain_id in sorted(chains):
        chain = chains[chain_id]
        residues = list(chain["residues"].values())
        sequence = "".join(THREE_TO_ONE.get(res["resname"], "X") for res in residues)
        ca_count = sum(res["has_ca"] for res in residues)
        chain_summaries.append({
            "chain_id": chain_id,
            "residue_count": len(residues),
            "sequence_one_letter": sequence,
            "ca_atom_residues": ca_count,
            "ca_coverage_fraction": ca_count / len(residues) if residues else 0.0,
            "atom_records": chain["atoms"],
            "mean_b_factor": sum(chain["bfactors"]) / len(chain["bfactors"]) if chain["bfactors"] else None,
            "mean_occupancy": sum(chain["occupancies"]) / len(chain["occupancies"]) if chain["occupancies"] else None,
        })
    return {
        "input_file": source.name,
        "total_atom_records": total_atoms,
        "chain_count": len(chain_summaries),
        "chains": chain_summaries,
        "scope_note": "Inventory only; not a complete structure-quality validation.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Protein PDB file")
    parser.add_argument("--out", required=True, help="Output JSON report")
    args = parser.parse_args()
    report = inspect_pdb(args.input)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote structure inventory to {output}")
    print(report["scope_note"])


if __name__ == "__main__":
    main()
