import pytest

from portfolio_tools.pdb_inventory import inspect_pdb


def test_toy_pdb_inventory():
    report = inspect_pdb("examples/toy_protein.pdb")
    assert report["chain_count"] == 2
    assert report["total_atom_records"] == 9
    assert report["chains"][0]["chain_id"] == "A"
    assert report["chains"][0]["sequence_one_letter"] == "AG"
    assert report["chains"][0]["ca_coverage_fraction"] == 1.0
    assert "not a complete" in report["scope_note"]


def test_no_atom_records_is_rejected(tmp_path):
    file = tmp_path / "empty.pdb"
    file.write_text("HEADER    TEST\nEND\n", encoding="utf-8")
    with pytest.raises(ValueError, match="No ATOM records"):
        inspect_pdb(file)
