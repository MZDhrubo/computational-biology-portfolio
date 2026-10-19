import pytest

from portfolio_tools.docking_descriptor_triage import triage_table


def test_triage_sorts_and_flags_synthetic_compounds():
    table = triage_table("examples/synthetic_docking_descriptors.csv")
    assert table.iloc[0]["compound_id"] == "SYN_D"
    assert table.iloc[0]["rank_by_input_docking_score"] == 1
    assert bool(table.loc[table.compound_id == "SYN_A", "passes_illustrative_descriptor_rules"].iloc[0])
    assert not bool(table.loc[table.compound_id == "SYN_B", "passes_illustrative_descriptor_rules"].iloc[0])


def test_missing_column_is_rejected(tmp_path):
    file = tmp_path / "bad.csv"
    file.write_text("compound_id,docking_score\nSYN_A,-7.0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Missing required columns"):
        triage_table(file)
