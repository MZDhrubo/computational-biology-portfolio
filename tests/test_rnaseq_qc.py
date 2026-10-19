import pandas as pd
import pytest

from portfolio_tools.rnaseq_qc import load_count_matrix, sample_qc, pca_coordinates


def test_load_and_qc_synthetic_counts():
    counts = load_count_matrix("examples/synthetic_gene_counts.csv")
    qc = sample_qc(counts)
    pca = pca_coordinates(counts)
    assert counts.shape == (8, 4)
    assert qc["sample"].tolist() == ["Sample_A", "Sample_B", "Sample_C", "Sample_D"]
    assert (qc["library_size"] > 0).all()
    assert pca.shape == (4, 3)
    assert pca[["PC1", "PC2"]].notna().all().all()


def test_negative_count_is_rejected(tmp_path):
    file = tmp_path / "bad.csv"
    file.write_text("gene,S1,S2\ng1,1,2\ng2,-1,3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="negative"):
        load_count_matrix(file)


def test_non_integer_count_is_rejected(tmp_path):
    file = tmp_path / "bad.csv"
    file.write_text("gene,S1,S2\ng1,1.5,2\ng2,3,4\n", encoding="utf-8")
    with pytest.raises(ValueError, match="whole numbers"):
        load_count_matrix(file)
