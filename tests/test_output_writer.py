"""Tests for src/output_writer.py."""

from __future__ import annotations

import pandas as pd
import pytest

from src.output_writer import OutputError, write_output


class TestOutputWriter:
    def test_parquet_roundtrip(self, tmp_path, gse1563_expr):
        de = pd.DataFrame({"logFC": [0.5], "pval": [0.01], "padj": [0.05]})
        files = write_output(gse1563_expr, de, {"gse_id": "GSE1563"}, str(tmp_path), "GSE1563", "parquet")
        loaded = pd.read_parquet(files["expression"])
        assert loaded.shape == gse1563_expr.shape

    def test_csv_roundtrip(self, tmp_path, gse1563_expr):
        de = pd.DataFrame({"logFC": [0.5], "pval": [0.01], "padj": [0.05]})
        files = write_output(gse1563_expr, de, {"gse_id": "GSE1563"}, str(tmp_path), "GSE1563", "csv")
        loaded = pd.read_csv(files["expression"], index_col=0)
        assert loaded.shape == gse1563_expr.shape

    def test_invalid_format_raises(self, tmp_path):
        with pytest.raises(OutputError, match="Invalid output format"):
            write_output(pd.DataFrame(), pd.DataFrame(), {}, str(tmp_path), "GSE0", "feather")

    def test_metadata_written(self, tmp_path):
        files = write_output(pd.DataFrame(), pd.DataFrame(), {"gse_id": "X"}, str(tmp_path), "GSEX")
        assert "metadata" in files
        assert files["metadata"].exists()
