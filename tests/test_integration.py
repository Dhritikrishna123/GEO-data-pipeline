"""Integration tests: full pipeline on real GEO data."""

from __future__ import annotations

import pytest

from src.config import PipelineConfig
from src.de_analyzer import compute_de
from src.expression_processor import build_expression_matrix, normalize_expression
from src.geo_fetcher import fetch_gse, get_expression_table, get_platform_table, get_sample_characteristics
from src.metadata_extractor import extract_sample_groups
from src.output_writer import write_output
from src.pipeline import process_single_gse


class TestIntegration:
    def test_full_gse1563(self, gse1563, gse1563_raw, gse1563_probe, gse1563_expr):
        chars = get_sample_characteristics(gse1563)
        control, disease = extract_sample_groups(chars, disease_term="kidney transplant")
        assert len(control.samples) >= 2 and len(disease.samples) >= 2
        de = compute_de(gse1563_expr, control, disease)
        assert len(de) > 0

    def test_output_integration(self, tmp_path, gse1563_expr):
        de = gse1563_expr.head(0)
        files = write_output(gse1563_expr, de, {"gse_id": "GSE1563"}, str(tmp_path), "GSE1563")
        assert all(f.exists() for f in files.values())

    def test_pipeline_process_single(self, config):
        result = process_single_gse("GSE1563", config, disease_term="kidney transplant")
        assert result["gse_id"] == "GSE1563"
        assert result["n_genes"] > 0
        assert result["disease"] == "kidney transplant"
