"""Tests for src/expression_processor.py."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.expression_processor import build_expression_matrix, normalize_expression
from src.geo_fetcher import get_expression_table, get_platform_table


class TestExpressionProcessor:
    def test_build_matrix_shape(self, gse1563_raw, gse1563_probe):
        expr = build_expression_matrix(gse1563_raw, gse1563_probe)
        assert expr.shape[0] > 0
        assert expr.shape[1] > 0

    def test_normalize_no_nans(self, gse1563_expr):
        assert not gse1563_expr.isna().any().any()

    def test_duplicate_genes_merged(self, gse1563_expr):
        assert not gse1563_expr.index.duplicated().any()

    def test_gene_symbols_valid(self, gse1563_expr):
        bad = [g for g in gse1563_expr.index if not isinstance(g, str) or len(g) < 1]
        assert len(bad) == 0

    def test_all_positive_values(self, gse1563_expr):
        assert (gse1563_expr.values >= 0).all()

    def test_different_platforms(self, gse45643):
        raw = get_expression_table(gse45643)
        probe = get_platform_table(gse45643)
        expr = build_expression_matrix(raw, probe)
        expr = normalize_expression(expr, probe)
        assert expr.shape[0] > 0

    def test_build_matrix_empty(self):
        expr = build_expression_matrix({})
        assert expr.empty

    def test_build_matrix_no_mapping(self, gse1563_raw):
        expr = build_expression_matrix(gse1563_raw, probe_to_gene=None)
        assert expr.shape[0] > 0

    def test_gene_count_gse68918(self, gse68918):
        raw = get_expression_table(gse68918)
        probe = get_platform_table(gse68918)
        expr = build_expression_matrix(raw, probe)
        assert len(expr) > 0
