"""Tests for src/de_analyzer.py."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.de_analyzer import compute_de
from src.metadata_extractor import SampleGroup, extract_sample_groups
from src.expression_processor import build_expression_matrix, normalize_expression
from src.geo_fetcher import get_expression_table, get_platform_table, get_sample_characteristics


class TestDEAnalyzer:
    def test_de_with_known_differential(self):
        np.random.seed(42)
        n_genes = 100
        n_a, n_b = 10, 10
        expr_a = np.random.randn(n_a, n_genes) + 5
        expr_b = np.random.randn(n_b, n_genes) + 6
        all_vals = np.vstack([expr_a, expr_b])
        all_vals = np.clip(all_vals, 0.001, None)
        all_vals = np.log2(all_vals)

        expr = pd.DataFrame(
            np.log2(np.vstack([all_vals[:n_a], all_vals[n_a:]]).T),
            index=[f"gene_{i}" for i in range(n_genes)],
            columns=[f"a_{i}" for i in range(n_a)] + [f"b_{i}" for i in range(n_b)],
        )
        group_a = SampleGroup(name="A", samples=[f"a_{i}" for i in range(n_a)])
        group_b = SampleGroup(name="B", samples=[f"b_{i}" for i in range(n_b)])

        de = compute_de(expr, group_a, group_b)
        assert len(de) == n_genes
        assert "logFC" in de.columns
        assert "pval" in de.columns
        assert "padj" in de.columns

    def test_de_all_zero_expression(self):
        expr = pd.DataFrame(
            np.zeros((5, 4)),
            index=[f"g{i}" for i in range(5)],
            columns=[f"s{i}" for i in range(4)],
        )
        group_a = SampleGroup(name="A", samples=["s0", "s1"])
        group_b = SampleGroup(name="B", samples=["s2", "s3"])
        de = compute_de(expr, group_a, group_b)
        assert len(de) == 0

    def test_de_insufficient_samples(self):
        expr = pd.DataFrame(
            np.random.randn(3, 2),
            index=[f"g{i}" for i in range(3)],
            columns=["s0", "s1"],
        )
        group_a = SampleGroup(name="A", samples=["s0"])
        group_b = SampleGroup(name="B", samples=["s1"])
        de = compute_de(expr, group_a, group_b)
        assert len(de) == 0

    def test_pvalues_between_0_and_1(self):
        np.random.seed(0)
        expr = pd.DataFrame(
            np.random.randn(50, 8) + 5,
            index=[f"g{i}" for i in range(50)],
            columns=[f"s{i}" for i in range(8)],
        )
        group_a = SampleGroup(name="A", samples=[f"s{i}" for i in range(4)])
        group_b = SampleGroup(name="B", samples=[f"s{i}" for i in range(4, 8)])
        de = compute_de(expr, group_a, group_b)
        assert (de["pval"] >= 0).all() and (de["pval"] <= 1).all()

    def test_logfc_finite(self):
        np.random.seed(1)
        expr = pd.DataFrame(
            np.random.randn(30, 6) + 5,
            index=[f"g{i}" for i in range(30)],
            columns=[f"s{i}" for i in range(6)],
        )
        group_a = SampleGroup(name="A", samples=["s0", "s1", "s2"])
        group_b = SampleGroup(name="B", samples=["s3", "s4", "s5"])
        de = compute_de(expr, group_a, group_b)
        assert np.isfinite(de["logFC"]).all()

    def test_de_with_real_gse(self, gse1563, gse1563_raw, gse1563_probe, gse1563_expr):
        chars = get_sample_characteristics(gse1563)
        control, disease = extract_sample_groups(chars)
        if len(control.samples) >= 2 and len(disease.samples) >= 2:
            de = compute_de(gse1563_expr, control, disease)
            assert len(de) > 0
            assert "logFC" in de.columns
            assert "padj" in de.columns
            assert (de["padj"] >= 0).all() and (de["padj"] <= 1).all()
