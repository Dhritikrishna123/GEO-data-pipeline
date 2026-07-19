"""Tests for src/disease_search.py."""

from __future__ import annotations

import pytest

from src.disease_search import search_gse_by_disease


class TestDiseaseSearch:
    def test_search_returns_results(self, config):
        results = search_gse_by_disease("breast cancer", config, max_results=5)
        assert len(results) > 0

    def test_search_gse_ids_valid(self, config):
        results = search_gse_by_disease("breast cancer", config, max_results=5)
        for r in results:
            assert r.gse_id.startswith("GSE")

    def test_search_with_limit(self, config):
        results = search_gse_by_disease("diabetes", config, max_results=3)
        assert len(results) <= 3

    def test_search_empty_raises(self, config):
        with pytest.raises(ValueError, match="non-empty"):
            search_gse_by_disease("", config)

    def test_search_whitespace_raises(self, config):
        with pytest.raises(ValueError, match="non-empty"):
            search_gse_by_disease("   ", config)

    def test_search_has_metadata(self, config):
        results = search_gse_by_disease("leukemia", config, max_results=2)
        for r in results:
            assert r.title != ""
            assert r.gse_id != ""
