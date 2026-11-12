"""Tests for disease search module."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from src.config import PipelineConfig
from src.disease_search import (
    DRG_SEARCH_QUERIES,
    GEOSearchResult,
    EntrezError,
    _parse_gds_summary,
    search_drg_studies,
    search_gse_by_disease,
)


@pytest.fixture
def config():
    return PipelineConfig(email="test@example.com")


# ---------------------------------------------------------------------------
# DRG_SEARCH_QUERIES constant
# ---------------------------------------------------------------------------

def test_drg_search_queries_has_three_entries():
    assert len(DRG_SEARCH_QUERIES) == 3


def test_drg_search_queries_contains_canonical_terms():
    combined = " ".join(DRG_SEARCH_QUERIES).lower()
    assert "dorsal root ganglion" in combined
    assert "neuropathic" in combined
    assert "injured" in combined


# ---------------------------------------------------------------------------
# _parse_gds_summary
# ---------------------------------------------------------------------------

def test_parse_gds_summary_valid():
    summary = [{"GSE": "12345", "title": "Test", "summary": "A study", "taxon": "Homo sapiens", "gdsType": "Expression profiling by high throughput sequencing", "n_samples": 10}]
    result = _parse_gds_summary(summary)
    assert result is not None
    assert result.gse_id == "GSE12345"
    assert result.n_samples == 10


def test_parse_gds_summary_empty():
    assert _parse_gds_summary([]) is None


def test_parse_gds_summary_adds_gse_prefix():
    summary = [{"GSE": "99", "title": "", "summary": "", "taxon": "", "gdsType": "", "n_samples": 0}]
    result = _parse_gds_summary(summary)
    assert result.gse_id == "GSE99"


# ---------------------------------------------------------------------------
# search_drg_studies
# ---------------------------------------------------------------------------

def _make_mock_result(gse_id, taxon="Homo sapiens"):
    return GEOSearchResult(
        gse_id=gse_id,
        title=f"DRG study {gse_id}",
        summary="Human DRG single nucleus neuropathic pain",
        taxon=taxon,
        gds_type="Expression profiling by high throughput sequencing",
        n_samples=50,
    )


def test_search_drg_studies_deduplicates(config):
    """Results appearing in multiple queries should only appear once."""
    shared = _make_mock_result("GSE11111")
    unique_q2 = _make_mock_result("GSE22222")

    call_count = {"n": 0}

    def mock_search(disease, cfg, max_results=None):
        call_count["n"] += 1
        if call_count["n"] == 1:
            return [shared]
        elif call_count["n"] == 2:
            return [shared, unique_q2]  # shared appears again
        else:
            return []

    with patch("src.disease_search.search_gse_by_disease", side_effect=mock_search):
        results = search_drg_studies(config)

    gse_ids = [r.gse_id for r in results]
    assert gse_ids.count("GSE11111") == 1, "Duplicate GSE should be deduplicated"
    assert "GSE22222" in gse_ids


def test_search_drg_studies_warns_if_fewer_than_min(config, caplog):
    """Warning is logged when fewer than min_human_studies are found."""
    import logging

    one_result = [_make_mock_result("GSE11111")]

    with patch("src.disease_search.search_gse_by_disease", return_value=one_result):
        with caplog.at_level(logging.WARNING, logger="src.disease_search"):
            results = search_drg_studies(config)

    warning_texts = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert any("human DRG studies" in w for w in warning_texts)


def test_search_drg_studies_uses_config_queries(config):
    """Uses queries from config.drg.search_queries."""
    config.drg.search_queries = ["test query one"]
    captured = []

    def mock_search(disease, cfg, max_results=None):
        captured.append(disease)
        return []

    with patch("src.disease_search.search_gse_by_disease", side_effect=mock_search):
        search_drg_studies(config)

    assert captured == ["test query one"]


def test_search_drg_studies_handles_entrez_error(config, caplog):
    """EntrezError on one query does not crash the whole search."""
    import logging

    call_count = {"n": 0}

    def mock_search(disease, cfg, max_results=None):
        call_count["n"] += 1
        if call_count["n"] == 1:
            raise EntrezError("Network error")
        return [_make_mock_result("GSE33333")]

    with patch("src.disease_search.search_gse_by_disease", side_effect=mock_search):
        with caplog.at_level(logging.WARNING, logger="src.disease_search"):
            results = search_drg_studies(config)

    # Should still return results from successful queries
    assert any(r.gse_id == "GSE33333" for r in results)


def test_search_drg_studies_empty_queries_returns_empty(config, caplog):
    """Empty query list returns [] without crashing."""
    import logging

    config.drg.search_queries = []
    with caplog.at_level(logging.WARNING, logger="src.disease_search"):
        results = search_drg_studies(config, queries=[])
    assert results == []


def test_search_gse_by_disease_raises_on_empty(config):
    with pytest.raises(ValueError, match="non-empty"):
        search_gse_by_disease("", config)
