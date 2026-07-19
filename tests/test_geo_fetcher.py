"""Tests for src/geo_fetcher.py."""

from __future__ import annotations

import pytest

from src.config import PipelineConfig
from src.geo_fetcher import (
    GEOFetchError,
    fetch_gse,
    get_expression_table,
    get_platform_table,
    get_sample_characteristics,
)


class TestGEOFetcher:
    def test_fetch_real_gse1563(self, gse1563):
        assert gse1563 is not None
        assert len(gse1563.gsms) > 0

    def test_fetch_real_gse68918(self, gse68918):
        assert gse68918 is not None
        assert len(gse68918.gsms) > 0

    def test_fetch_real_gse45643(self, gse45643):
        assert gse45643 is not None
        assert len(gse45643.gsms) > 0

    def test_invalid_gse_raises(self, config):
        with pytest.raises(ValueError, match="Invalid GSE accession"):
            fetch_gse("INVALID", config)

    def test_expression_table_not_empty(self, gse1563):
        samples = get_expression_table(gse1563)
        assert len(samples) > 0

    def test_platform_mapping_exists(self, gse1563):
        mapping = get_platform_table(gse1563)
        assert mapping is not None
        assert len(mapping) > 0

    def test_sample_characteristics_exist(self, gse1563):
        chars = get_sample_characteristics(gse1563)
        assert len(chars) > 0

    def test_fetch_handles_nonexistent_gse(self, config):
        with pytest.raises(GEOFetchError):
            fetch_gse("GSE99999999", config)
