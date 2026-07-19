"""Shared test fixtures for GEO pipeline tests."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

# Ensure project root is importable
_root = str(Path(__file__).resolve().parent.parent)
if _root not in sys.path:
    sys.path.insert(0, _root)

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=".*PyGEO.*")

from src.config import PipelineConfig
from src.geo_fetcher import fetch_gse, get_expression_table, get_platform_table, get_sample_characteristics
from src.expression_processor import build_expression_matrix, normalize_expression
from src.metadata_extractor import extract_sample_groups


@pytest.fixture(scope="session")
def config():
    """Pipeline configuration loaded from config.yaml."""
    return PipelineConfig.from_yaml("config.yaml")


@pytest.fixture(scope="session")
def gse1563(config):
    """Parsed GSE1563 dataset (session-scoped, fetched once)."""
    return fetch_gse("GSE1563", config)


@pytest.fixture(scope="session")
def gse68918(config):
    """Parsed GSE68918 dataset (session-scoped, fetched once)."""
    return fetch_gse("GSE68918", config)


@pytest.fixture(scope="session")
def gse45643(config):
    """Parsed GSE45643 dataset (session-scoped, fetched once)."""
    return fetch_gse("GSE45643", config)


@pytest.fixture(scope="session")
def gse1563_raw(gse1563):
    """Raw expression samples for GSE1563."""
    return get_expression_table(gse1563)


@pytest.fixture(scope="session")
def gse1563_probe(gse1563):
    """Probe-to-gene mapping for GSE1563."""
    return get_platform_table(gse1563)


@pytest.fixture(scope="session")
def gse1563_expr(gse1563_raw, gse1563_probe):
    """Normalized expression matrix for GSE1563."""
    expr = build_expression_matrix(gse1563_raw, gse1563_probe)
    return normalize_expression(expr, gse1563_probe)


@pytest.fixture(scope="session")
def gse1563_groups(gse1563):
    """Control/disease sample groups for GSE1563."""
    chars = get_sample_characteristics(gse1563)
    return extract_sample_groups(chars)
