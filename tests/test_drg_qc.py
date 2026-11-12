"""Unit tests for DRG quality control checker."""

from __future__ import annotations

import pytest

from src.drg_qc import DRGQualityChecker, QCResult, _count_groups, _detect_sequencing_type, _is_human


# ---------------------------------------------------------------------------
# _is_human
# ---------------------------------------------------------------------------

def test_is_human_homo_sapiens():
    assert _is_human("Homo sapiens") is True

def test_is_human_case_insensitive():
    assert _is_human("homo sapiens") is True

def test_is_human_human_keyword():
    assert _is_human("human") is True

def test_is_human_mouse():
    assert _is_human("Mus musculus") is False

def test_is_human_rat():
    assert _is_human("Rattus norvegicus") is False


# ---------------------------------------------------------------------------
# _count_groups
# ---------------------------------------------------------------------------

def test_count_groups_basic():
    chars = {
        "GSM1": ["condition: CCI", "tissue: DRG"],
        "GSM2": ["condition: sham", "tissue: DRG"],
        "GSM3": ["condition: injured", "tissue: DRG"],
        "GSM4": ["condition: naive", "tissue: DRG"],
    }
    n_ctrl, n_inj = _count_groups(chars)
    assert n_inj >= 2
    assert n_ctrl >= 2

def test_count_groups_no_labels():
    chars = {
        "GSM1": ["tissue: DRG"],
        "GSM2": ["tissue: DRG"],
    }
    n_ctrl, n_inj = _count_groups(chars)
    # Should not crash; values may both be 0
    assert n_ctrl + n_inj >= 0

def test_count_groups_all_injured():
    chars = {
        "GSM1": ["condition: CCI"],
        "GSM2": ["condition: SNL"],
    }
    n_ctrl, n_inj = _count_groups(chars)
    assert n_inj == 2
    assert n_ctrl == 0


# ---------------------------------------------------------------------------
# _detect_sequencing_type
# ---------------------------------------------------------------------------

def test_detect_snrna_seq():
    seq = _detect_sequencing_type("Human DRG single nucleus RNA-seq", "")
    assert seq == "snRNA-seq"

def test_detect_scrna_seq():
    seq = _detect_sequencing_type("DRG scRNA-seq neuropathic pain", "")
    assert seq == "scRNA-seq"

def test_detect_bulk():
    seq = _detect_sequencing_type("Bulk RNA-seq of DRG neurons", "")
    assert seq == "bulk RNA-seq"

def test_detect_microarray():
    seq = _detect_sequencing_type("Affymetrix microarray DRG", "")
    assert seq == "microarray"

def test_detect_unknown():
    seq = _detect_sequencing_type("DRG gene expression", "")
    assert seq == "unknown"


# ---------------------------------------------------------------------------
# DRGQualityChecker.check
# ---------------------------------------------------------------------------

@pytest.fixture
def default_checker():
    return DRGQualityChecker(require_human=True, require_labels=True, min_cells_per_group=None)


def test_qc_passes_clean_human_study(default_checker):
    chars = {
        "GSM1": ["condition: CCI injured"],
        "GSM2": ["condition: CCI injured"],
        "GSM3": ["condition: sham control"],
        "GSM4": ["condition: sham control"],
    }
    result = default_checker.check(
        gse_id="GSE99999",
        taxon="Homo sapiens",
        title="Human DRG single nucleus RNA-seq neuropathic pain",
        summary="",
        sample_characteristics=chars,
    )
    assert result.passed is True
    assert result.is_human is True
    assert result.has_labels is True
    assert result.n_injured >= 2
    assert result.n_control >= 2
    assert len(result.reasons) == 0


def test_qc_fails_non_human(default_checker):
    chars = {
        "GSM1": ["condition: CCI"],
        "GSM2": ["condition: sham"],
    }
    result = default_checker.check(
        gse_id="GSE99998",
        taxon="Mus musculus",
        title="Mouse DRG CCI model",
        summary="",
        sample_characteristics=chars,
    )
    assert result.passed is False
    assert result.is_human is False
    assert any("Non-human" in r for r in result.reasons)


def test_qc_fails_no_labels(default_checker):
    chars = {
        "GSM1": ["tissue: DRG"],
        "GSM2": ["tissue: DRG"],
    }
    result = default_checker.check(
        gse_id="GSE99997",
        taxon="Homo sapiens",
        title="Human DRG",
        summary="",
        sample_characteristics=chars,
    )
    assert result.passed is False
    assert result.has_labels is False


def test_qc_min_cells_threshold():
    checker = DRGQualityChecker(require_human=True, require_labels=True, min_cells_per_group=3)
    chars = {
        "GSM1": ["condition: injured"],
        "GSM2": ["condition: sham"],
    }
    result = checker.check(
        gse_id="GSE99996",
        taxon="Homo sapiens",
        title="Human DRG",
        summary="",
        sample_characteristics=chars,
    )
    assert result.passed is False
    assert any("samples" in r for r in result.reasons)


def test_qc_no_min_cells_no_threshold():
    checker = DRGQualityChecker(require_human=True, require_labels=True, min_cells_per_group=None)
    chars = {
        "GSM1": ["condition: injured"],
        "GSM2": ["condition: sham"],
    }
    result = checker.check(
        gse_id="GSE99995",
        taxon="Homo sapiens",
        title="Human DRG",
        summary="",
        sample_characteristics=chars,
    )
    # Should pass (no cell count threshold)
    assert result.passed is True


def test_qc_low_sample_count_generates_warning(default_checker):
    chars = {
        "GSM1": ["condition: CCI injured"],
        "GSM2": ["condition: sham"],
    }
    result = default_checker.check(
        gse_id="GSE99994",
        taxon="Homo sapiens",
        title="Human DRG",
        summary="",
        sample_characteristics=chars,
    )
    # Should have a warning about low total sample count
    assert any("low total sample count" in w.lower() for w in result.warnings)


def test_from_config():
    """DRGQualityChecker.from_config round-trips config settings."""
    from src.config import PipelineConfig, DRGConfig
    config = PipelineConfig(
        drg={"require_human": False, "require_labels": True, "min_cells_per_group": 50}
    )
    checker = DRGQualityChecker.from_config(config)
    assert checker.require_human is False
    assert checker.require_labels is True
    assert checker.min_cells_per_group == 50
