"""Tests for src/metadata_extractor.py."""

from __future__ import annotations

import pytest

from src.metadata_extractor import SampleGroup, extract_sample_groups


class TestMetadataExtractor:
    def test_classify_disease_vs_control(self):
        chars = {
            "sample1": ["tissue: normal"],
            "sample2": ["tissue: normal"],
            "sample3": ["tissue: breast cancer"],
            "sample4": ["tissue: breast cancer"],
        }
        control, disease = extract_sample_groups(chars)
        assert len(control.samples) == 2
        assert len(disease.samples) == 2

    def test_all_control_raises(self):
        chars = {
            "sample1": ["tissue: normal"],
            "sample2": ["tissue: normal"],
        }
        with pytest.raises(Exception):
            extract_sample_groups(chars)

    def test_insufficient_samples_raises(self):
        chars = {"sample1": ["tissue: normal"]}
        with pytest.raises(Exception):
            extract_sample_groups(chars)

    def test_empty_chars_raises(self):
        with pytest.raises(Exception):
            extract_sample_groups({})

    def test_ambiguous_grouping(self):
        chars = {
            "sample1": ["unknown characteristic"],
            "sample2": ["another unknown"],
            "sample3": ["more unknown"],
            "sample4": ["yet another"],
        }
        control, disease = extract_sample_groups(chars)
        assert len(control.samples) + len(disease.samples) == 4

    def test_mixed_characteristics(self):
        chars = {
            "s1": ["age: 50", "disease: breast cancer"],
            "s2": ["age: 45", "disease: breast cancer"],
            "s3": ["age: 52", "tissue: normal"],
            "s4": ["age: 48", "tissue: normal"],
        }
        control, disease = extract_sample_groups(chars)
        assert len(control.samples) >= 1
        assert len(disease.samples) >= 1

    def test_disease_term_boost(self):
        chars = {
            "s1": ["treatment: vehicle"],
            "s2": ["treatment: vehicle"],
            "s3": ["treatment: drug_a"],
            "s4": ["treatment: drug_a"],
        }
        control, disease = extract_sample_groups(chars, disease_term="drug_a")
        assert len(control.samples) >= 1
        assert len(disease.samples) >= 1

    def test_no_common_keywords(self):
        chars = {
            "s1": ["color: red"],
            "s2": ["color: blue"],
            "s3": ["color: green"],
            "s4": ["color: yellow"],
        }
        control, disease = extract_sample_groups(chars)
        assert len(control.samples) + len(disease.samples) == 4

    def test_large_dataset(self):
        chars = {}
        for i in range(50):
            chars[f"normal_{i}"] = ["condition: normal"]
        for i in range(50):
            chars[f"disease_{i}"] = ["condition: disease"]
        control, disease = extract_sample_groups(chars)
        assert len(control.samples) == 50
        assert len(disease.samples) == 50

    def test_three_groups(self):
        chars = {
            "s1": ["condition: normal"],
            "s2": ["condition: normal"],
            "s3": ["condition: disease"],
            "s4": ["condition: disease"],
            "s5": ["condition: treated"],
            "s6": ["condition: treated"],
        }
        control, disease = extract_sample_groups(chars)
        assert len(control.samples) + len(disease.samples) == 6
