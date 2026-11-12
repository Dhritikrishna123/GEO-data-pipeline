"""Extract metadata from GEO sample characteristics."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

# Keywords for detecting disease vs control groups
DISEASE_KEYWORDS = [
    "disease", "tumor", "tumour", "cancer", "leukemia", "leukaemia",
    "lymphoma", "melanoma", "carcinoma", "sarcoma", "glioblastoma",
    "diabetic", "alzheimer", "parkinson", "infection", "infected",
    "hepatitis", "cirrhosis", "fibrosis", "inflammation", "inflammatory",
    "autoimmune", "hiv", "covid", "sars-cov", "ebola", "malaria",
    "tuberculosis", "sepsis", "psoriasis", "asthma", "copd",
    "myocardial", "ischemia", "stroke", "epilepsy", "schizophrenia",
    "depression", "bipolar", "autism", "diabetes", "obesity",
    "nephritis", "lupus", "scleroderma", "arthritis", "gastritis",
    "colitis", "pancreatitis", "meningitis", "encephalitis",
    "bronchitis", "pneumonia", "sinusitis", "otitis", "dermatitis",
    "eczema", "vitiligo", "alopecia", "keratosis", "fibromyalgia",
    "sarcoidosis", "amyloidosis", "hemochromatosis", "porphyria",
    "down syndrome", "turner syndrome", "klinefelter",
    "fragile x", "marfan", "cystic fibrosis", "sickle cell",
    "thalassemia", "hemophilia", "retinoblastoma", "neuroblastoma",
    "medulloblastoma", "pheochromocytoma", "meningioma",
    "schwannoma", "neurofibroma", "liposarcoma", "rhabdomyosarcoma",
    "osteosarcoma", "chondrosarcoma", "ewing", "wilms",
    "hodgkin", "non-hodgkin", "burkitt", "mantle cell",
    "follicular", "diffuse", "large b-cell", "small lymphocytic",
    "chronic lymphocytic", "acute lymphoblastic", "acute myeloid",
    "chronic myeloid", "myelodysplastic", "myelofibrosis",
    "polycythemia", "thrombocythemia", "multiple myeloma",
    "blastic", "blast crisis", "refractory", "relapsed",
    "rejection", "dysfunction",
    # DRG / neuropathic pain — injured side
    "injured", "injury", "neuropathic", "neuropathy", "nerve injury",
    "cci", "snl", "sni", "axotomy", "spared nerve",
    "chronic constriction", "spinal nerve ligation",
    "sciatic nerve ligation", "paclitaxel", "oxaliplatin",
    "allodynia", "hyperalgesia", "pain",
]

CONTROL_KEYWORDS = [
    "control", "normal", "healthy", "wild type", "wild-type", "wt",
    "untreated", "vehicle", "sham", "placebo", "baseline",
    "non-diseased", "disease-free", "remission", "recovery",
    # DRG / neuropathic pain — uninjured / naive side
    "naive", "naïve", "uninjured", "intact", "contralateral",
    "ipsilateral naive", "no injury",
]

# Patterns to match sample types (tissue, cell line, etc.)
SAMPLE_TYPE_PATTERNS = [
    r"(?:tissue|tissues)\s*[:;]?\s*(.+)",
    r"(?:cell\s*line|lines)\s*[:;]?\s*(.+)",
    r"(?:sample\s*type)\s*[:;]?\s*(.+)",
    r"(?:material)\s*[:;]?\s*(.+)",
]


@dataclass
class SampleGroup:
    """A group of samples (disease or control)."""
    name: str
    samples: list[str] = field(default_factory=list)
    characteristics: list[str] = field(default_factory=list)
    confidence: float = 0.0


class MetadataExtractionError(Exception):
    """Raised when metadata extraction fails."""


def extract_sample_groups(
    characteristics: dict[str, list[str]],
    disease_term: str | None = None,
) -> tuple[SampleGroup, SampleGroup]:
    """Classify samples into disease and control groups from characteristics.

    Args:
        characteristics: Dict mapping sample_name -> list of characteristics strings.
        disease_term: Optional disease term to boost classification confidence.

    Returns:
        Tuple of (control_group, disease_group).

    Raises:
        MetadataExtractionError: If grouping fails completely.
    """
    if not characteristics:
        raise MetadataExtractionError("No sample characteristics provided")

    n_samples = len(characteristics)
    if n_samples < 2:
        raise MetadataExtractionError(f"Need at least 2 samples, got {n_samples}")

    disease_score = _score_samples_for_disease(characteristics, disease_term)

    control_samples = [s for s, sc in disease_score.items() if sc < 0]
    disease_samples = [s for s, sc in disease_score.items() if sc > 0]
    uncertain_samples = [s for s, sc in disease_score.items() if sc == 0]

    # Ensure both control and disease groups exist
    if not control_samples or not disease_samples:
        raise MetadataExtractionError(
            f"Cannot classify samples into control and disease groups "
            f"(found {len(control_samples)} control, {len(disease_samples)} disease)"
        )

    # Reassign uncertain samples based on characteristics similarity
    for sample in uncertain_samples:
        chars = characteristics[sample]
        best_group = _classify_uncertain_sample(chars, characteristics, control_samples, disease_samples)
        if best_group == "control":
            control_samples.append(sample)
        else:
            disease_samples.append(sample)

    control = SampleGroup(
        name="control",
        samples=control_samples,
        confidence=_group_confidence(control_samples, characteristics),
    )
    disease = SampleGroup(
        name="disease",
        samples=disease_samples,
        confidence=_group_confidence(disease_samples, characteristics),
    )

    logger.info(
        "Grouped %d samples: %d control, %d disease",
        n_samples, len(control_samples), len(disease_samples),
    )
    return control, disease


def _score_samples_for_disease(
    characteristics: dict[str, list[str]],
    disease_term: str | None = None,
) -> dict[str, int]:
    """Score each sample: +1 for disease keyword, -1 for control keyword, 0 uncertain."""
    scores = {}
    for sample, chars in characteristics.items():
        score = 0
        text = " ".join(chars).lower()

        for kw in DISEASE_KEYWORDS:
            if kw in text:
                score += 1
                break

        for kw in CONTROL_KEYWORDS:
            if kw in text:
                score -= 1
                break

        if disease_term:
            dt = disease_term.lower()
            if dt in text:
                score += 2

        scores[sample] = score
    return scores


def _classify_uncertain_sample(
    chars: list[str],
    all_chars: dict[str, list[str]],
    control_samples: list[str],
    disease_samples: list[str],
) -> str:
    """Classify an uncertain sample by majority vote of most similar classified samples."""
    text = " ".join(chars).lower()

    ctrl_overlap = 0
    for cs in control_samples:
        ctrl_text = " ".join(all_chars[cs]).lower()
        ctrl_overlap += sum(1 for w in ctrl_text.split() if w in text)

    dis_overlap = 0
    for ds in disease_samples:
        dis_text = " ".join(all_chars[ds]).lower()
        dis_overlap += sum(1 for w in dis_text.split() if w in text)

    return "control" if ctrl_overlap >= dis_overlap else "disease"


def _group_confidence(samples: list[str], all_chars: dict[str, list[str]]) -> float:
    """Estimate classification confidence by characteristic consistency."""
    all_words = []
    for s in samples:
        text = " ".join(all_chars[s]).lower()
        all_words.extend(text.split())
    if not all_words:
        return 0.0
    unique_ratio = len(set(all_words)) / len(all_words)
    return round(1.0 - unique_ratio, 3)
