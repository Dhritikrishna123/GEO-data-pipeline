"""DRG study quality control checks.

Applies the team's three QC gates before a study is used for training:
    1. Human check  — taxon must be Homo sapiens
    2. Label check  — must have injured vs. healthy/sham labels
    3. Cell count   — each group must meet the minimum cell threshold (if set)
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

# Keywords that indicate an injured / disease condition in DRG studies
_INJURY_KEYWORDS = {
    "injured", "injury", "neuropathic", "neuropathy",
    "cci", "snl", "sni", "chronic constriction",
    "spinal nerve ligation", "spared nerve",
    "sciatic nerve ligation", "axotomy",
    "paclitaxel", "oxaliplatin",
    "allodynia", "hyperalgesia",
    # broader injury model terms
    "nerve injury", "nerve ligation", "nerve crush",
    "pain model", "pain behavior", "pain sensitivity",
    "peripheral nerve", "dorsal root ganglion injury",
    "nociception", "nociceptive", "sensitization",
    "thermal hyperalgesia", "mechanical allodynia",
    "chemotherapy-induced", "diabetic neuropathy",
    "post-surgical pain", "inflammatory pain",
}

# Keywords that indicate a healthy / sham / naive control condition
_CONTROL_KEYWORDS = {
    "sham", "naive", "naïve", "uninjured", "intact",
    "contralateral", "control", "healthy", "normal",
    # broader control terms
    "vehicle", "untreated", "baseline", "wild-type",
    "non-injured", "no treatment", "saline",
    "ipsilateral naive", "pain-free",
}

# Sequencing type detection patterns (checked against title + summary)
_SEQ_TYPE_PATTERNS: list[tuple[str, re.Pattern]] = [
    ("snRNA-seq", re.compile(r"\bsn\s*rna[\s-]?seq\b|\bsingle[\s-]?nucleus\b", re.I)),
    ("scRNA-seq", re.compile(r"\bsc\s*rna[\s-]?seq\b|\bsingle[\s-]?cell\b", re.I)),
    ("bulk RNA-seq", re.compile(r"\bbulk[\s-]?rna[\s-]?seq\b|\brna[\s-]?seq\b", re.I)),
    ("microarray", re.compile(r"\bmicroarray\b|\baffymetrix\b|\billumina\s+bead\b", re.I)),
]


@dataclass
class QCResult:
    """Result of a DRG QC check for a single study."""

    # Whether the study passed all mandatory checks
    passed: bool

    # Study-level metadata populated during QC
    study_accession: str = ""
    is_human: bool = False
    has_labels: bool = False
    n_control: int = 0
    n_injured: int = 0
    sequencing_type: str = "unknown"

    # Failure reasons (non-empty only when passed=False)
    reasons: list[str] = field(default_factory=list)

    # Advisory warnings that do not block the study
    warnings: list[str] = field(default_factory=list)


class DRGQualityChecker:
    """Apply the three-gate DRG QC pipeline to a study.

    Parameters
    ----------
    require_human:
        Reject studies that are not from *Homo sapiens*.
    require_labels:
        Reject studies whose sample metadata does not contain clearly
        distinguishable injured vs. control labels.
    min_cells_per_group:
        Minimum number of cells/samples required in **each** group
        (injured and control).  ``None`` disables this check.
    """

    def __init__(
        self,
        require_human: bool = True,
        require_labels: bool = True,
        min_cells_per_group: int | None = None,
    ) -> None:
        self.require_human = require_human
        self.require_labels = require_labels
        self.min_cells_per_group = min_cells_per_group

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def check(
        self,
        gse_id: str,
        taxon: str,
        title: str,
        summary: str,
        sample_characteristics: dict[str, list[str]],
    ) -> QCResult:
        """Run all QC checks and return a :class:`QCResult`.

        Args:
            gse_id: GEO accession string (e.g. "GSE12345").
            taxon: Organism name from GEO metadata (e.g. "Homo sapiens").
            title: Dataset title from GEO.
            summary: Dataset summary / abstract from GEO.
            sample_characteristics: Dict mapping sample_name -> list of
                characteristics strings (from
                :func:`geo_fetcher.get_sample_characteristics`).

        Returns:
            QCResult with ``passed`` set appropriately.
        """
        reasons: list[str] = []
        warnings: list[str] = []

        # 1. Human check
        is_human = _is_human(taxon)
        if not is_human and self.require_human:
            reasons.append(f"Non-human organism: {taxon!r} (require_human=True)")

        # 2. Label check — first from sample characteristics, then from study text
        n_control, n_injured = _count_groups(sample_characteristics)
        has_labels = n_control > 0 and n_injured > 0

        if not has_labels:
            # Fallback: scan title + summary for injury model keywords
            study_text = f"{title} {summary}".lower()
            has_injury_text = any(kw in study_text for kw in _INJURY_KEYWORDS)
            has_control_text = any(kw in study_text for kw in _CONTROL_KEYWORDS)
            if has_injury_text and has_control_text:
                has_labels = True
                # Estimate group counts from n_samples if available
                n_total = len(sample_characteristics)
                n_injured = max(1, n_total // 2)
                n_control = n_total - n_injured
                logger.info(
                    "%s: labels detected from study text (not sample chars); "
                    "estimating %d injured, %d control",
                    gse_id, n_injured, n_control,
                )

        if not has_labels and self.require_labels:
            reasons.append(
                f"Could not detect injured vs. control labels "
                f"(found {n_injured} injured, {n_control} control samples; "
                f"title/summary also lacked injury+control keywords)"
            )

        # 3. Cell count check (only if we have labels)
        if has_labels and self.min_cells_per_group is not None:
            if n_control < self.min_cells_per_group:
                reasons.append(
                    f"Control group has {n_control} samples, "
                    f"need >= {self.min_cells_per_group}"
                )
            if n_injured < self.min_cells_per_group:
                reasons.append(
                    f"Injured group has {n_injured} samples, "
                    f"need >= {self.min_cells_per_group}"
                )

        # Advisory warnings
        if n_control + n_injured < 10:
            warnings.append(
                f"Very low total sample count ({n_control + n_injured}); "
                "DE results may be unreliable"
            )

        # Sequencing type detection
        seq_type = _detect_sequencing_type(title, summary)

        passed = len(reasons) == 0

        result = QCResult(
            passed=passed,
            study_accession=gse_id,
            is_human=is_human,
            has_labels=has_labels,
            n_control=n_control,
            n_injured=n_injured,
            sequencing_type=seq_type,
            reasons=reasons,
            warnings=warnings,
        )

        if passed:
            logger.info(
                "QC PASSED %s | human=%s labels=%s injured=%d control=%d seq=%s",
                gse_id, is_human, has_labels, n_injured, n_control, seq_type,
            )
        else:
            logger.warning(
                "QC FAILED %s | %s",
                gse_id, "; ".join(reasons),
            )

        return result

    @classmethod
    def from_config(cls, config) -> "DRGQualityChecker":
        """Construct from a :class:`~config.PipelineConfig` instance."""
        return cls(
            require_human=config.drg.require_human,
            require_labels=config.drg.require_labels,
            min_cells_per_group=config.drg.min_cells_per_group,
        )


# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

def _is_human(taxon: str) -> bool:
    """Return True if the organism string refers to Homo sapiens."""
    lower = taxon.lower()
    return "homo sapiens" in lower or "human" in lower


def _count_groups(
    sample_characteristics: dict[str, list[str]],
) -> tuple[int, int]:
    """Count how many samples look like controls vs. injured.

    Returns:
        (n_control, n_injured)
    """
    n_control = 0
    n_injured = 0

    for _sample, chars in sample_characteristics.items():
        text = " ".join(chars).lower()

        has_injury = any(kw in text for kw in _INJURY_KEYWORDS)
        has_control = any(kw in text for kw in _CONTROL_KEYWORDS)

        if has_injury and not has_control:
            n_injured += 1
        elif has_control and not has_injury:
            n_control += 1
        elif has_injury and has_control:
            # Ambiguous — lean toward injured (e.g. "CCI vs sham" in same field)
            n_injured += 1

    return n_control, n_injured


def _detect_sequencing_type(title: str, summary: str) -> str:
    """Detect the sequencing type from the study title and summary."""
    text = f"{title} {summary}"
    for seq_type, pattern in _SEQ_TYPE_PATTERNS:
        if pattern.search(text):
            return seq_type
    return "unknown"
