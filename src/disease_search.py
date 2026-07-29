"""Search NCBI GEO for datasets matching a disease term."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass

from Bio import Entrez

from .config import PipelineConfig

logger = logging.getLogger(__name__)

VALID_GSE_PATTERN = r"^GSE\d+$"


@dataclass
class GEOSearchResult:
    """A single GEO dataset found by search."""
    gse_id: str
    title: str
    summary: str
    taxon: str
    gds_type: str
    n_samples: int = 0


class EntrezError(Exception):
    """Raised when an NCBI Entrez request fails."""


def search_gse_by_disease(
    disease: str,
    config: PipelineConfig,
    max_results: int | None = None,
) -> list[GEOSearchResult]:
    """Search GEO DataSets for GSE accessions matching a disease term.

    Args:
        disease: Disease name to search for.
        config: Pipeline configuration with email, API key, rate limits.
        max_results: Override max results from config.

    Returns:
        List of GEOSearchResult with dataset metadata.

    Raises:
        ValueError: If disease is empty.
        EntrezError: If NCBI request fails after retries.
    """
    if not disease or not disease.strip():
        raise ValueError("Disease term must be a non-empty string")

    Entrez.email = config.email
    if config.ncbi_api_key:
        Entrez.api_key = config.ncbi_api_key

    retmax = max_results or config.search.max_results
    delay = config.rate_limit.entrez_delay

    query = (
        f'"{disease}"[Title/Abstract] '
        f"AND gds[Filter] "
        f"AND {config.search.organism}[Organism]"
    )

    logger.info("Searching GEO for: %s", disease)

    last_error = None
    for attempt in range(1, config.rate_limit.max_retries + 1):
        try:
            with Entrez.esearch(db="gds", term=query, retmax=retmax, usehistory="y") as handle:
                results = Entrez.read(handle)
            break
        except Exception as e:
            last_error = e
            logger.warning("Search attempt %d failed: %s", attempt, e)
            if attempt < config.rate_limit.max_retries:
                time.sleep(delay * attempt)
    else:
        raise EntrezError(f"Search failed after {config.rate_limit.max_retries} attempts: {last_error}")

    count = int(results["Count"])
    logger.info("Found %d datasets for '%s'", count, disease)

    if count == 0:
        return []

    id_list = results["IdList"]
    time.sleep(delay)

    # Fetch summaries with retries
    records = []
    for uid in id_list:
        summary = None
        for attempt in range(1, config.rate_limit.max_retries + 1):
            try:
                with Entrez.esummary(db="gds", id=uid) as handle:
                    summary = Entrez.read(handle)
                time.sleep(delay)
                break
            except Exception as e:
                logger.warning("Attempt %d to fetch summary for UID %s failed: %s", attempt, uid, e)
                if attempt < config.rate_limit.max_retries:
                    time.sleep(delay * attempt)

        if summary is not None:
            record = _parse_gds_summary(summary)
            if record:
                records.append(record)

    logger.info("Parsed %d GSE records for '%s'", len(records), disease)
    return records


def _parse_gds_summary(summary) -> GEOSearchResult | None:
    """Extract useful fields from a GDS eSummary record."""
    try:
        record = summary[0]
    except (IndexError, KeyError):
        return None

    gse_id = record.get("GSE", "")
    if not gse_id:
        accession = record.get("Accession", "")
        if accession.startswith("GSE"):
            gse_id = accession
        else:
            return None

    if not gse_id.startswith("GSE"):
        gse_id = f"GSE{gse_id}"

    n_samples = record.get("n_samples", 0)
    if hasattr(n_samples, "__int__"):
        n_samples = int(n_samples)

    return GEOSearchResult(
        gse_id=gse_id,
        title=record.get("title", ""),
        summary=record.get("summary", ""),
        taxon=record.get("taxon", ""),
        gds_type=record.get("gdsType", ""),
        n_samples=n_samples,
    )


def search_diseases(
    diseases: list[str],
    config: PipelineConfig,
    max_results: int | None = None,
) -> dict[str, list[GEOSearchResult]]:
    """Search for multiple diseases.

    Returns:
        Dict mapping disease name -> list of GEOSearchResult.
    """
    return {
        d: search_gse_by_disease(d, config, max_results)
        for d in diseases
    }
