"""Search NCBI GEO for datasets matching a disease term."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass

from Bio import Entrez

from .config import PipelineConfig

logger = logging.getLogger(__name__)

VALID_GSE_PATTERN = r"^GSE\d+$"

# Canonical DRG search queries — used by search_drg_studies()
DRG_SEARCH_QUERIES = [
    "human dorsal root ganglion single nucleus",
    "DRG scRNA-seq neuropathic pain",
    "injured vs control DRG",
]


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

    # GEO's gds database does not index Organism reliably, so we use free-text
    # search (all fields) which matches 200+ results vs 0-2 with [Organism].
    # Organism filtering is instead handled downstream by DRGQualityChecker.
    # We append "human" as a hint to bias toward Homo sapiens results.
    organism_hint = "human" if "homo sapiens" in config.search.organism.lower() else ""
    hint_part = f" {organism_hint}" if organism_hint and organism_hint not in disease.lower() else ""
    query = f"{disease}{hint_part}"

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


def search_drg_studies(
    config: PipelineConfig,
    max_results: int | None = None,
    queries: list[str] | None = None,
) -> list[GEOSearchResult]:
    """Search GEO using the canonical DRG queries and return deduplicated results.

    Runs each query in ``config.drg.search_queries`` (or the provided
    ``queries`` override) sequentially, deduplicates results by GSE ID, and
    warns if fewer than ``config.drg.min_human_studies`` studies are found.

    The three default queries are:
        - "human dorsal root ganglion single nucleus"
        - "DRG scRNA-seq neuropathic pain"
        - "injured vs control DRG"

    Args:
        config: Pipeline configuration.
        max_results: Maximum results per individual query. Defaults to
            ``config.search.max_results``.
        queries: Override the search query list. Defaults to
            ``config.drg.search_queries``.

    Returns:
        Deduplicated list of :class:`GEOSearchResult`, ordered by first-seen
        query (highest-priority query first).
    """
    active_queries = queries if queries is not None else config.drg.search_queries
    if not active_queries:
        logger.warning("No DRG search queries configured — returning empty results")
        return []

    seen: dict[str, GEOSearchResult] = {}  # gse_id -> result (preserves first-seen order)

    for query in active_queries:
        logger.info("DRG search query: %r", query)
        try:
            results = search_gse_by_disease(query, config, max_results)
        except EntrezError as e:
            logger.warning("DRG query %r failed: %s", query, e)
            continue

        for r in results:
            if r.gse_id not in seen:
                seen[r.gse_id] = r

    deduped = list(seen.values())
    n_human = sum(
        1 for r in deduped
        if "homo sapiens" in r.taxon.lower() or "human" in r.taxon.lower()
    )

    logger.info(
        "DRG search complete: %d unique studies (%d human) across %d queries",
        len(deduped), n_human, len(active_queries),
    )

    if n_human < config.drg.min_human_studies:
        logger.warning(
            "Found only %d human DRG studies; team target is >= %d. "
            "Consider broadening search queries or increasing max_results.",
            n_human, config.drg.min_human_studies,
        )

    return deduped


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

    # GDS summary returns all linked GSE IDs as "12345;67890;..." — take the first
    if isinstance(gse_id, str) and ";" in gse_id:
        gse_id = gse_id.split(";")[0].strip()
    gse_id = str(gse_id).strip()

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
