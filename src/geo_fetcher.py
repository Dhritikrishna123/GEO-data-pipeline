"""Download and parse GEO Series records with caching and retry."""

from __future__ import annotations

import logging
import re
import time
from pathlib import Path

import GEOparse

from .config import PipelineConfig

logger = logging.getLogger(__name__)

GSE_ID_PATTERN = re.compile(r"^GSE\d+$")

# Known column names for gene symbol mapping in GPL platforms
GENE_SYMBOL_COLUMNS = ["Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment", "ORF", "SPOT_ID"]

# Values that indicate missing/invalid gene symbols
MISSING_GENE_VALUES = frozenset({"nan", "---", "", "n/a"})


class GEOFetchError(Exception):
    """Raised when a GEO dataset cannot be fetched."""


def fetch_gse(
    gse_id: str,
    config: PipelineConfig,
) -> GEOparse.GSE:
    """Download and parse a GEO Series (GSE) record.

    Caches SOFT files in config.paths.raw_dir to avoid re-downloading.

    Args:
        gse_id: GEO accession (e.g. "GSE1563").
        config: Pipeline configuration.

    Returns:
        Parsed GSE object.

    Raises:
        ValueError: If gse_id is not a valid GSE accession.
        GEOFetchError: If download or parsing fails after retries.
    """
    if not GSE_ID_PATTERN.match(gse_id):
        raise ValueError(f"Invalid GSE accession: {gse_id!r} (expected format: GSE<number>)")

    raw_path = Path(config.paths.raw_dir)
    raw_path.mkdir(parents=True, exist_ok=True)

    # Try cached file first
    cached = list(raw_path.glob(f"{gse_id}_family.soft.gz"))
    if cached:
        logger.info("Loading cached %s from %s", gse_id, cached[0])
        try:
            return GEOparse.get_GEO(filepath=str(cached[0]))
        except Exception as e:
            logger.warning("Cached file corrupted for %s: %s", gse_id, e)

    # Download with retry
    delay = config.rate_limit.geo_delay
    last_error = None

    for attempt in range(1, config.rate_limit.max_retries + 1):
        try:
            logger.info("Downloading %s (attempt %d/%d)...", gse_id, attempt, config.rate_limit.max_retries)
            gse = GEOparse.get_GEO(geo=gse_id, destdir=str(raw_path), silent=True)
            time.sleep(delay)
            return gse
        except Exception as e:
            last_error = e
            logger.warning("Download attempt %d failed for %s: %s", attempt, gse_id, e)
            if attempt < config.rate_limit.max_retries:
                wait_time = delay * (2 ** (attempt - 1))
                logger.info("Retrying in %.1fs...", wait_time)
                time.sleep(wait_time)

    raise GEOFetchError(f"Failed to download {gse_id} after {config.rate_limit.max_retries} attempts: {last_error}")


def get_expression_table(gse: GEOparse.GSE) -> dict[str, list[tuple[str, float]]]:
    """Extract expression data from all samples in a GSE.

    Returns:
        Dict mapping sample_name -> list of (probe_id, value) tuples.
        Empty dict if no expression data found.
    """
    samples = {}
    for gsm_name, gsm in gse.gsms.items():
        if len(gsm.table) > 0 and "VALUE" in gsm.table.columns:
            table = gsm.table[["ID_REF", "VALUE"]].dropna()
            if len(table) > 0:
                samples[gsm_name] = list(zip(table["ID_REF"], table["VALUE"].astype(float)))
    return samples


def get_platform_table(gse: GEOparse.GSE) -> dict[str, str] | None:
    """Extract probe-to-gene mapping from the GPL platform.

    Returns:
        Dict mapping probe_id -> gene_symbol, or None if not found.
    """
    for gpl_name, gpl in gse.gpls.items():
        table = gpl.table

        for col in GENE_SYMBOL_COLUMNS:
            if col in table.columns:
                mapping = {}
                for _, row in table.iterrows():
                    probe = str(row["ID"])
                    gene = str(row[col]).strip()
                    if gene and gene.lower() not in MISSING_GENE_VALUES:
                        if "//" in gene:
                            gene = gene.split("//")[0].strip()
                        mapping[probe] = gene
                if mapping:
                    logger.debug("Mapped %d probes via column '%s' in %s", len(mapping), col, gpl_name)
                    return mapping

        logger.warning("No gene symbol column found in %s", gpl_name)
    return None


def get_sample_characteristics(gse: GEOparse.GSE) -> dict[str, list[str]]:
    """Extract sample characteristics for grouping (disease vs control).

    Returns:
        Dict mapping sample_name -> list of characteristics strings.
    """
    chars = {}
    for gsm_name, gsm in gse.gsms.items():
        sample_chars = []
        for key in sorted(gsm.metadata.keys()):
            if key.startswith("characteristics_ch"):
                for val in gsm.metadata[key]:
                    if val.strip():
                        sample_chars.append(val.strip())
        chars[gsm_name] = sample_chars
    return chars
