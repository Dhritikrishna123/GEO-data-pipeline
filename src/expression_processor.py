"""Expression matrix normalization and processing."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Minimum non-zero values per probe to keep it
MIN_NONZERO_RATIO = 0.5
# Maximum fraction of missing values per probe
MAX_MISSING_RATIO = 0.5


def normalize_expression(
    expression: pd.DataFrame,
    probe_to_gene: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Normalize expression matrix: log-transform, probe mapping, gene aggregation.

    Args:
        expression: Raw expression DataFrame with probes as index, samples as columns.
        probe_to_gene: Optional mapping of probe_id -> gene_symbol.

    Returns:
        Normalized expression DataFrame with gene symbols as index.
    """
    if expression.empty:
        return expression

    expr = expression.copy()

    # Log-transform if needed
    global_max = expr.values.max()
    if not np.isnan(global_max) and global_max > 30:
        expr = np.log2(expr + 1)

    # Map probes to genes
    if probe_to_gene:
        expr.index = expr.index.map(lambda p: probe_to_gene.get(str(p), str(p)))

    # Drop rows with no gene symbol
    expr = expr[expr.index != ""]

    # Aggregate duplicate genes (take mean)
    if expr.index.duplicated().any():
        expr = expr.groupby(expr.index).mean()

    # Filter low-expression probes
    expr = _filter_low_expression(expr)

    # Handle missing values
    expr = _handle_missing(expr)

    return expr


def build_expression_matrix(
    raw_samples: dict[str, list[tuple[str, float]]],
    probe_to_gene: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Build an expression DataFrame from raw sample data.

    Args:
        raw_samples: Dict mapping sample_name -> list of (probe_id, value) tuples.
        probe_to_gene: Optional probe-to-gene mapping.

    Returns:
        DataFrame with probes as index, samples as columns.
    """
    if not raw_samples:
        return pd.DataFrame()

    # Collect all probes
    all_probes = set()
    for vals in raw_samples.values():
        all_probes.update(p for p, _ in vals)

    # Build matrix
    matrix = {}
    for sample, vals in raw_samples.items():
        probe_dict = dict(vals)
        matrix[sample] = [probe_dict.get(p, np.nan) for p in all_probes]

    df = pd.DataFrame(matrix, index=sorted(all_probes))

    # Map probe IDs to gene symbols
    if probe_to_gene:
        df.index = df.index.map(lambda p: probe_to_gene.get(str(p), str(p)))

    # Drop rows with no gene symbol
    df = df[df.index != ""]

    # Aggregate duplicate genes
    if df.index.duplicated().any():
        df = df.groupby(df.index).mean()

    return df


def _filter_low_expression(
    expr: pd.DataFrame,
    min_nonzero_ratio: float = MIN_NONZERO_RATIO,
) -> pd.DataFrame:
    """Remove probes with too many zeros."""
    n_samples = len(expr.columns)
    if n_samples == 0:
        return expr

    min_nonzero = int(n_samples * min_nonzero_ratio)
    mask = (expr != 0).sum(axis=1) >= min_nonzero
    filtered = expr[mask]
    n_removed = len(expr) - len(filtered)
    if n_removed > 0:
        logger.info("Filtered %d low-expression probes", n_removed)
    return filtered


def _handle_missing(
    expr: pd.DataFrame,
    max_missing_ratio: float = MAX_MISSING_RATIO,
) -> pd.DataFrame:
    """Remove probes with too many missing values, fill remaining."""
    n_samples = len(expr.columns)
    if n_samples == 0:
        return expr

    max_missing = int(n_samples * max_missing_ratio)
    mask = expr.isna().sum(axis=1) <= max_missing
    filtered = expr[mask]
    filtered = filtered.fillna(0)

    n_removed = len(expr) - len(filtered)
    if n_removed > 0:
        logger.info("Removed %d probes with too many missing values", n_removed)
    return filtered
