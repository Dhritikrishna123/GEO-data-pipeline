"""Compute differential expression with Welch's t-test and BH FDR correction."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy import stats

from .metadata_extractor import SampleGroup

logger = logging.getLogger(__name__)

# A gene is considered "expressed" if at least this fraction of samples
# have non-zero values
EXPRESSION_THRESHOLD_RATIO = 0.5

# Thresholds for labeling DE results
LOGFC_THRESHOLD = 0.5
SIGNIFICANCE_THRESHOLD = 0.05


def compute_de(
    expression: pd.DataFrame,
    group_a: SampleGroup,
    group_b: SampleGroup,
) -> pd.DataFrame:
    """Compute differential expression between two sample groups.

    Uses Welch's t-test (unequal variances) and BH FDR correction.

    Args:
        expression: DataFrame with genes as index, samples as columns.
        group_a: Control/normal group.
        group_b: Disease/test group.

    Returns:
        DataFrame with columns: logFC, pval, padj, gene_symbol, mean_a, mean_b, label.
    """
    n_a = len(group_a.samples)
    n_b = len(group_b.samples)
    if n_a < 2 or n_b < 2:
        logger.warning(
            "Not enough samples for DE analysis: group_a=%d, group_b=%d",
            n_a, n_b,
        )
        return pd.DataFrame()

    if n_a + n_b < 10:
        logger.info(
            "Low total sample count (%d) - results may be unreliable",
            n_a + n_b,
        )

    all_samples = [s for s in group_a.samples if s in expression.columns]
    all_samples += [s for s in group_b.samples if s in expression.columns]
    expr = expression[all_samples].dropna()
    expr = expr[expr.index != ""]

    if len(expr) == 0:
        return pd.DataFrame()

    # Log-transform if raw counts detected
    expr = _ensure_log_scale(expr)

    # Filter low-expression genes
    expr = _filter_low_expression(expr, len(all_samples))

    if len(expr) == 0:
        return pd.DataFrame()

    # Compute t-test per gene
    vals_a = expr[group_a.samples].values
    vals_b = expr[group_b.samples].values

    t_stats, p_vals = stats.ttest_ind(vals_a, vals_b, axis=1, equal_var=False)

    t_stats = np.nan_to_num(t_stats, nan=0.0)
    p_vals = np.nan_to_num(p_vals, nan=1.0)

    # BH FDR correction
    padj = _benjamini_hochberg(p_vals)

    # logFC: already log-transformed, so use difference
    mean_a = vals_a.mean(axis=1)
    mean_b = vals_b.mean(axis=1)
    log_fc = mean_b - mean_a

    result = pd.DataFrame({
        "logFC": log_fc,
        "pval": p_vals,
        "padj": padj,
        "gene_symbol": expr.index,
        "mean_a": mean_a,
        "mean_b": mean_b,
    })

    result["label"] = "NS"
    result.loc[
        (result["padj"] < SIGNIFICANCE_THRESHOLD) & (result["logFC"].abs() >= LOGFC_THRESHOLD),
        "label",
    ] = "Significant"

    return result


def _ensure_log_scale(expr: pd.DataFrame) -> pd.DataFrame:
    """Detect and apply log2 transformation if needed.

    If global max < 30, assume data is already log-transformed.
    Otherwise, apply log2(x + 1).
    """
    global_max = expr.values.max()
    if global_max < 30:
        return expr
    return np.log2(expr + 1)


def _filter_low_expression(
    expr: pd.DataFrame,
    n_samples: int,
    threshold_ratio: float = EXPRESSION_THRESHOLD_RATIO,
) -> pd.DataFrame:
    """Remove genes with low expression across samples.

    A gene must have non-zero values in at least `threshold_ratio * n_samples` samples.
    """
    min_nonzero = int(n_samples * threshold_ratio)
    mask = (expr != 0).sum(axis=1) >= min_nonzero
    return expr[mask]


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """BH FDR correction for multiple testing.

    Returns:
        Array of adjusted p-values in original order.
    """
    n = len(p_values)
    ranked_order = np.argsort(p_values)
    ranked_pvals = p_values[ranked_order]
    adjusted = np.zeros(n)

    for i, p in enumerate(ranked_pvals):
        adjusted[i] = p * n / (i + 1)

    # Enforce monotonicity (step-up)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    # Map back to original order
    result = np.zeros(n)
    result[ranked_order] = adjusted
    return result
