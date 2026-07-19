"""Write results to Parquet or CSV with validation."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

VALID_FORMATS = frozenset({"parquet", "csv"})


class OutputError(Exception):
    """Raised when output writing fails."""


def write_output(
    expression: pd.DataFrame,
    de_results: pd.DataFrame,
    metadata: dict,
    output_dir: str,
    gse_id: str,
    fmt: str = "parquet",
) -> dict[str, Path]:
    """Write expression and DE results to disk.

    Args:
        expression: Normalized expression DataFrame.
        de_results: Differential expression results.
        metadata: Dataset metadata.
        output_dir: Output directory path.
        gse_id: GEO accession.
        fmt: Output format ("parquet" or "csv").

    Returns:
        Dict mapping file type to output path.

    Raises:
        OutputError: If format is invalid or write fails.
    """
    if fmt not in VALID_FORMATS:
        raise OutputError(
            f"Invalid output format {fmt!r}. Must be one of: {sorted(VALID_FORMATS)}"
        )

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    prefix = out_path / gse_id
    output_files = {}

    try:
        # Write expression
        expr_path = prefix.with_suffix(f".expression.{fmt}")
        _write_dataframe(expression, expr_path, fmt)
        output_files["expression"] = expr_path

        # Write DE results
        de_path = prefix.with_suffix(f".de_results.{fmt}")
        _write_dataframe(de_results, de_path, fmt)
        output_files["de_results"] = de_path

        # Write metadata as JSON
        meta_path = prefix.with_suffix(".metadata.json")
        _write_metadata(metadata, meta_path)
        output_files["metadata"] = meta_path

    except OutputError:
        raise
    except Exception as e:
        raise OutputError(f"Failed to write output: {e}") from e

    logger.info("Output files written to %s", out_path)
    return output_files


def _write_dataframe(df: pd.DataFrame, path: Path, fmt: str) -> None:
    """Write a DataFrame to disk."""
    if df is None or df.empty:
        logger.warning("DataFrame is empty, writing empty file: %s", path)
        df = pd.DataFrame()

    if fmt == "parquet":
        df.to_parquet(path, index=True)
    elif fmt == "csv":
        df.to_csv(path, index=True)


def _write_metadata(metadata: dict, path: Path) -> None:
    """Write metadata dict as JSON."""
    import json

    clean = {k: v for k, v in metadata.items() if v is not None}
    with open(path, "w") as f:
        json.dump(clean, f, indent=2, default=str)
