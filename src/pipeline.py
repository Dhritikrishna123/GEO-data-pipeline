"""Pipeline orchestration for GEO data processing."""

from __future__ import annotations

import argparse
import logging
import os
import sys

from .config import PipelineConfig
from .de_analyzer import compute_de
from .disease_search import search_gse_by_disease
from .expression_processor import build_expression_matrix, normalize_expression
from .geo_fetcher import fetch_gse, get_expression_table, get_platform_table, get_sample_characteristics
from .metadata_extractor import MetadataExtractionError, extract_sample_groups
from .output_writer import OutputError, write_output

logger = logging.getLogger(__name__)


def process_single_gse(
    gse_id: str,
    config: PipelineConfig,
    disease_term: str | None = None,
) -> dict:
    """Process a single GEO Series.

    Args:
        gse_id: GEO accession (e.g. "GSE1563").
        config: Pipeline configuration.
        disease_term: Disease term used for metadata extraction context.

    Returns:
        Dict with keys: gse_id, n_genes, n_sig, n_total, files, disease, source.
    """
    if not gse_id or not gse_id.strip():
        raise ValueError("gse_id must be a non-empty string")
    gse_id = gse_id.strip()

    logger.info("Processing %s", gse_id)

    gse = fetch_gse(gse_id, config)
    raw_samples = get_expression_table(gse)
    probe_to_gene = get_platform_table(gse)
    sample_chars = get_sample_characteristics(gse)

    if not raw_samples:
        logger.warning("No expression data for %s", gse_id)
        return _empty_result(gse_id, "no_expression_data")

    expression = build_expression_matrix(raw_samples, probe_to_gene)
    if expression.empty:
        logger.warning("Empty expression matrix for %s", gse_id)
        return _empty_result(gse_id, "empty_expression_matrix")

    control, disease = _extract_groups(sample_chars, disease_term, gse_id)
    de_results = _compute_de_if_possible(expression, control, disease, gse_id)

    expression = normalize_expression(expression, probe_to_gene)
    disease_label = _resolve_disease_label(disease_term, disease)

    output_files = write_output(
        expression=expression,
        de_results=de_results if de_results is not None else expression.head(0),
        metadata={
            "gse_id": gse_id,
            "disease": disease_label,
            "platform": list(gse.gpls.keys()),
            "n_samples": len(raw_samples),
            "n_genes": len(expression),
        },
        output_dir=config.paths.output_dir,
        gse_id=gse_id,
        fmt=config.output.format,
    )

    n_sig = (
        int((de_results["padj"] < 0.05).sum())
        if de_results is not None and not de_results.empty
        else 0
    )

    return {
        "gse_id": gse_id,
        "n_genes": len(expression),
        "n_sig": n_sig,
        "n_total": len(de_results) if de_results is not None else 0,
        "files": {k: str(v) for k, v in output_files.items()},
        "disease": disease_label,
        "source": "GEO",
    }


def search_and_process(
    disease: str,
    config: PipelineConfig,
    max_results: int = 10,
    dataset_type: str | None = None,
) -> list[dict]:
    """Search GEO by disease and process all matching datasets."""
    if not disease or not disease.strip():
        raise ValueError("Disease term must be a non-empty string")

    results = search_gse_by_disease(disease, config, max_results=max_results)

    if dataset_type:
        results = [r for r in results if dataset_type.lower() in r.gds_type.lower()]
        logger.info("Filtered to %d '%s' datasets", len(results), dataset_type)

    outputs = []
    for r in results:
        try:
            result = process_single_gse(r.gse_id, config, disease_term=disease)
            outputs.append(result)
        except Exception as e:
            logger.error("Failed to process %s: %s", r.gse_id, e)
            outputs.append(_empty_result(r.gse_id, str(e)))

    return outputs


def run_pipeline(
    diseases: list[str],
    config: PipelineConfig,
    max_results: int = 10,
    dataset_type: str | None = None,
) -> list[dict]:
    """Run full pipeline: search + process for multiple diseases."""
    all_results = []
    for disease in diseases:
        results = search_and_process(disease, config, max_results, dataset_type)
        all_results.extend(results)
    return all_results


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _extract_groups(sample_chars, disease_term, gse_id):
    """Classify samples into control/disease groups, returning (None, None) on failure."""
    try:
        return extract_sample_groups(sample_chars, disease_term or "")
    except MetadataExtractionError as e:
        logger.warning("Metadata extraction failed for %s: %s", gse_id, e)
        return None, None


def _compute_de_if_possible(expression, control, disease, gse_id):
    """Run DE analysis if both groups have enough samples."""
    if control and disease and len(control.samples) >= 2 and len(disease.samples) >= 2:
        return compute_de(expression, control, disease)
    logger.warning("Insufficient groups for DE in %s", gse_id)
    return None


def _resolve_disease_label(disease_term, disease):
    """Pick the best disease label from input term or metadata."""
    if disease_term:
        return disease_term
    if disease and disease.samples:
        return disease.name
    return ""


def _empty_result(gse_id: str, reason: str) -> dict:
    """Return a placeholder result dict."""
    return {
        "gse_id": gse_id,
        "n_genes": 0,
        "n_sig": 0,
        "n_total": 0,
        "files": {},
        "disease": "",
        "source": "GEO",
        "error": reason,
    }


def _print_result(result: dict) -> None:
    """Pretty-print a processing result."""
    print(f"\n{result['gse_id']}:")
    print(f"  Disease: {result['disease'] or 'N/A'}")
    print(f"  Genes: {result['n_genes']:,}")
    print(f"  Significant: {result['n_sig']:,} / {result['n_total']:,}")
    if result.get("error"):
        print(f"  Error: {result['error']}")
    for ftype, fpath in result.get("files", {}).items():
        print(f"  {ftype}: {fpath}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="GEO Disease Dataset Processor",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    search_p = subparsers.add_parser("search", help="Search GEO for datasets")
    search_p.add_argument("disease", help="Disease name to search")
    search_p.add_argument("--max-results", type=int, default=10)
    search_p.add_argument("--config", default="config.yaml")

    process_p = subparsers.add_parser("process", help="Process a single GSE")
    process_p.add_argument("gse_id", help="GSE accession (e.g. GSE1563)")
    process_p.add_argument("--disease", default=None, help="Disease term for metadata")
    process_p.add_argument("--config", default="config.yaml")

    pipeline_p = subparsers.add_parser("pipeline", help="Search and process")
    pipeline_p.add_argument("diseases", nargs="+", help="Disease terms")
    pipeline_p.add_argument("--max-results", type=int, default=10)
    pipeline_p.add_argument("--dataset-type", default=None)
    pipeline_p.add_argument("--config", default="config.yaml")

    for sub in [search_p, process_p, pipeline_p]:
        sub.add_argument("--format", choices=["parquet", "csv"], default=None)

    return parser


def _load_config(args) -> PipelineConfig:
    """Load config from YAML, overlay CLI flags and env vars."""
    config = PipelineConfig.from_yaml(args.config)
    if args.format:
        config.output.format = args.format
    if not config.ncbi_api_key:
        config.ncbi_api_key = os.environ.get("NCBI_API_KEY")
    if not config.email:
        config.email = os.environ.get("GEO_PIPELINE_EMAIL", "")
    return config


def _cmd_search(args, config: PipelineConfig) -> None:
    results = search_gse_by_disease(args.disease, config, args.max_results)
    for r in results:
        print(f"{r.gse_id}\t{r.title}\tSamples: {r.n_samples}\t{r.gds_type}")
    print(f"\nTotal: {len(results)} datasets")


def _cmd_process(args, config: PipelineConfig) -> None:
    result = process_single_gse(args.gse_id, config, args.disease)
    _print_result(result)


def _cmd_pipeline(args, config: PipelineConfig) -> None:
    results = run_pipeline(args.diseases, config, args.max_results, args.dataset_type)
    print(f"\nProcessed {len(results)} datasets")
    total_sig = sum(r.get("n_sig", 0) for r in results)
    print(f"Total significant: {total_sig}")


def main() -> None:
    """CLI entry point with top-level error handling."""
    parser = _build_parser()
    args = parser.parse_args()

    try:
        config = _load_config(args)
    except (ValueError, FileNotFoundError) as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        sys.exit(1)

    dispatch = {
        "search": _cmd_search,
        "process": _cmd_process,
        "pipeline": _cmd_pipeline,
    }

    try:
        dispatch[args.command](args, config)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        logger.error("Command '%s' failed: %s", args.command, e)
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
