"""Pipeline orchestration for GEO data processing."""

from __future__ import annotations

import argparse
import logging
import os
import sys

from .config import PipelineConfig
from .de_analyzer import compute_de
from .disease_search import search_gse_by_disease, search_drg_studies
from .drg_qc import DRGQualityChecker
from .expression_processor import build_expression_matrix, build_raw_counts_matrix, normalize_expression
from .geo_fetcher import fetch_gse, get_expression_table, get_platform_table, get_sample_characteristics
from .geo_supplementary import detect_supplementary_type, load_supplementary_counts
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
        Dict with keys: gse_id, n_genes, n_sig, n_total, files, disease, source,
        plus enriched fields: n_injured, n_control, is_human, sequencing_type,
        qc_passed, qc_reasons.
    """
    if not gse_id or not gse_id.strip():
        raise ValueError("gse_id must be a non-empty string")
    gse_id = gse_id.strip()

    logger.info("Processing %s", gse_id)

    gse = fetch_gse(gse_id, config)
    sample_chars = get_sample_characteristics(gse)

    # Determine taxon and study-level text for QC / sequencing type detection
    taxon = _get_gse_taxon(gse)
    title = _get_gse_field(gse, "title")
    summary = _get_gse_field(gse, "summary")

    # ------------------------------------------------------------------
    # Attempt to load supplementary counts (snRNA-seq / scRNA-seq)
    # ------------------------------------------------------------------
    supp_type = detect_supplementary_type(gse)
    supplementary_counts: "pd.DataFrame | None" = None  # noqa: F821
    if supp_type is not None:
        logger.info(
            "%s has %s supplementary files — attempting to load raw counts",
            gse_id, supp_type,
        )
        supplementary_counts = load_supplementary_counts(
            gse,
            cache_dir=config.paths.raw_dir,
            geo_delay=config.rate_limit.geo_delay,
            max_retries=config.rate_limit.max_retries,
        )

    # ------------------------------------------------------------------
    # Inline GSM expression table (bulk / microarray fallback)
    # ------------------------------------------------------------------
    raw_samples = get_expression_table(gse)
    probe_to_gene = get_platform_table(gse)

    # Choose primary expression source
    if supplementary_counts is not None and not supplementary_counts.empty:
        # Supplementary raw counts are authoritative
        import pandas as pd
        raw_counts_df = supplementary_counts.copy()
        expression = supplementary_counts.copy().astype(float)
        logger.info("Using supplementary counts (%s) for %s", supp_type, gse_id)
    elif raw_samples:
        # Fallback to inline GSM table
        raw_counts_df = build_raw_counts_matrix(raw_samples, probe_to_gene)
        expression = build_expression_matrix(raw_samples, probe_to_gene)
    else:
        logger.warning("No expression data for %s", gse_id)
        return _empty_result(gse_id, "no_expression_data")

    if expression.empty:
        logger.warning("Empty expression matrix for %s", gse_id)
        return _empty_result(gse_id, "empty_expression_matrix")

    control, disease = _extract_groups(sample_chars, disease_term, gse_id)
    de_results = _compute_de_if_possible(expression, control, disease, gse_id)

    expression = normalize_expression(expression, probe_to_gene if not supplementary_counts is not None else None)
    disease_label = _resolve_disease_label(disease_term, disease)

    n_control = len(control.samples) if control else 0
    n_injured = len(disease.samples) if disease else 0

    output_files = write_output(
        expression=expression,
        de_results=de_results if de_results is not None else expression.head(0),
        metadata={
            "gse_id": gse_id,
            "study_accession": gse_id,
            "disease": disease_label,
            "platform": list(gse.gpls.keys()),
            "n_samples": len(raw_samples) if raw_samples else expression.shape[1],
            "n_genes": len(expression),
            "n_control": n_control,
            "n_injured": n_injured,
            "is_human": "homo sapiens" in taxon.lower() or "human" in taxon.lower(),
            "taxon": taxon,
            "sequencing_type": supp_type or "bulk",
            "supplementary_format": supp_type,
        },
        output_dir=config.paths.output_dir,
        gse_id=gse_id,
        fmt=config.output.format,
        raw_counts=raw_counts_df if "raw_counts_df" in dir() else None,
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
        # Enriched fields
        "n_control": n_control,
        "n_injured": n_injured,
        "is_human": "homo sapiens" in taxon.lower() or "human" in taxon.lower(),
        "taxon": taxon,
        "sequencing_type": supp_type or "bulk",
        "study_accession": gse_id,
    }


def process_drg_studies(
    config: PipelineConfig,
    max_results: int | None = None,
) -> list[dict]:
    """Search, QC-filter, and process human DRG studies.

    This is the top-level entry point for the team's DRG evidence collection
    workflow.  It:

    1. Runs :func:`~disease_search.search_drg_studies` with the three
       canonical DRG queries from ``config.drg.search_queries``.
    2. Applies :class:`~drg_qc.DRGQualityChecker` to each candidate study
       — studies that fail QC are logged and skipped.
    3. Warns if fewer than ``config.drg.min_human_studies`` pass QC.
    4. Calls :func:`process_single_gse` for each passing study.
    5. Returns a list of result dicts enriched with QC metadata.

    Args:
        config: Pipeline configuration.
        max_results: Max results per search query (overrides config default).

    Returns:
        List of result dicts, one per study that passed QC.  Each dict
        includes all fields from :func:`process_single_gse` plus
        ``qc_passed``, ``qc_reasons``, and ``qc_warnings``.
    """
    logger.info(
        "Starting DRG study collection with queries: %s",
        config.drg.search_queries,
    )

    # Step 1: Search
    candidates = search_drg_studies(config, max_results=max_results)
    logger.info("DRG search returned %d unique candidate studies", len(candidates))

    if not candidates:
        logger.warning("No DRG candidates found — check search queries and NCBI connectivity")
        return []

    # Step 2: QC each candidate
    qc = DRGQualityChecker.from_config(config)
    passing = []

    for candidate in candidates:
        # Fetch minimal info for QC (we need sample characteristics)
        try:
            gse = fetch_gse(candidate.gse_id, config)
            sample_chars = get_sample_characteristics(gse)
        except Exception as e:
            logger.warning("Could not fetch %s for QC: %s — skipping", candidate.gse_id, e)
            continue

        qc_result = qc.check(
            gse_id=candidate.gse_id,
            taxon=candidate.taxon,
            title=candidate.title,
            summary=candidate.summary,
            sample_characteristics=sample_chars,
        )

        if qc_result.passed:
            passing.append((candidate, gse, qc_result))
        else:
            logger.warning(
                "Skipping %s (QC failed): %s",
                candidate.gse_id, "; ".join(qc_result.reasons),
            )

    n_passing = len(passing)
    logger.info("%d / %d DRG studies passed QC", n_passing, len(candidates))

    if n_passing < config.drg.min_human_studies:
        logger.warning(
            "Only %d studies passed QC; team target is >= %d. "
            "Results may be fragile — consider adding more search queries "
            "or relaxing QC thresholds.",
            n_passing, config.drg.min_human_studies,
        )

    if n_passing < 2:
        logger.warning(
            "Fewer than 2 studies available — cross-checking between studies "
            "is not possible. Proceed with caution."
        )

    # Step 3: Process each passing study
    outputs = []
    for candidate, _gse, qc_result in passing:
        try:
            result = process_single_gse(
                candidate.gse_id,
                config,
                disease_term="neuropathic pain",
            )
            result["qc_passed"] = True
            result["qc_reasons"] = []
            result["qc_warnings"] = qc_result.warnings
            result["n_control"] = qc_result.n_control
            result["n_injured"] = qc_result.n_injured
            result["sequencing_type"] = qc_result.sequencing_type
            outputs.append(result)
        except Exception as e:
            logger.error("Failed to process DRG study %s: %s", candidate.gse_id, e)
            err = _empty_result(candidate.gse_id, str(e))
            err["qc_passed"] = True
            err["qc_reasons"] = []
            err["qc_warnings"] = qc_result.warnings
            outputs.append(err)

    logger.info(
        "DRG pipeline complete: %d studies processed, %d with DE results",
        len(outputs),
        sum(1 for o in outputs if o.get("n_sig", 0) > 0),
    )
    return outputs


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

def _get_gse_taxon(gse) -> str:
    """Extract organism taxon from GSE metadata."""
    for key in ("organism_ch1", "organism", "taxon"):
        vals = gse.metadata.get(key, [])
        if vals:
            return vals[0].strip()
    # Fall back to checking GSMs
    for gsm in gse.gsms.values():
        for key in ("organism_ch1", "organism"):
            vals = gsm.metadata.get(key, [])
            if vals:
                return vals[0].strip()
    return ""


def _get_gse_field(gse, field: str) -> str:
    """Extract a string field from GSE metadata."""
    vals = gse.metadata.get(field, [])
    return " ".join(v.strip() for v in vals if v.strip())


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
        "n_control": 0,
        "n_injured": 0,
        "is_human": False,
        "taxon": "",
        "sequencing_type": "unknown",
        "study_accession": gse_id,
    }


def _print_result(result: dict) -> None:
    """Pretty-print a processing result."""
    print(f"\n{result['gse_id']}:")
    print(f"  Disease: {result['disease'] or 'N/A'}")
    print(f"  Genes: {result['n_genes']:,}")
    print(f"  Significant: {result['n_sig']:,} / {result['n_total']:,}")
    print(f"  Control samples: {result.get('n_control', 'N/A')}")
    print(f"  Injured samples: {result.get('n_injured', 'N/A')}")
    print(f"  Sequencing type: {result.get('sequencing_type', 'N/A')}")
    print(f"  Human: {result.get('is_human', 'N/A')}")
    if result.get("error"):
        print(f"  Error: {result['error']}")
    if result.get("qc_warnings"):
        for w in result["qc_warnings"]:
            print(f"  Warning: {w}")
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

    drg_p = subparsers.add_parser(
        "drg",
        help=(
            "Search GEO using the three canonical DRG queries, apply quality "
            "control, and process all passing human DRG studies"
        ),
    )
    drg_p.add_argument("--max-results", type=int, default=None,
                       help="Max results per query (default: config.search.max_results)")
    drg_p.add_argument("--config", default="config.yaml")

    for sub in [search_p, process_p, pipeline_p, drg_p]:
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


def _cmd_drg(args, config: PipelineConfig) -> None:
    """Run the DRG-specific pipeline."""
    print("\nRunning DRG study collection pipeline...")
    print(f"Search queries: {config.drg.search_queries}")
    print(f"QC: require_human={config.drg.require_human}, "
          f"require_labels={config.drg.require_labels}, "
          f"min_cells_per_group={config.drg.min_cells_per_group}")
    print()

    results = process_drg_studies(config, max_results=args.max_results)

    print(f"\n{'='*60}")
    print(f"DRG pipeline complete: {len(results)} studies processed")
    print(f"{'='*60}")

    for r in results:
        _print_result(r)

    # Summary table
    print(f"\n{'Study':<12} {'Human':<8} {'Seq type':<14} {'Control':<9} {'Injured':<9} {'Genes':<8} {'Sig DE'}")
    print("-" * 74)
    for r in results:
        print(
            f"{r['gse_id']:<12} "
            f"{str(r.get('is_human','?')):<8} "
            f"{r.get('sequencing_type','?'):<14} "
            f"{r.get('n_control','?'):<9} "
            f"{r.get('n_injured','?'):<9} "
            f"{r['n_genes']:<8} "
            f"{r['n_sig']}"
        )


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
        "drg": _cmd_drg,
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
