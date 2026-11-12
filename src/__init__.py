"""GEO Data Pipeline."""

from .config import DRGConfig, PipelineConfig
from .disease_search import (
    DRG_SEARCH_QUERIES,
    GEOSearchResult,
    search_drg_studies,
    search_gse_by_disease,
)
from .drg_qc import DRGQualityChecker, QCResult
from .expression_processor import (
    build_expression_matrix,
    build_raw_counts_matrix,
    normalize_expression,
)
from .geo_fetcher import GEOFetchError, fetch_gse
from .geo_supplementary import detect_supplementary_type, load_supplementary_counts
from .metadata_extractor import MetadataExtractionError, SampleGroup, extract_sample_groups
from .output_writer import OutputError, write_output, write_raw_counts
from .pipeline import process_drg_studies, process_single_gse, run_pipeline

__all__ = [
    # Config
    "DRGConfig",
    "PipelineConfig",
    # Search
    "DRG_SEARCH_QUERIES",
    "GEOSearchResult",
    "search_drg_studies",
    "search_gse_by_disease",
    # QC
    "DRGQualityChecker",
    "QCResult",
    # Expression
    "build_expression_matrix",
    "build_raw_counts_matrix",
    "normalize_expression",
    # Fetcher
    "GEOFetchError",
    "fetch_gse",
    # Supplementary
    "detect_supplementary_type",
    "load_supplementary_counts",
    # Metadata
    "MetadataExtractionError",
    "SampleGroup",
    "extract_sample_groups",
    # Output
    "OutputError",
    "write_output",
    "write_raw_counts",
    # Pipeline
    "process_drg_studies",
    "process_single_gse",
    "run_pipeline",
]
