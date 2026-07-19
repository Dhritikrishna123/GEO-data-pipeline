"""GEO Data Pipeline - Fetch gene expression data for target discovery.

Public API:
    from src import PipelineConfig, search_gse_by_disease, process_single_gse

Usage:
    config = PipelineConfig.from_yaml("config.yaml")
    results = search_gse_by_disease("breast cancer", config)
    result = process_single_gse("GSE1563", config)
"""

from .config import PipelineConfig, PathsConfig, RateLimitConfig, SearchConfig, OutputConfig
from .disease_search import (
    GEOSearchResult,
    EntrezError,
    search_gse_by_disease,
    search_diseases,
)
from .geo_fetcher import (
    GEOFetchError,
    fetch_gse,
    get_expression_table,
    get_platform_table,
    get_sample_characteristics,
)
from .metadata_extractor import (
    MetadataExtractionError,
    SampleGroup,
    extract_sample_groups,
)
from .expression_processor import (
    build_expression_matrix,
    normalize_expression,
)
from .de_analyzer import (
    compute_de,
)
from .output_writer import (
    OutputError,
    write_output,
)
from .pipeline import (
    process_single_gse,
    search_and_process,
    run_pipeline,
)

__all__ = [
    # Config
    "PipelineConfig",
    "PathsConfig",
    "RateLimitConfig",
    "SearchConfig",
    "OutputConfig",
    # Search
    "GEOSearchResult",
    "EntrezError",
    "search_gse_by_disease",
    "search_diseases",
    # Fetch
    "GEOFetchError",
    "fetch_gse",
    "get_expression_table",
    "get_platform_table",
    "get_sample_characteristics",
    # Metadata
    "MetadataExtractionError",
    "SampleGroup",
    "extract_sample_groups",
    # Expression
    "build_expression_matrix",
    "normalize_expression",
    # DE
    "compute_de",
    # Output
    "OutputError",
    "write_output",
    # Pipeline
    "process_single_gse",
    "search_and_process",
    "run_pipeline",
]
