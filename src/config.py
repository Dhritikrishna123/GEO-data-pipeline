"""Pipeline configuration with env var support and validation."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class PathsConfig:
    raw_dir: str = "data/raw"
    processed_dir: str = "data/processed"
    output_dir: str = "data/output"

    def ensure_dirs(self) -> None:
        """Create output directories if they don't exist."""
        for d in (self.raw_dir, self.processed_dir, self.output_dir):
            try:
                Path(d).mkdir(parents=True, exist_ok=True)
            except OSError as e:
                raise ValueError(f"Cannot create directory {d!r}: {e}") from e


@dataclass
class RateLimitConfig:
    entrez_delay: float = 0.34   # 3 req/s without API key
    geo_delay: float = 1.0
    max_retries: int = 3
    retry_backoff: float = 2.0

    def __post_init__(self) -> None:
        if self.entrez_delay < 0:
            raise ValueError(f"entrez_delay must be >= 0, got {self.entrez_delay}")
        if self.geo_delay < 0:
            raise ValueError(f"geo_delay must be >= 0, got {self.geo_delay}")
        if self.max_retries < 1:
            raise ValueError(f"max_retries must be >= 1, got {self.max_retries}")


@dataclass
class SearchConfig:
    max_results: int = 100
    organism: str = "Homo sapiens"


@dataclass
class OutputConfig:
    format: str = "parquet"  # "parquet" or "csv"

    def __post_init__(self) -> None:
        if self.format not in ("parquet", "csv"):
            raise ValueError(f"output format must be 'parquet' or 'csv', got {self.format!r}")


# Default DRG search queries — the three canonical evidence queries
_DEFAULT_DRG_QUERIES = [
    "human dorsal root ganglion single nucleus",
    "DRG scRNA-seq neuropathic pain",
    "injured vs control DRG",
]


@dataclass
class DRGConfig:
    """Configuration for DRG-specific disease evidence collection."""

    # Minimum number of human DRG studies we expect to find across all queries.
    # If fewer are found the pipeline warns but continues.
    min_human_studies: int = 3

    # Minimum number of cells/samples required per group (injured / control).
    # None means no hard threshold is applied.
    min_cells_per_group: int | None = None

    # Reject non-human studies during QC.
    require_human: bool = True

    # Reject studies that lack clear injured vs. healthy labels.
    require_labels: bool = True

    # The three canonical DRG search queries.
    search_queries: list[str] = field(default_factory=lambda: list(_DEFAULT_DRG_QUERIES))

    def __post_init__(self) -> None:
        if self.min_human_studies < 1:
            raise ValueError(f"min_human_studies must be >= 1, got {self.min_human_studies}")
        if self.min_cells_per_group is not None and self.min_cells_per_group < 1:
            raise ValueError(
                f"min_cells_per_group must be >= 1 or None, got {self.min_cells_per_group}"
            )
        if not self.search_queries:
            raise ValueError("drg.search_queries must not be empty")


@dataclass
class PipelineConfig:
    email: str = ""
    ncbi_api_key: str | None = None
    paths: PathsConfig = field(default_factory=PathsConfig)
    rate_limit: RateLimitConfig = field(default_factory=RateLimitConfig)
    search: SearchConfig = field(default_factory=SearchConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    drg: DRGConfig = field(default_factory=DRGConfig)

    def __post_init__(self) -> None:
        self.paths.ensure_dirs()

    def __init__(self, **kwargs: Any) -> None:
        """Allow dict values for nested config objects."""
        if isinstance(kwargs.get("paths"), dict):
            kwargs["paths"] = PathsConfig(**kwargs["paths"])
        if isinstance(kwargs.get("rate_limit"), dict):
            kwargs["rate_limit"] = RateLimitConfig(**kwargs["rate_limit"])
        if isinstance(kwargs.get("search"), dict):
            kwargs["search"] = SearchConfig(**kwargs["search"])
        if isinstance(kwargs.get("output"), dict):
            kwargs["output"] = OutputConfig(**kwargs["output"])
        if isinstance(kwargs.get("drg"), dict):
            drg_data = kwargs["drg"]
            kwargs["drg"] = DRGConfig(**drg_data)
        self.email = kwargs.get("email", "")
        self.ncbi_api_key = kwargs.get("ncbi_api_key")
        self.paths = kwargs.get("paths", PathsConfig())
        self.rate_limit = kwargs.get("rate_limit", RateLimitConfig())
        self.search = kwargs.get("search", SearchConfig())
        self.output = kwargs.get("output", OutputConfig())
        self.drg = kwargs.get("drg", DRGConfig())
        self.paths.ensure_dirs()

    @classmethod
    def from_yaml(cls, path: str | Path) -> "PipelineConfig":
        """Load config from a YAML file.

        Falls back to env vars and defaults if file not found.
        """
        path = Path(path)
        if not path.exists():
            logger_msg = f"Config file {path} not found, using defaults"
            import logging
            logging.getLogger(__name__).warning(logger_msg)
            return cls(
                email=os.environ.get("GEO_PIPELINE_EMAIL", ""),
                ncbi_api_key=os.environ.get("NCBI_API_KEY"),
            )

        try:
            with open(path) as f:
                data = yaml.safe_load(f) or {}
        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in {path}: {e}") from e

        return cls(
            email=data.get("email", os.environ.get("GEO_PIPELINE_EMAIL", "")),
            ncbi_api_key=data.get("ncbi_api_key") or os.environ.get("NCBI_API_KEY"),
            paths=data.get("paths", {}),
            rate_limit=data.get("rate_limit", {}),
            search=data.get("search", {}),
            output=data.get("output", {}),
            drg=data.get("drg", {}),
        )

    @classmethod
    def from_env(cls) -> "PipelineConfig":
        """Load config from environment variables only."""
        return cls(
            email=os.environ.get("GEO_PIPELINE_EMAIL", ""),
            ncbi_api_key=os.environ.get("NCBI_API_KEY"),
        )
