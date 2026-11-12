"""Download and parse supplementary count files for single-nucleus / single-cell studies.

Many GEO snRNA-seq / scRNA-seq datasets store raw counts as supplementary files
(e.g. 10x Genomics ``matrix.mtx.gz`` + ``barcodes.tsv.gz`` + ``features.tsv.gz``,
or ``*.h5`` / ``*.loom`` files) rather than as inline GSM expression tables.
This module detects and downloads those files and loads them into a dense
``pandas.DataFrame`` suitable for Geneformer training.

Supported formats
-----------------
- 10x MEX (``matrix.mtx.gz`` + ``barcodes.tsv.gz`` + ``features.tsv.gz``)
- 10x HDF5 (``*.h5`` — requires ``anndata`` / ``scanpy``)
- AnnData / loom (``*.loom`` — requires ``loompy``)

If none of the above are detected the function returns ``None`` and the caller
falls back to the inline GSM expression table.
"""

from __future__ import annotations

import gzip
import logging
import re
import time
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import requests

if TYPE_CHECKING:
    import GEOparse

logger = logging.getLogger(__name__)

# Patterns to detect supplementary file types
_MEX_MATRIX_PATTERN = re.compile(r"matrix\.mtx(\.gz)?$", re.I)
_MEX_BARCODES_PATTERN = re.compile(r"barcodes\.tsv(\.gz)?$", re.I)
_MEX_FEATURES_PATTERN = re.compile(r"(features|genes)\.tsv(\.gz)?$", re.I)
_H5_PATTERN = re.compile(r"\.(h5|hdf5|h5ad)$", re.I)
_LOOM_PATTERN = re.compile(r"\.loom$", re.I)
_CSV_COUNT_PATTERN = re.compile(r"\.(csv|txt)(\.gz)?$", re.I)
_TAR_PATTERN = re.compile(r"(_RAW\.tar|\.tar\.gz|\.tgz)$", re.I)

# Base URL for GEO supplementary FTP
_GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


class SupplementaryDownloadError(Exception):
    """Raised when a supplementary file download fails."""


def detect_supplementary_type(gse: "GEOparse.GSE") -> str | None:
    """Detect whether a GSE has supplementary count files.

    Checks both GSE-level and GSM-level supplementary_file metadata for known
    single-nucleus / single-cell file patterns, including .tar archives which
    GEO uses to bundle per-sample MEX directories.

    Args:
        gse: Parsed GEOparse GSE object.

    Returns:
        ``"mex"`` for 10x MEX format, ``"h5"`` for HDF5, ``"loom"`` for loom,
        ``"tar"`` for .tar bundles containing MEX/HDF5/loom files,
        or ``None`` if no supported supplementary files are detected.
    """
    # Collect from both GSE-level and per-GSM-level supplementary files
    supp_files = _collect_supplementary_urls(gse)
    # Also scan each GSM
    for gsm in gse.gsms.values():
        for key, vals in gsm.metadata.items():
            if "supplementary_file" in key.lower():
                for v in vals:
                    v = v.strip()
                    if v and v.lower() != "none":
                        supp_files.append(v)

    if not supp_files:
        return None

    # Check for .tar archives first (common GEO pattern for RAW bundles)
    if any(_TAR_PATTERN.search(u) for u in supp_files):
        return "tar"

    has_matrix = any(_MEX_MATRIX_PATTERN.search(u) for u in supp_files)
    has_barcodes = any(_MEX_BARCODES_PATTERN.search(u) for u in supp_files)
    has_features = any(_MEX_FEATURES_PATTERN.search(u) for u in supp_files)

    if has_matrix and has_barcodes and has_features:
        return "mex"

    if any(_H5_PATTERN.search(u) for u in supp_files):
        return "h5"

    if any(_LOOM_PATTERN.search(u) for u in supp_files):
        return "loom"

    return None


def load_supplementary_counts(
    gse: "GEOparse.GSE",
    cache_dir: str | Path,
    geo_delay: float = 1.0,
    max_retries: int = 3,
) -> pd.DataFrame | None:
    """Download and load raw count matrix from GEO supplementary files.

    Tries each supported format in order: tar → MEX → HDF5 → loom.
    Returns a dense DataFrame with genes as rows and cell barcodes as columns,
    or ``None`` if no supported supplementary files are found.

    Args:
        gse: Parsed GEOparse GSE object.
        cache_dir: Local directory for caching downloaded files.
        geo_delay: Seconds to wait between HTTP requests.
        max_retries: Number of download retries per file.

    Returns:
        Dense count DataFrame (genes × cells) or ``None``.
    """
    supp_type = detect_supplementary_type(gse)
    if supp_type is None:
        logger.debug("No supplementary count files detected for %s", gse.name)
        return None

    cache_path = Path(cache_dir) / gse.name / "supplementary"
    cache_path.mkdir(parents=True, exist_ok=True)

    # Collect all supplementary URLs (GSE + GSM level)
    supp_urls = _collect_supplementary_urls(gse)
    for gsm in gse.gsms.values():
        for key, vals in gsm.metadata.items():
            if "supplementary_file" in key.lower():
                for v in vals:
                    v = v.strip()
                    if v and v.lower() != "none":
                        supp_urls.append(v)

    if supp_type == "tar":
        return _load_tar(supp_urls, cache_path, geo_delay, max_retries)
    if supp_type == "mex":
        return _load_mex(supp_urls, cache_path, geo_delay, max_retries)
    if supp_type == "h5":
        return _load_h5(supp_urls, cache_path, geo_delay, max_retries)
    if supp_type == "loom":
        return _load_loom(supp_urls, cache_path, geo_delay, max_retries)

    return None


# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

def _collect_supplementary_urls(gse: "GEOparse.GSE") -> list[str]:
    """Extract all supplementary file URLs from GSE metadata."""
    urls: list[str] = []
    for key, vals in gse.metadata.items():
        if "supplementary_file" in key.lower():
            for v in vals:
                v = v.strip()
                if v and v.lower() != "none":
                    urls.append(v)
    return urls


def _ftp_to_https(url: str) -> str:
    """Convert NCBI FTP URLs to HTTPS equivalents.

    GEO stores supplementary files on NCBI FTP servers but all URLs are also
    accessible via HTTPS.  ``requests`` does not support FTP, so we rewrite.
    """
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov"):
        return url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov", 1)
    return url


def _download_file(
    url: str,
    dest: Path,
    geo_delay: float,
    max_retries: int,
) -> Path:
    """Download a URL to a local path with retry logic."""
    url = _ftp_to_https(url)  # Convert FTP → HTTPS for NCBI

    if dest.exists():
        logger.debug("Using cached supplementary file: %s", dest)
        return dest

    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            logger.info("Downloading supplementary file %s (attempt %d)", url, attempt)
            response = requests.get(url, timeout=120, stream=True)
            response.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in response.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
            time.sleep(geo_delay)
            return dest
        except Exception as e:
            last_error = e
            logger.warning("Download attempt %d failed for %s: %s", attempt, url, e)
            if attempt < max_retries:
                time.sleep(geo_delay * (2 ** (attempt - 1)))

    raise SupplementaryDownloadError(
        f"Failed to download {url} after {max_retries} attempts: {last_error}"
    )



def _load_tar(
    urls: list[str],
    cache_path: Path,
    geo_delay: float,
    max_retries: int,
) -> pd.DataFrame | None:
    """Download a GEO _RAW.tar archive and load count matrices from it.

    GEO commonly bundles per-sample MEX directories or count files into
    a single .tar archive (e.g. GSE168243_RAW.tar).  This function:
    1. Downloads and extracts the archive.
    2. Looks for MEX files (matrix.mtx.gz + barcodes + features) inside.
    3. Concatenates counts from all samples into one DataFrame.
    """
    import tarfile

    tar_url = next((u for u in urls if _TAR_PATTERN.search(u)), None)
    if tar_url is None:
        return None

    tar_dest = cache_path / Path(tar_url.split("/")[-1])
    extract_dir = cache_path / "extracted"
    extract_dir.mkdir(parents=True, exist_ok=True)

    try:
        tar_dest = _download_file(tar_url, tar_dest, geo_delay, max_retries)
    except SupplementaryDownloadError as e:
        logger.error("TAR download failed: %s", e)
        return None

    # Extract archive
    try:
        logger.info("Extracting %s → %s", tar_dest.name, extract_dir)
        with tarfile.open(str(tar_dest)) as tf:
            tf.extractall(str(extract_dir))
    except Exception as e:
        logger.error("Failed to extract TAR: %s", e)
        return None

    # Find all MEX / H5 files inside the extracted directory
    all_files = list(extract_dir.rglob("*"))
    all_file_strs = [str(f) for f in all_files]

    # Try MEX first
    matrices = [f for f in all_file_strs if _MEX_MATRIX_PATTERN.search(f)]
    barcodes_files = [f for f in all_file_strs if _MEX_BARCODES_PATTERN.search(f)]
    features_files = [f for f in all_file_strs if _MEX_FEATURES_PATTERN.search(f)]

    if matrices and barcodes_files and features_files:
        logger.info(
            "Found %d MEX matrix files in TAR — loading and concatenating",
            len(matrices),
        )
        dfs = []
        for mat_path in sorted(matrices):
            # Find matching barcodes/features in same directory
            mat_dir = str(Path(mat_path).parent)
            bc = next((f for f in barcodes_files if f.startswith(mat_dir)), None)
            feat = next((f for f in features_files if f.startswith(mat_dir)), None)
            if bc and feat:
                df = _load_mex([mat_path, bc, feat], cache_path, geo_delay, max_retries)
                if df is not None and not df.empty:
                    sample_name = Path(mat_path).parent.name
                    df.columns = [f"{sample_name}_{c}" for c in df.columns]
                    dfs.append(df)

        if dfs:
            combined = pd.concat(dfs, axis=1, join="outer").fillna(0)
            logger.info(
                "TAR: combined %d samples → %d genes x %d cells",
                len(dfs), combined.shape[0], combined.shape[1],
            )
            return combined

    # Try H5 files in the archive
    h5_files = [f for f in all_file_strs if _H5_PATTERN.search(f)]
    if h5_files:
        logger.info("Found %d H5 files in TAR", len(h5_files))
        dfs = []
        for h5_path in sorted(h5_files):
            df = _load_h5([h5_path], cache_path, geo_delay, max_retries)
            if df is not None and not df.empty:
                sample_name = Path(h5_path).stem
                df.columns = [f"{sample_name}_{c}" for c in df.columns]
                dfs.append(df)
        if dfs:
            combined = pd.concat(dfs, axis=1, join="outer").fillna(0)
            return combined

    # Try CSV/TXT count files (genes x cells, first column = gene names)
    csv_files = [
        f for f in all_file_strs
        if _CSV_COUNT_PATTERN.search(f)
        and not any(skip in f.lower() for skip in ["readme", "metadata", "annotation", "cluster", "umap", "pca", "barcode"])
    ]
    if csv_files:
        logger.info("Found %d CSV/TXT count files in TAR", len(csv_files))
        dfs = []
        for csv_path in sorted(csv_files):
            df = _load_csv_counts(csv_path)
            if df is not None and not df.empty:
                sample_name = Path(csv_path).stem.replace(".csv", "").replace(".txt", "")
                # Prefix cell barcodes with sample name to keep them unique
                df.columns = [f"{sample_name}_{c}" for c in df.columns]
                dfs.append(df)
        if dfs:
            combined = pd.concat(dfs, axis=1, join="outer").fillna(0).astype("int64")
            logger.info(
                "TAR (CSV): combined %d samples → %d genes x %d cells",
                len(dfs), combined.shape[0], combined.shape[1],
            )
            return combined

    logger.warning("TAR extracted but no supported count files found inside")
    return None


def _load_csv_counts(path: str) -> "pd.DataFrame | None":
    """Load a genes × cells count matrix from a (possibly gzipped) CSV/TXT file.

    Expected format:
    - First column: gene symbols (row names)
    - Remaining columns: cell barcodes
    - Values: integer read counts
    """
    import gzip as _gz

    try:
        df = pd.read_csv(
            path,
            index_col=0,
            compression="gzip" if path.endswith(".gz") else None,
        )
        # Validate — should have numeric values
        if df.shape[0] == 0 or df.shape[1] == 0:
            logger.warning("Empty CSV count matrix: %s", path)
            return None
        # Round to integer
        df = df.round().astype("int64")
        logger.info(
            "Loaded CSV counts from %s: %d genes x %d cells",
            Path(path).name, df.shape[0], df.shape[1],
        )
        return df
    except Exception as e:
        logger.warning("Failed to load CSV count file %s: %s", path, e)
        return None


def _load_mex(
    urls: list[str],
    cache_path: Path,
    geo_delay: float,
    max_retries: int,
) -> pd.DataFrame | None:
    """Load a 10x MEX format count matrix."""
    try:
        from scipy.io import mmread
        from scipy.sparse import issparse
    except ImportError:
        logger.warning("scipy not available; cannot load MEX format")
        return None

    matrix_url = next((u for u in urls if _MEX_MATRIX_PATTERN.search(u)), None)
    barcodes_url = next((u for u in urls if _MEX_BARCODES_PATTERN.search(u)), None)
    features_url = next((u for u in urls if _MEX_FEATURES_PATTERN.search(u)), None)

    if not all([matrix_url, barcodes_url, features_url]):
        logger.warning("Incomplete MEX file set — cannot load matrix")
        return None

    try:
        matrix_path = _download_file(
            matrix_url, cache_path / Path(matrix_url).name, geo_delay, max_retries
        )
        barcodes_path = _download_file(
            barcodes_url, cache_path / Path(barcodes_url).name, geo_delay, max_retries
        )
        features_path = _download_file(
            features_url, cache_path / Path(features_url).name, geo_delay, max_retries
        )
    except SupplementaryDownloadError as e:
        logger.error("MEX download failed: %s", e)
        return None

    try:
        # Load sparse matrix
        mat = mmread(str(matrix_path) if not str(matrix_path).endswith(".gz") else _gunzip_to_tmp(matrix_path))
        if issparse(mat):
            mat = mat.toarray()

        # Load barcodes
        barcodes = _read_tsv_column(barcodes_path, col=0)

        # Load features (gene names in column 1, gene IDs in column 0)
        features = _read_tsv_column(features_path, col=1) or _read_tsv_column(features_path, col=0)

        if len(features) != mat.shape[0] or len(barcodes) != mat.shape[1]:
            logger.warning(
                "MEX dimension mismatch: matrix %s, features %d, barcodes %d",
                mat.shape, len(features), len(barcodes),
            )
            return None

        df = pd.DataFrame(mat, index=features, columns=barcodes)
        logger.info("Loaded MEX counts: %d genes x %d cells", *df.shape)
        return df

    except Exception as e:
        logger.error("Failed to parse MEX files: %s", e)
        return None


def _load_h5(
    urls: list[str],
    cache_path: Path,
    geo_delay: float,
    max_retries: int,
) -> pd.DataFrame | None:
    """Load a 10x HDF5 or AnnData h5ad count matrix."""
    h5_url = next((u for u in urls if _H5_PATTERN.search(u)), None)
    if h5_url is None:
        return None

    try:
        h5_path = _download_file(
            h5_url, cache_path / Path(h5_url).name, geo_delay, max_retries
        )
    except SupplementaryDownloadError as e:
        logger.error("HDF5 download failed: %s", e)
        return None

    # Try anndata first (h5ad)
    if str(h5_path).endswith(".h5ad"):
        try:
            import anndata
            adata = anndata.read_h5ad(str(h5_path))
            df = pd.DataFrame(
                adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
                index=adata.obs_names,
                columns=adata.var_names,
            ).T  # genes × cells
            logger.info("Loaded h5ad counts: %d genes x %d cells", *df.shape)
            return df
        except ImportError:
            logger.warning("anndata not available; cannot load h5ad")
        except Exception as e:
            logger.error("Failed to parse h5ad: %s", e)

    # Try 10x HDF5
    try:
        import tables  # noqa: F401  (PyTables)
        import scipy.sparse as sp
        import h5py
        with h5py.File(str(h5_path), "r") as f:
            # 10x v3 format
            if "matrix" in f:
                grp = f["matrix"]
                data = grp["data"][:]
                indices = grp["indices"][:]
                indptr = grp["indptr"][:]
                shape = tuple(grp["shape"][:])
                barcodes = [b.decode() for b in grp["barcodes"][:]]
                features = [g.decode() for g in grp["features"]["name"][:]]
                mat = sp.csc_matrix((data, indices, indptr), shape=shape).toarray()
                df = pd.DataFrame(mat, index=features, columns=barcodes)
                logger.info("Loaded 10x HDF5 counts: %d genes x %d cells", *df.shape)
                return df
    except ImportError:
        logger.warning("h5py not available; cannot load 10x HDF5")
    except Exception as e:
        logger.error("Failed to parse 10x HDF5: %s", e)

    return None


def _load_loom(
    urls: list[str],
    cache_path: Path,
    geo_delay: float,
    max_retries: int,
) -> pd.DataFrame | None:
    """Load a loom file count matrix."""
    loom_url = next((u for u in urls if _LOOM_PATTERN.search(u)), None)
    if loom_url is None:
        return None

    try:
        loom_path = _download_file(
            loom_url, cache_path / Path(loom_url).name, geo_delay, max_retries
        )
    except SupplementaryDownloadError as e:
        logger.error("Loom download failed: %s", e)
        return None

    try:
        import loompy
        with loompy.connect(str(loom_path), mode="r") as ds:
            df = pd.DataFrame(ds[:, :], index=ds.ra["Gene"], columns=ds.ca["CellID"])
        logger.info("Loaded loom counts: %d genes x %d cells", *df.shape)
        return df
    except ImportError:
        logger.warning("loompy not available; cannot load loom file")
    except Exception as e:
        logger.error("Failed to parse loom file: %s", e)

    return None


def _read_tsv_column(path: Path, col: int) -> list[str]:
    """Read a single column from a (possibly gzipped) TSV file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(str(path), "rt", encoding="utf-8") as f:
            return [line.strip().split("\t")[col] for line in f if line.strip()]
    except Exception as e:
        logger.warning("Failed to read column %d from %s: %s", col, path, e)
        return []


def _gunzip_to_tmp(path: Path) -> str:
    """Decompress a .gz file to a temp path and return the path."""
    out_path = path.with_suffix("")
    if out_path.exists():
        return str(out_path)
    with gzip.open(str(path), "rb") as f_in, open(str(out_path), "wb") as f_out:
        f_out.write(f_in.read())
    return str(out_path)
