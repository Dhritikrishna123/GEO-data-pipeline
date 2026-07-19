# GEO Data Pipeline

Fetch gene expression data from NCBI GEO for pharmaceutical target discovery.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Configuration

Set your email (required by NCBI):

```bash
export GEO_PIPELINE_EMAIL="you@example.com"
# Optional: increases rate limit from 3 to 10 req/s
export NCBI_API_KEY="your_api_key"
```

Or edit `config.yaml`:

```yaml
email: "you@example.com"
ncbi_api_key: null
```

## Usage

### Search for datasets

```bash
geo-pipeline search "breast cancer" --max-results 10
```

### Process a specific GSE

```bash
geo-pipeline process GSE1563 --disease "kidney transplant"
```

### Full pipeline (search + process)

```bash
geo-pipeline pipeline "breast cancer" "lung cancer" --max-results 10
```

### Programmatic use

```python
from src import PipelineConfig, search_gse_by_disease, process_single_gse

config = PipelineConfig.from_yaml("config.yaml")

# Search
results = search_gse_by_disease("breast cancer", config, max_results=5)
for r in results:
    print(f"{r.gse_id}: {r.title}")

# Process a specific dataset
result = process_single_gse("GSE1563", config, disease_term="kidney transplant")
print(result)
```

## Output

Each processed dataset produces three files in `data/output/`:

| File | Description |
|------|-------------|
| `GSE*.expression.parquet` | Normalized expression matrix (genes x samples) |
| `GSE*.de_results.parquet` | Differential expression results |
| `GSE*.metadata.json` | Dataset metadata |

### DE Results Schema

| Column | Description |
|--------|-------------|
| `logFC` | Log2 fold change (disease vs control) |
| `pval` | Raw p-value (Welch's t-test) |
| `padj` | BH-corrected p-value |
| `gene_symbol` | Gene symbol |
| `mean_a` | Mean expression in control group |
| `mean_b` | Mean expression in disease group |
| `label` | "Significant" or "NS" |

## Running Tests

```bash
pip install -e ".[test]"
pytest tests/
```

## Rate Limits

- Without API key: 3 requests/second to NCBI
- With API key: 10 requests/second
- SOFT files are cached in `data/raw/` to avoid re-downloading
