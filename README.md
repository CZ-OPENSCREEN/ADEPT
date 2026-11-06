# ADEPT: Automated Differential Expression Pipeline for high-throughput Transcriptomics screens

Docker-containerized pipeline for differential expression analysis of DRUG-seq data and related plate-based perturbation screens (DRUG-seq, PLATE-seq, etc.). 
With a single command, it takes quantified count matrices from Alevin/Alevin-Fry or STARsolo and produces publication-ready DE results, volcano plots, heatmaps, UpSet plots, 
and UMAP/PCA projections.

## Quick start

### 1. Get the image (first run only)

Choose **one** of the two options. Both give you an image named `adept:latest`, so all commands below work unchanged.

**Option A: from Docker Hub**

```bash
docker pull bojke91/adept:latest
docker tag bojke91/adept:latest adept:latest
```

The second line gives the downloaded image the short name `adept` used throughout this README. (Alternatively, skip it and write `bojke91/adept` instead of `adept` in every command.)

**Option B: from a file**

If you downloaded the image as a file from Zenodo [23161389](https://zenodo.org/records/23161389):

```bash
docker load -i adept_1.0_2026-10-20.tar.gz
```

This restores the image as both `adept:latest` and `adept:1.0.0`; no extra tagging is needed.

Check that it worked with `docker images adept`.

### 2. Configure and run

```bash
# Generate a config file interactively
docker run -it --rm -v $(pwd)/my_experiment:/data adept --configure

# Run the pipeline
docker run --rm -v $(pwd)/my_experiment:/data adept
```

The mounted directory (`my_experiment/`) must contain your `config.yml`, metadata file, and quantification output directory. Results are written to the output path specified in the config.

On Windows use `${PWD}` (PowerShell) or `%cd%` (cmd.exe) instead of `$(pwd)`.

## Running the example dataset

An example dataset is available on Zenodo: [23161389](https://zenodo.org/records/23161389). The count matrices are already quantified with Salmon alevin-fry (`cr-like` UMI resolution), so you can run ADEPT on them directly without any upstream processing.

### 1. Download the files

From the Zenodo record, download:

- `2026_01_23_fry_quant_cr_like.zip`: per-plate alevin-fry quantification output

From the this repository, download:

- `metadata_adept.csv`: well-level metadata
- `config.yml`: a ready-made configuration for this dataset

### 2. Set up the data directory

```bash
mkdir adept_example
unzip 2026_01_23_fry_quant_cr_like.zip -d adept_example/
mv metadata_adept.csv config.yml adept_example/
```

The directory should now look like this:

```
adept_example/
├── config.yml
├── metadata_adept.csv
└── 2026_01_23_fry_quant_cr_like/
    ├── VH02001612/
    ├── VH02001614/
    └── ...
```

### 3. Choose a config

Pick **one** of the following:

- **Use the provided `config.yml`** (recommended). It is already set up for this dataset and reproduces the published analysis.
- **Build one in the browser.** Open [`adept_configurator.html`](adept_documentation/adept_configurator.html) in any web browser (no installation needed), fill in the fields, and download the resulting `config.yml` into `adept_example/`. For the example data, set:
  - Metadata file: `metadata_adept.csv`
  - Matrix directory: `2026_01_23_fry_quant_cr_like`
  - Quantification tool: `alevin`
- **Use the command-line wizard:** `docker run -it --rm -v $(pwd)/adept_example:/data adept --configure`, and enter the same values as above.

A config built with the configurator or the wizard uses the default analysis settings, which can differ from the published run. To reproduce the publication results exactly, use the provided `config.yml` together with the versioned image `adept:1.0.0` (see [Reproducibility](#reproducibility)).

### 4. Run the pipeline

```bash
docker run --rm -v $(pwd)/adept_example:/data adept
```

When the run finishes, open `results/padj_<X>_LFC_<Y>/final_report_*.html` in the output directory. It links to every table and figure the pipeline produced.

> [!NOTE]
> The example `config.yml` sets `min_count` to 1, while the pipeline default is 5:
>
> ```yaml
> analysis:
>   de:
>     min_count: 1
> ```
>
> This permissive gene filter was used only for the cross-pipeline benchmark in the ADEPT publication, so that ADEPT tests every gene the comparator pipelines test. For your own data, keep the default of 5 (omit the field or set it to 5).
>
> The configurator and the `--configure` wizard do not set this field, so a config built with them uses 5. To reproduce the published example results, use the provided `config.yml`.
>
> The provided `config.yml` also sets `analysis.dmso.top_n: 3`, the value used in the publication. The pipeline default is now `"all"` (every non-flagged DMSO well, at least 6 per plate), which configs built with the configurator or the wizard use.

### Environment variables

| Variable | Purpose | Default |
|----------|---------|---------|
| `DRUGSEQ_CONFIG` | Path to `config.yml` inside the container | `/data/config.yml` |
| `DRUGSEQ_WORKERS` | Number of parallel `crew` workers | `max(2, min(8, cores - 1))` |

## What goes in the data directory

```
my_experiment/
├── config.yml                    # pipeline configuration (generated by --configure)
├── metadata.csv                  # well-level annotations (see format below)
└── quant/                        # quantification output
    ├── PLATE_001/                # one subdirectory per plate, named to match plate_barcode
    ├── PLATE_002/
    └── ...
```

Plate subdirectory names must match the `plate_barcode` values in the metadata. The expected contents of each plate subdirectory depend on `analysis.quant_tool` (see below).

## Metadata format

The metadata file (CSV, TSV, XLSX, or XLS) must contain one row per well with these columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `plate_barcode` | string | Unique plate identifier | `PLATE_001` |
| `treatment` | string | Compound name or vehicle label | `Brusatol`, `DMSO` |
| `treatment_dose` | numeric | Dose value | `0.1`, `1.0`, `10` |
| `treatment_dose_unit` | string | Dose unit | `uM`, `nM` |
| `timepoint` | string | Treatment duration | `24hr`, `48hr` |
| `plate_well` | string | Well position (letter + number) | `A01`, `H12` |
| `well_index` | integer | Numeric well index (1-based) | `1`, `384` |
| `cell_line_name` | string | Cell line identifier | `A549`, `HeLa` |
| `batch` | integer | Batch number (plates processed together) | `1`, `2` |
| `row` | string | Plate row letter | `A`, `P` |
| `col` | integer | Plate column number | `1`, `24` |

`row`, `col`, `plate_well` and `well_index` are validated when the file is read; the remaining columns are validated immediately afterwards. `plate_barcode`, `batch` and `well_index` must additionally be non-blank in every row. Merged or blank cells in an XLSX layout are the usual cause of a failure here, and the pipeline stops with an explicit message rather than failing further downstream.

Control wells (e.g. DMSO) must use the same label in the `treatment` column consistently across all plates and batches. This label must match the `control_compound` field in `config.yml`.

## Pipeline outputs

After a successful run, the output directory contains:

```
annotations/
└── gene_annotation.rds                   # cached BioMart annotation (shared across runs)

results/padj_01_LFC_1/                    # directory name derives from the configured thresholds
├── final_report_padj_01_LFC_1.html       # single-page index linking all outputs
│
├── DE_genes/                             # differential expression results (all genes, unfiltered)
│   ├── DE_batch_1_plate_PLATE_001.csv    # per-plate DE
│   ├── DE_batch_1_collapsed_dupCor.csv   # collapsed DE with plate correction (primary result)
│   ├── DE_batch_1_dose_trend.csv         # per-compound log-dose slope test
│   ├── *_skipped.csv                     # contrasts not fitted, with the reason
│   └── HTML_reports/                     # interactive filterable tables (significant genes only)
│       └── DE_batch_1_plate_*.html
│
├── Volcano_plots/
│   ├── <contrast>.svg                    # e.g. DrugA_1uM_24hr.svg, one per selected contrast
│   ├── <contrast>.html                   # interactive twin, with live padj / |LFC| sliders
│   └── volcano_DE_<treatment>_batch<b>.svg   # static multi-dose overview per compound
│
├── Heatmaps_DE/                          # top DE genes × doses expression heatmaps
│   ├── heatmap_<treatment>_<aggregate_by>.html
│   └── heatmap_<treatment>_<aggregate_by>.png
│
├── upset/                                # UpSet plots showing DE gene overlap across contrasts
│   └── <treatments>_upset.svg (+ .html)
│
├── UMAP_plots/                           # dimensionality reduction (wells and samples)
│   ├── umap_by_treatment_batch.html
│   ├── umap_by_leiden_batch.html
│   ├── umap_dmso_plate_qc_batch.html
│   ├── pca_by_treatment_batch.html
│   ├── pca_dmso_plate_qc_batch.html
│   └── sample_umap_DE_gene_basis.html    # samples projected on DE genes only
│
├── heatmaps/                             # plate-level QC heatmaps (UMI counts per well)
│   └── plate_<barcode>.html (+ .png)
│
├── de_count_heatmaps/                    # plate grids coloured by DE gene count per well
│   └── plate_de_<barcode>.html (+ .png)
│
├── knee_plots/                           # UMI knee plots per plate + well QC table
│   ├── <barcode>_knee_plot.png
│   └── well_qc.csv                       # per-well UMIs, genes, thresholds, keep flag, reason
│
├── dmso_wells/                           # DMSO QC diagnostics
│   ├── best_DMSO_keep_pca.csv            # wells selected as controls
│   ├── bad_DMSO_remove_pca.csv           # wells flagged as outliers
│   ├── dmso_qc_batch_*.csv               # per-well scoring details
│   └── dmso_dose_summary_batch_*.csv     # per-dose quality comparison
│
└── correlation_plots/                    # pairwise log2FC scatter matrices across doses
    └── pairs_<treatment>_batch<b>.svg
```

The **collapsed dupCor CSV** (`DE_batch_*_collapsed_dupCor.csv`) is the primary DE result. It pools wells across plates within a batch and uses `duplicateCorrelation` to account for within-plate technical correlation. The per-plate CSVs are provided for QC and diagnostics, and for experiments with only one plate per batch.

Each DE CSV contains these columns:

| Column | Description |
|--------|-------------|
| `ENSEMBLID` | Ensembl gene ID |
| `LFC` | Log2 fold change (treatment vs DMSO) |
| `AveExpr` | Average log2 expression across contrast wells |
| `t` | Moderated t-statistic |
| `pval` | Raw p-value |
| `padj` | Benjamini-Hochberg adjusted p-value (within this contrast) |
| `B` | Log-odds of differential expression |
| `n_genes_tested` | Number of genes in this contrast's gene universe |
| `Sample` | Contrast label (compound_dose_timepoint) |
| `treatment` | Compound name |
| `treatment_dose` | Dose value |
| `treatment_dose_unit` | Dose unit |
| `timepoint` | Treatment duration |
| `batch` | Batch number |
| `SYMBOL` | Gene symbol |
| `GENENAME` | Full gene name |

The collapsed CSV additionally carries `plate_well` (the plate/well pairs pooled into the contrast); the per-plate CSV carries `plate_barcode` instead.

## Supported quantification tools

| Tool | Config value | Expected directory structure |
|------|-------------|------------------------------|
| Alevin-Fry (read with `fishpond::loadFry`) | `alevin` | `plate_dir/` — a complete alevin-fry output directory |
| Classic Alevin (read with `tximeta`) | `alevin_classic` | `plate_dir/alevin/quants_mat.gz` |
| Bare alevin-fry MatrixMarket export | `alevin_mtx` | `plate_dir/alevin/quants_mat.mtx` + `quants_mat_rows.txt` + `quants_mat_cols.txt` |
| STARsolo | `starsolo` | `plate_dir/Solo.out/Gene/filtered/` |

Every reader checks that the matrix dimensions match the accompanying gene and barcode label files and stops on any mismatch.

## Key methodological choices

**DMSO QC:** Control wells are scored using a combined robust-z + PCA outlier metric. The best DMSO dose is selected automatically per batch, wells above median + 2 MAD of the score are flagged, and by default (`top_n: "all"`) every remaining well is used as a control, with flagged wells readmitted in rank order if a plate would otherwise keep fewer than `min_keep` (6). A numeric `top_n` keeps that many best-ranked wells per plate instead. This is deterministic, fast, and avoids the computational cost of Monte Carlo permutation approaches.

**Well QC (optional):** `qc.well` drops failed or near-empty wells before DMSO selection and DE, on absolute floors (`min_umis`, `min_genes`) and, if configured, a plate-relative floor (`min_frac_median`). Dropped wells never enter a contrast but remain visible in the QC views. A config with no `qc.well` block is treated as QC disabled.

**Gene filtering:** A two-stage presence-aware strategy: (1) a global `edgeR::filterByExpr()` pass grouped by condition (`min.count = 5` by default), which (unlike a mean-based filter) keeps a gene that is expressed only in one condition; (2) a per-contrast `filterByExpr()` on just the wells being compared, so a gene cannot survive on the strength of counts in an unrelated treatment.

**Differential expression:** limma with `eBayes(robust = TRUE)`, with the model chosen from the design (`analysis.de.method: auto`, the default). Batches spanning two or more plates use limma-trend on log-CPM with `duplicateCorrelation` to model within-plate technical correlation, giving honest p-values that reflect the true effective sample size across plates. Single-plate batches and the per-plate results use limma-voom. Both choices were made on DMSO-versus-DMSO null contrasts; set `method: trend` or `method: voom` to force one model. Contrasts with insufficient replication or too few surviving genes are skipped and logged to `*_skipped.csv`.

**Dose-trend test:** Alongside the per-dose contrasts, each compound's full dose series is tested for a monotonic log-dose response on a deliberately more permissive gene universe (`min.count = 1`). This recovers low-expression, dose-coherent responders that no single per-dose contrast calls significant, while rejecting dose-independent noise.

**UMAP:** Plate-level batch effects are removed with `limma::removeBatchEffect()` (preserving treatment signal) before PCA → Leiden clustering → UMAP layout. Matching PCA plots are produced as a sanity check that the UMAP amplifies real structure rather than inventing it.

## Configuration

Use the interactive configurator (`--configure`) or create `config.yml` manually. The `paths` block (`data_dir`, `output_dir`, `cache_dir`) is mandatory! The pipeline stops immediately if any is missing. `experiment`, `analysis`, `input_files` and `filters` blocks are also mandatory. The `qc`, `selection` and `umap` blocks are optional.

See the [documentation page](adept_documentation/adept_documentation.html) for a full reference of all config fields.

## Requirements

Docker only. All R packages, Python dependencies, and system libraries are bundled in the image. The pipeline uses `renv.lock` for exact version pinning of all R packages.

For local development without Docker: R 4.6.0, the packages listed in `renv.lock`, Python 3 with `numpy`, `umap-learn`, and `scikit-learn`, plus Chromium and pandoc for the HTML/PNG report output.

## Reproducibility

The pipeline is fully reproducible via Docker + renv + targets:

- **Docker image** pins the OS, system libraries, R version (4.6.0), and Python packages. `latest` always points to the newest release; to reproduce a specific analysis exactly (e.g. the results in the ADEPT publication), use the versioned image `bojke91/adept:1.0.0` (Docker Hub) or `adept:1.0.0` (after `docker load`).
- **renv.lock** pins all R package versions (CRAN and Bioconductor).
- **targets** tracks the dependency graph so re-runs only recompute what changed.
- **config.yml** serves as a complete record of the analysis parameters.
