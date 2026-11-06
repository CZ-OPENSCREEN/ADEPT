library(yaml)

# readline() doesn't work under Rscript (non-interactive).
# readLines("stdin") works in both interactive and non-interactive modes,
# but returns character(0) (not "") when stdin is closed / no TTY attached —
# typically because the user ran the container without `-it`. Catch that and
# bail out with a clear message rather than a cryptic downstream R error.
ask <- function(prompt) {
  cat(prompt)
  input <- readLines(con = "stdin", n = 1)
  if (length(input) == 0) {
    cat("\n\n", strrep("=", 60), "\n", sep = "")
    cat("ERROR: no input received from stdin.\n\n")
    cat("The configuration wizard needs an interactive terminal.\n")
    cat("Re-run the container with the -it Docker flag:\n\n")
    cat("  docker run -it --rm -v $(pwd):/data adept:latest --configure\n\n")
    cat(strrep("=", 60), "\n", sep = "")
    quit(status = 1, save = "no")
  }
  input
}

cat("\n=== ADEPT Configuration Wizard ===\n\n")

# Template lives inside image
config <- yaml::read_yaml("/app/config_template.yml")

# ---- Experiment ----
cat("── Experiment ──\n")

title <- ask("Experiment title: ")
if (title != "") config$experiment$title <- title

cell_line <- ask("Cell line: ")
if (cell_line != "") config$experiment$cell_line <- cell_line

citation <- ask("Citation (optional): ")
if (citation != "") config$experiment$citation <- citation

# ---- Input files ----
cat("\n── Input Files ──\n")

metadata <- ask("Metadata file name: ")
if (metadata != "") config$input_files$metadata <- metadata

alevin_dir <- ask("Counts matrix file location: ")
if (alevin_dir != "") config$input_files$matrix_directory <- alevin_dir

# ---- Analysis ----
cat("\n── Analysis ──\n")

species <- ask("Biomart dataset (default hsapiens_gene_ensembl): ")
if (species != "") config$analysis$biomart_dataset <- species

control <- ask("Control compound (default DMSO): ")
if (control != "") config$analysis$control_compound <- control

plate_fmt <- ask("Plate format [96 / 384 / 1536] (default 384): ")
if (plate_fmt != "") config$analysis$plate_format <- as.integer(plate_fmt)

quant_tool <- ask("Quantification tool [alevin / alevin_classic / starsolo / alevin_mtx] (default alevin): ")
if (quant_tool != "") config$analysis$quant_tool <- quant_tool

# ---- DMSO QC ----
cat("\n── DMSO QC ──\n")

z_thresh <- ask("DMSO z-score threshold (default 3): ")
if (z_thresh != "") config$analysis$dmso$z_thresh <- as.numeric(z_thresh)

dmso_top_n <- ask("DMSO wells kept per plate ['all' = every non-flagged well, or a number] (default all): ")
if (dmso_top_n != "") {
  config$analysis$dmso$top_n <-
    if (identical(tolower(dmso_top_n), "all")) "all" else as.integer(dmso_top_n)
}

dmso_min_keep <- ask("Minimum DMSO controls per plate, used with 'all' (default 6): ")
if (dmso_min_keep != "") config$analysis$dmso$min_keep <- as.integer(dmso_min_keep)

dmso_npc <- ask("DMSO PCs used for outlier scoring (default 10, or 'auto'): ")
if (dmso_npc != "") {
  config$analysis$dmso$npc <-
    if (identical(dmso_npc, "auto")) "auto" else as.integer(dmso_npc)
}

dmso_pca_weight <- ask("Weight of the PCA component in the DMSO score (default 1.0): ")
if (dmso_pca_weight != "") config$analysis$dmso$pca_weight <- as.numeric(dmso_pca_weight)

# ---- DE model ----
cat("\n── DE model ──\n")

de_method <- ask("DE model [auto / trend / voom] (default auto): ")
if (de_method != "") config$analysis$de$method <- de_method

# ---- Well QC ----
cat("\n── Well QC (drop failed / near-empty wells before DE) ──\n")

well_qc_enabled <- ask("Enable well QC filter? [TRUE / FALSE] (default TRUE): ")
if (well_qc_enabled != "") config$qc$well$enabled <- as.logical(well_qc_enabled)

well_min_umis <- ask("Min UMIs per well (default 1000): ")
if (well_min_umis != "") config$qc$well$min_umis <- as.integer(well_min_umis)

well_min_genes <- ask("Min genes detected per well (default 200): ")
if (well_min_genes != "") config$qc$well$min_genes <- as.integer(well_min_genes)

well_frac <- ask("Drop wells below this fraction of plate-median UMIs [0 = off] (default 0): ")
if (well_frac != "") config$qc$well$min_frac_median <- as.numeric(well_frac)

well_method <- ask("UMI cut method [fixed / inflection / knee] (default fixed; inflection/knee only after inspecting the knee plots): ")
if (well_method != "") config$qc$well$method <- well_method

# ---- Filters ----
cat("\n── DE Filters ──\n")

lfc <- ask("LFC threshold (default 1): ")
if (lfc != "") config$filters$lfc_threshold <- as.numeric(lfc)

padj <- ask("padj threshold (default 0.1): ")
if (padj != "") config$filters$padj_threshold <- as.numeric(padj)

# ---- Selection ----
cat("\n── Selection (leave blank to use all) ──\n")

treatments_input <- ask("Treatments to select (comma-separated, e.g. CA-90-VK59,JC-43-ZO95): ")
if (treatments_input != "") {
  config$selection$treatments <- trimws(strsplit(treatments_input, ",")[[1]])
}

contrasts_input <- ask("Contrasts to plot (comma-separated, e.g. CA-90-VK59_3.16456uM_24hr): ")
if (contrasts_input != "") {
  config$selection$contrasts <- trimws(strsplit(contrasts_input, ",")[[1]])
}

batches_input <- ask("Batches to use (comma-separated, e.g. 1,2,3 or blank for all): ")
if (batches_input != "") {
  config$selection$batches <- as.integer(trimws(strsplit(batches_input, ",")[[1]]))
}

doses_input <- ask("Doses to select (comma-separated, e.g. 10.0,3.16456 or blank for all): ")
if (doses_input != "") {
  config$selection$doses <- as.numeric(trimws(strsplit(doses_input, ",")[[1]]))
}

# ---- Heatmaps ----
cat("\n── DE Heatmaps ──\n")

htop <- ask("Top N genes per heatmap (default 25): ")
if (htop != "") config$selection$heatmaps$top_n <- as.integer(htop)

agg <- ask("Aggregate by [none / concentration / replicate / timepoint] (default none): ")
if (agg != "") config$selection$heatmaps$aggregate_by <- agg

# ---- UMAP ----
cat("\n── UMAP ──\n")

n_pca <- ask("Number of PCs (default 20): ")
if (n_pca != "") config$umap$n_pca <- as.integer(n_pca)

n_neighbors <- ask("UMAP n_neighbors (default 20): ")
if (n_neighbors != "") config$umap$n_neighbors <- as.integer(n_neighbors)

leiden_k <- ask("Leiden k neighbors (default 15): ")
if (leiden_k != "") config$umap$leiden_k <- as.integer(leiden_k)

leiden_res <- ask("Leiden resolution (default 0.5): ")
if (leiden_res != "") config$umap$leiden_resolution <- as.numeric(leiden_res)

remove_batch <- ask("Remove batch effect for UMAP? [TRUE / FALSE] (default TRUE): ")
if (remove_batch != "") config$umap$remove_batch <- as.logical(remove_batch)

# ---- Save ----
write_yaml(config, "/data/config.yml")
cat("\n✔ config.yml successfully created in /data\n\n")