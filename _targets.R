# _targets.R

# Load core packages ------------------------------------------------------
library(targets)
library(tarchetypes)
library(future)
library(crew)


# Set pipeline options ---------------------------------------------------

n_workers <- as.integer(
  Sys.getenv(
    "DRUGSEQ_WORKERS",
    max(2, min(8, parallelly::availableCores() - 1))
  )
)

tar_option_set(
  packages = c(
    "DT", "tidyverse", "foreach", "limma", "edgeR",
    "umap", "data.table", "GGally", "plotly", "ggrepel", "viridis",
    "RANN", "igraph", "RColorBrewer", "patchwork", "UpSetR", "EnhancedVolcano",
    "biomaRt", "htmlwidgets", "readxl", "DropletUtils",
    "glue", "leiden", "svglite"
  ),
  format = "rds",
  error = "continue",
  controller = crew::crew_controller_local(
    workers            = n_workers,
    seconds_idle       = 60,
    garbage_collection = TRUE
  )
)

config_path <- Sys.getenv("DRUGSEQ_CONFIG", "/data/config.yml")

if (!file.exists(config_path)) {
  stop(
    "Config file not found: ", config_path, "\n\n",
    "Run configuration wizard:\n",
    "docker run -it --rm -v $(pwd):/data adept --configure\n"
  )
}

# NOTE: the targets store location is taken from config$paths$cache_dir by
# run_pipeline.R and passed to tar_make(store = ...). It is not set here,
# because tar_config_set() would write a _targets.yaml side-effect into /app
# every time this script is parsed and would then be overridden anyway.

# Source custom functions ------------------------------------------------
# Inside _targets.R
tar_source("R")

list(
  # ===============================================================
  # CONFIGURATION
  # ===============================================================
  tar_target(config_file, config_path, format = "file"),
  tar_target(config, read_config(config_file)),
  tar_target(
    meta_file,
    file.path(config$paths$data_dir, config$input_files$metadata),
    format = "file"
  ),
  tar_target(
    matrix_directory,
    file.path(config$paths$data_dir, config$input_files$matrix_directory)
  ),
  
  # Track the actual quantification files (matrix.mtx / features.tsv /
  # barcodes.tsv) as a file target so that regenerating them IN PLACE
  # invalidates raw_quant_list. Without this, matrix_directory is only a path
  # string: targets never sees on-disk edits and silently replays a stale
  # import - which reproduces the "identical(rownames(...)) is not TRUE"
  # crash from old cached matrices even after the files have been fixed.
  tar_target(
    quant_files,
    list.files(
      matrix_directory,
      pattern    = "\\.(mtx|tsv)$|quants_mat_(cols|rows)\\.txt$",
      recursive  = TRUE,
      full.names = TRUE
    ),
    format = "file"
  ),
  
  # ===============================================================
  # INPUT PATHS FROM YAML
  # ===============================================================
  
  tar_target(padj_val, config$filters$padj_threshold),
  tar_target(LFC_val, config$filters$lfc_threshold),
  tar_target(heatmap_top_n, config$selection$heatmaps$top_n),
  tar_target(heatmap_aggregate_by, config$selection$heatmaps$aggregate_by),
  tar_target(
    run_id,
    paste0(
      "padj_", gsub("\\.", "", as.character(padj_val)),
      "_LFC_", gsub("\\.", "", as.character(LFC_val))
    )
  ),
  # ===============================================================
  # OUTPUT DIRECTORIES
  # ===============================================================
  
  tar_target(results_dir, file.path(config$paths$output_dir, "results", run_id)),
  tar_target(path_VP, file.path(results_dir, "Volcano_plots")),
  tar_target(path_UPSET, file.path(results_dir, "upset")),
  tar_target(path_DE, file.path(results_dir, "DE_genes")),
  tar_target(path_UMAP, file.path(results_dir, "UMAP_plots")),
  tar_target(path_HTML, file.path(results_dir, "DE_genes/HTML_reports")),
  tar_target(path_heatmaps_DE, file.path(results_dir, "Heatmaps_DE")),
  tar_target(path_KNEE, file.path(results_dir, "knee_plots")),
  tar_target(path_HEATMAP, file.path(results_dir, "heatmaps")),
  tar_target(path_HEATMAP_DE_COUNT, file.path(results_dir, "de_count_heatmaps")),
  tar_target(path_DMSO, file.path(results_dir, "dmso_wells")),
  tar_target(path_CORR, file.path(results_dir, "correlation_plots")),
  tar_target(
    create_dirs,
    {
      dirs <- c(
        path_VP, path_UPSET, path_DE, path_UMAP,
        path_heatmaps_DE, path_CORR, path_KNEE,
        path_HEATMAP, path_DMSO, path_HTML, path_HEATMAP_DE_COUNT
      )
      purrr::walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
      dirs
    }
  ),
  
  # ===============================================================
  # METADATA
  # ===============================================================
  
  tar_target(meta_experiment, load_experiment_meta(config)),
  tar_target(
    plate_layout_raw,
    read_plate_layout(
      filename      = meta_file,
      sheet         = 1,
      required_cols = c("row", "col", "plate_well", "well_index")
    )
  ),
  tar_target(plate_layout, process_plate_metadata(plate_layout_raw)),
  
  # ===============================================================
  # QUANTIFICATION IMPORT
  # ===============================================================
  
  tar_target(
    raw_quant_list,
    {
      quant_files # dependency only: forces re-import when any quant file changes
      read_in_quant_files(
        dir          = matrix_directory,
        tool         = config$analysis$quant_tool,
        plate_layout = plate_layout
      )
    }
  ),
  tar_target(combined_counts, combine_metadata_counts(raw_quant_list, plate_layout), deployment = "main", memory = "transient"),
  # ---------------------------------------------------------------
  # WELL-LEVEL QC - drop failed / near-empty wells before DE.
  # `min_umis` is always enforced; `min_frac_median` (relative-to-plate-median)
  # is the robust lever for near-dead cytotoxic wells. Method "fixed" (default)
  # uses only those floors. "inflection"/"knee" additionally raise the cut to
  # the DropletUtils::barcodeRanks threshold (the curve the knee plots draw);
  # that threshold was designed for droplet data with thousands of empty
  # barcodes, so on a plate where every well holds cells it can land inside the
  # real-well distribution - opt in only after inspecting the knee plots.
  #
  # A config WITHOUT a qc$well block is treated as QC disabled, so configs
  # written before well QC existed still run unchanged. The published U2OS
  # benchmark config enables well QC explicitly (qc: well: enabled: yes).
  # ---------------------------------------------------------------
  tar_target(
    well_qc,
    compute_well_qc(
      combined_counts,
      min_umis        = config$qc$well$min_umis %||% 1000,
      min_genes       = config$qc$well$min_genes %||% 200,
      min_frac_median = config$qc$well$min_frac_median %||% 0,
      method          = config$qc$well$method %||% "fixed",
      out_csv         = file.path(path_KNEE, "well_qc.csv")
    ),
    deployment = "main"
  ),
  tar_target(
    plate_layout_qc,
    if (!isTRUE(config$qc$well$enabled %||% FALSE)) {
      plate_layout
    } else {
      apply_well_qc_to_layout(plate_layout, well_qc)
    },
    deployment = "main"
  ),
  tar_target(
    all_knee_plots,
    make_all_knee_plots(combined_counts, well_qc, path_KNEE),
    format = "file"
  ),
  tar_target(
    all_heatmaps_plots,
    make_all_heatmaps(
      combined_counts, plate_layout, meta_experiment, path_HEATMAP,
      plate_format = config$analysis$plate_format
    ),
    format = "file"
  ),
  # experiment_object is built from the UNFILTERED plate_layout on purpose: it
  # feeds the QC views (plate heatmaps, sample UMAP/PCA) and the DE compound
  # heatmaps, where seeing a QC-dropped well is informative. DE and DMSO
  # selection receive plate_layout_qc instead, so dropped wells never enter a
  # contrast. If a QC-dropped well appears in a DE heatmap, that is why.
  tar_target(experiment_object, make_experiment_list(plate_layout, combined_counts), deployment = "main", memory = "transient"),
  
  # ===============================================================
  # GENE ANNOTATION (CACHED)
  # ===============================================================
  
  tar_target(
    all_gene_ids,
    unique(unlist(lapply(experiment_object, function(batch) {
      unique(unlist(lapply(batch, function(plate) plate$UMI.counts$gene)))
    })))
  ),
  tar_target(
    gene_annotation_file,
    create_gene_annotation_file(
      gene_ids = all_gene_ids,
      species  = config$analysis$biomart_dataset,
      host     = config$analysis$biomart_host,
      outfile  = file.path(config$paths$output_dir, "annotations", "gene_annotation.rds")
    ),
    format = "file",
    deployment = "main" # run in main process so mirror-attempt messages are
    # visible live; also correct since this is network I/O,
    # not CPU work that benefits from a worker subprocess
  ),
  tar_target(
    gene_annotation,
    {
      # The cached annotation file may be either:
      #   * a plain data.frame (legacy format, before version-tracking), or
      #   * a list bundle with $annotation + provenance metadata (new format
      #     produced by build_gene_annotation in 02_input_functions.R).
      # Handle both transparently so old caches keep working and downstream
      # consumers always receive a data.frame.
      bundle <- readRDS(gene_annotation_file)
      if (is.data.frame(bundle)) {
        bundle # legacy format
      } else if (is.list(bundle) && is.data.frame(bundle$annotation)) {
        bundle$annotation # new bundle format
      } else {
        stop(
          "Unrecognised gene annotation file format. ",
          "Expected a data.frame (legacy) or list with $annotation ",
          "(new bundle format). Got: ", class(bundle)[1]
        )
      }
    }
  ),
  
  # ===============================================================
  # DMSO SELECTION
  # ===============================================================
  
  tar_target(
    best_dmso,
    {
      create_dirs # force dependency edge; value is unused
      find_best_DMSO_by_z_score(
        experiment_object,
        plate_layout_qc,
        config$analysis$control_compound,
        z_thresh        = config$analysis$dmso$z_thresh,
        per_plate_top_n = config$analysis$dmso$top_n %||% "all",
        min_keep_per_plate = config$analysis$dmso$min_keep %||% 6,
        npc             = config$analysis$dmso$npc,
        pca_weight      = config$analysis$dmso$pca_weight %||% 1.0,
        out_data_dir    = path_DMSO
      )
    },
    deployment = "main"
  ),
  # ===============================================================
  # DIFFERENTIAL EXPRESSION
  # ===============================================================
  
  # Two DE passes, each a single target run in the main process (they
  # serialise experiment_object, which OOMs crew workers):
  #   de_results          per-plate fits (limma-voom under method "auto");
  #                       QC / plate-consistency view
  #   de_results_collapse batches pooled with duplicateCorrelation across
  #                       plates (single-plate batches fit within-plate);
  #                       the PRIMARY result every downstream target consumes
  #
  # min.count controls the presence-aware gene filter (edgeR::filterByExpr).
  # Two deliberate settings exist:
  #   * 5 (default) — standard analyses (e.g. GDPx4, production runs).
  #   * 1 — the U2OS cross-pipeline benchmark ONLY, set via config
  #     (analysis: de: min_count: 1). Permissive so ADEPT's tested universe
  #     subsumes each comparator's and comparator-unique genes keep an AveExpr;
  #     at min.count = 5 ~93% of comparator-unique genes are censored from the
  #     expression comparison, biasing the primary test (Supplementary Methods).
  # Change either only with a re-run of the corresponding validation.
  #
  # method (analysis: de: method) = "auto" (default), "trend" or "voom".
  # auto = voom for single-plate fits, trend for batches pooled across plates
  # with duplicateCorrelation; chosen by the DMSO-vs-DMSO nulls of each path
  # (GDPx4 single-plate: trend 40% vs voom 30% FWER; Novartis cross-plate:
  # trend 4.2% vs voom 56.7%). "trend" reproduces runs before Oct 2026.
  
  tar_target(
    de_results,
    DE_limma_zQC(
      Exp = experiment_object,
      CTRL = config$analysis$control_compound,
      meta = plate_layout_qc,
      best_DMSO_wells = best_dmso$best_DMSO,
      out_dir = path_DE,
      annotated_genes = gene_annotation,
      html_reports = TRUE,
      html_dir = path_HTML,
      P_val = padj_val,
      LFC_value = LFC_val,
      max_html_rows = 2000,
      min.count = config$analysis$de$min_count %||% 5,
      method = config$analysis$de$method %||% "auto"
    ),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    de_results_collapse,
    DE_limma_zQC_collapse_dupCor(
      Exp = experiment_object,
      CTRL = config$analysis$control_compound,
      meta = plate_layout_qc,
      best_DMSO_wells = best_dmso$best_DMSO,
      out_dir = path_DE,
      annotated_genes = gene_annotation,
      min.count = config$analysis$de$min_count %||% 5,
      method = config$analysis$de$method %||% "auto"
    ),
    deployment = "main"
  ),
  
  # ===============================================================
  # SELECTION LOGIC FROM YAML
  # ===============================================================
  
  tar_target(
    selected_treatments,
    {
      trt <- config$selection$treatments
      if (is.null(trt)) unique(plate_layout$treatment) else trt
    }
  ),
  tar_target(
    selected_batches,
    {
      b <- config$selection$batches
      if (is.null(b)) sort(unique(plate_layout$batch)) else b
    }
  ),
  tar_target(
    selected_doses,
    {
      d <- config$selection$doses
      # same type as the treatment_dose column make_lt() filters on, so the
      # fallback is exactly "all doses" (previously read a nonexistent column)
      if (is.null(d)) unique(plate_layout$treatment_dose) else d
    }
  ),
  tar_target(
    selected_contrasts,
    {
      ctr <- config$selection$contrasts
      if (is.null(ctr)) unique(de_results_collapse$Sample) else ctr
    }
  ),
  
  # ===============================================================
  # VOLCANO PLOTS
  # ===============================================================
  
  # Generate volcano for each selected treatment
  tar_target(
    volcano_plots,
    plot_volcano(
      de_result       = de_results_collapse,
      contrast        = selected_contrasts,
      output_dir      = path_VP,
      PADJ_THRESHOLD  = padj_val,
      LFC_THRESHOLD   = LFC_val,
      meta_experiment = meta_experiment
    ),
    pattern = map(selected_contrasts),
    format = "file"
  ),
  tar_target(
    volcano_plots_multiple,
    volcano_plot_static(
      this_contrast_toptabel = de_results_collapse,
      this_arms              = selected_treatments,
      batch_select           = selected_batches,
      p_val                  = padj_val,
      lfc_val                = LFC_val,
      save_dir               = path_VP,
      file_prefix            = "volcano_DE"
    ),
    pattern = map(selected_treatments),
    format = "file"
  ),
  
  
  # ===============================================================
  # HEATMAPS
  # ===============================================================
  
  tar_target(
    heatmaps_DE,
    plot_drugseq_compound_heatmap(
      compound_name     = selected_treatments,
      de_results        = de_results_collapse,
      experiment_object = experiment_object,
      top_n             = heatmap_top_n,
      out_dir           = path_heatmaps_DE,
      aggregate_by      = heatmap_aggregate_by,
      meta_experiment   = meta_experiment
    ),
    pattern = map(selected_treatments),
    format = "file"
  ),
  
  # DE-gene-count plate heatmaps: one plate grid per plate, each well coloured
  # by how many DE genes its compound x dose contrast produced. Uses the
  # collapsed result and the full (unfiltered) layout so every well is drawn.
  tar_target(
    de_count_heatmaps,
    make_all_de_count_heatmaps(
      de_results      = de_results_collapse,
      plate_layout    = plate_layout,
      meta_experiment = meta_experiment,
      outdir          = path_HEATMAP_DE_COUNT,
      plate_format    = config$analysis$plate_format,
      lfc_cutoff      = LFC_val,
      padj_cutoff     = padj_val,
      CTRL            = config$analysis$control_compound
    ),
    format = "file"
  ),
  
  # ===============================================================
  # UPSET PLOTS
  # ===============================================================
  
  tar_target(
    lt_sets,
    make_lt(
      de_results_collapse,
      logFC_cutoff = LFC_val,
      padj_cutoff  = padj_val,
      treatments   = selected_treatments,
      doses        = selected_doses
    )
  ),
  tar_target(
    upset_title,
    build_upset_title(
      treatments  = selected_treatments,
      padj_cutoff = padj_val,
      lfc_cutoff  = LFC_val
    )
  ),
  tar_target(
    upset_plot,
    save_upset_plot(
      lt_list         = lt_sets,
      output_dir      = path_UPSET,
      title           = upset_title,
      filename        = paste0(paste(selected_treatments, collapse = "_vs_"), "_upset.svg"),
      meta_experiment = meta_experiment
    ),
    format = "file"
  ),
  
  # ===============================================================
  # CORRELATION PLOTS
  # ===============================================================
  
  tar_target(
    ggpairs_plots_multiple,
    {
      this_arms <- de_results_collapse %>%
        dplyr::filter(treatment == selected_treatments, batch %in% selected_batches) %>%
        dplyr::arrange(treatment_dose) %>%
        dplyr::distinct(Sample) %>%
        dplyr::pull(Sample)
      
      ggpairs_custom_tidy(
        toptable_meta   = de_results_collapse,
        this_arms       = this_arms,
        padj_threshold  = padj_val,
        logFC_threshold = LFC_val,
        corr_size       = 3,
        batch_s         = selected_batches,
        strip_size      = 8,
        save_dir        = path_CORR,
        meta_experiment = meta_experiment,
        file_prefix     = "pairs",
        treatment_name  = selected_treatments
      )
    },
    pattern = map(selected_treatments),
    format = "file"
  ),
  
  # ===============================================================
  # 6. Global UMAP + Leiden clustering
  #
  # Gene filtering uses filter_and_normalize() (filterByExpr with a
  # group-aware min.count = 3); see the note at umap_dge below.
  #
  # Flow:
  #   umap_UMI_all          raw genes × wells matrix (all plates)
  #   umap_meta             metadata aligned to raw matrix columns
  #   umap_meta_annotated   + plate spatial position annotation
  #   umap_dge              filter_and_normalize → list(logCPM, dge)
  #                         wells dropped here if colSums == 0
  #   umap_meta_filtered    metadata re-aligned to SURVIVING wells
  #                         (must happen AFTER filtering, before batch correction)
  #   umap_logCPM           genes × surviving-wells logCPM
  #   umap_logCPM_bc        batch-corrected logCPM
  #   umap_pca / layout / clusters / df  → plots
  # ===============================================================
  
  tar_target(umap_n_pca, config$umap$n_pca),
  tar_target(umap_n_neighbors, config$umap$n_neighbors),
  tar_target(leiden_k, config$umap$leiden_k),
  tar_target(leiden_resolution, config$umap$leiden_resolution %||% 0.5),
  tar_target(remove_batch_umap, config$umap$remove_batch),
  
  # -- Raw UMI matrix across all plates ---------------------------------------
  tar_target(
    umap_UMI_all,
    build_global_umi_matrix(experiment_object),
    memory = "transient"
  ),
  
  # -- Align metadata to raw matrix columns -----------------------------------
  tar_target(
    umap_meta,
    align_meta_to_umi(umap_UMI_all, plate_layout)
  ),
  
  # -- Annotate well spatial position -----------------------------------------
  tar_target(
    umap_meta_annotated,
    annotate_plate_position(umap_meta)
  ),
  
  # -- Filter genes + normalise ------------------------------------------------
  # UMAP uses filter_and_normalize(): a single filterByExpr pass with
  # min.count = 3, grouped by treatment. UMAP has no multiple-testing penalty,
  # so it keeps every gene with adequate expression in at least one treatment
  # group. The DE path uses the stricter filter_and_normalize_hybrid
  # (filterByExpr, min.count = 5, then per-contrast re-filtering). There is
  # no quantile-based gene filter anywhere in the pipeline.
  tar_target(
    umap_dge,
    filter_and_normalize(
      counts    = umap_UMI_all,
      groups    = umap_meta_annotated$treatment,
      min.count = 3
    ),
    deployment = "main", memory = "transient"
  ),
  tar_target(
    umap_logCPM,
    umap_dge$logCPM
  ),
  
  # -- Re-align metadata to surviving wells after filtering -------------------
  # filter_and_normalize may drop empty wells (colSums == 0).
  # umap_meta_filtered must match umap_logCPM columns exactly
  # before batch correction and all downstream steps.
  tar_target(
    umap_meta_filtered,
    umap_meta_annotated[
      match(colnames(umap_logCPM), umap_meta_annotated$plate_id), ,
      drop = FALSE
    ]
  ),
  
  # -- Batch correction -------------------------------------------------------
  tar_target(
    umap_logCPM_bc,
    remove_umap_batch_effect(
      umap_logCPM,
      meta         = umap_meta_filtered,
      CTRL         = config$analysis$control_compound,
      batch_var    = "plate_barcode",
      remove_batch = remove_batch_umap
    ),
    deployment = "main", memory = "transient"
  ),
  
  # -- PCA → Leiden → UMAP layout --------------------------------------------
  tar_target(umap_pca, run_umap_pca(umap_logCPM_bc), deployment = "main", memory = "transient"),
  tar_target(
    leiden_clusters,
    run_leiden_from_pca(umap_pca$x,
                        n_pc = umap_n_pca, k = leiden_k,
                        resolution = leiden_resolution
    )
  ),
  tar_target(
    umap_layout,
    run_umap_from_pca(umap_pca$x, n_pca = umap_n_pca, n_neighbors = umap_n_neighbors)
  ),
  
  # -- Assemble plot table (layout + filtered meta + clusters) ----------------
  tar_target(
    umap_df,
    assemble_umap_table(umap_layout, umap_meta_filtered, leiden_clusters)
  ),
  
  # -- Assemble PCA plot table (PC1/PC2 + filtered meta + clusters) -----------
  # Mirrors umap_df but uses the first two principal components directly.
  # Useful as a sanity check on the UMAP layout — if treatment clusters are
  # visible in the raw PCA, the UMAP is amplifying real structure rather
  # than inventing it.
  tar_target(
    pca_df,
    assemble_pca_table(umap_pca, umap_meta_filtered, leiden_clusters)
  ),
  
  
  # -- Plots ------------------------------------------------------------------
  tar_target(
    umap_plot_treatment_html,
    {
      p <- plotly::plot_ly(
        umap_df,
        x = ~UMAP1, y = ~UMAP2, color = ~treatment,
        type = "scatter", mode = "markers",
        text = ~ paste(
          "Plate:", plate_id, "<br>Treatment:", treatment,
          "<br>Dose:", treatment_dose, "<br>Dose_unit:", treatment_dose_unit,
          "<br>Timepoint:", timepoint, "<br>Batch:", batch
        ),
        hoverinfo = "text"
      ) %>% plotly::layout(title = "UMAP colored by treatment")
      
      save_plotly_html(
        plot = p, filename = "umap_by_treatment_batch.html",
        out_dir = path_UMAP, meta_experiment = meta_experiment
      )
    },
    format = "file"
  ),
  tar_target(
    umap_plot_dmso_qc_html,
    {
      df_dmso <- umap_df %>%
        dplyr::filter(treatment == config$analysis$control_compound)
      
      p <- plotly::plot_ly(
        df_dmso,
        x = ~UMAP1, y = ~UMAP2, color = ~well_category,
        type = "scatter", mode = "markers",
        text = ~ paste(
          "Plate:", plate_id, "<br>Batch:", batch,
          "<br>Well category:", well_category
        ),
        hoverinfo = "text"
      ) %>% plotly::layout(title = "UMAP — DMSO wells (plate QC)")
      
      save_plotly_html(
        plot = p, filename = "umap_dmso_plate_qc_batch.html",
        out_dir = path_UMAP, meta_experiment = meta_experiment
      )
    },
    format = "file"
  ),
  tar_target(
    umap_plot_leiden_html,
    {
      p <- plotly::plot_ly(
        umap_df,
        x = ~UMAP1, y = ~UMAP2, color = ~cluster_leiden,
        type = "scatter", mode = "markers",
        text = ~ paste(
          "Plate:", plate_id, "<br>Leiden cluster:", cluster_leiden,
          "<br>Treatment:", treatment, "<br>Dose:", treatment_dose,
          "<br>Timepoint:", timepoint, "<br>Batch:", batch
        ),
        hoverinfo = "text"
      ) %>% plotly::layout(title = "UMAP colored by Leiden cluster")
      
      save_plotly_html(
        plot = p, filename = "umap_by_leiden_batch.html",
        out_dir = path_UMAP, meta_experiment = meta_experiment
      )
    },
    format = "file"
  ),
  
  # -- PCA plots (parallel to the UMAP plots above) ---------------------------
  # These use the same input table layout as the UMAP plots, but with PC1/PC2
  # axes labeled by their variance-explained percentages. Useful as a sanity
  # check that the UMAP is amplifying real structure rather than inventing it.
  tar_target(
    pca_plot_treatment_html,
    {
      pc1_pct <- attr(pca_df, "pc1_var_pct")
      pc2_pct <- attr(pca_df, "pc2_var_pct")
      p <- plotly::plot_ly(
        pca_df,
        x = ~PC1, y = ~PC2, color = ~treatment,
        type = "scatter", mode = "markers",
        text = ~ paste(
          "Plate:", plate_id, "<br>Treatment:", treatment,
          "<br>Dose:", treatment_dose, "<br>Dose_unit:", treatment_dose_unit,
          "<br>Timepoint:", timepoint, "<br>Batch:", batch
        ),
        hoverinfo = "text"
      ) %>% plotly::layout(
        title = "PCA colored by treatment (all wells)",
        xaxis = list(title = sprintf("PC1 (%.1f%%)", pc1_pct)),
        yaxis = list(title = sprintf("PC2 (%.1f%%)", pc2_pct))
      )
      
      save_plotly_html(
        plot = p, filename = "pca_by_treatment_batch.html",
        out_dir = path_UMAP, meta_experiment = meta_experiment
      )
    },
    format = "file"
  ),
  tar_target(
    pca_plot_dmso_qc_html,
    {
      df_dmso <- pca_df %>%
        dplyr::filter(treatment == config$analysis$control_compound)
      pc1_pct <- attr(pca_df, "pc1_var_pct")
      pc2_pct <- attr(pca_df, "pc2_var_pct")
      
      p <- plotly::plot_ly(
        df_dmso,
        x = ~PC1, y = ~PC2, color = ~well_category,
        type = "scatter", mode = "markers",
        text = ~ paste(
          "Plate:", plate_id, "<br>Batch:", batch,
          "<br>Well category:", well_category
        ),
        hoverinfo = "text"
      ) %>% plotly::layout(
        title = "PCA - DMSO wells (plate QC)",
        xaxis = list(title = sprintf("PC1 (%.1f%%)", pc1_pct)),
        yaxis = list(title = sprintf("PC2 (%.1f%%)", pc2_pct))
      )
      
      save_plotly_html(
        plot = p, filename = "pca_dmso_plate_qc_batch.html",
        out_dir = path_UMAP, meta_experiment = meta_experiment
      )
    },
    format = "file"
  ),
  
  
  # ===============================================================
  # 7. Gene-level UMAP (DE gene programs)
  # ===============================================================
  
  tar_target(gene_umap_lfc_cutoff, LFC_val),
  tar_target(gene_umap_padj_cutoff, padj_val),
  tar_target(gene_umap_n_pca, config$umap$n_pca),
  tar_target(gene_umap_n_neighbors, config$umap$n_neighbors),
  tar_target(gene_leiden_k, config$umap$leiden_k),
  tar_target(gene_leiden_resolution, config$umap$leiden_resolution %||% 0.5),
  tar_target(
    gene_umap_genes,
    define_de_gene_universe(
      de_results_collapse,
      lfc_cutoff  = gene_umap_lfc_cutoff,
      padj_cutoff = gene_umap_padj_cutoff
    )
  ),
  tar_target(
    gene_umap_matrix_raw,
    build_gene_contrast_matrix(
      de_results_collapse,
      genes     = gene_umap_genes,
      value_col = "LFC"
    ),
    deployment = "main",
    memory = "transient"
  ),
  tar_target(gene_umap_matrix, scale_gene_matrix(gene_umap_matrix_raw)),
  tar_target(
    sample_umap_pca,
    prcomp(t(gene_umap_matrix), scale. = FALSE)
  ),
  tar_target(
    sample_leiden_clusters,
    run_leiden_from_pca(
      sample_umap_pca$x,
      n_pc = gene_umap_n_pca, k = gene_leiden_k,
      resolution = gene_leiden_resolution
    )
  ),
  tar_target(
    sample_umap_layout,
    run_umap_from_pca(
      sample_umap_pca$x,
      n_pca       = gene_umap_n_pca,
      n_neighbors = gene_umap_n_neighbors
    )
  ),
  tar_target(
    sample_umap_df,
    {
      df <- sample_umap_layout
      df$Sample <- rownames(sample_umap_pca$x)
      df$cluster_leiden <- factor(sample_leiden_clusters)
      
      dplyr::left_join(
        df,
        dplyr::distinct(
          de_results_collapse,
          Sample, treatment, treatment_dose, treatment_dose_unit, timepoint, batch
        ),
        by = "Sample"
      )
    }
  ),
  tar_target(
    sample_umap_plot_html,
    {
      p <- plotly::plot_ly(
        sample_umap_df,
        x = ~UMAP1, y = ~UMAP2, color = ~treatment, symbol = ~treatment,
        type = "scatter", mode = "markers",
        text = ~ paste(
          "Sample:", Sample, "<br>Treatment:", treatment,
          "<br>Dose:", treatment_dose, "<br>Time:", timepoint,
          "<br>Leiden:", cluster_leiden
        ),
        hoverinfo = "text"
      ) %>% plotly::layout(title = "Sample UMAP (PCA on DE genes only)")
      
      save_plotly_html(
        plot = p, filename = "sample_umap_DE_gene_basis.html",
        out_dir = path_UMAP, meta_experiment = meta_experiment
      )
    },
    format = "file"
  ),
  
  # ===============================================================
  # FINAL REPORT — one HTML index linking every output
  # ===============================================================
  
  tar_target(
    final_report,
    {
      # Referencing these symbols makes them dependencies of final_report,
      # so this target re-runs whenever any upstream output changes.
      .deps <- list(
        volcano_plots, volcano_plots_multiple,
        heatmaps_DE, upset_plot, ggpairs_plots_multiple,
        all_knee_plots, all_heatmaps_plots, de_count_heatmaps,
        umap_plot_treatment_html, umap_plot_dmso_qc_html,
        umap_plot_leiden_html, sample_umap_plot_html,
        pca_plot_treatment_html, pca_plot_dmso_qc_html,
        de_results, de_results_collapse, best_dmso,
        plate_layout
      )
      
      build_final_report(
        results_dir = results_dir,
        config = config,
        run_id = run_id,
        padj = padj_val,
        lfc = LFC_val,
        meta_experiment = meta_experiment,
        plate_layout = plate_layout,
        de_results = de_results_collapse,
        best_dmso = best_dmso,
        paths = list(
          "Volcano plots" = path_VP,
          "UpSet plot" = path_UPSET,
          "DE results" = path_DE,
          "DE heatmaps" = path_heatmaps_DE,
          "UMAP plots" = path_UMAP,
          "Correlation plots" = path_CORR,
          "Knee plots (QC)" = path_KNEE,
          "Plate heatmaps (QC)" = path_HEATMAP,
          "DE-gene count heatmaps" = path_HEATMAP_DE_COUNT,
          "DMSO QC" = list(dir = path_DMSO, list_only = TRUE)
        )
      )
    },
    format = "file"
  )
)