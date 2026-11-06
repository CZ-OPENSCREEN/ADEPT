# =============================================================================
# DE_limma_percontrast.R
#
# Changes:
#   1. filter_and_normalize() → filter_and_normalize_hybrid()
#        Global pre-filter + TMM normalisation, returning the dge object so
#        per-contrast subsetting can reuse it. Default (use_filterByExpr = TRUE)
#        is a presence-aware edgeR::filterByExpr(min.count = 5) grouped by
#        condition; the older top-25%-by-mean-CPM quantile filter is kept
#        behind use_filterByExpr = FALSE for A/B comparison only.
#        standard analyses use 5; 
#        the cross-pipeline benchmark comparison runs min.count = 1 
#
#   2. run_single_contrast()
#        Filters genes specifically for one compound × dose × timepoint vs DMSO,
#        fits a two-group model, and returns topTable. Accepts optional
#        block + correlation arguments for dupCor use.
#
#   3. DE_limma_zQC()
#        Replaces single lmFit-over-all-contrasts with a per-contrast loop.
#        Each compound × dose × timepoint gets its own gene universe.
#
#   4. DE_limma_zQC_collapse_dupCor()
#        corfit estimated once globally (needs full matrix for stable estimate),
#        then passed into per-contrast run_single_contrast calls.
#        Same per-contrast gene filtering as above.
#
#   5. method = "auto" | "trend" | "voom"  (config: analysis: de: method)
#        "auto" (recommended): voom where plates are fitted on their own (the
#        per-plate pass, and single-plate batches in the collapsed pass), trend
#        where a batch spans >= 2 plates and duplicateCorrelation is used.
#        Chosen by the DMSO-vs-DMSO nulls of each path:
#          single-plate (GDPx4, 280 contrasts): trend 40%, voom 30% FWER
#          cross-plate  (Novartis, 120):        trend 4.2%, voom 56.7% FWER
#        (voom's consensus plate correlation was ~0.45 vs ~0.01 on logCPM,
#        which made voom + dupCor anti-conservative.)
#        "trend": limma-trend on the batch logCPM (prior.count 1) -> original.
#        "voom":  limma-voom precision weights on the contrast's counts (same
#                 TMM factors and library sizes as the logCPM), eBayes(robust).
#                 With dupCor the global correlation is estimated on voom
#                 values and passed to the per-contrast voom + lmFit (standard
#                 voom -> duplicateCorrelation -> voom). AveExpr is reported
#                 from the batch logCPM in both modes, so it stays comparable.
#        Why: GDPx4 DMSO-vs-DMSO null (2 fake vs ~6 DMSO wells, 280 per-plate
#        contrasts): trend FWER 33%, voom 20% (macpie, also voom, 16%).
#        Trend's excess false positives were low-count, large-|LFC| genes in
#        low-depth wells. arrayWeights / whole-plate fits made it worse.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SHARED HELPER: global pre-filter + normalization
# -----------------------------------------------------------------------------

#' Global pre-filter and normalization
#'
#' Drops empty wells, TMM-normalises, applies a global gene filter, then
#' re-normalises. Returns both logCPM (for modelling) and the dge object
#' (needed so run_single_contrast can call filterByExpr on per-contrast
#' subsets). Every contrast reads logCPM from this ONE batch-level
#' normalisation, so AveExpr is on a common scale across contrasts.
#'
#' @param counts  Raw count matrix (genes × wells)
#' @param groups  Condition label per well (used by filterByExpr)
#' @param use_filterByExpr  TRUE (default in both DE callers): presence-aware
#'                edgeR::filterByExpr grouped by condition. FALSE: legacy
#'                top-`pct` by mean CPM (kept for A/B only).
#' @param min.count  Passed to filterByExpr. Both DE callers pass 5.
#' @param pct     Fraction kept by the legacy quantile filter (default 0.25)
#' @return list(logCPM, dge)
#' 
filter_and_normalize_hybrid <- function(counts, pct = 0.25, groups = NULL,
                                        use_filterByExpr = FALSE, min.count = 5) {
  
  dge <- edgeR::DGEList(counts = counts)
  
  # Drop empty wells, keeping `groups` aligned to the surviving columns.
  keep_cols <- colSums(dge$counts) > 0
  dge <- dge[, keep_cols, keep.lib.sizes = FALSE]
  if (!is.null(groups)) groups <- groups[keep_cols]
  dge <- edgeR::calcNormFactors(dge)
  
  if (use_filterByExpr) {
    # Presence-aware global filter (recommended). filterByExpr keeps a gene if
    # it clears a count threshold in at least the smallest group's worth of
    # wells, so a gene induced only at the top dose(s) survives because it is
    # present in those wells. Unlike the mean-CPM quantile below, this is NOT
    # biased against condition-specific (e.g. high-dose-induced) genes, which
    # is exactly the class of real low-expression DE genes the quantile drops.
    keep <- edgeR::filterByExpr(dge, group = groups, min.count = min.count)
    message("      Global filter (filterByExpr, min.count = ", min.count,
            "): keeping ", sum(keep), " / ", length(keep), " genes")
  } else {
    # A gene that is
    # only expressed at high dose has low mean CPM by construction and is
    # dropped here before any contrast is fit. Kept for backward-compatible A/B.
    mean_cpm <- rowMeans(edgeR::cpm(dge, log = FALSE))
    keep     <- mean_cpm >= quantile(mean_cpm, 1 - pct)
    message("      Global filter (top ", round(pct * 100),
            "% mean CPM): keeping ", sum(keep), " / ", length(keep), " genes")
  }
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)   # re-normalize after gene removal
  
  list(
    logCPM = edgeR::cpm(dge, log = TRUE, prior.count = 1),
    dge    = dge
  )
}


# -----------------------------------------------------------------------------
# 2. per-contrast filter + fit
# -----------------------------------------------------------------------------

# Design note: variance estimation
#
# Each contrast is fitted independently (treatment + DMSO wells only).
# This means gene-wise variance is estimated from ~6 samples rather than
# the full plate. Trade-off:
#   Pro:  Variance is specific to this contrast. Heterogeneous treatments
#         don't contaminate each other's variance estimates.
#   Con:  Fewer residual d.f. for variance estimation.
#   Mitigation: eBayes(robust=TRUE) shrinks gene-wise variances toward
#         a global trend, borrowing information across genes.
#
# The collapsed (dupCor) path pools all plates within a batch, giving
# more d.f. and explicit plate-effect modelling. Use collapsed results
# as the primary analysis; per-plate results for QC/diagnostics.

#' Filter genes and fit model for one contrast (compound × dose × timepoint)
#'
#' @param logCPM_global  Full logCPM matrix from filter_and_normalize_hybrid
#' @param dge_global     Full dge object from filter_and_normalize_hybrid
#' @param treat_wells    Character vector of well IDs for the treatment group
#' @param dmso_wells     Character vector of well IDs for matching DMSO
#' @param treat_label    sample_clean label for treatment (e.g. "DrugA_1uM_4hr")
#' @param ctrl_label     sample_clean label for DMSO (e.g. "DMSO_0uM_4hr")
#' @param well_to_group  Named character vector: well_id -> sample_clean
#' @param min.count      Passed to filterByExpr
#' @param block          Named vector: well_id -> plate_barcode (dupCor only)
#' @param correlation    corfit$consensus from duplicateCorrelation (dupCor only)
#' @param min_resid_df   Minimum residual df (n_wells - 2) required to fit the
#'                       contrast. With residual df = 1 eBayes inference is
#'                       dominated by the prior; with 0 it fails outright.
#'                       NOTE: this is a TOTAL-df check; 1 treatment + 3 DMSO
#'                       wells gives df = 2 and passes it. Per-group
#'                       replication is enforced separately by min_group_n.
#' @param method         "trend" (limma-trend on logCPM_global) or "voom"
#'                       (voom weights on the contrast counts). See header.
#' @param min_group_n    Minimum wells required in EACH of the two groups
#'                       (default 2). A contrast with a single treatment well
#'                       (e.g. one replicate lost to well QC) is skipped with
#'                       reason "insufficient_group_replication".
#' @return list with two slots:
#'           $results   data.frame of topTable results, or NULL if skipped
#'           $skip_info NULL if the contrast ran, otherwise a 1-row data.frame
#'                       (treat_label, ctrl_label, n_treat, n_dmso,
#'                        resid_df, n_genes_kept, reason)
run_single_contrast <- function(logCPM_global,
                                dge_global,
                                treat_wells,
                                dmso_wells,
                                treat_label,
                                ctrl_label,
                                well_to_group,
                                min.count    = 1,
                                min_resid_df = 2,
                                min_group_n  = 2,
                                block        = NULL,
                                correlation  = NULL,
                                method       = c("trend", "voom")) {
  
  method         <- match.arg(method)
  contrast_wells <- c(treat_wells, dmso_wells)
  groups_sub     <- well_to_group[contrast_wells]
  
  # Helper: build a structured skip record so callers can write a skip log
  make_skip <- function(reason, resid_df = NA_integer_, n_genes_kept = NA_integer_) {
    data.frame(
      treat_label  = treat_label,
      ctrl_label   = ctrl_label,
      n_treat      = length(treat_wells),
      n_dmso       = length(dmso_wells),
      resid_df     = resid_df,
      n_genes_kept = n_genes_kept,
      reason       = reason,
      stringsAsFactors = FALSE
    )
  }
  
  # ---- Replication guards --------------------------------------------------
  # (a) Each group needs >= min_group_n wells. Without this, a contrast with
  #     1 treatment + 3 DMSO wells passes the residual-df check below and is
  #     fitted against a single unreplicated well.
  n_min_group <- min(length(treat_wells), length(dmso_wells))
  if (n_min_group < min_group_n) {
    message("      Skipping [", treat_label,
            "] — ", length(treat_wells), " treat + ",
            length(dmso_wells), " DMSO wells (need >= ", min_group_n,
            " in each group)")
    return(list(
      results   = NULL,
      skip_info = make_skip("insufficient_group_replication",
                            resid_df = length(contrast_wells) - 2L)
    ))
  }
  
  # (b) 2-group design uses 2 df, so residual df = n_wells - 2.
  #     eBayes needs >= 1 residual df to estimate gene-wise variances at all.
  resid_df <- length(contrast_wells) - length(unique(groups_sub))
  if (resid_df < min_resid_df) {
    message("      Skipping [", treat_label,
            "] — ", length(treat_wells), " treat + ",
            length(dmso_wells), " DMSO wells (residual df = ",
            resid_df, ", need >= ", min_resid_df, ")")
    return(list(
      results   = NULL,
      skip_info = make_skip("insufficient_replication", resid_df = resid_df)
    ))
  }
  
  # ---- Per-contrast gene filter --------------------------------------------
  # filterByExpr sees only the wells relevant to this contrast,
  # so a gene expressed only in other groups cannot inflate the keep set.
  dge_sub <- dge_global[, contrast_wells, keep.lib.sizes = FALSE]
  keep    <- edgeR::filterByExpr(dge_sub, group = groups_sub,
                                 min.count = min.count)
  
  if (sum(keep) < 10) {
    message("      Only ", sum(keep), " genes pass filter for [",
            treat_label, "] — skipping")
    return(list(
      results   = NULL,
      skip_info = make_skip("too_few_genes_after_filter",
                            resid_df = resid_df, n_genes_kept = sum(keep))
    ))
  }
  
  logCPM_sub <- logCPM_global[keep, contrast_wells, drop = FALSE]
  
  # ---- Two-group design ----------------------------------------------------
  design           <- model.matrix(~ 0 + groups_sub)
  colnames(design) <- sub("^groups_sub", "", colnames(design))
  
  if (!all(c(treat_label, ctrl_label) %in% colnames(design))) {
    message("      Design levels missing for [", treat_label, "] — skipping")
    return(list(
      results   = NULL,
      skip_info = make_skip("design_levels_missing",
                            resid_df = resid_df, n_genes_kept = sum(keep))
    ))
  }
  
  # Contrast vector: treatment - DMSO
  contr_vec <- ifelse(colnames(design) == treat_label,  1L,
                      ifelse(colnames(design) == ctrl_label,  -1L, 0L))
  contr_mat <- matrix(contr_vec, ncol = 1,
                      dimnames = list(colnames(design), treat_label))
  
  # ---- Fit (trend or voom; plain or dupCor) --------------------------------
  # dupCor path: block and consensus correlation come from the global
  # duplicateCorrelation() call made once before the contrast loop
  use_block <- !is.null(block) && !is.null(correlation)
  blk       <- if (use_block) block[contrast_wells] else NULL
  
  if (method == "voom") {
    # Same genes, TMM factors and library sizes as logCPM_sub
    y_sub <- dge_sub[keep, , keep.lib.sizes = TRUE]
    obj   <- if (use_block) {
      limma::voom(y_sub, design, block = blk, correlation = correlation)
    } else {
      limma::voom(y_sub, design)
    }
  } else {
    obj <- logCPM_sub
  }
  
  fit <- if (use_block) {
    limma::lmFit(obj, design, block = blk, correlation = correlation)
  } else {
    limma::lmFit(obj, design)
  }
  
  fit2 <- limma::contrasts.fit(fit, contr_mat)
  # robust = TRUE guards against the occasional near-zero-variance gene
  # that survives filtering when replicate count is very low.
  # voom carries the mean-variance relationship in its weights -> trend = FALSE
  fit3 <- limma::eBayes(fit2, trend = (method == "trend"), robust = TRUE)
  # AveExpr from the batch logCPM in both modes (comparable across contrasts
  # and with earlier trend runs); voom's own E uses prior.count 0.5
  fit3$Amean <- rowMeans(logCPM_sub)
  fit3$genes <- data.frame(GeneId = rownames(logCPM_sub))
  
  tt <- limma::topTable(fit3, coef = 1, number = Inf, sort.by = "none") %>%
    dplyr::rename(LFC = logFC, padj = adj.P.Val, pval = P.Value) %>%
    dplyr::mutate(
      sample_clean   = treat_label,
      n_genes_tested = sum(keep)     # useful QC: track gene universe per contrast
    )
  
  list(results = tt, skip_info = NULL)
  
}

# -----------------------------------------------------------------------------
# 2b. PER-COMPOUND DOSE-TREND TEST
# -----------------------------------------------------------------------------

#' Test each gene for a monotonic log-dose response across a compound's whole
#' dose series (complementary to the per-dose contrasts).
#'
#' Rationale: a lowly-expressed gene that is only induced/repressed at the top
#' doses may not clear significance in any single per-dose contrast, yet shows a
#' clear monotonic trend across the series. Modelling log10(dose) as a continuous
#' covariate recovers these genes while still rejecting dose-independent noise
#' (which has no trend). DMSO wells anchor the low-dose end of the axis.
#'
#' @param logCPM_global  Full logCPM matrix from filter_and_normalize_hybrid
#' @param dge_global     Full dge object from filter_and_normalize_hybrid
#' @param compound_wells Named numeric vector: well_id -> dose (uM), treated wells
#' @param dmso_wells     Character vector of matching DMSO well IDs
#' @param treat_label    Label for this compound x timepoint
#' @param min.count      Passed to filterByExpr
#' @param min_doses      Minimum distinct doses required to fit a trend
#' @param block,correlation  dupCor block + consensus correlation (optional)
#' @return list(results = topTable on the log-dose slope, or NULL; skip_info)
run_dose_trend <- function(logCPM_global,
                           dge_global,
                           compound_wells,
                           dmso_wells,
                           treat_label,
                           min.count   = 1,
                           min_doses   = 4,
                           block       = NULL,
                           correlation = NULL) {
  
  treated_ids <- names(compound_wells)
  all_wells   <- c(treated_ids, dmso_wells)
  doses_treat <- as.numeric(compound_wells)
  
  # Guard: keep only wells present in the supplied (possibly permissive) matrix,
  # so a well dropped during normalization can't crash the fit.
  avail          <- colnames(logCPM_global)
  compound_wells <- compound_wells[names(compound_wells) %in% avail]
  dmso_wells     <- dmso_wells[dmso_wells %in% avail]
  treated_ids    <- names(compound_wells)
  all_wells      <- c(treated_ids, dmso_wells)
  doses_treat    <- as.numeric(compound_wells)
  
  make_skip <- function(reason, n = NA_integer_) {
    data.frame(treat_label = treat_label, reason = reason,
               n_doses = length(unique(doses_treat)), n_genes = n,
               stringsAsFactors = FALSE)
  }
  
  # ---- Need enough distinct doses to define a trend ------------------------
  if (length(unique(doses_treat)) < min_doses) {
    message("      Skipping dose-trend [", treat_label, "] — only ",
            length(unique(doses_treat)), " distinct dose(s) (need >= ",
            min_doses, ")")
    return(list(results = NULL, skip_info = make_skip("too_few_doses")))
  }
  
  # log10-dose covariate; DMSO anchored one log unit below the lowest dose
  floor_dose <- min(doses_treat) / 10
  log_dose   <- log10(c(doses_treat, rep(floor_dose, length(dmso_wells))))
  
  # ---- Presence-aware per-compound gene filter -----------------------------
  # group = dose level, so a gene detected only at high dose still survives.
  dge_sub <- dge_global[, all_wells, keep.lib.sizes = FALSE]
  keep    <- edgeR::filterByExpr(dge_sub, group = factor(round(log_dose, 3)),
                                 min.count = min.count)
  if (sum(keep) < 10) {
    message("      Only ", sum(keep), " genes pass filter for dose-trend [",
            treat_label, "] — skipping")
    return(list(results = NULL, skip_info = make_skip("too_few_genes", sum(keep))))
  }
  
  logCPM_sub <- logCPM_global[keep, all_wells, drop = FALSE]
  
  # ---- Continuous log-dose design: intercept + slope -----------------------
  design <- model.matrix(~ log_dose)
  
  if (!is.null(block) && !is.null(correlation)) {
    fit <- limma::lmFit(logCPM_sub, design,
                        block = block[all_wells], correlation = correlation)
  } else {
    fit <- limma::lmFit(logCPM_sub, design)
  }
  
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  fit$genes <- data.frame(GeneId = rownames(logCPM_sub))
  
  # coef "log_dose" = dose-response slope (logCPM change per log10 dose unit).
  # A significant slope == monotonic dose response, which is the property that
  # distinguishes real low-count responders from dose-independent noise.
  tt <- limma::topTable(fit, coef = "log_dose", number = Inf, sort.by = "none") %>%
    dplyr::rename(dose_slope = logFC, padj = adj.P.Val, pval = P.Value) %>%
    dplyr::mutate(treat_label = treat_label, n_genes_tested = sum(keep))
  
  list(results = tt, skip_info = NULL)
}

# --------------------------------------------------------------
# FIND BEST DMSO control z-threshold

# DMSO well quality control — robust-z + PCA-aware scoring
#
# Identifies the best DMSO wells per plate using:
#   Part 1 — gene-wise robust z-score outlier counts
#   Part 2 — PCA-score robust z-score outlier counts
#   Part 3 — combined scaled score (lower = better well)
#
# Returns a list:
#   $best_DMSO  data.frame of per-plate top-n wells to KEEP
#   $bad_DMSO   data.frame of wells whose combined score exceeds
#               median + 2*MAD (computed within the selected DMSO dose)

# --------------------------------------------------------------

find_best_DMSO_by_z_score <- function(
    Exp,
    meta,
    CTRL            = "DMSO",
    z_thresh        = 3,
    per_plate_top_n = "all",
    min_keep_per_plate = 6,
    npc             = 10,
    pca_weight      = 1.0,
    scale_pca       = FALSE,
    out_data_dir    = "output"
) {
  # -------------------------------------------------------------------------
  # Argument notes (added in response to reviewer feedback):
  #
  # per_plate_top_n
  #   "all" (default) retains every DMSO well that does not exceed the
  #   bad-well threshold (median + 2xMAD of the combined score), subject to
  #   min_keep_per_plate below. A number N instead keeps the N best-scoring
  #   wells per plate (the published U2OS benchmark used 3). In the
  #   DMSO-vs-DMSO null tests "all" was at least as well calibrated as 3 on
  #   both datasets (GSE176150 cross-plate FWER 4.2% vs 5.0%; GDPx4 per-plate
  #   limma-trend 39.6% vs 44.3%), so it is the default.
  #
  # min_keep_per_plate  (only with per_plate_top_n = "all")
  #   Flagged wells are dropped only while the plate keeps at least this many
  #   DMSO controls; below that, the best-scoring flagged wells are kept. On
  #   plates with few DMSO wells, losing a control costs more than including a
  #   mildly deviant one: GDPx4 per-plate null (8 DMSO/plate, 2 used as the fake
  #   compound) FWER 51% / 39% / 19.5% with 4 / 5 / 6 controls kept. Plates
  #   with many DMSO wells (e.g. 24 on Novartis) never reach the floor.
  #
  # scale_pca
  #   Whether to z-scale each gene before PCA (i.e. prcomp(..., scale. = TRUE)).
  #   FALSE (default, original behaviour) lets gene variance contribute its
  #   biological weight: highly expressed, highly variable genes dominate the
  #   first PCs. TRUE gives each gene unit variance, so the PCA reflects
  #   well-level structure across the full transcriptome rather than being
  #   driven by a few high-variance genes. The scaled version is more
  #   defensible for outlier WELL detection (the QC question), the unscaled
  #   version is more defensible for biology-aware QC. Both are reasonable;
  #   the default is unchanged for backward compatibility.
  # -------------------------------------------------------------------------
  
  stopifnot(
    is.numeric(per_plate_top_n) && per_plate_top_n > 0 ||
      identical(per_plate_top_n, "all")
  )
  
  message("Identifying best DMSO wells using robust-z + PCA-aware scoring")
  message("  scale_pca        = ", scale_pca)
  message("  per_plate_top_n  = ", per_plate_top_n)
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
  
  best_DMSO_all <- list()
  bad_DMSO_all <- list()
  
  batch_keys <- names(Exp)
  
  for (batch_key in batch_keys) {
    message("\n Processing batch: ", batch_key)
    
    # ---- Metadata ----
    meta_sub <- meta %>%
      dplyr::filter(batch == batch_key) %>%
      dplyr::mutate(
        barcode_well = paste0(plate_well, "_", well_index),
        plate_barcode_short = substr(
          plate_barcode,
          nchar(plate_barcode) - 2,
          nchar(plate_barcode)
        )
      )
    stopifnot("plate_id" %in% colnames(meta_sub))
    
    UMI_list <- Exp[[batch_key]]
    
    # ---- Build UMI matrix ----
    # Build column names explicitly as <plate_barcode>_<plate_well>_<well_index>
    # BEFORE cbind, so we don't rely on cbind's automatic "<name>." prefixing
    # and a regex to undo it. The regex approach assumes alphanumeric-only
    # plate barcodes, which silently fails for any barcode containing
    # underscores, hyphens, dots, or lowercase letters. This explicit approach
    # matches what DE_limma_zQC_collapse_dupCor() does later in this file.
    UMI_batch <- lapply(names(UMI_list), function(plate_name) {
      m <- as.data.frame(UMI_list[[plate_name]]$UMI.counts)
      rownames(m) <- m$gene
      m <- m[, !colnames(m) %in% c("gene", "sum"), drop = FALSE]
      m <- m[order(rownames(m)), , drop = FALSE]
      colnames(m) <- paste0(plate_name, "_", colnames(m))
      m
    })
    
    stopifnot(all(sapply(
      UMI_batch,
      function(m) identical(rownames(m), rownames(UMI_batch[[1]]))
    )))
    
    UMI_batch <- do.call(cbind, UMI_batch)
    UMI_batch <- as.matrix(UMI_batch)
    
    valid_ids <- intersect(colnames(UMI_batch), meta_sub$plate_id)
    UMI_batch <- UMI_batch[, valid_ids, drop = FALSE]
    meta_sub <- meta_sub[match(valid_ids, meta_sub$plate_id), , drop = FALSE]
    
    stopifnot(identical(colnames(UMI_batch), meta_sub$plate_id))
    
    # ---- Gene filtering ----
    dge_all  <- edgeR::DGEList(counts = UMI_batch)
    keep_genes <- edgeR::filterByExpr(dge_all, min.count = 5)
    
    UMI_f <- UMI_batch[keep_genes, , drop = FALSE]
    UMI_f <- UMI_f[, colSums(UMI_f) > 0, drop = FALSE]
    
    d <- edgeR::DGEList(counts = UMI_f)
    d <- edgeR::calcNormFactors(d, method = "TMM")
    
    X <- edgeR::cpm(d, log = TRUE, prior.count = 1)
    
    # Re-align meta after potential empty-well removal
    meta_sub <- meta_sub[meta_sub$plate_id %in% colnames(X), , drop = FALSE]
    meta_sub <- meta_sub[match(colnames(X), meta_sub$plate_id), , drop = FALSE]
    
    # ---- DMSO subset ----
    dmso_meta <- meta_sub[meta_sub$treatment == CTRL, , drop = FALSE]
    dmso_ids <- dmso_meta$plate_id
    dmso_X <- X[, dmso_ids, drop = FALSE]
    
    # ============================================================
    # PART 1: Gene-wise robust z-score outliers
    # ============================================================
    
    g_med <- matrixStats::rowMedians(dmso_X, na.rm = TRUE)
    g_mad <- matrixStats::rowMads(dmso_X, constant = 1.4826, na.rm = TRUE)
    g_mad[g_mad <= 0 | !is.finite(g_mad)] <- 1e-6
    
    rz <- sweep(dmso_X, 1, g_med, "-")
    rz <- sweep(rz, 1, g_mad, "/")
    
    gene_outlier_counts <- colSums(abs(rz) > z_thresh, na.rm = TRUE)
    
    # ============================================================
    # PART 2: PCA-aware outlier detection
    # ============================================================
    
    # PCA on (wells × genes); scale_pca controls per-gene unit-variance scaling.
    # See argument notes at the top of this function for trade-offs.
    pca <- prcomp(t(dmso_X), center = TRUE, scale. = scale_pca)
    
    if (is.character(npc) && npc == "auto") {
      var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
      npc_auto      <- which(var_explained >= 0.80)[1]
      npc_use       <- max(3, min(npc_auto, ncol(dmso_X) - 1, na.rm = TRUE))
      message("   PCs auto-selected: ", npc_use,
              " (", round(var_explained[npc_use] * 100, 1), "% variance)")
    } else {
      npc_use <- min(as.integer(npc), ncol(dmso_X) - 1)
      var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
      message("   PCs used: ", npc_use,
              " (", round(var_explained[npc_use] * 100, 1), "% variance)")
    }
    
    pc_scores <- pca$x[, seq_len(npc_use), drop = FALSE]
    
    pc_med <- apply(pc_scores, 2, median, na.rm = TRUE)
    pc_mad <- apply(pc_scores, 2, mad, constant = 1.4826, na.rm = TRUE)
    pc_mad[pc_mad <= 0 | !is.finite(pc_mad)] <- 1e-6
    
    pc_rz <- sweep(pc_scores, 2, pc_med, "-")
    pc_rz <- sweep(pc_rz, 2, pc_mad, "/")
    
    pc_outlier_counts <- rowSums(abs(pc_rz) > z_thresh, na.rm = TRUE)
    
    # ============================================================
    # PART 3: Combine gene + PC information
    # ============================================================
    
    cor_components <- cor(gene_outlier_counts, pc_outlier_counts,
                          use = "complete.obs", method = "spearman")
    message("   Gene-PCA outlier correlation (Spearman): ",
            round(cor_components, 3))
    
    # scale() divides by SD, which is 0 when an outlier-count vector is constant
    # (clean DMSO wells all share the same count, or a batch has very few wells).
    # That yields NaN scores -> NA median/threshold -> the scalar if() below
    # crashes with "missing value where TRUE/FALSE needed". safe_scale returns 0
    # for a zero-variance vector, so identical-looking wells simply tie: all
    # kept, none flagged.
    safe_scale <- function(x) {
      s <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(s) || s == 0) return(rep(0, length(x)))
      as.numeric(scale(x))
    }
    combined_score <- as.numeric(
      safe_scale(gene_outlier_counts) +
        pca_weight * safe_scale(pc_outlier_counts)
    )
    
    # ---- Store metrics ----
    dmso_meta$num_outlier_genes <- gene_outlier_counts[dmso_meta$plate_id]
    dmso_meta$num_outlier_pcs <- pc_outlier_counts[dmso_meta$plate_id]
    names(combined_score) <- rownames(pc_scores)
    dmso_meta$combined_score <- combined_score[dmso_meta$plate_id]
    
    
    # ============================================================
    # PART 4: Select best DMSO dose
    # ============================================================
    
    dose_summary <- dmso_meta %>%
      dplyr::group_by(treatment_dose) %>%
      dplyr::summarise(
        median_score = median(combined_score, na.rm = TRUE),
        n_wells      = dplyr::n(),
        .groups      = "drop"
      ) %>%
      dplyr::arrange(median_score)
    
    message("   Dose scores:")
    for (i in seq_len(nrow(dose_summary))) {
      message("     dose ", dose_summary$treatment_dose[i],
              ": median_score = ", round(dose_summary$median_score[i], 2),
              " (n = ", dose_summary$n_wells[i], " wells)")
    }
    
    utils::write.csv(
      dose_summary,
      file.path(out_data_dir, paste0("dmso_dose_summary_batch_", batch_key, ".csv")),
      row.names = FALSE
    )
    
    best_dose <- dose_summary$treatment_dose[1]
    
    best_unit <- dmso_meta %>%
      dplyr::filter(treatment_dose == best_dose) %>%
      dplyr::pull(treatment_dose_unit) %>%
      unique() %>%
      paste(collapse = ", ")
    
    message(
      " Best DMSO dose in batch ", batch_key, ": ",
      best_dose, " ", best_unit
    )
    
    dmso_meta <- dmso_meta %>% dplyr::filter(treatment_dose == best_dose)
    
    # ============================================================
    # PART 5: Select best / worst wells
    # ============================================================
    
    score_median  <- median(dmso_meta$combined_score, na.rm = TRUE)
    score_mad     <- mad(dmso_meta$combined_score, constant = 1.4826, na.rm = TRUE)
    bad_threshold <- score_median + 2 * score_mad
    
    if (identical(per_plate_top_n, "all")) {
      # Keep every well that is NOT flagged as bad — the most permissive
      # setting, useful when DMSO replicates are abundant (e.g. 16+ on a
      # 384-well plate).
      # ... but never fewer than min_keep_per_plate per plate (best scores first)
      best_per_plate <- dmso_meta %>%
        dplyr::group_by(plate_barcode) %>%
        dplyr::arrange(combined_score, .by_group = TRUE) %>%
        dplyr::filter(combined_score <= bad_threshold |
                        dplyr::row_number() <= min_keep_per_plate) %>%
        dplyr::ungroup()
      n_rescued <- sum(best_per_plate$combined_score > bad_threshold, na.rm = TRUE)
      if (n_rescued > 0)
        message("   Kept ", n_rescued, " flagged DMSO well(s) so that every plate keeps >= ",
                min_keep_per_plate, " controls")
    } else {
      best_per_plate <- dmso_meta %>%
        dplyr::group_by(plate_barcode) %>%
        dplyr::arrange(combined_score) %>%
        dplyr::slice_head(n = per_plate_top_n) %>%
        dplyr::ungroup()
    }
    
    bad_batch <- dmso_meta %>%
      dplyr::filter(combined_score > bad_threshold,
                    !plate_id %in% best_per_plate$plate_id)   # not the wells kept by min_keep_per_plate
    
    if (nrow(bad_batch) == 0 && nrow(dmso_meta) > 3) {
      worst <- dmso_meta %>%
        dplyr::arrange(dplyr::desc(combined_score)) %>%
        dplyr::slice_head(n = 1)
      if (isTRUE(worst$combined_score > score_median + score_mad)) {
        bad_batch <- worst
      }
    }
    
    message("   Bad DMSO wells flagged: ", nrow(bad_batch),
            " (threshold: ", round(bad_threshold, 2), ")")
    message("   Best DMSO wells kept:   ", nrow(best_per_plate))
    
    best_DMSO_all[[batch_key]] <- best_per_plate %>%
      dplyr::mutate(batch = batch_key, best_dose = best_dose)
    
    bad_DMSO_all[[batch_key]] <- bad_batch %>%
      dplyr::mutate(batch = batch_key)
    
    # ---- Optional diagnostics ----
    utils::write.csv(
      dmso_meta,
      file.path(out_data_dir, paste0("dmso_qc_batch_", batch_key, ".csv")),
      row.names = FALSE
    )
  }
  
  # ---- Output ----
  best_DMSO_df <- dplyr::bind_rows(best_DMSO_all)
  bad_DMSO_df <- dplyr::bind_rows(bad_DMSO_all)
  
  utils::write.csv(
    best_DMSO_df,
    file.path(out_data_dir, "best_DMSO_keep_pca.csv"),
    row.names = FALSE
  )
  
  utils::write.csv(
    bad_DMSO_df,
    file.path(out_data_dir, "bad_DMSO_remove_pca.csv"),
    row.names = FALSE
  )
  
  message("Finished. Results written to: ", out_data_dir)
  
  list(
    best_DMSO = best_DMSO_df,
    bad_DMSO  = bad_DMSO_df
  )
}

# --------------------------------------------------------------
# Clean DE table
# --------------------------------------------------------------

# Perform plate-aware differential expression analysis using limma with
# QC-filtered DMSO controls and per-timepoint contrasts.

#' Differential expression is computed for DRUG-seq data using a plate-aware
#' limma workflow that preserves experimental structure and minimizes technical
#' confounding. Analysis is performed independently for each batch and each
#' plate within a batch.
#'
#' For each plate, raw UMI counts are filtered with a presence-aware
#' edgeR::filterByExpr (min.count = 5, grouped by condition), normalized using
#' edgeR TMM normalization, and transformed to log2 CPM. Each contrast then
#' re-applies filterByExpr to its own treatment + DMSO wells.
#'
#' For each treatment–dose–timepoint, a two-group design (treatment vs DMSO)
#' is fitted on that contrast's wells only. Within each plate and timepoint,
#' treatment samples are contrasted against the corresponding DMSO control
#' measured on the same plate, ensuring strictly plate-local comparisons.

#'
#' Optionally, DMSO wells may be pre-filtered using an external QC procedure
#' (e.g. robust z-score–based stability scoring), such that only high-quality
#' DMSO controls contribute to model fitting and contrast definition.
#'
#' Linear models are fit using limma, followed by empirical Bayes moderation
#' with a mean–variance trend. Differential expression results from all valid
#' contrasts are collected and annotated with experimental metadata.
#'
#' Output tables are standardized by renaming key statistics, rounding
#' numerical values for readability, and joining external gene annotations.
#' Differential expression results are written per plate and combined into a
#' single data frame for downstream analysis.


clean_de_table <- function(de_df) {
  de_df %>%
    dplyr::rename(
      ENSEMBLID = GeneId,
    ) %>%
    dplyr::mutate(
      pval = signif(pval, 2),
      padj = signif(padj, 2)
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(LFC, AveExpr, t, B),
        ~ round(.x, 2)
      )
    )
}

annotate_drugseq_de <- function(de_results, gene_annotation) {
  stopifnot("ENSEMBLID" %in% colnames(de_results))
  stopifnot("ENSEMBLID" %in% colnames(gene_annotation))
  
  dplyr::left_join(
    de_results,
    gene_annotation,
    by = "ENSEMBLID"
  )
}

# =============================================================================
# DE limma helpers + two analysis functions (per-plate and collapsed/dupCor)
# Refactored for consistency, correctness, and robustness
# =============================================================================

# -----------------------------------------------------------------------------
# SHARED HELPER: filter, normalize → logCPM matrix
# -----------------------------------------------------------------------------
# filter_and_normalize is used by the UMAP path. 
# It performs a single filterByExpr pass with min.count = 3 — appropriate because
# UMAP needs all reasonably-expressed genes (no multiple-testing penalty) but does
# want to drop near-empty rows that would just add noise.
# The DE path does not use this function. It uses filter_and_normalize_hybrid
# (filterByExpr, min.count = 5, grouped by condition) globally, followed by a
# per-contrast filterByExpr with the SAME min.count inside run_single_contrast.
# The dose-trend path re-normalises the same wells with trend_min_count = 1 to
# build its own, more permissive gene universe.
#
#' @param counts  raw integer count matrix (genes × wells)
#' @param groups  character vector of group labels (length == ncol(counts))
#' @param min.count  passed to filterByExpr; default 3 is appropriate for
#'   the full-matrix UMAP filter (see file-level comment above for context).
#' 
filter_and_normalize <- function(counts, groups, min.count = 3) {
  dge <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, group = groups, min.count = min.count)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- dge[, colSums(dge$counts) > 0, keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  list(
    logCPM = edgeR::cpm(dge, log = TRUE, prior.count = 1),
    dge    = dge
  )
}

# -----------------------------------------------------------------------------
# 3. PER-PLATE FUNCTION
# -----------------------------------------------------------------------------

DE_limma_zQC <- function(
    Exp,
    CTRL            = "DMSO",
    meta,
    best_DMSO_wells = NULL,
    out_dir         = "out_data",
    annotated_genes,
    html_reports    = TRUE,
    html_dir        = "html_reports",
    P_val           = 0.05,
    LFC_value       = 1,
    max_html_rows   = 2000,
    min.count       = 5,
    top_pct          = 0.25,          # legacy mean-CPM pre-filter fraction
    use_filterByExpr = TRUE,          # NEW: presence-aware global filter (recommended)
    min_resid_df     = 2,             # replication threshold per contrast
    method           = c("auto", "trend", "voom")   # see header; auto = voom here
) {
  method <- match.arg(method)
  if (method == "auto") method <- "voom"     # per-plate fits are single-plate by construction
  message(" DE method (per-plate): limma-", method)
  dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
  if (html_reports) dir.create(html_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- HTML helper -----------------------------------------------------------
  make_html <- function(df_table, contrast_name, outfile) {
    if (nrow(df_table) == 0) return(NULL)
    
    has_col  <- function(nm) nm %in% colnames(df_table)
    ensg_idx <- which(colnames(df_table) == "ENSEMBLID") - 1L   # 0-based for JS
    
    # Concentrations: coerce to numeric (if it parses cleanly) so the column
    # sorts numerically and can be rounded for display below.
    if (has_col("treatment_dose") && !is.numeric(df_table$treatment_dose)) {
      td <- suppressWarnings(as.numeric(as.character(df_table$treatment_dose)))
      if (!any(is.na(td) & !is.na(df_table$treatment_dose)))
        df_table$treatment_dose <- td
    }
    
    # Modern self-contained CSS injected before the widget HTML
    page_css <- htmltools::tags$style(htmltools::HTML('
      body {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", sans-serif;
        background: #f9fafb;
        margin: 0;
        padding: 24px;
        color: #111827;
        font-size: 13px;
      }
      .dataTables_wrapper {
        background: #fff;
        border-radius: 12px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.04);
        padding: 20px 24px 16px;
      }
      table.dataTable { border-collapse: collapse !important; }

      /* Header */
      table.dataTable thead th {
        background: #f8fafc !important;
        color: #4b5563 !important;
        font-size: 11px !important;
        font-weight: 700 !important;
        text-transform: uppercase !important;
        letter-spacing: 0.06em !important;
        border-top: none !important;
        border-bottom: 2px solid #e2e8f0 !important;
        padding: 10px 14px !important;
        white-space: nowrap;
      }
      /* Header tooltips: dotted underline + help cursor hint that a
         hover explanation is available. The hdr-tip class is added to each
         relevant <th> by the headerCallback (see below) after render. */
      table.dataTable thead th.hdr-tip {
        cursor: help;
        text-decoration: underline dotted #9ca3af;
        text-underline-offset: 3px;
      }
      /* Column filter row */
      table.dataTable thead tr:last-child th {
        background: #fff !important;
        padding: 4px 8px !important;
        border-bottom: 1px solid #e2e8f0 !important;
      }
      table.dataTable thead tr:last-child input {
        width: 100%;
        border: 1px solid #e2e8f0;
        border-radius: 4px;
        padding: 3px 6px;
        font-size: 11px;
        color: #374151;
        outline: none;
      }
      table.dataTable thead tr:last-child input:focus {
        border-color: #0077c8;
        box-shadow: 0 0 0 2px rgba(0,119,200,0.1);
      }

      /* Body */
      table.dataTable tbody td {
        padding: 8px 14px !important;
        font-size: 13px !important;
        border-bottom: 1px solid #f1f5f9 !important;
        border-top: none !important;
        border-left: none !important;
        border-right: none !important;
        vertical-align: middle;
      }
      table.dataTable tbody tr:hover td { background: #f8fafc !important; }
      table.dataTable tbody tr:last-child td { border-bottom: none !important; }

      /* ENSEMBL links */
      a.ensg-link {
        color: #0077c8;
        text-decoration: none;
        font-family: ui-monospace, "IBM Plex Mono", "Fira Code", monospace;
        font-size: 12px;
        letter-spacing: 0.01em;
      }
      a.ensg-link:hover { text-decoration: underline; color: #005fa3; }

      /* Caption */
      caption {
        caption-side: top;
        text-align: left;
        font-size: 15px;
        font-weight: 600;
        color: #111827;
        padding: 0 0 14px;
      }

      /* Buttons */
      .dt-buttons { margin-bottom: 14px; }
      .dt-button {
        background: #fff !important;
        border: 1px solid #d1d5db !important;
        border-radius: 6px !important;
        color: #374151 !important;
        font-size: 12px !important;
        font-weight: 500 !important;
        padding: 6px 14px !important;
        margin-right: 6px !important;
        box-shadow: 0 1px 2px rgba(0,0,0,0.04) !important;
        cursor: pointer !important;
      }
      .dt-button:hover {
        background: #f3f4f6 !important;
        border-color: #9ca3af !important;
        color: #111827 !important;
      }
      .dt-button-collection {
        border-radius: 8px !important;
        box-shadow: 0 4px 12px rgba(0,0,0,0.12) !important;
        padding: 6px !important;
      }
      .dt-button-collection .dt-button {
        display: block !important;
        width: 100% !important;
        margin: 2px 0 !important;
      }

      /* Search */
      .dataTables_filter label { font-size: 13px; color: #6b7280; }
      .dataTables_filter input {
        border: 1px solid #d1d5db !important;
        border-radius: 6px !important;
        padding: 5px 10px !important;
        font-size: 13px !important;
        margin-left: 8px;
        outline: none;
        transition: border-color 0.15s, box-shadow 0.15s;
      }
      .dataTables_filter input:focus {
        border-color: #0077c8 !important;
        box-shadow: 0 0 0 3px rgba(0,119,200,0.1);
      }

      /* Length */
      .dataTables_length label { font-size: 13px; color: #6b7280; }
      .dataTables_length select {
        border: 1px solid #d1d5db !important;
        border-radius: 6px !important;
        padding: 4px 8px !important;
        font-size: 13px !important;
        margin: 0 6px;
      }

      /* Info + pagination */
      .dataTables_info { font-size: 12px; color: #9ca3af; padding-top: 10px; }
      .dataTables_paginate { padding-top: 10px; }
      .paginate_button {
        border-radius: 6px !important;
        font-size: 12px !important;
        padding: 4px 10px !important;
        margin: 0 2px !important;
        border: 1px solid transparent !important;
        color: #374151 !important;
      }
      .paginate_button.current, .paginate_button.current:hover {
        background: #0077c8 !important;
        border-color: #0077c8 !important;
        color: #fff !important;
      }
      .paginate_button:hover {
        background: #f3f4f6 !important;
        border-color: #d1d5db !important;
        color: #111827 !important;
      }
      .paginate_button.disabled, .paginate_button.disabled:hover {
        color: #d1d5db !important;
        cursor: default !important;
        background: none !important;
      }
    '))
    
    # Column definitions
    col_defs <- if (length(ensg_idx)) {
      list(list(
        targets = ensg_idx,
        render  = DT::JS("
          function(data, type, row, meta) {
            if (type === 'display')
              return '<a class=\"ensg-link\" href=\"https://www.ensembl.org/Gene/Summary?g=' +
                     data + '\" target=\"_blank\">' + data + '</a>';
            return data;
          }")
      ))
    } else list()
    
    # ---- Header tooltips -----------------------------------------------------
    # One short explanation per column, shown as a native hover tooltip. We do
    # NOT rename the columns (DT's named `colnames` rewrites the column keys,
    # which then breaks formatRound/formatStyle that reference them by name).
    # Instead a headerCallback runs after the header renders and, for each <th>
    # whose text matches a known column, sets a title attribute + a dotted
    # underline class. Column names stay intact, so the format* calls below work.
    col_tips <- c(
      ENSEMBLID = "Ensembl gene ID. Click a value to open the gene on Ensembl.",
      SYMBOL    = "Gene symbol (HGNC).",
      GENENAME  = "Full gene name.",
      LFC       = "log2 fold change, treatment vs DMSO. +1 = 2x higher in treatment, -1 = 2x lower. Sign gives direction, magnitude gives effect size.",
      AveExpr   = "Average expression of the gene across all wells in this contrast, on the log2-CPM scale. Low values flag genes near the detection floor.",
      t         = "Moderated t-statistic (limma empirical Bayes): the LFC divided by its shrunken standard error. Larger magnitude = stronger, more reliable signal.",
      pval      = "Raw (unadjusted) p-value. Not corrected for testing thousands of genes - do not use it on its own to call significance.",
      padj      = "Benjamini-Hochberg adjusted p-value (FDR). This is the column to threshold on; e.g. padj < 0.1 keeps the expected false-discovery rate at 10%.",
      B         = "Empirical Bayes log-odds that the gene is differentially expressed. B = 0 means 50:50 odds; B ~ 3 is about 95% probability. It is a monotonic transform of t, so it ranks genes the same way padj does and adds little beyond it.",
      n_genes_tested = "Number of genes that passed expression filtering and were tested in this contrast.",
      Sample    = "Contrast label (treatment + dose + timepoint vs control)."
    )
    tips_json <- as.character(jsonlite::toJSON(as.list(col_tips), auto_unbox = TRUE))
    header_cb <- DT::JS(sprintf(
      "function(thead, data, start, end, display) {
         var tips = %s;
         var ths = thead.getElementsByTagName('th');
         for (var i = 0; i < ths.length; i++) {
           var th = ths[i];
           var key = (th.textContent || '').trim();
           if (Object.prototype.hasOwnProperty.call(tips, key)) {
             th.setAttribute('title', tips[key]);
             if (th.className.indexOf('hdr-tip') === -1) th.className += ' hdr-tip';
           }
         }
       }", tips_json))
    
    dt <- DT::datatable(
      df_table,
      filter     = "top",
      rownames   = FALSE,
      escape     = FALSE,
      width      = "100%",
      class      = "display compact",
      extensions = c("Buttons", "FixedHeader"),
      caption    = paste("DEGs \u2014", contrast_name),
      options    = list(
        pageLength  = 25,
        lengthMenu  = c(10, 25, 50, 100),
        autoWidth   = FALSE,
        scrollX     = TRUE,
        dom         = "Bipfrtpl",
        fixedHeader = TRUE,
        buttons     = list(
          list(extend = "copy",  text = "Copy"),
          list(extend = "print", text = "Print"),
          list(extend  = "collection",
               text    = "Download \u25be",
               buttons = c("csv", "excel", "pdf"))
        ),
        columnDefs = col_defs,
        headerCallback = header_cb
      )
    )
    
    # Display formatting (filter still operates on raw values)
    round_cols  <- intersect(c("LFC", "AveExpr", "t", "B"), colnames(df_table))
    signif_cols <- intersect(c("pval", "padj"),              colnames(df_table))
    dose_cols   <- intersect("treatment_dose", colnames(df_table))
    dose_cols   <- dose_cols[vapply(dose_cols,
                                    function(c) is.numeric(df_table[[c]]),
                                    logical(1))]
    if (length(round_cols))  dt <- DT::formatRound(dt,  columns = round_cols,  digits = 3)
    if (length(signif_cols)) dt <- DT::formatSignif(dt, columns = signif_cols, digits = 3)
    if (length(dose_cols))   dt <- DT::formatSignif(dt, columns = dose_cols,   digits = 3)
    
    # LFC coloring: upregulated = red (#fee2e2 bg, #b91c1c text)
    #               downregulated = blue (#dbeafe bg, #1d4ed8 text)
    if (has_col("LFC")) {
      dt <- dt %>% DT::formatStyle(
        "LFC",
        backgroundColor = DT::styleInterval(
          c(-LFC_value, LFC_value),
          c("#dbeafe", "transparent", "#fee2e2")
        ),
        color = DT::styleInterval(
          c(-LFC_value, LFC_value),
          c("#1d4ed8", "#374151", "#b91c1c")
        ),
        fontWeight = DT::styleInterval(
          c(-LFC_value, LFC_value),
          c("600", "400", "600")
        )
      )
    }
    
    # padj coloring: violet < 0.001, purple < 0.01, amber < 0.05, grey otherwise
    if (has_col("padj")) {
      dt <- dt %>% DT::formatStyle(
        "padj",
        color = DT::styleInterval(
          c(0.001, 0.01, 0.05),
          c("#6d28d9", "#7c3aed", "#d97706", "#9ca3af")
        ),
        fontWeight = DT::styleInterval(
          c(0.001, 0.05),
          c("600", "500", "400")
        )
      )
    }
    
    dt <- dt %>% htmlwidgets::prependContent(page_css)
    htmlwidgets::saveWidget(dt, outfile, selfcontained = TRUE,
                            title = paste("DEGs |", contrast_name))
  }
  
  # ---- Main loop -----------------------------------------------------------
  for (b in names(Exp)) {
    message("\n Processing batch ", b)
    
    meta_b <- meta %>%
      dplyr::filter(batch == b) %>%
      dplyr::mutate(
        barcode_well = paste0(plate_well, "_", well_index),
        tp_clean     = paste0(stringr::str_extract(timepoint, "\\d+"), "hr"),
        dose_clean   = paste0(treatment_dose, treatment_dose_unit),
        sample_label = paste0(treatment, "_", dose_clean, "_", tp_clean),
        sample_clean = make.names(sample_label)
      )
    
    if (nrow(meta_b) == 0) next
    
    for (plate_name in names(Exp[[b]])) {
      message(" Plate ", plate_name)
      
      meta_p <- meta_b %>% dplyr::filter(plate_barcode == plate_name)
      if (!nrow(meta_p)) next
      
      # -- Best DMSO filtering -----------------------------------
      if (!is.null(best_DMSO_wells)) {
        dmso_keep <- best_DMSO_wells$barcode_well[
          best_DMSO_wells$batch         == b &
            best_DMSO_wells$plate_barcode == plate_name
        ]
        meta_p <- meta_p %>%
          dplyr::filter(treatment != CTRL | barcode_well %in% dmso_keep)
      }
      
      # -- Build count matrix ------------------------------------
      mat           <- as.data.frame(Exp[[b]][[plate_name]]$UMI.counts)
      rownames(mat) <- mat$gene
      mat           <- mat %>% dplyr::select(-gene, -sum)
      
      meta_p$well_id <- paste0(meta_p$plate_well, "_", meta_p$well_index)
      keep_cols      <- intersect(colnames(mat), meta_p$well_id)
      if (!length(keep_cols)) next
      
      mat    <- mat[, keep_cols, drop = FALSE]
      meta_p <- meta_p[match(keep_cols, meta_p$well_id), ]
      
      # global pre-filter (presence-aware filterByExpr by default) ----------
      norm <- filter_and_normalize_hybrid(
        mat,
        pct              = top_pct,
        groups           = meta_p$sample_clean,
        use_filterByExpr = use_filterByExpr,
        min.count        = min.count
      )
      logCPM_global <- norm$logCPM
      dge_global    <- norm$dge
      
      # Re-align meta to wells that survived the filter
      meta_p <- meta_p[meta_p$well_id %in% colnames(logCPM_global), ]
      if (!nrow(meta_p)) next
      
      # Named vector for quick well -> group lookup inside the contrast loop
      well_to_group <- setNames(meta_p$sample_clean, meta_p$well_id)
      
      # per-contrast loop (compound × dose × timepoint) ----------
      per_contrast_results <- list()
      skipped_contrasts    <- list()
      
      for (tp in unique(meta_p$tp_clean)) {
        
        m_tp     <- dplyr::filter(meta_p, tp_clean == tp)
        dmso_w   <- m_tp$well_id[m_tp$treatment == CTRL]
        ctrl_lbl <- unique(m_tp$sample_clean[m_tp$treatment == CTRL])
        
        if (!length(dmso_w) || !length(ctrl_lbl)) next
        
        # One group per unique compound × dose at this timepoint
        treat_groups <- m_tp %>%
          dplyr::filter(treatment != CTRL) %>%
          dplyr::group_by(sample_clean) %>%
          dplyr::summarise(wells = list(well_id), .groups = "drop")
        
        for (i in seq_len(nrow(treat_groups))) {
          lbl     <- treat_groups$sample_clean[i]
          treat_w <- treat_groups$wells[[i]]
          
          message("      Contrast: ", lbl, " vs ", ctrl_lbl,
                  "  (", length(treat_w), " treat + ",
                  length(dmso_w), " DMSO wells)")
          
          res_i <- run_single_contrast(
            logCPM_global = logCPM_global,
            dge_global    = dge_global,
            treat_wells   = treat_w,
            dmso_wells    = dmso_w,
            treat_label   = lbl,
            ctrl_label    = ctrl_lbl,
            well_to_group = well_to_group,
            min.count     = min.count,
            min_resid_df  = min_resid_df,
            method        = method
          )
          
          if (!is.null(res_i$results))
            per_contrast_results[[paste(tp, lbl, sep = "_")]] <- res_i$results
          if (!is.null(res_i$skip_info))
            skipped_contrasts[[length(skipped_contrasts) + 1]] <- res_i$skip_info
        }
      }
      
      # -- Write skip log (even if empty so QC always finds the file) --------
      skip_path <- file.path(
        out_dir,
        paste0("DE_batch_", b, "_plate_", plate_name, "_skipped.csv")
      )
      if (length(skipped_contrasts)) {
        skip_df <- dplyr::bind_rows(skipped_contrasts) %>%
          dplyr::mutate(batch = b, plate_barcode = plate_name, .before = 1)
        write.csv(skip_df, skip_path, row.names = FALSE)
        message("      Skipped ", nrow(skip_df), " contrast(s); log: ",
                basename(skip_path))
      } else if (file.exists(skip_path)) {
        file.remove(skip_path)   # stale log from a previous run
      }
      
      if (!length(per_contrast_results)) next
      res <- dplyr::bind_rows(per_contrast_results)
      
      # -- Annotation (same as original) -------------------------------------
      ann <- meta_p %>%
        dplyr::filter(treatment != CTRL) %>%
        dplyr::group_by(sample_clean) %>%
        dplyr::summarise(
          Sample              = dplyr::first(sample_label),
          treatment           = dplyr::first(treatment),
          treatment_dose      = dplyr::first(treatment_dose),
          treatment_dose_unit = dplyr::first(treatment_dose_unit),
          plate_barcode       = dplyr::first(plate_barcode),
          batch               = dplyr::first(batch),
          timepoint           = dplyr::first(timepoint),
          .groups = "drop"
        )
      
      DE <- res %>%
        dplyr::left_join(ann, by = "sample_clean") %>%
        clean_de_table() %>%
        annotate_drugseq_de(annotated_genes) %>%
        dplyr::select(-sample_clean)
      
      # -- Write CSV ---------------------------------------------------------
      write.csv(
        DE,
        file.path(out_dir,
                  paste0("DE_batch_", b, "_plate_", plate_name, ".csv")),
        row.names = FALSE
      )
      
      # -- Write HTML --------------------------------------------------------
      if (html_reports) {
        html_df <- DE %>%
          dplyr::filter(padj < P_val & abs(LFC) > LFC_value) %>%
          dplyr::arrange(padj, dplyr::desc(abs(LFC))) %>%
          dplyr::select(-dplyr::any_of("B"))   # drop B from HTML; still in CSV
        
        if (nrow(html_df) > max_html_rows)
          html_df <- dplyr::slice_head(html_df, n = max_html_rows)
        
        make_html(
          df_table      = html_df,
          contrast_name = paste(b, plate_name),
          outfile       = file.path(
            html_dir,
            paste0("DE_batch_", b, "_plate_", plate_name, ".html"))
        )
      }
    } 
  }
  
  gc()
  files <- list.files(out_dir, "DE_batch_.*_plate_.*\\.csv", full.names = TRUE)
  invisible(files)
}


# -----------------------------------------------------------------------------
# 4. COLLAPSED DUPCOR FUNCTION
# -----------------------------------------------------------------------------

DE_limma_zQC_collapse_dupCor <- function(
    Exp,
    CTRL            = "DMSO",
    meta,
    best_DMSO_wells = NULL,
    out_dir         = "out_data_collapse",
    annotated_genes,
    min.count       = 5,
    top_pct          = 0.25,
    use_filterByExpr = TRUE,          # presence-aware global filter (recommended)
    run_trend        = TRUE,          # also run per-compound dose-trend test
    trend_min_count  = 1,             # permissive filter for the dose-trend path
    #      (dose-coherence controls noise, not an
    #       expression cutoff) — keep this LOW so the
    #       trend can reach low-expression genes the
    #       strict per-dose filter (min.count) removes
    min_resid_df     = 2,             # replication threshold per contrast
    method           = c("auto", "trend", "voom")   # see header; auto decided per batch
) {
  method <- match.arg(method)
  message(" DE method (collapsed): ", method)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_batches <- list()
  
  for (b in names(Exp)) {
    message("\n Processing batch ", b, " (collapsed with duplicateCorrelation)")
    
    meta_sub <- meta[as.character(meta$batch) == as.character(b), , drop = FALSE]
    if (nrow(meta_sub) == 0) next
    
    meta_sub$barcode_well <- paste0(meta_sub$plate_well, "_", meta_sub$well_index)
    meta_sub$plate_barcode_short <- substr(
      meta_sub$plate_barcode,
      nchar(meta_sub$plate_barcode) - 2,
      nchar(meta_sub$plate_barcode)
    )
    
    # -- Best DMSO filtering -------------------------------------
    if (!is.null(best_DMSO_wells)) {
      dmso_keep <- best_DMSO_wells$barcode_well[
        as.character(best_DMSO_wells$batch) == as.character(b)
      ]
      meta_sub <- subset(meta_sub, treatment != CTRL | barcode_well %in% dmso_keep)
    }
    
    meta_sub$tp_clean     <- paste0(stringr::str_extract(meta_sub$timepoint, "\\d+"), "hr")
    meta_sub$dose_clean   <- paste0(meta_sub$treatment_dose, meta_sub$treatment_dose_unit)
    meta_sub$sample_label <- paste0(
      meta_sub$treatment, "_", meta_sub$dose_clean, "_", meta_sub$tp_clean
    )
    meta_sub$sample_clean <- make.names(meta_sub$sample_label)
    
    # -- Build UMI matrix across plates --------------------------
    UMI_list <- lapply(names(Exp[[b]]), function(plate_name) {
      mat           <- as.data.frame(Exp[[b]][[plate_name]]$UMI.counts)
      rownames(mat) <- mat$gene
      mat           <- mat[order(rownames(mat)), , drop = FALSE]
      colnames(mat) <- paste0(plate_name, "_", colnames(mat))
      mat
    })
    
    UMI <- do.call(cbind, UMI_list)
    UMI <- UMI[, !grepl("_gene$|_sum$", colnames(UMI)), drop = FALSE]
    
    meta_sub$plate_id_prefixed <- paste0(
      meta_sub$plate_barcode, "_",
      meta_sub$plate_well,    "_",
      meta_sub$well_index
    )
    
    common   <- intersect(colnames(UMI), meta_sub$plate_id_prefixed)
    if (!length(common)) next
    
    UMI      <- UMI[, common, drop = FALSE]
    meta_sub <- meta_sub[match(common, meta_sub$plate_id_prefixed), , drop = FALSE]
    
    # -- Degenerate-batch guard ------------------------------------------------
    # The collapse path needs >= 2 groups (DMSO + >= 1 treatment) to build any
    # contrast; a batch with a single group cannot and is skipped with a reason.
    # A SINGLE-PLATE batch is NOT skipped: it is fitted within-plate (plain lmFit,
    # no cross-plate dupCor), because with one plate there is no cross-plate
    # correlation to model. duplicateCorrelation is used only when a batch spans
    # >= 2 plates (see use_dupcor below). This lets fully independent-plate designs
    # (each compound on one plate) produce collapse results instead of nothing.
    n_groups_batch <- dplyr::n_distinct(meta_sub$sample_clean)
    if (n_groups_batch < 2) {
      message(" Skipping batch ", b, " for collapse (single_group_batch: ",
              n_groups_batch, " group). Need DMSO + >= 1 treatment group.")
      skip_path <- file.path(
        out_dir, paste0("DE_batch_", b, "_collapsed_dupCor_skipped.csv")
      )
      write.csv(
        data.frame(batch = b, reason = "single_group_batch",
                   n_plates = dplyr::n_distinct(meta_sub$plate_barcode),
                   n_groups = n_groups_batch),
        skip_path, row.names = FALSE
      )
      next
    }
    
    # global pre-filter (presence-aware filterByExpr by default) ------------
    norm <- filter_and_normalize_hybrid(
      UMI,
      pct              = top_pct,
      groups           = meta_sub$sample_clean,
      use_filterByExpr = use_filterByExpr,
      min.count        = min.count
    )
    logCPM_global <- norm$logCPM
    dge_global    <- norm$dge
    
    # Re-align meta to wells that survived the filter
    meta_sub <- meta_sub[meta_sub$plate_id_prefixed %in% colnames(logCPM_global), ]
    if (!nrow(meta_sub)) next
    
    # Re-check AFTER the empty-well drop: filter_and_normalize_hybrid removes
    # zero-count wells, which can eliminate a whole group. Still only skip on
    # < 2 groups (design cannot be built). n_plates decides whether cross-plate
    # dupCor is applicable for this batch.
    n_plates_batch <- dplyr::n_distinct(meta_sub$plate_barcode)
    n_groups_batch <- dplyr::n_distinct(meta_sub$sample_clean)
    if (n_groups_batch < 2) {
      message(" Skipping batch ", b, " for collapse (single_group_after_filter: ",
              n_groups_batch, " group survived the empty-well filter).")
      skip_path <- file.path(
        out_dir, paste0("DE_batch_", b, "_collapsed_dupCor_skipped.csv")
      )
      write.csv(
        data.frame(batch = b, reason = "single_group_after_filter",
                   n_plates = n_plates_batch, n_groups = n_groups_batch),
        skip_path, row.names = FALSE
      )
      next
    }
    
    # Cross-plate dupCor only applies when the batch spans >= 2 plates.
    use_dupcor <- n_plates_batch >= 2
    
    # auto: trend when plates are pooled with dupCor, voom for single-plate batches
    method_b <- if (method == "auto") { if (use_dupcor) "trend" else "voom" } else method
    message(" Batch ", b, ": limma-", method_b,
            if (use_dupcor) " with duplicateCorrelation" else " (single plate)")
    
    # Named vectors for the contrast loop
    well_to_group <- setNames(meta_sub$sample_clean, meta_sub$plate_id_prefixed)
    well_to_plate <- setNames(meta_sub$plate_barcode, meta_sub$plate_id_prefixed)
    
    # Estimate corfit ONCE globally before contrast loop (only when the batch
    # spans >= 2 plates). For a single-plate batch there is no cross-plate
    # correlation to estimate, so we fit within-plate (block/correlation NULL).
    if (use_dupcor) {
      design_global           <- model.matrix(~ 0 + sample_clean, data = meta_sub)
      colnames(design_global) <- sub("^sample_clean", "", colnames(design_global))
      
      # voom: first voom pass (no block) -> correlation on voom values; the
      # per-contrast voom + lmFit then use this consensus (second pass)
      dc_obj <- if (method_b == "voom") limma::voom(dge_global, design_global) else logCPM_global
      corfit <- limma::duplicateCorrelation(
        dc_obj, design_global,
        block = meta_sub$plate_barcode
      )
      message(" Consensus plate correlation: ", round(corfit$consensus, 3))
      
      # Flag compounds present on only one plate — dupCor cannot help them
      n_plates_per_compound <- meta_sub %>%
        dplyr::filter(treatment != CTRL) %>%
        dplyr::group_by(sample_clean) %>%
        dplyr::summarise(n_plates = dplyr::n_distinct(plate_barcode), .groups = "drop")
      
      single_plate <- n_plates_per_compound$sample_clean[
        n_plates_per_compound$n_plates == 1
      ]
      if (length(single_plate))
        message(" Single-plate compounds (dupCor provides no cross-plate benefit): ",
                paste(single_plate, collapse = ", "))
    } else {
      corfit <- list(consensus = NA_real_)
      message(" Batch ", b, " is single-plate — fitting within-plate ",
              "(no cross-plate dupCor borrowing).")
    }
    
    # per-contrast loop (compound × dose × timepoint) ------------
    per_contrast_results <- list()
    skipped_contrasts    <- list()
    
    for (tp in unique(meta_sub$tp_clean)) {
      
      m_tp     <- dplyr::filter(meta_sub, tp_clean == tp)
      dmso_w   <- m_tp$plate_id_prefixed[m_tp$treatment == CTRL]
      ctrl_lbl <- unique(m_tp$sample_clean[m_tp$treatment == CTRL])
      
      if (!length(dmso_w) || !length(ctrl_lbl)) next
      
      # One group per unique compound × dose at this timepoint
      treat_groups <- m_tp %>%
        dplyr::filter(treatment != CTRL) %>%
        dplyr::group_by(sample_clean) %>%
        dplyr::summarise(wells = list(plate_id_prefixed), .groups = "drop")
      
      for (i in seq_len(nrow(treat_groups))) {
        lbl     <- treat_groups$sample_clean[i]
        treat_w <- treat_groups$wells[[i]]
        
        message(" Contrast: ", lbl, " vs ", ctrl_lbl,
                "  (", length(treat_w), " treat + ",
                length(dmso_w), " DMSO wells)")
        
        res_i <- run_single_contrast(
          logCPM_global = logCPM_global,
          dge_global    = dge_global,
          treat_wells   = treat_w,
          dmso_wells    = dmso_w,
          treat_label   = lbl,
          ctrl_label    = ctrl_lbl,
          well_to_group = well_to_group,
          min.count     = min.count,
          min_resid_df  = min_resid_df,
          block         = if (use_dupcor) well_to_plate else NULL,
          correlation   = if (use_dupcor) corfit$consensus else NULL,
          method        = method_b
        )
        
        if (!is.null(res_i$results))
          per_contrast_results[[paste(tp, lbl, sep = "_")]] <- res_i$results
        if (!is.null(res_i$skip_info))
          skipped_contrasts[[length(skipped_contrasts) + 1]] <- res_i$skip_info
      }
    }
    
    # -- Write skip log (even if empty so QC always finds the file) ----------
    skip_path <- file.path(
      out_dir, paste0("DE_batch_", b, "_collapsed_dupCor_skipped.csv")
    )
    if (length(skipped_contrasts)) {
      skip_df <- dplyr::bind_rows(skipped_contrasts) %>%
        dplyr::mutate(batch = b, .before = 1)
      write.csv(skip_df, skip_path, row.names = FALSE)
      message(" Skipped ", nrow(skip_df), " contrast(s); log: ",
              basename(skip_path))
    } else if (file.exists(skip_path)) {
      file.remove(skip_path)
    }
    
    if (!length(per_contrast_results)) next
    res <- dplyr::bind_rows(per_contrast_results)
    
    # -- Annotation (same as original) ----------------------------------------
    ann <- meta_sub %>%
      dplyr::filter(treatment != CTRL) %>%
      dplyr::group_by(sample_clean) %>%
      dplyr::summarise(
        Sample              = dplyr::first(sample_label),
        treatment           = dplyr::first(treatment),
        treatment_dose      = dplyr::first(treatment_dose),
        treatment_dose_unit = dplyr::first(treatment_dose_unit),
        timepoint           = dplyr::first(timepoint),
        batch               = dplyr::first(batch),
        plate_well          = paste(
          unique(paste(plate_barcode_short, plate_well)), collapse = ","
        ),
        .groups = "drop"
      )
    
    DE_batch <- res %>%
      dplyr::left_join(ann, by = "sample_clean") %>%
      clean_de_table() %>%
      annotate_drugseq_de(annotated_genes) %>%
      dplyr::select(-sample_clean)
    
    write.csv(
      DE_batch,
      file.path(out_dir, paste0("DE_batch_", b, "_collapsed_dupCor.csv")),
      row.names = FALSE
    )
    
    all_batches[[b]] <- DE_batch
    
    # ---- DOSE-TREND analysis (per compound across the dose series) --------
    # Complementary to the per-dose contrasts above. Recovers low-expression,
    # dose-responsive genes that no single per-dose contrast calls significant,
    # while rejecting dose-independent noise (no trend). Written to its own CSV
    # so the per-dose output is unchanged; join on (treatment, ENSEMBLID)
    # downstream to flag "trend-only" rescued genes. Reuses the global corfit.
    if (run_trend) {
      # ---- SEPARATE permissive universe for the trend test -------------------
      # The per-dose path above uses the strict global filter (min.count) for
      # specificity, which floors expression around AveExpr ~4.4 at DRUG-seq
      # depth. The trend test relies on dose-coherence (not an expression
      # cutoff) to control noise, so it must reach BELOW that floor. Re-normalize
      # the same wells with a low min.count just for this path; the strict
      # per-dose output is untouched. corfit$consensus (a single scalar plate
      # correlation) is reused — standard limma practice, transferable across
      # gene sets.
      norm_t <- filter_and_normalize_hybrid(
        UMI[, meta_sub$plate_id_prefixed, drop = FALSE],
        pct              = top_pct,
        groups           = meta_sub$sample_clean,
        use_filterByExpr = TRUE,
        min.count        = trend_min_count
      )
      logCPM_trend <- norm_t$logCPM
      dge_trend    <- norm_t$dge
      message(sprintf(
        " Dose-trend universe (batch %s): %d genes (permissive min.count=%g) "
        , b, nrow(logCPM_trend), trend_min_count),
        sprintf("vs %d strict — floor min AveExpr %.2f",
                nrow(logCPM_global),
                suppressWarnings(min(rowMeans(logCPM_trend)))))
      
      dose_trend_results <- list()
      dose_trend_skipped <- list()
      
      for (tp in unique(meta_sub$tp_clean)) {
        m_tp   <- dplyr::filter(meta_sub, tp_clean == tp)
        dmso_w <- m_tp$plate_id_prefixed[m_tp$treatment == CTRL]
        if (!length(dmso_w)) next
        
        treated <- dplyr::filter(m_tp, treatment != CTRL)
        for (cmpd in unique(treated$treatment)) {
          cw  <- dplyr::filter(treated, treatment == cmpd)
          compound_wells <- setNames(as.numeric(cw$treatment_dose),
                                     cw$plate_id_prefixed)
          lbl <- paste0(cmpd, "_", tp)
          
          message(" Dose-trend: ", lbl, " (",
                  length(unique(compound_wells)), " doses + ",
                  length(dmso_w), " DMSO wells)")
          
          res_i <- run_dose_trend(
            logCPM_global  = logCPM_trend,
            dge_global     = dge_trend,
            compound_wells = compound_wells,
            dmso_wells     = dmso_w,
            treat_label    = lbl,
            min.count      = trend_min_count,
            block          = if (use_dupcor) well_to_plate else NULL,
            correlation    = if (use_dupcor) corfit$consensus else NULL
          )
          
          if (!is.null(res_i$results))
            dose_trend_results[[lbl]] <- res_i$results %>%
            dplyr::mutate(treatment = cmpd, timepoint_clean = tp, batch = b)
          if (!is.null(res_i$skip_info))
            dose_trend_skipped[[length(dose_trend_skipped) + 1]] <- res_i$skip_info
        }
      }
      
      if (length(dose_trend_results)) {
        dt <- dplyr::bind_rows(dose_trend_results) %>%
          dplyr::rename(ENSEMBLID = GeneId) %>%
          dplyr::mutate(
            dplyr::across(c(dose_slope, AveExpr, t, B), ~ round(.x, 2)),
            pval = signif(pval, 2),
            padj = signif(padj, 2)
          ) %>%
          annotate_drugseq_de(annotated_genes)
        
        write.csv(
          dt,
          file.path(out_dir, paste0("DE_batch_", b, "_dose_trend.csv")),
          row.names = FALSE
        )
        message(" Dose-trend: wrote ", length(dose_trend_results),
                " compound x timepoint trend table(s)")
      }
      if (length(dose_trend_skipped)) {
        write.csv(
          dplyr::bind_rows(dose_trend_skipped),
          file.path(out_dir, paste0("DE_batch_", b, "_dose_trend_skipped.csv")),
          row.names = FALSE
        )
      }
    }
  }
  
  if (!length(all_batches)) return(NULL)
  dplyr::bind_rows(all_batches)
}