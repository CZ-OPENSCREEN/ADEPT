# --------------------------------------------------------------
# Read in count matrices
# --------------------------------------------------------------

read_starsolo_simple <- function(sample_dir) {
  
  solo_dir <- file.path(sample_dir, "Solo.out", "Gene", "filtered")
  
  if (!dir.exists(solo_dir)) {
    stop("STARsolo output directory not found: ", solo_dir)
  }
  
  mtx_file      <- file.path(solo_dir, "matrix.mtx")
  barcodes_file <- file.path(solo_dir, "barcodes.tsv")
  features_file <- file.path(solo_dir, "features.tsv")
  
  # Check all files exist
  missing <- c(mtx_file, barcodes_file, features_file)[
    !file.exists(c(mtx_file, barcodes_file, features_file))
  ]
  if (length(missing) > 0) {
    stop("Missing files: ", paste(missing, collapse = ", "))
  }
  
  mat      <- Matrix::readMM(mtx_file)
  barcodes <- read.delim(barcodes_file, header = FALSE)
  features <- read.delim(features_file, header = FALSE)
  
  # Fail loudly on a dimension mismatch. R's `rownames(mat) <- features$V1`
  # silently recycles/truncates when lengths differ, which relabels rows with
  # the WRONG gene IDs (and shifts every label after the first gap) rather than
  # erroring. That produces plausible-looking but corrupt DE results. These
  # checks turn that silent trap into an explicit, directory-named error.
  if (nrow(mat) != nrow(features)) {
    stop("Dimension mismatch in ", solo_dir, ": matrix.mtx has ", nrow(mat),
         " gene rows but features.tsv has ", nrow(features), " entries. ",
         "The matrix and features file are out of sync - regenerate them ",
         "together (do not edit features.tsv on its own).")
  }
  if (ncol(mat) != nrow(barcodes)) {
    stop("Dimension mismatch in ", solo_dir, ": matrix.mtx has ", ncol(mat),
         " cell columns but barcodes.tsv has ", nrow(barcodes), " entries.")
  }
  
  # STARsolo matrix.mtx is genes x barcodes
  rownames(mat) <- features$V1   # Ensembl ID
  colnames(mat) <- barcodes$V1
  
  return(mat)
}

read_alevin_mtx_simple <- function(sample_dir) {
  
  # Bare alevin-fry mtx export: a finished MatrixMarket triplet with no
  # quant.json / meta_info.json (so fishpond::loadFry cannot be used) and no
  # EDS quants_mat.gz (so tximeta cannot be used). Read it directly.
  alevin_dir <- file.path(sample_dir, "alevin")
  if (!dir.exists(alevin_dir)) {
    stop("No alevin/ subdirectory found in: ", sample_dir)
  }
  
  mtx_file  <- file.path(alevin_dir, "quants_mat.mtx")
  cols_file <- file.path(alevin_dir, "quants_mat_cols.txt")  # genes
  rows_file <- file.path(alevin_dir, "quants_mat_rows.txt")  # barcodes
  
  missing <- c(mtx_file, cols_file, rows_file)[
    !file.exists(c(mtx_file, cols_file, rows_file))
  ]
  if (length(missing) > 0) {
    stop("Missing alevin-fry mtx files: ", paste(missing, collapse = ", "))
  }
  
  mat   <- Matrix::readMM(mtx_file)   # alevin-fry: barcodes (rows) x genes (cols)
  genes <- readLines(cols_file)
  cells <- readLines(rows_file)
  
  # Fail loudly on a dimension mismatch. As in read_starsolo_simple, assigning
  # dimnames of the wrong length silently recycles/truncates and mislabels every
  # gene/barcode, producing plausible-looking but corrupt DE results. alevin-fry
  # writes cells as ROWS and genes as COLUMNS - the transpose of STARsolo - so
  # these checks also guard against the matrix ever arriving in the other
  # orientation (which would otherwise be silently transposed the wrong way).
  if (nrow(mat) != length(cells)) {
    stop("Dimension mismatch in ", alevin_dir, ": quants_mat.mtx has ", nrow(mat),
         " barcode rows but quants_mat_rows.txt has ", length(cells), " entries. ",
         "The matrix and row-label file are out of sync - regenerate them together.")
  }
  if (ncol(mat) != length(genes)) {
    stop("Dimension mismatch in ", alevin_dir, ": quants_mat.mtx has ", ncol(mat),
         " gene columns but quants_mat_cols.txt has ", length(genes), " entries. ",
         "The matrix and column-label file are out of sync - regenerate them together.")
  }
  
  rownames(mat) <- cells
  colnames(mat) <- genes
  
  # Return genes x barcodes to match read_starsolo_simple / loadFry output.
  Matrix::t(mat)
}

read_in_quant_files <- function(dir,
                                tool = c("alevin", "alevin_classic", "starsolo", "alevin_mtx"),
                                plate_layout = NULL) {
  tool <- match.arg(tool)
  
  samples <- list.dirs(dir, full.names = FALSE, recursive = FALSE)
  samples <- samples[nzchar(samples)]
  if (length(samples) == 0) stop("No sample directories found in: ", dir)
  
  result_list <- vector("list", length(samples))
  names(result_list) <- samples
  
  for (i in seq_along(samples)) {
    sample_dir <- file.path(dir, samples[i])
    
    if (tool == "alevin") {
      # Alevin-Fry path (unchanged)
      sce <- fishpond::loadFry(
        fryDir    = sample_dir,
        outputFormat = "all"
      )
      mat <- assay(sce, "counts")
      
    } else if (tool == "alevin_classic") {
      # Classic Alevin path — use tximeta::loadFry or alevinQC::alevinQCPathsFromDir
      # Most reliable: tximport
      alevin_dir <- file.path(sample_dir, "alevin")
      
      if (!dir.exists(alevin_dir)) {
        stop("No alevin/ subdirectory found in: ", sample_dir)
      }
      
      txi <- tximeta::tximeta(
        coldata    = data.frame(
          names   = samples[i],
          files   = file.path(alevin_dir, "quants_mat.gz"),
          stringsAsFactors = FALSE
        ),
        type       = "alevin",
        dropInfReps = TRUE
      )
      mat <- assay(txi, "counts")
      # genes x barcodes
      
    } else if (tool == "alevin_mtx") {
      # Bare alevin-fry mtx export (finished triplet, no quant.json)
      mat <- read_alevin_mtx_simple(sample_dir)
      mat <- as.matrix(mat)
      
    } else {
      # STARsolo path (unchanged)
      mat <- read_starsolo_simple(sample_dir)
      mat <- as.matrix(mat)
    }
    
    # ---- Shared post-processing ----
    # Barcode → well renaming
    if (!is.null(plate_layout) &&
        all(c("barcode", "well") %in% names(plate_layout))) {
      m <- match(colnames(mat), plate_layout$barcode)
      if (any(is.na(m))) {
        stop(
          "Barcodes not found in plate layout: ",
          paste(head(colnames(mat)[is.na(m)], 10), collapse = ", ")
        )
      }
      colnames(mat) <- plate_layout$well[m]
    }
    
    df        <- as.data.frame(as.matrix(mat))
    df        <- tibble::rownames_to_column(df, "gene")
    df$sum    <- rowSums(df[, -1, drop = FALSE])
    df        <- df[order(-df$sum), ]
    
    result_list[[i]] <- df
  }
  
  return(result_list)
}


# --------------------------------------------------------------
# Read Plate layout file .xlsx, .xlx, .csv, .tsv
# --------------------------------------------------------------

read_plate_layout <- function(filename,
                              sheet = 1,
                              required_cols = c("row", "col", "plate_well", "well_index")) {
  # 1. Detect extension based on filename only
  ext <- tolower(tools::file_ext(filename))
  
  # 2. Read file based on extension
  df <- switch(ext,
               "csv"  = read.csv(filename, stringsAsFactors = FALSE),
               "tsv"  = read.delim(filename, stringsAsFactors = FALSE),
               "xlsx" = readxl::read_excel(filename, sheet = sheet),
               "xls"  = readxl::read_excel(filename, sheet = sheet),
               stop("Unsupported file type: .", ext)
  )
  
  # 3. Convert to tibble
  df <- tibble::as_tibble(df)
  
  # 4. Check required columns FIRST, so a missing `col` produces this message
  #    rather than an opaque error from the sprintf() below.
  if (!all(required_cols %in% names(df))) {
    missing <- setdiff(required_cols, names(df))
    stop(
      "Missing required columns in plate layout: ",
      paste(missing, collapse = ", ")
    )
  }
  
  # 5. Zero-pad the column index so wells sort A01 < A02 < ... < A24
  df$col <- sprintf("%02d", as.integer(df$col))
  
  df
}

process_plate_metadata <- function(df) {
  # --- Clean column names ---
  df <- janitor::clean_names(df)
  
  # --- Check for essential columns (raw input must have at least these) ---
  required <- c(
    "plate_barcode", "treatment", "treatment_dose", "treatment_dose_unit",
    "timepoint", "plate_well", "well_index", "cell_line_name","batch"
  )
  
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required metadata columns: ", paste(missing, collapse = ", "))
  }
  
  # --- Normalise key column types -------------------------------------------
  # plate_barcode and well_index are identifiers, not numbers. When the sheet
  # provides them as integers (Ginkgo uses a 1..N per-plate well index), they
  # must be character so they (a) match the count-matrix column names and
  # (b) join cleanly against values derived from those names (always character).
  # This also resolves the well_index <chr> vs <int> left_join error in the
  # heatmap builder.
  df$plate_barcode <- as.character(df$plate_barcode)
  df$well_index    <- as.character(df$well_index)
  
  # --- Fail loudly on missing KEY values ------------------------------------
  # plate_barcode / batch / well_index key the per-plate/per-batch experiment
  # list downstream: make_experiment_list() loops over unique(batch) and indexes
  # combined_counts[[plate_barcode]]. A single NA here - commonly from blank or
  # merged cells in the layout sheet, where readxl fills only the first row of a
  # merged block - propagates to combined_counts[[NA]], which dies with an opaque
  # "subscript out of bounds" several targets later (experiment_object). Catch it
  # here, at the source, with an actionable message.
  is_blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
  key_cols <- c("plate_barcode", "batch", "well_index")
  na_counts <- vapply(key_cols, function(cn) sum(is_blank(df[[cn]])), integer(1))
  if (any(na_counts > 0)) {
    bad <- key_cols[na_counts > 0]
    stop(
      "Plate layout has missing/blank values in key column(s): ",
      paste(sprintf("%s (%d wells)", bad, na_counts[bad]), collapse = ", "),
      ". Every well needs a non-missing plate_barcode, batch and well_index. ",
      "If the sheet uses merged/blank cells per plate block, fill them down ",
      "before running, e.g. tidyr::fill(df, plate_barcode, batch, .direction = \"down\")."
    )
  }
  
  # --- Replicate numbering -----------------------------------------------
  # Numbered within batch x treatment x dose x unit x timepoint, i.e. per
  # CONDITION within a batch, NOT per plate: if the same condition sits on two
  # plates of one batch its wells are rep1..rep4, ordered by plate, row, col.
  # `replicate` is only ever used together with `plate_barcode` downstream
  # (heatmap aggregate_by = "replicate"), so this is safe, but do not read
  # plate_replicate as a within-plate index.
  df <- df %>%
    dplyr::group_by(
      batch,
      treatment,
      treatment_dose,
      treatment_dose_unit,
      timepoint
    ) %>%
    dplyr::arrange(plate_barcode,
                   row,
                   as.numeric(col),
                   .by_group = TRUE) %>%
    dplyr::mutate(
      plate_replicate = row_number(),
      replicate = paste0("rep", plate_replicate)
    ) %>%
    dplyr::ungroup()
  
  # --- Compute replicate label (rep1, rep2, rep3) ---
  df <- df %>%
    dplyr::mutate(
      replicate = paste0("rep", plate_replicate),
      timepoint = paste0(stringr::str_extract(timepoint, "\\d+"), "hr"),
      sample = paste0(treatment, "_", treatment_dose,treatment_dose_unit, "_", timepoint)
    )
  
  # --- Compute unique well ID ---
  df <- df %>%
    dplyr::mutate(
      plate_id = paste0(plate_barcode, "_", plate_well, "_", well_index)
    )
  
  # --- Reorder columns to desired layout ---
  column_order <- c(
    "sample", "plate_barcode", "treatment", "treatment_dose",
    "treatment_dose_unit", "cell_line_name", "timepoint",
    "plate_well", "batch", "well_index", "plate_id", "plate_replicate",
    "replicate", "col", "row"
  )
  
  # Keep only those that exist
  existing <- intersect(column_order, names(df))
  df <- df[, existing]
  
  return(df)
}


# --------------------------------------------------------------
# Combine meta and count matrices
# --------------------------------------------------------------

combine_metadata_counts <- function(df_list, meta) {
  
  meta$plate_barcode <- as.character(meta$plate_barcode)
  meta$well_index    <- as.character(meta$well_index)
  
  plate_ids <- unique(meta$plate_barcode)
  
  # --- Map each count table to its plate ---------------------------
  # well_index is only unique WITHIN a plate (Ginkgo: 1..N per plate), so the
  # rename map must be built per plate, keyed on that plate's (plate_barcode,
  # well_index). To do that we first need to know which plate each count table
  # belongs to. Prefer matching by name (df_list is named by its source folder,
  # which equals plate_barcode); fall back to positional order (table k <->
  # plate_ids[k]) when names do not line up, which reproduces the previous
  # behaviour exactly. For datasets where well_index is globally unique (e.g.
  # the Novartis/primary benchmark) a per-plate map yields output identical to
  # the old global map, so this does not disturb frozen results.
  if (!is.null(names(df_list)) && all(names(df_list) %in% plate_ids)) {
    plate_of <- names(df_list)
  } else {
    if (length(df_list) != length(plate_ids)) {
      stop("Cannot map count tables to plates: ", length(df_list),
           " count tables but ", length(plate_ids), " plate_barcode(s), and the ",
           "count-table names do not match plate_barcode. Name each count table ",
           "by its plate_barcode, or supply them in plate order.")
    }
    plate_of <- plate_ids
  }
  
  split_names <- function(x) {
    well <- sub("_.*", "", x)
    
    tibble::tibble(
      name = x,
      letter = stringr::str_extract(well, "^[A-Za-z]+"),
      number = as.integer(stringr::str_extract(well, "[0-9]+"))
    )
  }
  
  fixed_cols  <- c("gene", "sum")
  df_list_out <- vector("list", length(df_list))
  
  for (k in seq_along(df_list)) {
    plate <- df_list[[k]]
    pb    <- plate_of[[k]]
    pmeta <- meta[meta$plate_barcode == pb, , drop = FALSE]
    
    if (nrow(pmeta) == 0) {
      stop("No metadata rows for plate_barcode '", pb,
           "' while combining count table ", k, ".")
    }
    
    # Per-plate rename map: this plate's well_index -> <plate_well>_<well_index>
    rename_map  <- stats::setNames(
      paste0(pmeta$plate_well, "_", pmeta$well_index),
      pmeta$well_index
    )
    valid_wells <- unname(rename_map)
    
    plate <- plate %>%
      dplyr::rename_with(~ dplyr::recode(.x, !!!rename_map, .default = .x))
    
    missing_wells <- setdiff(valid_wells, colnames(plate))
    if (length(missing_wells) > 0) {
      message(
        "⚠ Plate ", pb, ": missing wells (", length(missing_wells), "): ",
        paste(utils::head(missing_wells, 10), collapse = ", "),
        if (length(missing_wells) > 10) " ..." else ""
      )
    }
    
    plate <- plate %>%
      dplyr::relocate(gene, sum)
    
    wells        <- intersect(colnames(plate), valid_wells)
    well_split   <- split_names(wells)
    wells_sorted <- well_split %>%
      dplyr::arrange(letter, number) %>%
      dplyr::pull(name)
    
    df_list_out[[k]] <- plate[, c(fixed_cols, wells_sorted), drop = FALSE]
  }
  
  names(df_list_out) <- plate_of
  
  df_list_out
}

# =============================================================================
# build_gene_annotation.R
#
# Drop-in replacement for the existing build_gene_annotation() and
# create_gene_annotation_file().
#
# Changes:
#   * Records Ensembl provenance (release version, query date, mirror used,
#     dataset, gene-coverage stats) alongside the annotation table.
#   * Saves both as a single bundle: list(annotation, ensembl_version,
#     query_date, mirror_used, dataset, n_genes_queried, n_genes_returned).
#   * No change to the connection logic or chunked-query retries.
# =============================================================================

build_gene_annotation <- function(
    gene_ids,
    biomart_dataset,
    biomart_host = "www.ensembl.org",
    chunk_size = 500
) {
  
  stopifnot(requireNamespace("biomaRt"))
  stopifnot(requireNamespace("dplyr"))
  stopifnot(requireNamespace("stringr"))
  
  options(timeout = 300)
  
  primary_mirror <- sub("\\..*", "", sub("^https?://", "", biomart_host))
  all_mirrors    <- unique(c(primary_mirror, "www", "uswest", "useast", "asia"))
  valid_mirrors  <- all_mirrors[all_mirrors %in% c("www", "uswest", "useast", "asia")]
  
  mart        <- NULL
  mirror_used <- NA_character_
  
  for (m in valid_mirrors) {
    message(sprintf(">> Attempting to connect to mirror: %s", m))
    for (attempt in 1:2) {
      mart <- tryCatch({
        biomaRt::useEnsembl(
          biomart = "genes",
          dataset = biomart_dataset,
          mirror  = m
        )
      }, error = function(e) {
        message(sprintf("   Attempt %d on '%s' failed: %s",
                        attempt, m, conditionMessage(e)))
        NULL
      })
      if (!is.null(mart)) break
      Sys.sleep(2)
    }
    if (!is.null(mart)) {
      message(sprintf("Successfully connected to mirror: %s", m))
      mirror_used <- m
      break
    }
  }
  
  if (is.null(mart)) {
    stop("CRITICAL: All Ensembl mirrors are unreachable. ",
         "Check your internet or Ensembl status.")
  }
  
  # ---------------------------------------------------------------------------
  # Capture Ensembl release version
  #
  # biomaRt's mart object stores the release version in the
  # @host slot's URL after useEnsembl() resolves it (e.g. "https://nov2024
  # .archive.ensembl.org" → 113), or in listEnsemblArchives() against the
  # active release date. Three strategies, tried in order:
  #
  #   1. listEnsemblArchives() and pick the row marked current_release == "*"
  #   2. Parse the version from the resolved mart's host URL
  #   3. Fall back to "unknown" — never let provenance recording break the
  #      pipeline
  # ---------------------------------------------------------------------------
  ensembl_version <- tryCatch({
    archives <- biomaRt::listEnsemblArchives()
    current_row <- archives[archives$current_release == "*", , drop = FALSE]
    if (nrow(current_row) >= 1 && nzchar(current_row$version[1])) {
      as.character(current_row$version[1])
    } else {
      # Strategy 2: try to parse from the mart's resolved URL
      host_str  <- tryCatch(mart@host, error = function(e) "")
      ver_match <- regmatches(host_str, regexpr("[0-9]{2,3}", host_str))
      if (length(ver_match) && nzchar(ver_match)) ver_match else "unknown"
    }
  }, error = function(e) "unknown")
  
  message(sprintf("Ensembl release detected: %s", ensembl_version))
  
  # ---- Querying in Chunks with Internal Retries (unchanged) -----------------
  gene_ids <- unique(gene_ids)
  chunks   <- split(gene_ids, ceiling(seq_along(gene_ids) / chunk_size))
  ann_list <- vector("list", length(chunks))
  
  for (i in seq_along(chunks)) {
    message(sprintf("Querying chunk %d/%d...", i, length(chunks)))
    success <- FALSE
    attempt_bm <- 1
    
    while (!success && attempt_bm <= 3) {
      ann_list[[i]] <- tryCatch({
        res <- biomaRt::getBM(
          attributes = c("ensembl_gene_id", "external_gene_name",
                         "entrezgene_id", "description"),
          filters = "ensembl_gene_id",
          values  = chunks[[i]],
          mart    = mart
        )
        success <- TRUE
        res
      }, error = function(e) {
        message(sprintf("   Query failed (attempt %d). Retrying in 5s...",
                        attempt_bm))
        Sys.sleep(5)
        attempt_bm <<- attempt_bm + 1
        NULL
      })
    }
    
    if (!success) {
      warning(sprintf(
        "Failed to retrieve data for chunk %d after 3 attempts.", i))
    }
  }
  
  # ---- Data Consolidation and Cleaning (unchanged) --------------------------
  # Normalise column types before binding. getBM() normally returns character
  # columns, but a chunk that comes back empty (e.g. a retired/unmapped Ensembl
  # ID, or IDs absent from this mirror's release) can yield a 0-row frame whose
  # columns default to logical/integer. dplyr::bind_rows() then aborts with
  # "Can't combine `..1$external_gene_name` <character> and ...<integer>".
  # Coercing each non-NULL chunk to a fixed character schema makes the bind
  # type-safe regardless of what any individual chunk returned.
  expected_cols <- c("ensembl_gene_id", "external_gene_name",
                     "entrezgene_id", "description")
  
  ann_list <- lapply(ann_list, function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    for (col in expected_cols) {
      if (col %in% names(df)) df[[col]] <- as.character(df[[col]])
    }
    df
  })
  
  ann <- dplyr::bind_rows(ann_list)
  if (nrow(ann) == 0) {
    stop("No gene annotations retrieved from Ensembl.")
  }
  
  ann_clean <- ann %>%
    dplyr::mutate(
      external_gene_name = dplyr::na_if(external_gene_name, ""),
      description        = dplyr::na_if(description, "")
    ) %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarise(
      SYMBOL   = dplyr::first(stats::na.omit(external_gene_name)),
      GENENAME = dplyr::first(stats::na.omit(description)),
      GENEALIASES = paste(
        unique(stats::na.omit(external_gene_name)),
        collapse = ";"
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      GENENAME = stringr::str_remove(GENENAME, "\\s*\\[.*?\\]")
    ) %>%
    dplyr::rename(
      ENSEMBLID = ensembl_gene_id
    )
  
  message("Annotation completed: ", nrow(ann_clean), " genes annotated.")
  
  # ---------------------------------------------------------------------------
  # Bundle the annotation with provenance metadata
  # ---------------------------------------------------------------------------
  list(
    annotation       = ann_clean,
    ensembl_version  = ensembl_version,
    query_date       = format(Sys.Date()),
    mirror_used      = mirror_used,
    dataset          = biomart_dataset,
    n_genes_queried  = length(gene_ids),
    n_genes_returned = nrow(ann_clean)
  )
}


# =============================================================================
# create_gene_annotation_file()
#
# Backward-compatible wrapper. The cache RDS now contains the bundle
# (list with $annotation + provenance), but downstream code that does
# readRDS(file)$annotation, OR readRDS(file) when the file is a legacy
# data.frame, both keep working — see the matching change in _targets.R.
# =============================================================================

create_gene_annotation_file <- function(gene_ids, species, host, outfile) {
  
  # Skip Ensembl entirely if the annotation file is already on disk.
  # Delete the file manually to force a fresh query
  # (e.g. after a major Ensembl release).
  if (file.exists(outfile)) {
    message("Annotation file already exists, skipping Ensembl query: ", outfile)
    return(outfile)
  }
  
  bundle <- build_gene_annotation(
    gene_ids        = gene_ids,
    biomart_dataset = species,
    biomart_host    = host
  )
  
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bundle, outfile)
  
  message(sprintf(
    "Saved annotation file: %s\n  Ensembl release %s, queried %s, %d / %d genes",
    outfile,
    bundle$ensembl_version,
    bundle$query_date,
    bundle$n_genes_returned,
    bundle$n_genes_queried
  ))
  
  outfile
}

# --------------------------------------------------------------
# Knee plot
# --------------------------------------------------------------

# =============================================================================
# WELL-LEVEL QC  (barcode-rank / knee based)
# =============================================================================

#' Per-well QC metrics + keep/drop decision for every plate.
#'
#' For each plate in `combined_counts` (columns are "<plate_well>_<well_index>",
#' plus the fixed "gene"/"sum" columns) this computes total UMIs and the number
#' of genes detected per well, then flags wells to drop. The UMI cut is an
#' absolute floor (`min_umis`) that a data-driven barcode-rank threshold can
#' only make STRICTER, never weaker -- so a poorly fit spline on a plate's small,
#' fixed set of wells can never silently drop everything. `min_genes` is a second
#' independent floor. This is the same barcodeRanks() object the knee plots draw,
#' so the plots visually justify the cut.
#'
#' @param combined_counts named list of per-plate count data.frames
#' @param min_umis  absolute minimum total UMIs to keep a well
#' @param min_genes absolute minimum genes detected to keep a well
#' @param min_frac_median  relative floor: drop wells below this fraction of the
#'                  plate's median UMIs. 0 (default) disables it. This is the
#'                  robust lever for near-dead wells (e.g. a cytotoxic compound
#'                  that leaves a well at a few % of plate depth) that sit above
#'                  a fixed UMI floor. Applied per plate.
#' @param method    "fixed" (default): apply only the absolute + relative
#'                  floors. "inflection" / "knee": additionally raise the cut to
#'                  the DropletUtils::barcodeRanks threshold when it exceeds the
#'                  floors. barcodeRanks was designed for droplet data, where
#'                  thousands of empty barcodes sit below the real cells; on a
#'                  DRUG-seq plate every well contains cells, so on a clean
#'                  plate the inflection can land INSIDE the real-well
#'                  distribution and drop good wells. Because the cut is the
#'                  strictest active floor, that mistake is in the harmful
#'                  direction. Inspect the knee plots before opting in.
#' @param out_csv   optional path to write the full per-well QC table
#' @return tibble: plate_barcode, well, plate_well, well_index, umis, n_genes,
#'                 knee, inflection, umi_cut, gene_cut, keep, reason
compute_well_qc <- function(combined_counts,
                            min_umis        = 1000,
                            min_genes       = 200,
                            min_frac_median = 0,
                            method          = c("fixed", "inflection", "knee"),
                            out_csv         = NULL) {
  method     <- match.arg(method)
  fixed_cols <- c("gene", "sum")
  
  per_plate <- purrr::imap(combined_counts, function(plate_df, plate_name) {
    mat   <- as.matrix(plate_df[, !(colnames(plate_df) %in% fixed_cols), drop = FALSE])
    umis  <- colSums(mat)
    ngene <- colSums(mat > 0)
    
    # Barcode-rank thresholds (guarded: a plate has few, fixed wells)
    knee_thr <- NA_real_; infl_thr <- NA_real_
    br <- tryCatch(DropletUtils::barcodeRanks(mat), error = function(e) NULL)
    if (!is.null(br)) {
      knee_thr <- suppressWarnings(as.numeric(metadata(br)$knee))
      infl_thr <- suppressWarnings(as.numeric(metadata(br)$inflection))
    }
    data_thr <- switch(method,
                       inflection = infl_thr,
                       knee       = knee_thr,
                       fixed      = NA_real_)
    rel_thr <- if (min_frac_median > 0)
      min_frac_median * stats::median(umis) else -Inf
    # Effective cut = the strictest active floor (absolute, relative, knee).
    umi_cut <- max(c(min_umis, rel_thr,
                     if (is.finite(data_thr)) data_thr else -Inf))
    
    tibble::tibble(
      plate_barcode = plate_name,
      well          = names(umis),
      plate_well    = sub("_.*", "", names(umis)),
      well_index    = sub("^[^_]*_", "", names(umis)),
      umis          = as.numeric(umis),
      n_genes       = as.integer(ngene),
      knee          = knee_thr,
      inflection    = infl_thr,
      umi_cut       = umi_cut,
      gene_cut      = min_genes
    ) %>%
      dplyr::mutate(
        keep   = .data$umis >= .data$umi_cut & .data$n_genes >= .data$gene_cut,
        reason = dplyr::case_when(
          .data$keep                       ~ "pass",
          .data$umis    < .data$umi_cut    ~ "low_umis",
          .data$n_genes < .data$gene_cut   ~ "low_genes",
          TRUE                             ~ "fail"
        )
      )
  })
  
  qc     <- dplyr::bind_rows(per_plate)
  n_drop <- sum(!qc$keep)
  message(sprintf(
    "Well QC (method = %s): dropped %d / %d wells (%.1f%%)  [min_umis = %s, min_genes = %s]",
    method, n_drop, nrow(qc), 100 * n_drop / max(nrow(qc), 1),
    format(min_umis, big.mark = ","), format(min_genes, big.mark = ",")))
  if (n_drop) {
    drp <- qc[!qc$keep, ]
    message("  Dropped wells:\n",
            paste0("    ", drp$plate_barcode, "  ", drp$plate_well, "  ",
                   drp$well_index, "  umis=", round(drp$umis),
                   "  genes=", drp$n_genes, "  (", drp$reason, ")",
                   collapse = "\n"))
  }
  if (!is.null(out_csv)) {
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(qc, out_csv, row.names = FALSE)
  }
  qc
}


#' Drop QC-failed wells from a plate layout.
#'
#' Joins the per-well keep flag from compute_well_qc() onto the plate layout by
#' (plate_barcode, plate_well, well_index) and returns only the passing rows.
#' Wells present in the layout but absent from the QC table (i.e. no counts at
#' all) are dropped. Downstream DE selects wells via this layout, so filtering
#' here removes failed wells from every contrast without touching the DE code.
apply_well_qc_to_layout <- function(plate_layout, well_qc) {
  key_layout <- paste(plate_layout$plate_barcode, plate_layout$plate_well,
                      plate_layout$well_index, sep = "|")
  key_qc     <- paste(well_qc$plate_barcode, well_qc$plate_well,
                      well_qc$well_index, sep = "|")
  keep_lookup <- stats::setNames(well_qc$keep, key_qc)
  keep_flag   <- keep_lookup[key_layout]
  keep_flag[is.na(keep_flag)] <- FALSE
  out <- plate_layout[keep_flag, , drop = FALSE]
  message(sprintf("plate_layout_qc: %d / %d wells retained after QC",
                  nrow(out), nrow(plate_layout)))
  out
}


#' Barcode-rank ("knee") plot for one plate, annotated with the QC decision.
#'
#' Draws total UMIs vs well rank (log-y). When `well_qc` is supplied, dropped
#' wells are coloured red and the applied UMI cut is drawn as a solid line, so
#' the plot documents exactly which wells the pipeline removed and why.
make_knee_plot_single <- function(plate_name, plate_df, outdir, well_qc = NULL) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outfile      <- file.path(outdir, paste0(plate_name, "_knee_plot.png"))
  count_matrix <- as.matrix(plate_df[, !(colnames(plate_df) %in% c("gene", "sum")), drop = FALSE])
  br           <- DropletUtils::barcodeRanks(count_matrix)
  
  keep_col <- rep("black", length(br$total))
  umi_cut  <- NA_real_
  n_drop   <- 0L
  if (!is.null(well_qc)) {
    qcp  <- well_qc[well_qc$plate_barcode == plate_name, ]
    kmap <- stats::setNames(qcp$keep, qcp$well)
    kept <- kmap[rownames(br)]
    keep_col <- ifelse(is.na(kept) | kept, "black", "firebrick")
    umi_cut  <- suppressWarnings(max(qcp$umi_cut, na.rm = TRUE))
    n_drop   <- sum(!qcp$keep)
  }
  
  png(outfile, width = 800, height = 600)
  plot(br$rank, br$total, log = "y", pch = 16, col = keep_col,
       xlab = "Well rank (by total UMIs)", ylab = "Total UMIs", main = plate_name,
       sub = if (!is.null(well_qc))
         sprintf("QC: %d well(s) dropped below %s UMIs (red)",
                 n_drop, format(round(umi_cut), big.mark = ",")) else "")
  o <- order(br$rank)
  lines(br$rank[o], br$fitted[o], col = "grey40")
  abline(h = metadata(br)$knee,       col = "dodgerblue",  lty = 2)
  abline(h = metadata(br)$inflection, col = "forestgreen", lty = 2)
  if (is.finite(umi_cut)) abline(h = umi_cut, col = "firebrick", lty = 1, lwd = 2)
  legend("bottomleft", bty = "n",
         lty = c(2, 2, 1), lwd = c(1, 1, 2),
         col = c("dodgerblue", "forestgreen", "firebrick"),
         legend = c("knee", "inflection", "applied QC cut"))
  dev.off()
  outfile
}


make_all_knee_plots <- function(combined_counts, well_qc = NULL, outdir) {
  purrr::imap_chr(
    .x = combined_counts,
    .f = function(plate_df, plate_name) {
      make_knee_plot_single(
        plate_name = plate_name,
        plate_df   = plate_df,
        outdir     = outdir,
        well_qc    = well_qc
      )
      file.path(outdir, paste0(plate_name, "_knee_plot.png"))
    }
  )
}


# --------------------------------------------------------------
# Complete heatmaps
# --------------------------------------------------------------

#' Render a static plate-heatmap thumbnail without a browser.
#'
#' The interactive plate heatmap loads Plotly from a remote CDN, so the old
#' webshot2 screenshot produced a blank/failed PNG whenever the rendering
#' machine could not reach cdn.plot.ly (e.g. an air-gapped container) — which
#' is why thumbnails stopped appearing. This draws the same count grid as a
#' plain ggplot2 tile map straight to PNG: no headless Chrome, no network, no
#' render delay, and the result is identical on every machine.
render_plate_thumbnail <- function(counts_matrix, row_levels, col_levels,
                                   title, out_png,
                                   width_px = 900, height_px = 600,
                                   palette  = c("#ffffff", "#ffc8c8", "#ff7878", "#b40000"),
                                   legend_name = "Counts") {
  mat <- counts_matrix
  dimnames(mat) <- list(row = as.character(row_levels),
                        col = as.character(col_levels))
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("row", "col", "counts")
  df$counts <- suppressWarnings(as.numeric(df$counts))
  # Rows top-to-bottom (A at the top), columns left-to-right.
  df$row <- factor(df$row, levels = rev(as.character(row_levels)))
  df$col <- factor(df$col, levels = as.character(col_levels))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = counts)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradientn(
      colours  = palette,
      na.value = "#f0f0f0",
      name     = legend_name
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text   = ggplot2::element_text(face = "bold"),
      axis.ticks  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "right"
    )
  
  ggplot2::ggsave(
    out_png, p,
    width  = width_px  / 100,
    height = height_px / 100,
    dpi    = 100,
    bg     = "white"
  )
  invisible(out_png)
}

make_plate_heatmap <- function(plate_name, df_plate, meta,
                               save_heatmap     = NULL,
                               meta_experiment  = NULL,
                               plate_format     = 96,
                               write_thumbnail  = TRUE,
                               thumbnail_width  = 900,
                               thumbnail_height = 600) {
  if (is.null(df_plate)) stop("Plate not found in df_list: ", plate_name)
  
  header_html <- if (!is.null(meta_experiment)) {
    build_simple_header(meta_experiment)
  } else {
    ""
  }
  
  # --- Compute per-well stats ---
  df_plate <- df_plate %>% dplyr::select(-c(gene, sum))
  df_plate <- as.matrix(df_plate)
  
  counts_per_well <- round(colSums(df_plate)) # total read counts
  n_detected <- colSums(df_plate > 0) # number of genes with counts > 0
  
  # --- Subset meta ---
  meta_subset <- meta %>% dplyr::filter(plate_barcode == plate_name)
  
  # --- Parse column names ---
  col_info <- tibble::tibble(colname = colnames(df_plate)) %>%
    tidyr::separate(colname, into = c("plate_well", "well_index"), 
                    sep = "_") %>%
    dplyr::mutate(well_index = as.character(well_index))
  
  # --- Join meta ---
  well_info <- col_info %>%
    dplyr::left_join(meta_subset,
                     by = c("well_index", "plate_well")
    ) %>%
    dplyr::distinct(plate_well, .keep_all = TRUE) %>%   # remove duplicates
    dplyr::mutate(
      counts = counts_per_well[colnames(df_plate)],
      n_genes = n_detected[colnames(df_plate)],
      row_id = stringr::str_extract(plate_well, "^[A-Z]+"),
      column_id = as.numeric(stringr::str_extract(plate_well, "\\d+"))
    )
  
  # --- Determine plate geometry ---
  plate_dims <- switch(
    as.character(plate_format),
    "96"   = list(rows = LETTERS[1:8],  cols = 1:12),
    "384"  = list(rows = LETTERS[1:16], cols = 1:24),
    "1536" = list(rows = LETTERS[1:32], cols = 1:48),
    stop("Unsupported plate format: ", plate_format)
  )
  
  row_levels <- plate_dims$rows
  col_levels <- plate_dims$cols
  
  n_rows <- length(row_levels)
  n_cols <- length(col_levels)
  
  cmp_matrix <- matrix("",
                       nrow = n_rows,
                       ncol = n_cols,
                       dimnames = list(row_levels, col_levels)
  )
  
  counts_matrix <- matrix(NA,
                          nrow = n_rows,
                          ncol = n_cols,
                          dimnames = list(row_levels, col_levels)
  )
  
  hover_matrix <- matrix("",
                         nrow = n_rows,
                         ncol = n_cols,
                         dimnames = list(row_levels, col_levels)
  )
  
  dose_lab_matrix <- matrix("",       nrow = n_rows, ncol = n_cols,
                            dimnames = list(row_levels, col_levels))
  dose_key_matrix <- matrix(NA_real_, nrow = n_rows, ncol = n_cols,
                            dimnames = list(row_levels, col_levels))
  wi_dose_lab <- plate_dose_label(well_info$treatment_dose, well_info$treatment_dose_unit)
  wi_dose_key <- plate_dose_sort_key(well_info$treatment_dose, well_info$treatment_dose_unit)
  
  for (i in seq_len(nrow(well_info))) {
    r <- match(well_info$row_id[i], row_levels)
    c <- match(well_info$column_id[i], col_levels)
    
    cmp_matrix[r, c] <- well_info$treatment[i]
    dose_lab_matrix[r, c] <- wi_dose_lab[i]
    dose_key_matrix[r, c] <- wi_dose_key[i]
    counts_matrix[r, c] <- well_info$counts[i]
    
    hover_matrix[r, c] <- paste0(
      "Sample: ", well_info$well_index[i], "<br>",
      "Compound: ", well_info$treatment[i], "<br>",
      "Conc: ", well_info$treatment_dose[i], " ",
      well_info$treatment_dose_unit[i], "<br>",
      "Counts: ", well_info$counts[i], "<br>",
      "Detected genes: ", well_info$n_genes[i]
    )
  }
  
  # --- JSON encode ---
  cmp_json <- jsonlite::toJSON(cmp_matrix, auto_unbox = TRUE)
  counts_json <- jsonlite::toJSON(counts_matrix, auto_unbox = TRUE, na = "null")
  hover_json <- jsonlite::toJSON(hover_matrix, auto_unbox = TRUE)
  dose_lab_json <- jsonlite::toJSON(dose_lab_matrix, auto_unbox = TRUE)
  dose_key_json <- jsonlite::toJSON(dose_key_matrix, auto_unbox = TRUE, na = "null", digits = NA)
  dose_controls <- plate_dose_controls_html()     # defined in 04_plot_functions.R
  dose_js       <- plate_dose_highlight_js("heatmap", show_labels = n_rows <= 16)
  
  # --- HTML output ---
  html <- glue::glue('
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>{plate_format}-well plate heatmap — {plate_name}</title>
  <!-- Pinned: "plotly-latest" is frozen at 1.58 and unversioned. Viewing this
       page needs network access to the CDN; the PNG twin next to it does not. -->
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
</head>
<body>

<div style="margin:20px;">

  {header_html}

  <div id="controls" style="margin-bottom:10px;">
    <button id="clearBtn">Clear selection</button>
    {dose_controls}
  </div>

  <div id="heatmap" style="width:100%;height:80vh;min-height:600px;"></div>

</div>


<script>
const compoundMatrix = {cmp_json};
const countsMatrix   = {counts_json};
const hoverMatrix    = {hover_json};
const doseLabelMatrix = {dose_lab_json};
const doseKeyMatrix   = {dose_key_json};

const rows = compoundMatrix.length;
const cols = compoundMatrix[0].length;


const xVals = Array.from({{length: cols}}, (_, i) => i + 1);
const yVals = Array.from({{length: rows}}, (_, i) => i + 1);
const xTickText = xVals.map(v => String(v));
const yTickText = Array.from({{length: rows}}, (_,i) => String.fromCharCode(65+i));

const baseTrace = {{
  z: countsMatrix,
  x: xVals,
  y: yVals,
  text: hoverMatrix,
  customdata: compoundMatrix,
  hovertemplate: "%{{text}}<extra></extra>",
  type: "heatmap",
  colorscale: [
   [0.00, "rgb(255,255,255)"],
   [0.33, "rgb(255,200,200)"],
   [0.66, "rgb(255,120,120)"],
   [1.00, "rgb(180,0,0)"]
  ],
  zauto: true,
  showscale: true,
  colorbar: {{ title: "Counts" }}
}};

const layout = {{
  title: "{meta_experiment$title} – Plate {plate_name}",
  font: {{ family: "Helvetica Neue, Helvetica, Arial, sans-serif" }},
  xaxis: {{
    title: "",
    tickmode: "array",
    tickvals: xVals,
    ticktext: xTickText,
    side: "top"
  }},
  yaxis: {{
    title: "",
    tickmode: "array",
    tickvals: yVals,
    ticktext: yTickText,
    autorange: "reversed"
  }},
  shapes: []
}};

Plotly.newPlot("heatmap", [baseTrace], layout);

{dose_js}
</script>
</body>
</html>
')
  
  # --- Save HTML ---
  if (!is.null(save_heatmap)) {
    file_path <- file.path(save_heatmap, glue::glue("plate_{plate_name}.html"))
    writeLines(html, file_path)
    message("✅ Heatmap saved to ", file_path, ". Open it in your browser.")
    
    # --- PNG thumbnail (offline, no browser) -----------------------------
    # The report's companion-image lookup checks for <stem>.png alongside
    # <stem>.html, so the thumbnail lands at plate_<plate_name>.png in the
    # same directory. The interactive HTML above loads Plotly from a remote
    # CDN (cdn.plot.ly); screenshotting it with webshot2 produced a blank PNG
    # whenever Chrome could not reach that CDN. We instead redraw the same
    # count grid directly in R, which always works and is machine-independent.
    # Failure is non-fatal — the report falls back to its generic HTML icon.
    if (isTRUE(write_thumbnail)) {
      png_path <- file.path(
        save_heatmap, glue::glue("plate_{plate_name}.png")
      )
      tryCatch({
        render_plate_thumbnail(
          counts_matrix = counts_matrix,
          row_levels    = row_levels,
          col_levels    = col_levels,
          title         = if (!is.null(meta_experiment))
            glue::glue("{meta_experiment$title} \u2013 Plate {plate_name}")
          else
            glue::glue("Plate {plate_name}"),
          out_png       = png_path,
          width_px      = thumbnail_width,
          height_px     = thumbnail_height
        )
        size_kb <- if (file.exists(png_path)) file.info(png_path)$size / 1024 else 0
        message(sprintf(
          "  Plate heatmap thumbnail saved: %s (%.0f KB)",
          basename(png_path), size_kb
        ))
      }, error = function(e) {
        message("  Plate heatmap thumbnail failed (",
                conditionMessage(e), ") - skipping PNG for ", plate_name)
      })
    }
  }
}

make_all_heatmaps <- function(combined_counts, plate_layout, meta_experiment, outdir, plate_format) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE) 
  purrr::imap_chr(combined_counts, function(df_plate, plate_name) {
    make_plate_heatmap(
      plate_name      = plate_name,
      df_plate        = df_plate,
      meta            = plate_layout,
      meta_experiment = meta_experiment,
      plate_format = plate_format,
      save_heatmap    = outdir
    )
    
    file.path(outdir, glue::glue("plate_{plate_name}.html"))
  })
}