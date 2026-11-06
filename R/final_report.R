#' 
# R/final_report.R
# -----------------------------------------------------------------------------
# Build a single-page HTML index (results_dir/final_report_<run_id>.html)
# linking to every output produced by the pipeline. Runs as the last target so
# it only fires after all upstream outputs are ready.
#
# Design goals:
#   * Self-contained: relies only on base R + fs + htmltools (already in deps).
#   * Portable: all links are RELATIVE to results_dir, so you can zip the
#     folder and send it, or mount it through nginx unchanged.
#   * Tolerant: missing directories are silently skipped.
#   * Cheap: images are referenced by <img src=> (browser lazy-loads them),
#     not inlined as base64 — the index stays tiny even with hundreds of SVGs.
# -----------------------------------------------------------------------------

build_final_report <- function(results_dir,
                               config,
                               run_id,
                               padj,
                               lfc,
                               meta_experiment,
                               plate_layout,
                               de_results = NULL,
                               best_dmso  = NULL,
                               paths,
                               ...) {
  
  stopifnot(dir.exists(results_dir))
  
  # -- 1. Scan each output directory ------------------------------------------
  
  sections <- lapply(names(paths), function(title) {
    entry <- paths[[title]]
    # Accept either a plain directory string or list(dir = "...", list_only = TRUE, exclude_ext = "html")
    if (is.list(entry)) {
      dir         <- entry$dir
      list_only   <- isTRUE(entry$list_only)
      exclude_ext <- if (!is.null(entry$exclude_ext)) tolower(entry$exclude_ext) else character(0)
    } else {
      dir         <- entry
      list_only   <- FALSE
      exclude_ext <- character(0)
    }
    if (is.null(dir) || !dir.exists(dir)) return(NULL)
    
    files <- list.files(
      dir, recursive = TRUE, full.names = TRUE,
      include.dirs = FALSE, no.. = TRUE
    )
    files <- files[!grepl("(^|/)[^/]+_files(/|$)", files)]
    files <- files[!grepl("^final_report_.*\\.html$", basename(files))]
    files <- files[!grepl("^widget.*\\.html$", basename(files))]
    
    # Hide companion thumbnails: a PNG/SVG that exists alongside an HTML
    # of the same stem is rendered as the HTML's preview by render_card(),
    # so showing it as a separate card would be redundant clutter. The
    # Volcano section uses a dedicated pair-aware renderer that needs both
    # files explicitly, so it's exempt.
    if (!identical(title, "Volcano plots")) {
      stems <- tools::file_path_sans_ext(files)
      exts  <- tolower(tools::file_ext(files))
      html_stems <- stems[exts %in% c("html", "htm")]
      keep <- !(stems %in% html_stems & exts %in% c("png", "svg", "jpg", "jpeg"))
      files <- files[keep]
    }
    
    if (!length(files)) return(NULL)
    
    info <- file.info(files)
    df <- data.frame(
      name     = basename(files),
      rel_path = fs::path_rel(files, results_dir),
      ext      = tolower(tools::file_ext(files)),
      size     = info$size,
      mtime    = info$mtime,
      stringsAsFactors = FALSE
    )
    df <- df[order(df$name), , drop = FALSE]
    
    # Exclude specific extensions (e.g. html files that have their own section)
    if (length(exclude_ext)) {
      df <- df[!(df$ext %in% exclude_ext), , drop = FALSE]
      if (!nrow(df)) return(NULL)
    }
    
    list(
      title     = title,
      rel_dir   = fs::path_rel(dir, results_dir),
      files     = df,
      list_only = list_only
    )
  })
  sections <- sections[!vapply(sections, is.null, logical(1))]
  
  # -- 2. Helpers -------------------------------------------------------------
  
  fmt_size <- function(bytes) {
    if (is.na(bytes)) return("\u2014")
    if (bytes < 1024)       sprintf("%d B",   as.integer(bytes))
    else if (bytes < 1024^2) sprintf("%.1f KB", bytes / 1024)
    else if (bytes < 1024^3) sprintf("%.1f MB", bytes / 1024^2)
    else                     sprintf("%.2f GB", bytes / 1024^3)
  }
  
  safe_id <- function(s) gsub("[^A-Za-z0-9]+", "_", s)
  esc     <- htmltools::htmlEscape
  urlenc  <- function(p) utils::URLencode(p, reserved = FALSE)
  
  # Small inline download glyph reused by every card's download control.
  DL_ICON <- paste0(
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" ',
    'stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">',
    '<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>',
    '<polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>'
  )
  
  IMG_EXT  <- c("svg", "png", "jpg", "jpeg", "gif", "webp")
  HTML_EXT <- c("html", "htm")
  DATA_EXT <- c("csv", "tsv", "rds", "rdata", "xlsx", "parquet")
  PDF_EXT  <- c("pdf")
  
  # -- 3. Render a visual card (images, interactive HTML, PDFs) ----------------
  
  render_card <- function(fr) {
    ext  <- fr$ext
    href <- urlenc(fr$rel_path)
    
    thumb <- if (ext %in% IMG_EXT) {
      # -- Image file: show directly as thumbnail ---
      sprintf('<img data-src="%s" loading="lazy" alt=""/>', href)
      
    } else if (ext %in% HTML_EXT) {
      # -- HTML file: look for a companion PNG/SVG thumbnail ---
      base_stem <- tools::file_path_sans_ext(fr$rel_path)
      companion <- NULL
      for (img_ext in c("png", "svg", "jpg", "jpeg")) {
        candidate <- paste0(base_stem, ".", img_ext)
        if (file.exists(file.path(results_dir, candidate))) {
          companion <- candidate
          break
        }
      }
      
      if (!is.null(companion)) {
        # Use the companion image as a clickable thumbnail
        sprintf('<img data-src="%s" loading="lazy" alt=""/>', urlenc(companion))
      } else {
        # No companion image — styled HTML icon (no iframe)
        sprintf(
          '<div class="html-thumb">
             <svg class="html-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
               <path d="M13 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V9z"/>
               <polyline points="13 2 13 9 20 9"/>
             </svg>
             <span class="html-label">Interactive</span>
             <span class="html-hint">Click to open</span>
           </div>'
        )
      }
      
    } else if (ext %in% PDF_EXT) {
      '<div class="icon">\U0001F4C4<div class="icon-label">PDF</div></div>'
    } else {
      sprintf(
        '<div class="icon">\U0001F4E6<div class="icon-label">%s</div></div>',
        toupper(ext)
      )
    }
    
    sprintf(
      '<div class="card"%s>
         <a class="card-body" href="%s" target="_blank" rel="noopener">
           <div class="thumb">%s</div>
           <div class="meta">
             <div class="name" title="%s">%s</div>
             <div class="sub">%s &middot; %s</div>
           </div>
         </a>
         <a class="card-dl" href="%s" download="%s" title="Download %s" aria-label="Download %s">%s</a>
       </div>',
      compound_attr(fr$name),
      href, thumb, esc(fr$name), esc(fr$name),
      toupper(ext), fmt_size(fr$size),
      href, esc(fr$name), esc(fr$name), esc(fr$name), DL_ICON
    )
  }
  
  # -- 4. Render a compact data-file link row ---------------------------------
  
  render_data_link <- function(fr) {
    href <- urlenc(fr$rel_path)
    ext  <- toupper(fr$ext)
    sprintf(
      '<a class="data-link"%s href="%s" target="_blank" rel="noopener">
         <span class="data-ext">%s</span>
         <span class="data-name" title="%s">%s</span>
         <span class="data-size">%s</span>
       </a>',
      compound_attr(fr$name),
      href, ext, esc(fr$name), esc(fr$name), fmt_size(fr$size)
    )
  }
  
  # -- 4b. Volcano-card renderer (paired SVG + HTML, in-place toggle) ---------
  #
  # For the Volcano plots section the SVG and HTML for one contrast are
  # twins: same stem, different extensions. We render a single card per
  # contrast that:
  #   * shows the SVG inline (large, scrollable) as the default view,
  #   * has a toggle button that swaps the SVG for an <iframe> of the HTML
  #     in-place — the interactive plotly version replaces the static SVG
  #     without leaving the report.
  # The toggle uses tiny inline JavaScript (no external dependencies) and
  # falls back gracefully: if JS is disabled, both an SVG <img> link and
  # a direct HTML link remain clickable.
  
  # -- 4b. Volcano renderers ---------------------------------------------------
  #
  # Two distinct sources write to the Volcano plots directory:
  #
  #   plot_volcano()        → {contrast}.svg  +  {contrast}.html  (always paired)
  #   volcano_plot_static() → volcano_{treatment}_batch{b}.svg    (SVG only, prefixed)
  #
  # The first group gets the original SVG ↔ Interactive toggle card (one card
  # per contrast, the viewer swaps between the static SVG and a plotly iframe).
  # The second group is rendered in a separate "Static overview" row below,
  # each card showing the SVG with a "Full size" link.
  #
  # Distinguishing the two: files from volcano_plot_static() always have a
  # basename starting with "volcano_".
  
  # Static card for plot_volcano() output. The SVG is shown as a compact,
  # static thumbnail; the interactive plotly version is reachable via a link
  # that opens in a new tab (no heavy in-report iframe). A download control
  # grabs the SVG directly.
  # Scroll performance: a volcano SVG holds one vector path per gene (tens of
  # thousands per plot), and the browser re-parses / re-rasterises that on the
  # main thread whenever the card is painted, which made the whole report stutter
  # while scrolling. The card therefore shows a lightweight raster instead:
  #   1. <stem>.png written next to the SVG by plot_volcano() /
  #      volcano_plot_static() is used when present (decoded off-thread, cheap);
  #   2. otherwise the SVG is loaded once and flattened into a <canvas> bitmap
  #      by rasterizeHeavyImages() in the page script (data-raster="1").
  # "Full size" and the download button always point at the vector SVG.
  volcano_preview <- function(svg_rel) {
    png_rel <- sub("\\.svg$", ".png", svg_rel, ignore.case = TRUE)
    if (!identical(png_rel, svg_rel) && file.exists(file.path(results_dir, png_rel)))
      list(src = urlenc(png_rel), attr = "")
    else
      list(src = urlenc(svg_rel), attr = ' data-raster="1"')
  }
  
  render_volcano_pair <- function(svg_fr, html_fr, card_id) {
    svg_href  <- urlenc(svg_fr$rel_path)
    prev      <- volcano_preview(svg_fr$rel_path)
    html_href <- if (!is.null(html_fr)) urlenc(html_fr$rel_path) else NA_character_
    label     <- tools::file_path_sans_ext(svg_fr$name)   # no "volcano_" prefix here
    
    links <- paste0(
      sprintf('<a class="vc-btn vc-btn-link" href="%s" download="%s" title="Download SVG" aria-label="Download SVG">%s</a>',
              svg_href, esc(basename(svg_fr$rel_path)), DL_ICON),
      sprintf('<a class="vc-btn vc-btn-link" href="%s" target="_blank" rel="noopener" title="Open SVG at full size">Full size \u2197</a>',
              svg_href),
      if (!is.null(html_fr))
        sprintf('<a class="vc-btn vc-btn-link" href="%s" target="_blank" rel="noopener" title="Open interactive version">Interactive \u2197</a>',
                html_href)
      else ""
    )
    
    sprintf(
      '<div class="vc-card" id="%s"%s>
         <div class="vc-toolbar">
           <span class="vc-label" title="%s">%s</span>
           <div class="vc-actions-links">%s</div>
         </div>
         <div class="vc-viewport">
           <img class="vc-svg" data-src="%s"%s loading="lazy" decoding="async" alt="Volcano plot for %s"/>
         </div>
       </div>',
      card_id, compound_attr(svg_fr$name), esc(label), esc(label), links, prev$src, prev$attr, esc(label)
    )
  }
  
  # Plain static card for volcano_plot_static() output (SVG only) -------------
  render_volcano_static_card <- function(svg_fr) {
    svg_href <- urlenc(svg_fr$rel_path)
    prev     <- volcano_preview(svg_fr$rel_path)
    label    <- sub("^volcano_", "", tools::file_path_sans_ext(svg_fr$name))
    sprintf(
      '<div class="vc-card"%s>
         <div class="vc-toolbar">
           <span class="vc-label" title="%s">%s</span>
           <div class="vc-actions-links">
             <a class="vc-btn vc-btn-link" href="%s" download="%s" title="Download SVG" aria-label="Download SVG">%s</a>
             <a class="vc-btn vc-btn-link" href="%s" target="_blank"
                rel="noopener" title="Open SVG at full resolution">Full size \u2197</a>
           </div>
         </div>
         <div class="vc-viewport">
           <img class="vc-svg" data-src="%s"%s loading="lazy" decoding="async" alt="%s"/>
         </div>
       </div>',
      compound_attr(svg_fr$name),
      esc(label), esc(label),
      svg_href, esc(basename(svg_fr$rel_path)), DL_ICON,
      svg_href, prev$src, prev$attr, esc(label)
    )
  }
  
  render_volcano_section <- function(df_visual) {
    # Split by filename prefix: "volcano_*" = volcano_plot_static() output
    is_static_only <- grepl("^volcano_", df_visual$name)
    df_toggle <- df_visual[!is_static_only, , drop = FALSE]
    df_static  <- df_visual[ is_static_only, , drop = FALSE]
    
    # -- Toggle cards (plot_volcano paired SVG + HTML) ------------------------
    toggle_cards <- character(0)
    if (nrow(df_toggle)) {
      df_toggle$stem <- tools::file_path_sans_ext(df_toggle$rel_path)
      for (i in seq_along(unique(df_toggle$stem))) {
        s     <- unique(df_toggle$stem)[i]
        group <- df_toggle[df_toggle$stem == s, , drop = FALSE]
        svg_row  <- group[group$ext == "svg",  , drop = FALSE]
        html_row <- group[group$ext == "html", , drop = FALSE]
        if (!nrow(svg_row)) {
          if (nrow(html_row)) toggle_cards <- c(toggle_cards, render_card(html_row[1, ]))
          next
        }
        toggle_cards <- c(toggle_cards, render_volcano_pair(
          svg_fr   = svg_row[1, ],
          html_fr  = if (nrow(html_row)) html_row[1, ] else NULL,
          card_id  = sprintf("vc-%d", i)
        ))
      }
    }
    
    # -- Static overview cards (volcano_plot_static SVG only) -----------------
    static_cards <- character(0)
    if (nrow(df_static)) {
      static_svg <- df_static[df_static$ext == "svg", , drop = FALSE]
      for (i in seq_len(nrow(static_svg)))
        static_cards <- c(static_cards, render_volcano_static_card(static_svg[i, ]))
    }
    
    # Helper: wrap cards in a labelled vc-grid, capped to the first 3 with a
    # "Show all" toggle (the rest stay unrendered until revealed).
    make_vc_grid <- function(cards, row_label) {
      if (!length(cards)) return("")
      sprintf('<div class="vc-row-label">%s</div>%s',
              row_label, cap_grid(cards, "vc-grid", n_show = 3))
    }
    
    paste(
      make_vc_grid(toggle_cards, "Interactive overview"),
      make_vc_grid(static_cards, "Static overview"),
      sep = "\n"
    )
  }
  
  # -- 4c. DE-results table renderer (merged CSV + HTML twin per row) --------
  #
  # The DE pipeline writes:
  #   <path_DE>/DE_batch_<b>_plate_<p>.csv                — full DE table
  #   <path_DE>/HTML_reports/DE_batch_<b>_plate_<p>.html  — interactive view
  #   <path_DE>/DE_batch_<b>_collapsed_dupCor.csv         — collapsed result
  #
  # The merged renderer pairs each CSV with its HTML twin (if any), parses
  # batch + plate barcode from the filename, joins plate_layout for cell
  # line + treatments on that plate, and renders one row per CSV.
  
  parse_de_filename <- function(name) {
    base <- tools::file_path_sans_ext(basename(name))
    if (grepl("^DE_batch_(.+)_collapsed_dupCor_skipped$", base)) {
      m <- regmatches(base, regexec("^DE_batch_(.+)_collapsed_dupCor_skipped$", base))[[1]]
      list(batch = m[2], plate = NA_character_, kind = "collapsed_skipped")
    } else if (grepl("^DE_batch_(.+)_collapsed_dupCor$", base)) {
      m <- regmatches(base, regexec("^DE_batch_(.+)_collapsed_dupCor$", base))[[1]]
      list(batch = m[2], plate = NA_character_, kind = "collapsed")
    } else if (grepl("^DE_batch_(.+)_dose_trend_skipped$", base)) {
      m <- regmatches(base, regexec("^DE_batch_(.+)_dose_trend_skipped$", base))[[1]]
      list(batch = m[2], plate = NA_character_, kind = "dose_trend_skipped")
    } else if (grepl("^DE_batch_(.+)_dose_trend$", base)) {
      m <- regmatches(base, regexec("^DE_batch_(.+)_dose_trend$", base))[[1]]
      list(batch = m[2], plate = NA_character_, kind = "dose_trend")
    } else if (grepl("^DE_batch_(.+)_plate_(.+)_skipped$", base)) {
      m <- regmatches(base, regexec("^DE_batch_(.+)_plate_(.+)_skipped$", base))[[1]]
      list(batch = m[2], plate = m[3], kind = "perplate_skipped")
    } else if (grepl("^DE_batch_(.+)_plate_(.+)$", base)) {
      m <- regmatches(base, regexec("^DE_batch_(.+)_plate_(.+)$", base))[[1]]
      list(batch = m[2], plate = m[3], kind = "perplate")
    } else {
      list(batch = NA_character_, plate = NA_character_, kind = "other")
    }
  }
  
  summarise_plate_treatments <- function(plate_layout, plate_barcode_arg,
                                         batch_arg, ctrl) {
    if (is.null(plate_layout) || !is.data.frame(plate_layout)) {
      return(list(cell_line = "", treatments = character(0)))
    }
    needed <- c("plate_barcode", "batch", "treatment", "cell_line_name")
    if (!all(needed %in% names(plate_layout))) {
      return(list(cell_line = "", treatments = character(0)))
    }
    rows <- plate_layout[
      !is.na(plate_layout$plate_barcode) &
        plate_layout$plate_barcode == plate_barcode_arg &
        as.character(plate_layout$batch) == as.character(batch_arg), ,
      drop = FALSE
    ]
    if (!nrow(rows)) {
      return(list(cell_line = "", treatments = character(0)))
    }
    # Drop the control compound from the treatments list — uninformative
    treatments <- unique(rows$treatment[!is.na(rows$treatment) &
                                          rows$treatment != ctrl])
    list(
      cell_line  = unique(rows$cell_line_name)[1],
      treatments = sort(treatments)
    )
  }
  
  summarise_batch_treatments <- function(plate_layout, batch_arg, ctrl) {
    if (is.null(plate_layout) || !is.data.frame(plate_layout)) {
      return(list(cell_line = "", treatments = character(0), n_plates = 0))
    }
    needed <- c("plate_barcode", "batch", "treatment", "cell_line_name")
    if (!all(needed %in% names(plate_layout))) {
      return(list(cell_line = "", treatments = character(0), n_plates = 0))
    }
    rows <- plate_layout[
      !is.na(plate_layout$batch) &
        as.character(plate_layout$batch) == as.character(batch_arg), ,
      drop = FALSE
    ]
    if (!nrow(rows)) return(list(cell_line = "", treatments = character(0),
                                 n_plates = 0))
    treatments <- unique(rows$treatment[!is.na(rows$treatment) &
                                          rows$treatment != ctrl])
    list(
      cell_line  = unique(rows$cell_line_name)[1],
      treatments = sort(treatments),
      n_plates   = length(unique(rows$plate_barcode))
    )
  }
  
  # Per-compound stats for the tag tooltips: DE-gene count (at the report's
  # thresholds), genes tested, doses, and well/replicate count. Built once from
  # the collapsed DE table + plate_layout; NULL-safe if either is missing.
  compound_stats <- local({
    stats <- new.env(parent = emptyenv())
    if (!is.null(de_results) && is.data.frame(de_results) &&
        all(c("Sample", "treatment", "LFC", "padj") %in% names(de_results))) {
      sig <- de_results[!is.na(de_results$padj) & de_results$padj < padj &
                          !is.na(de_results$LFC) & abs(de_results$LFC) >= lfc, ]
      de_by_cmp    <- tapply(seq_len(nrow(sig)), sig$treatment, length)
      genes_by_cmp <- if ("n_genes_tested" %in% names(de_results))
        tapply(de_results$n_genes_tested, de_results$treatment,
               function(x) max(x, na.rm = TRUE)) else NULL
      doses_by_cmp <- tapply(de_results$Sample, de_results$treatment,
                             function(x) length(unique(x)))
      for (cmp in unique(de_results$treatment)) {
        n_wells <- NA_integer_
        if (!is.null(plate_layout) && "treatment" %in% names(plate_layout))
          n_wells <- sum(plate_layout$treatment == cmp, na.rm = TRUE)
        assign(cmp, list(
          n_de    = if (!is.null(de_by_cmp) && cmp %in% names(de_by_cmp)) as.integer(de_by_cmp[[cmp]]) else 0L,
          n_genes = if (!is.null(genes_by_cmp) && cmp %in% names(genes_by_cmp)) as.integer(genes_by_cmp[[cmp]]) else NA_integer_,
          n_doses = if (cmp %in% names(doses_by_cmp)) as.integer(doses_by_cmp[[cmp]]) else NA_integer_,
          n_wells = n_wells
        ), envir = stats)
      }
    }
    stats
  })
  
  # Per-well counts for a compound, optionally scoped to a batch and/or plate.
  # DE-gene / dose / genes-tested stats stay compound-wide (they come from the
  # collapsed fit, which pools plates); only the well/replicate count is scoped,
  # because that is the thing that differs between a dupCor batch row and a
  # single-plate row.
  wells_for <- function(cmp, batch = NULL, plate = NULL) {
    if (is.null(plate_layout) || !("treatment" %in% names(plate_layout)))
      return(list(n = NA_integer_, per_plate = NULL))
    keep <- plate_layout$treatment == cmp & !is.na(plate_layout$treatment)
    if (!is.null(batch) && "batch" %in% names(plate_layout))
      keep <- keep & as.character(plate_layout$batch) == as.character(batch)
    if (!is.null(plate) && "plate_barcode" %in% names(plate_layout))
      keep <- keep & plate_layout$plate_barcode == plate
    sub <- plate_layout[keep, , drop = FALSE]
    per_plate <- if (!is.null(plate) || !("plate_barcode" %in% names(sub)))
      NULL else as.integer(table(sub$plate_barcode))
    list(n = nrow(sub), per_plate = per_plate)
  }
  
  # Tooltip text for one compound pill (plain text, used in a title= attribute).
  # `batch`/`plate` scope the well count; without them it is experiment-wide.
  compound_tip <- function(cmp, batch = NULL, plate = NULL) {
    s <- if (exists(cmp, envir = compound_stats, inherits = FALSE))
      get(cmp, envir = compound_stats, inherits = FALSE) else
        list(n_de = NA, n_doses = NA, n_genes = NA)
    w <- wells_for(cmp, batch = batch, plate = plate)
    # Well/replicate phrasing: on a per-plate row just "N wells (on this plate)";
    # on a batch row "N wells (a/b/c per plate)" so the reader sees the split.
    wells_txt <- if (is.na(w$n) || w$n == 0) "" else if (!is.null(plate)) {
      paste0(" \u00b7 ", w$n, " wells on this plate")
    } else if (!is.null(w$per_plate) && length(w$per_plate) > 1) {
      paste0(" \u00b7 ", w$n, " wells (",
             paste(w$per_plate, collapse = "/"), " per plate)")
    } else {
      paste0(" \u00b7 ", w$n, " wells")
    }
    paste0(
      cmp, " \u2014 ",
      if (!is.na(s$n_de)) paste0(s$n_de, " DE genes") else "0 DE genes",
      if (!is.na(s$n_doses)) paste0(" \u00b7 ", s$n_doses, " doses") else "",
      wells_txt,
      if (!is.na(s$n_genes)) paste0(" \u00b7 ", s$n_genes, " genes tested") else ""
    )
  }
  
  # Render a treatment pill list that shows the first `n_show` pills and hides
  # the rest behind a "+N more" toggle. Each pill carries a stats tooltip scoped
  # to the row's batch/plate so per-plate rows report per-plate replicates.
  render_treatment_pills <- function(treatments, n_show = 4,
                                     batch = NULL, plate = NULL) {
    if (!length(treatments)) return("")
    pill <- function(t)
      sprintf('<span class="de-pill de-pill-trt" title="%s">%s</span>',
              esc(compound_tip(t, batch = batch, plate = plate)), esc(t))
    if (length(treatments) <= n_show) {
      return(paste(vapply(treatments, pill, character(1)), collapse = " "))
    }
    shown  <- paste(vapply(treatments[seq_len(n_show)], pill, character(1)), collapse = " ")
    hidden <- paste(vapply(treatments[(n_show + 1):length(treatments)], pill, character(1)), collapse = " ")
    n_more <- length(treatments) - n_show
    sprintf(
      '<span class="pill-group">%s<span class="pill-hidden" hidden>%s</span><button type="button" class="pill-more" data-open="0">+%d more</button></span>',
      shown, hidden, n_more
    )
  }
  
  render_de_results_section <- function(df_files, plate_layout) {
    df_files$stem <- tools::file_path_sans_ext(basename(df_files$rel_path))
    csv_rows  <- df_files[df_files$ext == "csv",  , drop = FALSE]
    html_rows <- df_files[df_files$ext == "html", , drop = FALSE]
    
    if (!nrow(csv_rows)) return("")
    
    ctrl <- tryCatch(config$analysis$control_compound, error = function(e) "DMSO")
    if (is.null(ctrl) || !nzchar(ctrl)) ctrl <- "DMSO"
    
    # The HTML view lists this contrast's significant hits (padj < padj,
    # |LFC| > lfc), sorted by significance, capped at MAX_HTML_ROWS. Count them
    # from the CSV so each HTML link can say how many genes it shows by default.
    # Only the padj/LFC columns are read (a 350 MB per-plate CSV stays cheap),
    # results are cached, and any failure yields NA (the count is simply omitted).
    MAX_HTML_ROWS <- 2000L
    .hg_cache <- new.env(parent = emptyenv())
    html_gene_count <- function(csv_rel) {
      if (!is.null(.hg_cache[[csv_rel]])) return(.hg_cache[[csv_rel]])
      p <- file.path(results_dir, csv_rel)
      n <- tryCatch({
        d <- if (requireNamespace("data.table", quietly = TRUE)) {
          data.table::fread(p, select = c("padj", "LFC"), showProgress = FALSE)
        } else {
          readr::read_csv(p, progress = FALSE,
                          col_types = readr::cols_only(padj = readr::col_double(),
                                                       LFC  = readr::col_double()))
        }
        sum(!is.na(d$padj) & d$padj < padj & !is.na(d$LFC) & abs(d$LFC) > lfc)
      }, error = function(e) NA_integer_)
      val <- if (is.na(n)) NA_integer_ else min(as.integer(n), MAX_HTML_ROWS)
      .hg_cache[[csv_rel]] <- val
      val
    }
    
    # Build one record per CSV: its analysis kind, a display label, batch (for
    # ordering) and the CSV/HTML cells. Rows are then grouped by kind under a
    # header row so dupCor, dose-trend and per-plate outputs are visually
    # separated instead of interleaved by filename.
    rec <- vector("list", nrow(csv_rows))
    for (i in seq_len(nrow(csv_rows))) {
      csv  <- csv_rows[i, ]
      twin <- html_rows[html_rows$stem == csv$stem, , drop = FALSE]
      info <- parse_de_filename(csv$name)
      
      desc_extra <- ""
      if (info$kind == "perplate" && !is.na(info$plate)) {
        meta_info <- summarise_plate_treatments(plate_layout, info$plate,
                                                info$batch, ctrl)
        desc_main <- sprintf("Batch %s &middot; Plate %s",
                             esc(info$batch), esc(info$plate))
        desc_extra <- paste(
          if (nzchar(meta_info$cell_line))
            sprintf('<span class="de-pill">%s</span>', esc(meta_info$cell_line))
          else "",
          render_treatment_pills(meta_info$treatments,
                                 batch = info$batch, plate = info$plate),
          collapse = ""
        )
      } else if (info$kind == "collapsed") {
        meta_info <- summarise_batch_treatments(plate_layout, info$batch, ctrl)
        desc_main <- sprintf(
          "Batch %s &middot; <strong>Collapsed dupCor</strong> (%d plates)",
          esc(info$batch), meta_info$n_plates
        )
        desc_extra <- paste(
          if (nzchar(meta_info$cell_line))
            sprintf('<span class="de-pill">%s</span>', esc(meta_info$cell_line))
          else "",
          render_treatment_pills(meta_info$treatments, batch = info$batch),
          collapse = ""
        )
      } else if (info$kind == "dose_trend") {
        desc_main <- sprintf("Batch %s &middot; dose-trend (per-compound log-dose slope)",
                             esc(info$batch))
      } else if (info$kind == "dose_trend_skipped") {
        desc_main <- sprintf("Batch %s &middot; dose-trend skipped contrasts",
                             esc(info$batch))
      } else if (info$kind == "collapsed_skipped") {
        desc_main <- sprintf("Batch %s &middot; collapsed dupCor skipped contrasts",
                             esc(info$batch))
      } else if (info$kind == "perplate_skipped" && !is.na(info$plate)) {
        desc_main <- sprintf("Batch %s &middot; Plate %s &middot; skipped contrasts",
                             esc(info$batch), esc(info$plate))
      } else {
        desc_main  <- esc(csv$name)
      }
      
      html_cell <- if (nrow(twin)) {
        ng      <- html_gene_count(csv$rel_path)
        ng_span <- if (is.na(ng)) "" else sprintf(
          '<span class="de-size">%s%s gene%s</span>',
          format(ng, big.mark = ","),
          if (ng >= MAX_HTML_ROWS) "+" else "",
          if (ng == 1L) "" else "s"
        )
        ng_title <- if (is.na(ng)) "interactive table (significant hits)"
        else sprintf("interactive table: %s significant gene%s shown (padj < %s, |LFC| > %s, max %s rows)",
                     format(ng, big.mark = ","), if (ng == 1L) "" else "s",
                     padj, lfc, format(MAX_HTML_ROWS, big.mark = ","))
        sprintf(
          '<a class="de-link" href="%s" target="_blank" rel="noopener" title="%s">
             <span class="de-ext de-ext-html">HTML</span>
             %s
             <span class="de-size">%s</span>
           </a>',
          urlenc(twin$rel_path[1]), esc(ng_title), ng_span, fmt_size(twin$size[1])
        )
      } else {
        '<span class="de-link de-link-empty">&mdash;</span>'
      }
      
      csv_cell <- sprintf(
        '<a class="de-link" href="%s" target="_blank" rel="noopener" title="%s">
           <span class="de-ext de-ext-csv">CSV</span>
           <span class="de-size">%s</span>
         </a>',
        urlenc(csv$rel_path), esc(csv$name), fmt_size(csv$size)
      )
      
      rec[[i]] <- list(
        kind = info$kind, batch = info$batch, plate = info$plate,
        row = sprintf(
          '<tr>
             <td class="de-desc">
               <div class="de-desc-main">%s</div>
               <div class="de-desc-pills">%s</div>
             </td>
             <td class="de-cell">%s</td>
             <td class="de-cell">%s</td>
           </tr>',
          desc_main, desc_extra, csv_cell, html_cell)
      )
    }
    
    # Group order + human labels. Skipped-contrast logs sit with their analysis.
    group_order <- c(
      collapsed          = "Collapsed dupCor (primary result)",
      collapsed_skipped  = "Collapsed dupCor (primary result)",
      dose_trend         = "Dose-trend",
      dose_trend_skipped = "Dose-trend",
      perplate           = "Per-plate (QC / consistency)",
      perplate_skipped   = "Per-plate (QC / consistency)",
      other              = "Other"
    )
    kinds   <- vapply(rec, function(r) r$kind, character(1))
    batches <- vapply(rec, function(r) as.character(r$batch), character(1))
    # stable sort: by group, then NUMERIC batch (so 2 sorts before 10, not after),
    # then original order. Non-numeric batches (should not occur) sort last.
    batch_num <- suppressWarnings(as.integer(batches))
    batch_num[is.na(batch_num)] <- .Machine$integer.max
    grp     <- match(kinds, names(group_order))
    ord     <- order(grp, batch_num, seq_along(rec))
    rec     <- rec[ord]
    
    rows_html <- character(0)
    last_label <- ""
    for (r in rec) {
      lbl <- group_order[[r$kind]]
      if (!identical(lbl, last_label)) {
        rows_html <- c(rows_html, sprintf(
          '<tr class="de-group"><td colspan="3">%s</td></tr>', esc(lbl)))
        last_label <- lbl
      }
      rows_html <- c(rows_html, r$row)
    }
    
    sprintf(
      '<table class="de-table">
         <thead>
           <tr>
             <th>Contrast group</th>
             <th>CSV (full results)</th>
             <th>Interactive HTML</th>
           </tr>
         </thead>
         <tbody>%s</tbody>
       </table>',
      paste(rows_html, collapse = "\n")
    )
  }
  
  # -- 5. Section descriptions ------------------------------------------------
  
  section_desc <- list(
    "Volcano plots" = paste(
      "Differential expression volcano plots. Each point is a gene:",
      "x-axis is log2 fold change, y-axis is -log10(adjusted p-value).",
      "Interactive per-contrast plots (EnhancedVolcano): red = passes both the",
      "padj and |LFC| thresholds, blue = padj only, green = |LFC| only, grey = neither.",
      "Static per-compound overview plots (one panel per dose): red = up,",
      "blue = down (both thresholds), black = not significant."
    ),
    "UpSet plot" = paste(
      "UpSet plots show the overlap of DE gene sets across treatments and doses.",
      "Each vertical bar represents a unique intersection of conditions sharing DE genes.",
      "Useful for identifying shared vs. compound-specific transcriptional responses."
    ),
    "DE results" = paste(
      "Differential expression results, one row per analysis.",
      "Each row links to the full results CSV (all genes tested) and the",
      "interactive HTML table (filtered to significant hits).",
      "Per-plate rows contain a single plate's contrasts; collapsed dupCor",
      "rows pool all plates in a batch and are the recommended primary result."
    ),
    "DE heatmaps" = paste(
      "Heatmaps of the top DE genes per compound.",
      "Rows are genes (ranked by adjusted p-value), columns are wells or dose groups",
      "(see aggregate_by), ordered by dose. Colour is the per-gene z-score of",
      "log2(UMI + 1) expression across the columns, not log2 fold change;",
      "with two or fewer columns the raw log2(UMI + 1) is shown instead.",
      "PNG for static view, HTML for interactive zoom."
    ),
    "UMAP plots" = paste(
      "UMAP dimensionality reduction of all wells based on gene expression profiles.",
      "Each point is a well. Colored by treatment, Leiden cluster, or DMSO QC category.",
      "Useful for identifying batch effects, treatment-driven clustering, and outlier wells."
    ),
    "Correlation plots" = paste(
      "Pairwise correlation plots (GGally) of log2 fold change across doses within a treatment.",
      "Lower triangle shows scatter plots, diagonal shows density, upper triangle shows Pearson correlation.",
      "Strong correlations across doses indicate consistent dose-response transcriptional programs."
    ),
    "Knee plots (QC)" = paste(
      "Barcode-rank (knee) plot per plate: every point is one well, ranked by",
      "total UMIs (log scale). Dashed lines mark the DropletUtils knee and",
      "inflection for reference. When well QC is enabled the applied UMI cut is",
      "drawn as a solid red line and dropped wells are coloured red; the full",
      "per-well table with the reason for each decision is well_qc.csv.",
      "There are no droplets in DRUG-seq - a low point is a shallow or dead well,",
      "not an empty barcode."
    ),
    "Plate heatmaps (QC)" = paste(
      "Plate-level heatmaps showing total UMI counts per well in their physical plate layout.",
      "Helps identify spatial artefacts, edge effects, or failed wells on each plate."
    ),
    "DE-gene count heatmaps" = paste(
      "Plate-level heatmaps coloured by the number of DE genes each well's",
      "compound &times; dose contrast produced (padj and |log2FC| at the report",
      "thresholds). DE is defined per contrast, so every well of a condition",
      "shares its value; control wells are blank. Reveals which plate positions",
      "and doses drove the response and highlights inactive or failed conditions.",
      "Click a well to highlight all wells of the same compound."
    ),
    "DMSO QC" = paste(
      "Selection of DMSO control wells. Each well is scored by counting genes",
      "and principal components on which it is a robust-z outlier; the two",
      "counts are standardised and summed (lower = more stable). The DMSO dose",
      "with the lowest median score is chosen. With top_n = \"all\" (default) every",
      "well below median + 2 MAD is kept as a control, readmitting the best-ranked",
      "flagged wells until each plate has min_keep (default 6); with a numeric top_n",
      "the top_n lowest-scoring wells per plate are kept. Wells above median + 2 MAD",
      "are listed as flagged.",
      "Only kept wells enter DE contrasts."
    )
  )
  
  # -- 6. Render a section (split visual vs data files) -----------------------
  
  # Cap a grid to the first `n_show` cards; a "Show all" button reveals the
  # rest. Every card stays in ONE grid container and the overflow cards carry
  # .grid-over (display:none) instead of living in a separate hidden block.
  # That is what lets the cap be recomputed over an arbitrary subset of the
  # cards - e.g. only those passing the compound filter, which can sit at any
  # position in the grid. A container-level collapse cannot express that.
  # Hidden cards are display:none, so their data-src images are never fetched
  # until revealed - this is what keeps the initial report light regardless of
  # card count. Once revealed, the IntersectionObserver still lazy-loads them
  # as they scroll into view.
  cap_grid <- function(cards, grid_class = "grid", n_show = 3) {
    n <- length(cards)
    if (n == 0) return("")
    single_cls <- if (n == 1L) " single-item" else ""
    if (n <= n_show) {
      return(sprintf('<div class="%s%s">%s</div>',
                     grid_class, single_cls, paste(cards, collapse = "\n")))
    }
    # Every card string opens with the outer <div class="card"/"vc-card" ...>,
    # so the first class attribute is the card's own (sub() = first match).
    over <- (n_show + 1):n
    cards[over] <- sub('class="', 'class="grid-over ', cards[over], fixed = TRUE)
    sprintf(
      '<div class="grid-capped" data-nshow="%d">
         <div class="%s">%s</div>
         <button type="button" class="grid-more" data-open="0" data-n="%d">Show all %d \u25be</button>
       </div>',
      n_show, grid_class, paste(cards, collapse = "\n"), n, n
    )
  }
  
  render_section <- function(sec) {
    list_only <- isTRUE(sec$list_only)
    
    # Merged DE results: one row per CSV with twin HTML + plate metadata.
    # Detected by section title — sec$files contains both CSVs and HTMLs
    # because the section's directory is path_DE (recursive scan picks up
    # path_DE/HTML_reports/*.html as well).
    if (identical(sec$title, "DE results")) {
      desc <- section_desc[[sec$title]]
      desc_html <- if (!is.null(desc)) {
        sprintf('<p class="section-desc">%s</p>', desc)
      } else ""
      
      table_html <- render_de_results_section(sec$files, plate_layout)
      
      return(sprintf(
        '<section id="%s">
           <div class="section-head">
             <h2>%s</h2>
             <span class="count">%d file%s</span>
             <code class="dir">%s/</code>
           </div>
           %s
           %s
         </section>',
        safe_id(sec$title), esc(sec$title),
        nrow(sec$files), if (nrow(sec$files) == 1) "" else "s",
        esc(sec$rel_dir),
        desc_html, table_html
      ))
    }
    
    if (list_only) {
      df_visual <- sec$files[0, , drop = FALSE]
      df_data   <- sec$files
    } else {
      df_visual <- sec$files[!(sec$files$ext %in% DATA_EXT), , drop = FALSE]
      df_data   <- sec$files[sec$files$ext %in% DATA_EXT, , drop = FALSE]
    }
    
    visual_html <- ""
    if (nrow(df_visual) > 0) {
      if (identical(sec$title, "Volcano plots")) {
        # render_volcano_section now produces its own sub-grids (Static /
        # Interactive rows), each already tagged with single-item when needed.
        visual_html <- render_volcano_section(df_visual)
      } else {
        cards      <- vapply(seq_len(nrow(df_visual)),
                             function(i) render_card(df_visual[i, ]),
                             character(1))
        # Small square thumbnails pack ~4 per row, so previewing a full first
        # row looks better than a 3-card partial one; wider cards stay at 3.
        n_preview <- if (sec$title %in% c("Plate heatmaps (QC)",
                                          "DE-gene count heatmaps",
                                          "Knee plots (QC)")) 4L else 3L
        visual_html <- cap_grid(cards, "grid", n_show = n_preview)
      }
    }
    
    data_html <- ""
    if (nrow(df_data) > 0) {
      links <- vapply(seq_len(nrow(df_data)),
                      function(i) render_data_link(df_data[i, ]),
                      character(1))
      data_html <- sprintf(
        '<div class="data-list-container">
           <div class="data-list-header">Data files</div>
           <div class="data-list">%s</div>
         </div>',
        paste(links, collapse = "\n")
      )
    }
    
    # Section description (if defined)
    desc <- section_desc[[sec$title]]
    desc_html <- if (!is.null(desc)) {
      sprintf('<p class="section-desc">%s</p>', desc)
    } else ""
    
    sprintf(
      '<section id="%s">
         <div class="section-head">
           <h2>%s</h2>
           <span class="count">%d file%s</span>
           <code class="dir">%s/</code>
         </div>
         %s
         %s
         %s
       </section>',
      safe_id(sec$title), esc(sec$title),
      nrow(sec$files), if (nrow(sec$files) == 1) "" else "s",
      esc(sec$rel_dir),
      desc_html, visual_html, data_html
    )
  }
  
  render_nav <- function(secs) {
    paste(vapply(secs, function(sec) {
      sprintf('<a href="#%s">%s<span class="n">%d</span></a>',
              safe_id(sec$title), esc(sec$title), nrow(sec$files))
    }, character(1)), collapse = "\n")
  }
  
  # -- 6. Derive experiment metadata for header -------------------------------
  
  title     <- if (!is.null(meta_experiment$title))       meta_experiment$title       else "ADEPT run"
  cell_line <- if (!is.null(meta_experiment$cell_line))   meta_experiment$cell_line   else ""
  citation  <- if (!is.null(meta_experiment$citation))    meta_experiment$citation    else ""
  description <- if (!is.null(meta_experiment$description)) trimws(meta_experiment$description) else ""
  
  # Extract details from plate_layout (the metadata table)
  treatments   <- character(0)
  compounds    <- character(0)
  batches      <- character(0)
  plates       <- character(0)
  doses        <- character(0)
  dose_unit    <- ""
  timepoints   <- character(0)
  n_wells      <- 0L
  
  if (!is.null(plate_layout) && is.data.frame(plate_layout)) {
    if ("treatment" %in% names(plate_layout))
      treatments <- sort(unique(plate_layout$treatment))
    if ("compound" %in% names(plate_layout))
      compounds  <- sort(unique(plate_layout$compound[!is.na(plate_layout$compound) & plate_layout$compound != ""]))
    if ("batch" %in% names(plate_layout))
      batches    <- sort(unique(as.character(plate_layout$batch)))
    if ("plate_barcode" %in% names(plate_layout))
      plates     <- sort(unique(plate_layout$plate_barcode))
    if ("treatment_dose" %in% names(plate_layout))
      doses      <- sort(unique(as.numeric(plate_layout$treatment_dose)))
    if ("treatment_dose_unit" %in% names(plate_layout))
      dose_unit  <- unique(plate_layout$treatment_dose_unit)[1]
    if ("timepoint" %in% names(plate_layout))
      timepoints <- sort(unique(as.character(plate_layout$timepoint)))
    n_wells <- nrow(plate_layout)
  }
  
  # Remove control compound from treatment/compound lists for display
  ctrl <- if (!is.null(config$analysis$control_compound)) config$analysis$control_compound else "DMSO"
  display_treatments <- setdiff(treatments, ctrl)
  
  # Build treatment -> compound mapping for display
  trt_compound_map <- character(0)
  if (!is.null(plate_layout) && all(c("treatment", "compound") %in% names(plate_layout))) {
    tc <- unique(plate_layout[, c("treatment", "compound"), drop = FALSE])
    tc <- tc[!is.na(tc$compound) & tc$compound != "" & tc$treatment != ctrl, , drop = FALSE]
    tc <- tc[order(tc$treatment), , drop = FALSE]
    trt_compound_map <- setNames(tc$compound, tc$treatment)
  }
  
  # -- 6b. Compound filter ----------------------------------------------------
  # Tag every compound-specific output with the treatment(s) whose name appears
  # in its filename, so a header control can show only the chosen compounds.
  # Treatment names are matched as whole tokens (\Q..\E quotes any regex-special
  # characters; the boundaries stop "Cmp_33" matching inside "Cmp_334").
  filter_treatments <- sort(display_treatments)
  
  trt_token <- function(name) {
    if (!length(filter_treatments) || is.null(name)) return(character(0))
    nm <- basename(as.character(name))
    if (is.na(nm)) return(character(0))
    keep <- vapply(filter_treatments, function(t) {
      pat <- paste0("(^|[^A-Za-z0-9])\\Q", t, "\\E([^A-Za-z0-9]|$)")
      grepl(pat, nm, perl = TRUE)
    }, logical(1))
    filter_treatments[keep]
  }
  
  # ' data-compounds="id1 id2"' attribute for a filename ("" if no match).
  compound_attr <- function(name) {
    tags <- trt_token(name)
    if (!length(tags)) return("")
    sprintf(' data-compounds="%s"',
            paste(vapply(tags, safe_id, character(1)), collapse = " "))
  }
  
  # Treatments that actually have >= 1 filterable file in this run.
  present_treatments <- {
    all_names <- unlist(lapply(sections, function(s) s$files$name),
                        use.names = FALSE)
    if (length(all_names))
      sort(unique(unlist(lapply(all_names, trt_token), use.names = FALSE)))
    else character(0)
  }
  
  # Filter bar: one chip per treatment in the experiment. Treatments with no
  # results this run are shown disabled, so the full roster is still visible.
  compound_filter_html <- ""
  if (length(present_treatments)) {
    chips <- vapply(filter_treatments, function(t) {
      has  <- t %in% present_treatments
      cls  <- if (has) "cmp-chip" else "cmp-chip cmp-chip-off"
      attr <- if (has) sprintf(' data-cmp="%s"', safe_id(t))
      else ' aria-disabled="true" title="No results for this compound in this run"'
      sub  <- if (t %in% names(trt_compound_map) &&
                  nzchar(trt_compound_map[[t]]))
        sprintf('<span class="cmp-chip-sub">%s</span>',
                esc(trt_compound_map[[t]]))
      else ""
      sprintf('<button type="button" class="%s"%s>%s%s</button>',
              cls, attr, esc(t), sub)
    }, character(1), USE.NAMES = FALSE)
    
    compound_filter_html <- sprintf(
      '<div class="cmp-filter" id="cmpFilter">
         <span class="cmp-filter-label">Filter by compound</span>
         <input type="text" class="cmp-search" id="cmpSearch" placeholder="Search %d compounds\u2026" autocomplete="off" spellcheck="false" aria-label="Search compounds">
         <button type="button" class="cmp-chip cmp-all is-active" data-cmp="__all__">All</button>
         %s
         <button type="button" class="cmp-more" id="cmpMore" hidden></button>
         <span class="cmp-search-note" id="cmpSearchNote" hidden></span>
         <span class="cmp-filter-status" id="cmpFilterStatus"></span>
       </div>',
      length(filter_treatments),
      paste(chips, collapse = "\n")
    )
  }
  
  # -- 7. Summary stats -------------------------------------------------------
  
  n_files <- sum(vapply(sections, function(s) nrow(s$files), integer(1)))
  n_bytes <- sum(vapply(sections, function(s) sum(s$files$size, na.rm = TRUE),
                        numeric(1)))
  n_sections <- length(sections)
  
  # -- 7b. DE & QC summary stats -----------------------------------------------
  
  summary_stats <- list()
  
  if (!is.null(de_results) && is.data.frame(de_results) && nrow(de_results) > 0) {
    n_contrasts   <- length(unique(de_results$Sample))
    n_genes_total <- length(unique(de_results$ENSEMBLID))
    sig           <- de_results[!is.na(de_results$padj) & de_results$padj < padj &
                                  !is.na(de_results$LFC)  & abs(de_results$LFC) >= lfc, ]
    n_de_genes    <- length(unique(sig$ENSEMBLID))
    n_de_hits     <- nrow(sig)
    
    genes_per_contrast <- tapply(de_results$n_genes_tested, de_results$Sample,
                                 function(x) x[1])
    median_genes <- median(genes_per_contrast, na.rm = TRUE)
    
    if ("treatment" %in% names(sig) && nrow(sig) > 0) {
      trt_counts <- sort(table(sig$treatment), decreasing = TRUE)
      top_trt    <- names(trt_counts)[1]
      top_trt_n  <- unname(trt_counts[1])
    } else {
      top_trt <- "\u2014"; top_trt_n <- 0
    }
    
    summary_stats <- list(
      contrasts       = n_contrasts,
      genes_tested    = n_genes_total,
      median_genes    = round(median_genes),
      de_genes        = n_de_genes,
      de_hits         = n_de_hits,
      top_treatment   = top_trt,
      top_treatment_n = top_trt_n
    )
  }
  
  dmso_stats <- list()
  if (!is.null(best_dmso)) {
    n_kept    <- if (!is.null(best_dmso$best_DMSO)) nrow(best_dmso$best_DMSO) else 0
    n_flagged <- if (!is.null(best_dmso$bad_DMSO))  nrow(best_dmso$bad_DMSO)  else 0
    dmso_stats <- list(kept = n_kept, flagged = n_flagged)
  }
  
  # -- 8. Build header detail rows --------------------------------------------
  
  make_detail_row <- function(label, values, mono = TRUE) {
    if (!length(values) || all(values == "")) return("")
    tag <- if (mono) "code" else "span"
    pills <- paste(sprintf('<%s class="pill">%s</%s>', tag, esc(values), tag), collapse = " ")
    sprintf(
      '<div class="detail-row"><span class="detail-label">%s</span><div class="detail-values">%s</div></div>',
      esc(label), pills
    )
  }
  
  # Treatments with compound suffix
  TRT_COLLAPSE_N <- 10L   # show this many before "… show more"
  
  trt_pills <- character(0)
  if (length(display_treatments)) {
    trt_pills <- vapply(display_treatments, function(t) {
      cmp <- if (t %in% names(trt_compound_map)) paste0(" \u2014 ", trt_compound_map[t]) else ""
      sprintf('<code class="pill">%s%s</code>', esc(t), esc(cmp))
    }, character(1), USE.NAMES = FALSE)
  }
  trt_row <- if (length(trt_pills)) {
    n_trt <- length(trt_pills)
    if (n_trt <= TRT_COLLAPSE_N) {
      # All pills fit comfortably — render inline as before
      sprintf(
        '<div class="detail-row"><span class="detail-label">Treatments (%d)</span><div class="detail-values">%s</div></div>',
        n_trt, paste(trt_pills, collapse = " ")
      )
    } else {
      # Too many to show at once — render first N + a toggle button
      visible_pills  <- paste(trt_pills[seq_len(TRT_COLLAPSE_N)], collapse = " ")
      hidden_pills   <- paste(trt_pills[seq(TRT_COLLAPSE_N + 1L, n_trt)], collapse = " ")
      n_hidden <- n_trt - TRT_COLLAPSE_N
      sprintf(
        '<div class="detail-row trt-collapsible">
           <span class="detail-label">Treatments (%d)</span>
           <div class="detail-values">
             %s
             <span class="trt-overflow">%s</span>
             <button type="button" class="trt-toggle" data-less="show less"
                     data-more="+ %d more">+ %d more</button>
           </div>
         </div>',
        n_trt,
        visible_pills,
        hidden_pills,
        n_hidden, n_hidden
      )
    }
  } else ""
  
  batch_row <- make_detail_row(sprintf("Batches (%d)", length(batches)), batches)
  plate_row <- make_detail_row(sprintf("Plates (%d)", length(plates)), plates)
  
  # Doses: show a range pill (min – max) rather than one pill per dose so the
  # header stays compact even when there are many dilution steps.
  dose_row <- if (length(doses) >= 2) {
    fmt_d <- function(d) format(d, scientific = FALSE, trim = TRUE)
    u     <- if (nzchar(dose_unit)) paste0("\u00a0", dose_unit) else ""
    range_str <- sprintf("%s\u2013%s%s", fmt_d(min(doses)), fmt_d(max(doses)), u)
    sprintf(
      '<div class="detail-row"><span class="detail-label">Doses (%d)</span><div class="detail-values"><code class="pill">%s</code></div></div>',
      length(doses), esc(range_str)
    )
  } else if (length(doses) == 1) {
    lbl <- paste0(format(doses, scientific = FALSE, trim = TRUE),
                  if (nzchar(dose_unit)) paste0("\u00a0", dose_unit) else "")
    make_detail_row("Doses (1)", lbl)
  } else ""
  tp_row    <- make_detail_row("Timepoints", paste0(timepoints, " hr"))
  ctrl_row  <- make_detail_row("Control", ctrl)
  
  detail_block <- paste(trt_row, batch_row, plate_row, dose_row, tp_row, ctrl_row, sep = "\n")
  detail_block <- trimws(detail_block)
  
  # -- 8b. Build summary panel HTML -------------------------------------------
  
  summary_panel_html <- ""
  if (length(summary_stats) > 0 || length(dmso_stats) > 0) {
    rows <- c()
    if (length(summary_stats) > 0) {
      rows <- c(rows,
                sprintf('<tr><td class="sl">Contrasts tested</td><td class="sv">%s</td></tr>',
                        format(summary_stats$contrasts, big.mark = ",")),
                sprintf('<tr><td class="sl">Genes tested</td><td class="sv">%s</td></tr>',
                        format(summary_stats$genes_tested, big.mark = ",")),
                sprintf('<tr><td class="sl">Median genes / contrast</td><td class="sv">%s</td></tr>',
                        format(summary_stats$median_genes, big.mark = ",")),
                sprintf('<tr><td class="sl">DE genes (unique)</td><td class="sv sv-pass">%s</td></tr>',
                        format(summary_stats$de_genes, big.mark = ",")),
                sprintf('<tr><td class="sl">DE hits (total)</td><td class="sv">%s</td></tr>',
                        format(summary_stats$de_hits, big.mark = ",")),
                sprintf('<tr><td class="sl">Top treatment</td><td class="sv">%s <span class="sv-sub">(%s hits)</span></td></tr>',
                        esc(summary_stats$top_treatment), format(summary_stats$top_treatment_n, big.mark = ","))
      )
    }
    if (length(dmso_stats) > 0) {
      rows <- c(rows,
                sprintf('<tr><td class="sl">DMSO wells kept</td><td class="sv sv-pass">%d</td></tr>',
                        dmso_stats$kept),
                sprintf('<tr><td class="sl">DMSO wells flagged</td><td class="sv%s">%d</td></tr>',
                        if (dmso_stats$flagged > 0) " sv-warn" else "", dmso_stats$flagged)
      )
    }
    summary_panel_html <- sprintf(
      '<div class="summary-panel"><h3>Analysis Summary</h3><table class="summary-table">%s</table></div>',
      paste(rows, collapse = "\n")
    )
  }
  
  
  # -- 9b. Reproducibility block (added in response to reviewer feedback) -----
  # Embeds the resolved config and the R session info directly in the report
  # so the final HTML is a fully self-contained reproducibility record.
  # Also writes both as sidecar files for programmatic access.
  
  # Resolved config — print as YAML if the package is available, else as
  # str() output (always works, no extra dependencies)
  config_text <- tryCatch(
    paste(yaml::as.yaml(config), collapse = "\n"),
    error = function(e) paste(capture.output(utils::str(config, max.level = 4)),
                              collapse = "\n")
  )
  
  # sessionInfo — full version: package versions, locale, platform
  sinfo_text <- paste(capture.output(print(utils::sessionInfo())), collapse = "\n")
  
  # Ensembl version — read from gene_annotation if it's been recorded there
  # (see 02_input_functions.R / create_gene_annotation_file changes).
  # Falls back to "not recorded" if the annotation predates that change.
  ensembl_version_text <- tryCatch({
    ann_path <- file.path(config$paths$output_dir, "annotations",
                          "gene_annotation.rds")
    if (file.exists(ann_path)) {
      ann <- readRDS(ann_path)
      # Bundle is the new format (named list with $annotation +
      # provenance); a data.frame is the legacy format. Note that
      # data.frame also passes is.list(), so check !is.data.frame first.
      if (!is.data.frame(ann) && is.list(ann) &&
          !is.null(ann$ensembl_version))
        sprintf("Ensembl release %s (queried %s)",
                ann$ensembl_version,
                if (!is.null(ann$query_date)) ann$query_date else "unknown date")
      else
        "Ensembl release: not recorded (annotation file predates version-tracking)"
    } else "Ensembl release: annotation file not found"
  }, error = function(e) paste("Ensembl release: lookup failed —", conditionMessage(e)))
  
  # Sidecar files
  writeLines(config_text, file.path(results_dir, "run_config.yml"))
  writeLines(sinfo_text,  file.path(results_dir, "session_info.txt"))
  
  reproducibility_html <- sprintf('
<section class="repro-panel">
  <h3>Reproducibility</h3>
  <p class="repro-note">%s. The blocks below capture the exact configuration
  and software environment used to produce this report. Both are also saved
  as sidecar files (<code>run_config.yml</code>, <code>session_info.txt</code>)
  alongside this HTML.</p>
  <details>
    <summary>Resolved configuration</summary>
    <pre class="repro-block">%s</pre>
  </details>
  <details>
    <summary>R session info (package versions, platform, locale)</summary>
    <pre class="repro-block">%s</pre>
  </details>
</section>',
                                  esc(ensembl_version_text),
                                  esc(config_text),
                                  esc(sinfo_text)
  )
  
  # -- 9c. Assemble final HTML -----------------------------------------------
  
  body_html <- paste(vapply(sections, render_section, character(1)), collapse = "\n")
  body_html <- paste(body_html, reproducibility_html, sep = "\n")
  nav_html  <- render_nav(sections)
  
  css <- '
  *{box-sizing:border-box;margin:0;padding:0}
  body{font-family:"Helvetica Neue",Helvetica,Arial,sans-serif;background:#f5f5f5;color:#1f2937;min-height:100vh}
  code,.mono{font-family:"IBM Plex Mono",ui-monospace,monospace}

  header{background:#3d3d3d;color:#fff;padding:26px 36px}
  header h1{font-size:22px;font-weight:700;letter-spacing:-0.02em}
  header h1 span{color:#7ab648}
  header .subtitle{margin-top:4px;font-size:13px;color:#b0b0b0;max-width:700px;line-height:1.5}
  header .meta-line{margin-top:10px;font-size:13px;color:#b0b0b0;display:flex;gap:18px;flex-wrap:wrap}
  header .meta-line b{color:#e2e8f0;font-weight:600}
  header .chips{margin-top:12px;display:flex;gap:8px;flex-wrap:wrap}
  header .chip{background:#4a4a4a;color:#7ab648;padding:3px 10px;border-radius:12px;font-size:11.5px;font-family:"IBM Plex Mono",monospace;font-weight:500}
  header .chip b{color:#fff;margin-right:4px}

  .details-panel{background:#2d2d2d;padding:18px 36px 20px;border-top:1px solid #4a4a4a}
  .detail-row{display:flex;align-items:baseline;gap:10px;margin-bottom:8px;flex-wrap:wrap}
  .detail-label{font-size:12px;font-weight:600;color:#999;min-width:140px;flex-shrink:0;text-transform:uppercase;letter-spacing:0.04em}
  .detail-values{display:flex;flex-wrap:wrap;gap:6px}
  .pill{display:inline-block;background:#4a4a4a;color:#cbd5e1;padding:3px 10px;border-radius:6px;font-size:12px;font-family:"IBM Plex Mono",monospace;line-height:1.4;border:1px solid #5a5a5a}

  .summary-panel{background:#fff;border:1px solid #e0e0e0;border-radius:8px;max-width:1240px;margin:24px auto 0;padding:20px 28px}
  .summary-panel h3{font-size:13px;font-weight:700;color:#3d3d3d;margin-bottom:12px;text-transform:uppercase;letter-spacing:0.04em}
  .summary-table{width:100%;max-width:480px;border-collapse:collapse;font-size:13px}
  .summary-table tr{border-bottom:1px solid #f0f0f0}
  .summary-table tr:last-child{border-bottom:none}
  .sl{padding:6px 16px 6px 0;color:#666;font-weight:500}
  .sv{padding:6px 0;font-family:"IBM Plex Mono",monospace;font-weight:700;color:#3d3d3d;text-align:right}
  .sv-pass{color:#7ab648}
  .sv-warn{color:#d97706}
  .sv-sub{font-weight:400;font-size:11px;color:#999}

  main{max-width:1240px;margin:0 auto;display:grid;grid-template-columns:220px 1fr;gap:32px;padding:28px 24px}
  @media (max-width:800px){main{grid-template-columns:1fr}}

  nav.toc{position:sticky;top:20px;align-self:start;background:#fff;border:1px solid #e0e0e0;border-radius:8px;padding:12px;max-height:calc(100vh - 40px);overflow:auto}
  nav.toc a{display:flex;align-items:center;gap:8px;padding:7px 10px;text-decoration:none;color:#374151;font-size:13px;border-radius:5px}
  nav.toc a:hover{background:#e8f5e9;color:#2e7d32}
  nav.toc a .n{margin-left:auto;font-family:"IBM Plex Mono",monospace;font-size:11px;color:#9ca3af}

  section{margin-bottom:40px;scroll-margin-top:20px}
  .section-head{display:flex;align-items:baseline;gap:12px;margin-bottom:14px;border-bottom:2px solid #3d3d3d;padding-bottom:8px;flex-wrap:wrap}
  .section-head h2{font-size:17px;font-weight:700;color:#3d3d3d}
  .section-head .count{font-size:12px;color:#6b7280;font-family:"IBM Plex Mono",monospace}
  .section-head .dir{margin-left:auto;font-size:11.5px;color:#9ca3af;background:#f3f4f6;padding:2px 8px;border-radius:4px}
  .section-desc{margin:-4px 0 14px;font-size:13px;color:#6b7280;line-height:1.55;max-width:720px}

  .grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(180px,1fr));gap:14px}
  .card{position:relative;display:flex;flex-direction:column;background:#fff;border:1px solid #e0e0e0;border-radius:8px;overflow:hidden;transition:all .15s}
  .card:hover{border-color:#7ab648;transform:translateY(-1px);box-shadow:0 6px 18px -10px rgba(122,182,72,0.4)}
  .card-body{display:flex;flex-direction:column;flex:1;text-decoration:none;color:inherit}
  .card-dl{position:absolute;top:6px;right:6px;width:28px;height:28px;display:flex;align-items:center;justify-content:center;background:rgba(255,255,255,0.94);border:1px solid #e0e0e0;border-radius:6px;color:#374151;opacity:0;transition:opacity .12s,color .12s,border-color .12s;z-index:3}
  .card-dl svg{width:15px;height:15px}
  .card:hover .card-dl{opacity:1}
  .card-dl:hover{color:#0077c8;border-color:#0077c8;background:#eff6ff}
  .thumb{aspect-ratio:4/3;background:#fafafa;display:flex;align-items:center;justify-content:center;overflow:hidden;position:relative}
  .thumb img{max-width:100%;max-height:100%;object-fit:contain}
  .html-thumb{display:flex;flex-direction:column;align-items:center;justify-content:center;gap:6px;width:100%;height:100%;background:linear-gradient(135deg,#e8f5e9 0%,#c8e6c9 50%,#e8f5e9 100%)}
  .html-icon{width:40px;height:40px;color:#0077c8}
  .html-label{font-size:11px;font-weight:700;color:#0077c8;letter-spacing:0.03em;text-transform:uppercase;font-family:"IBM Plex Mono",monospace}
  .html-hint{font-size:10px;color:#7ab648;font-weight:500}
  .card:hover .html-thumb{background:linear-gradient(135deg,#c8e6c9 0%,#a5d6a7 50%,#c8e6c9 100%)}
  .card:hover .html-icon{color:#005a9e}
  .thumb .icon{font-size:40px;display:flex;flex-direction:column;align-items:center;gap:4px}
  .thumb .icon-label{font-size:10px;font-family:"IBM Plex Mono",monospace;font-weight:600;color:#6b7280;letter-spacing:0.05em}
  .meta{padding:8px 10px;border-top:1px solid #f1f5f9}
  .meta .name{font-size:12px;font-weight:500;color:#1f2937;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;font-family:"IBM Plex Mono",monospace}
  .meta .sub{font-size:11px;color:#9ca3af;margin-top:2px;font-family:"IBM Plex Mono",monospace}

  .data-list-container{margin-top:14px;border:1px solid #e0e0e0;border-radius:8px;background:#fff;overflow:hidden}
  .data-list-header{padding:8px 14px;font-size:12px;font-weight:600;color:#6b7280;background:#fafafa;border-bottom:1px solid #e0e0e0;text-transform:uppercase;letter-spacing:0.04em}
  .data-list{display:flex;flex-direction:column}
  .data-link{display:flex;align-items:center;gap:12px;padding:8px 14px;text-decoration:none;color:#374151;border-bottom:1px solid #f3f4f6;transition:background .1s;font-size:13px}
  .data-link:last-child{border-bottom:none}
  .data-link:hover{background:#e8f5e9}
  .data-ext{display:inline-block;width:38px;text-align:center;font-size:10px;font-weight:700;font-family:"IBM Plex Mono",monospace;color:#fff;background:#3d3d3d;padding:2px 6px;border-radius:4px;flex-shrink:0;letter-spacing:0.03em}
  .data-name{flex:1;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;font-family:"IBM Plex Mono",monospace;font-size:12.5px}
  .data-size{font-family:"IBM Plex Mono",monospace;font-size:11px;color:#9ca3af;flex-shrink:0}

  .repro-panel{background:#fff;border:1px solid #e0e0e0;border-radius:8px;max-width:1240px;margin:24px auto 0;padding:20px 28px}
  .repro-panel h3{font-size:13px;font-weight:700;color:#3d3d3d;margin-bottom:8px;text-transform:uppercase;letter-spacing:0.04em}
  .repro-note{font-size:12.5px;color:#666;margin-bottom:14px;line-height:1.6}
  .repro-panel details{border-top:1px solid #f0f0f0;padding:8px 0}
  .repro-panel details summary{cursor:pointer;font-size:13px;font-weight:600;color:#3d3d3d;padding:6px 0;outline:none}
  .repro-panel details summary:hover{color:#2e7d32}
  .repro-block{background:#1e1e2e;color:#cdd6f4;font-family:"IBM Plex Mono",monospace;font-size:11.5px;line-height:1.55;padding:14px 16px;margin:8px 0 4px;border-radius:6px;overflow-x:auto;white-space:pre-wrap;word-break:break-word;max-height:480px;overflow-y:auto}

  /* ------- Volcano-pair cards (SVG ↔ HTML toggle) ------------------------ */
  /* ------- DE results table (merged CSV + HTML twin per row) ------------- */
  .de-table{width:100%;border-collapse:collapse;margin-top:4px;background:#fff;border:1px solid #e0e0e0;border-radius:8px;overflow:hidden}
  .de-table thead th{background:#f3f4f6;color:#3d3d3d;font-size:11px;font-weight:700;letter-spacing:0.04em;text-transform:uppercase;text-align:left;padding:10px 14px;border-bottom:1px solid #e0e0e0}
  .de-table thead th:nth-child(2),.de-table thead th:nth-child(3){text-align:center;width:160px}
  .de-table tbody tr{border-bottom:1px solid #f0f0f0;transition:background .12s}
  .de-table tbody tr:hover{background:#f9fafb}
  .de-table tbody tr:last-child{border-bottom:none}
  .de-table tr.de-group td{background:#eef2f7;color:#334155;font-size:11px;font-weight:700;letter-spacing:0.05em;text-transform:uppercase;padding:8px 14px;border-top:1px solid #e0e0e0}
  .de-table tr.de-group:hover td{background:#eef2f7}
  .de-desc{padding:10px 14px;vertical-align:top}
  .de-desc-main{font-size:13px;font-weight:600;color:#1f2937;letter-spacing:-0.005em;line-height:1.4}
  .de-desc-pills{margin-top:4px;display:flex;flex-wrap:wrap;gap:4px}
  .de-pill{display:inline-block;background:#eef2f7;color:#475569;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:10.5px;padding:2px 8px;border-radius:10px;line-height:1.4;font-weight:500}
  .de-pill-trt{background:#e0f2fe;color:#0c4a6e}
  .pill-group{display:inline}
  .pill-more{background:#f1f5f9;border:1px solid #cbd5e1;color:#334155;font-size:10.5px;font-weight:600;padding:2px 9px;border-radius:10px;cursor:pointer;line-height:1.4;margin-left:2px}
  .pill-more:hover{background:#e2e8f0}
  img[data-src]{background:#f1f5f9;min-height:80px}
  img.lazy-loaded{animation:lazyfade .25s ease-in}
  @keyframes lazyfade{from{opacity:0}to{opacity:1}}
  .grid-capped{display:flex;flex-direction:column;gap:12px}
  /* Cards beyond the cap. !important so it beats .card/.vc-card display:flex. */
  .grid-over{display:none !important}
  .pill-hidden[hidden]{display:none !important}
  .grid-more{align-self:flex-start;background:#0077c8;color:#fff;border:none;border-radius:8px;font-size:13px;font-weight:600;padding:8px 18px;cursor:pointer;box-shadow:0 1px 2px rgba(0,0,0,.08)}
  .grid-more:hover{background:#005fa3}
  .de-cell{padding:10px 14px;text-align:center;vertical-align:middle}
  .de-link{display:inline-flex;align-items:center;gap:8px;text-decoration:none;color:#1f2937;background:#fff;border:1px solid #d1d5db;padding:5px 10px;border-radius:6px;font-size:11.5px;transition:border-color .12s, color .12s, background .12s}
  .de-link:hover{border-color:#0077c8;background:#eff6ff;color:#0077c8}
  .de-link-empty{color:#cbd5e1;border-style:dashed;cursor:default}
  .de-link-empty:hover{border-color:#cbd5e1;background:#fff;color:#cbd5e1}
  .de-ext{display:inline-block;font-family:"IBM Plex Mono",monospace;font-size:9.5px;font-weight:700;letter-spacing:0.04em;color:#fff;padding:2px 6px;border-radius:3px;line-height:1}
  .de-ext-csv{background:#3d3d3d}
  .de-ext-html{background:#0077c8}
  .de-size{font-family:"IBM Plex Mono",monospace;font-size:11px;color:#9ca3af}

  .vc-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:16px;margin-top:4px}
  .vc-grid.single-item{grid-template-columns:minmax(0,520px);justify-content:center}
  .vc-row-label{font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:0.06em;color:#9ca3af;margin:18px 0 6px;padding-bottom:6px;border-bottom:1px solid #f0f0f0}
  .vc-card{background:#fff;border:1px solid #e0e0e0;border-radius:8px;overflow:hidden;display:flex;flex-direction:column;transition:border-color .15s, box-shadow .15s}
  .vc-card:hover{border-color:#cbd5e1}
  .vc-toolbar{display:flex;align-items:center;gap:10px;padding:8px 12px;background:#fafafa;border-bottom:1px solid #e0e0e0;flex-wrap:wrap}
  .vc-label{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;font-weight:600;color:#3d3d3d;flex:1;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;letter-spacing:0.01em}
  .vc-actions-links{display:inline-flex;gap:6px;flex-shrink:0;flex-wrap:wrap}
  .vc-btn-link{display:inline-flex;align-items:center;justify-content:center;gap:5px;font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:11px;font-weight:600;padding:4px 8px;color:#6b7280;border:1px solid #d1d5db;border-radius:6px;background:#fff;text-decoration:none;white-space:nowrap;line-height:1}
  .vc-btn-link svg{width:14px;height:14px}
  .vc-btn-link:hover{color:#0077c8;border-color:#0077c8;background:#eff6ff}
  .vc-viewport{position:relative;background:#fff;display:flex;align-items:center;justify-content:center;min-height:160px;overflow:hidden}
  .vc-svg{display:block;width:100%;height:auto;max-height:300px;object-fit:contain}

  /* single-item centering for regular grids */
  .grid.single-item{grid-template-columns:minmax(0,360px);justify-content:center}

  /* treatment collapse */
  .trt-overflow{display:none}
  .trt-overflow.expanded{display:contents}
  .trt-toggle{background:none;border:1px solid #5a5a5a;color:#9ca3af;font-size:11px;font-family:"IBM Plex Mono",monospace;padding:3px 9px;border-radius:6px;cursor:pointer;flex-shrink:0;transition:color .12s,border-color .12s}
  .trt-toggle:hover{color:#7ab648;border-color:#7ab648}

  /* ------- Compound filter bar ------------------------------------------ */
  .cmp-filter{position:sticky;top:0;z-index:30;display:flex;align-items:center;gap:8px;flex-wrap:wrap;padding:10px 36px;background:#fff;border-bottom:1px solid #e0e0e0;box-shadow:0 2px 8px -5px rgba(0,0,0,0.3)}
  @media (max-width:800px){.cmp-filter{padding:10px 16px}}
  .cmp-filter-label{font-size:12px;font-weight:700;color:#3d3d3d;text-transform:uppercase;letter-spacing:0.04em;margin-right:4px}
  .cmp-chip{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;color:#374151;background:#f3f4f6;border:1px solid #d1d5db;border-radius:14px;padding:4px 12px;cursor:pointer;display:inline-flex;align-items:center;gap:6px;line-height:1.3;transition:background .12s,border-color .12s,color .12s}
  .cmp-chip:hover{border-color:#7ab648;color:#2e7d32}
  .cmp-chip.is-active{background:#7ab648;border-color:#7ab648;color:#fff}
  .cmp-chip .cmp-chip-sub{font-size:10.5px;opacity:0.7}
  .cmp-chip.is-active .cmp-chip-sub{opacity:0.95}
  .cmp-chip-off{opacity:0.4;cursor:not-allowed;text-decoration:line-through}
  .cmp-chip-off:hover{border-color:#d1d5db;color:#374151}
  .cmp-filter-status{font-size:12px;color:#6b7280;margin-left:auto;font-family:"IBM Plex Mono",ui-monospace,monospace}
  .cmp-search{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;color:#374151;background:#fff;border:1px solid #d1d5db;border-radius:14px;padding:4px 12px;min-width:190px;outline:none;transition:border-color .12s,box-shadow .12s}
  .cmp-search:focus{border-color:#7ab648;box-shadow:0 0 0 2px rgba(122,182,72,0.2)}
  .cmp-search::placeholder{color:#9ca3af}
  .cmp-more{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:12px;color:#2e7d32;background:transparent;border:1px dashed #b6d99a;border-radius:14px;padding:4px 12px;cursor:pointer;transition:background .12s,border-color .12s}
  .cmp-more:hover{background:#f0f7e9;border-color:#7ab648}
  .cmp-more[hidden]{display:none}
  .cmp-search-note{font-size:12px;color:#9ca3af;font-family:"IBM Plex Mono",ui-monospace,monospace}
  .cmp-search-note[hidden]{display:none}
  .cmp-chip-collapsed{display:none !important}
  .cmp-hidden{display:none !important}

  /* ------- Font uniformity safety net ----------------------------------- */
  /* Tables, code, and pre elements should follow the explicit stack for
     consistency. Without this, browsers fall back to UA defaults (Times
     for table content in some configs; monospace for pre that does not
     match Plex Mono). */
  table{font-family:"Helvetica Neue",Helvetica,Arial,sans-serif}
  pre{font-family:"IBM Plex Mono",ui-monospace,monospace}
  button,input,select{font-family:inherit}

  footer{max-width:1240px;margin:0 auto;padding:20px 24px 40px;color:#9ca3af;font-size:12px;border-top:1px solid #e0e0e0;margin-top:20px}
  '
  
  html <- paste0(
    '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>',
    esc(title),
    ' — Pipeline Report</title>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&amp;display=swap" rel="stylesheet"/>
<style>',
    css,
    '</style>
</head>
<body>

<header>
  <h1><span>ADEPT</span> ',
    esc(title),
    '</h1>
  ',
    if (nzchar(description)) sprintf('<div class="subtitle">%s</div>', esc(description)) else "",
    '
  <div class="meta-line">
    <span>Run <b>',
    esc(run_id),
    '</b></span>
    ',
    if (nzchar(cell_line)) sprintf("<span>Cell line <b>%s</b></span>", esc(cell_line)) else "",
    '
    ',
    if (nzchar(citation))  sprintf("<span>%s</span>", esc(citation)) else "",
    '
    <span>padj &le; <b>',
    esc(format(padj)),
    '</b></span>
    <span>|log2FC| &ge; <b>',
    esc(format(lfc)),
    '</b></span>
    <span>Generated <b>',
    format(Sys.time(), "%Y-%m-%d %H:%M"),
    '</b></span>
    <span>Wells <b>',
    format(n_wells, big.mark = ","),
    '</b></span>
  </div>
  <div class="chips">
    <span class="chip"><b>',
    n_sections,
    '</b>sections</span>
    <span class="chip"><b>',
    n_files,
    '</b>files</span>
    <span class="chip"><b>',
    fmt_size(n_bytes),
    '</b>total</span>
  </div>
</header>

<div class="details-panel">',
    detail_block,
    '</div>

',
    paste0(compound_filter_html, "\n", summary_panel_html),
    '

<main>
  <nav class="toc">',
    nav_html,
    '</nav>
  <div>',
    body_html,
    '</div>
</main>

<footer>Generated by ADEPT (Automated Differential Expression Pipeline for high-throughput Transcriptomics screens), <code>build_final_report()</code>. All links are relative to this file &mdash; zip <code>',
    esc(basename(results_dir)),
    '/</code> to share the full report.</footer>

<script>
// ---- Treatment pills collapse/expand ------------------------------------
(function() {
  document.addEventListener("click", function(ev) {
    var btn = ev.target.closest(".trt-toggle");
    if (!btn) return;
    var row      = btn.closest(".trt-collapsible");
    if (!row) return;
    var overflow = row.querySelector(".trt-overflow");
    if (!overflow) return;
    var expanded = overflow.classList.toggle("expanded");
    btn.textContent = expanded ? btn.dataset.less : btn.dataset.more;
  });
})();

// ---- Grid caps ("Show all N") -------------------------------------------
// One controller for every capped grid, shared by the "Show all" button and
// the compound filter. A cap is recomputed over the cards that currently pass
// the filter (those without .cmp-hidden): the first data-nshow of them stay
// visible, the remainder get .grid-over, and the button label counts only the
// passing cards. Lives outside the filter block so it works in reports that
// have no filter bar.
(function() {
  function gridCards(wrap) {
    var grid = wrap.querySelector(".grid, .vc-grid");
    if (!grid) return [];
    return Array.prototype.filter.call(grid.children, function(ch) {
      return ch.nodeType === 1;
    });
  }

  // Hand newly-visible images to the lazy-loader (absent before boot, in
  // which case initLazy() picks them up itself).
  function reveal(el) {
    if (!window.__lazyObserve) return;
    Array.prototype.forEach.call(el.querySelectorAll("img[data-src]"),
                                 window.__lazyObserve);
  }

  function recap(wrap) {
    var btn = wrap.querySelector(".grid-more");
    if (!btn) return;
    var nShow = parseInt(wrap.getAttribute("data-nshow"), 10) || 3;
    var open  = btn.getAttribute("data-open") === "1";
    var shown = 0, total = 0;
    gridCards(wrap).forEach(function(card) {
      // Filtered-out cards are hidden by .cmp-hidden and never count towards
      // the cap; park them in .grid-over so clearing the filter is coherent.
      if (card.classList.contains("cmp-hidden")) {
        card.classList.add("grid-over");
        return;
      }
      total++;
      var over = !open && shown >= nShow;
      if (!over) { shown++; reveal(card); }
      card.classList.toggle("grid-over", over);
    });
    if (total <= nShow) {                 // nothing left to reveal
      btn.setAttribute("data-open", "0");
      btn.style.display = "none";
      open = false;
    } else {
      btn.style.display = "";
    }
    btn.textContent = open
      ? "Show fewer " + String.fromCharCode(0x25B4)
      : "Show all " + total + " " + String.fromCharCode(0x25BE);
  }

  // collapse = true re-closes every grid first, so a new compound selection
  // always starts from the cap instead of inheriting an earlier "Show all".
  window.__recapGrids = function(collapse) {
    Array.prototype.forEach.call(
      document.querySelectorAll(".grid-capped"),
      function(wrap) {
        if (collapse) {
          var btn = wrap.querySelector(".grid-more");
          if (btn) btn.setAttribute("data-open", "0");
        }
        recap(wrap);
      });
  };

  window.__toggleGrid = function(btn) {
    var wrap = btn.closest(".grid-capped");
    if (!wrap) return;
    btn.setAttribute("data-open",
                     btn.getAttribute("data-open") === "1" ? "0" : "1");
    recap(wrap);
  };
})();

// ---- Compound filter ----------------------------------------------------
// Toggle chips to show only the chosen compounds. No chips selected = show
// everything. Cards/links carry data-compounds="id1 id2"; items without that
// attribute (global QC, PCA, UMAP, DE tables, ...) are never hidden. Grids,
// row labels and whole sections collapse when all their compound-specific
// items are hidden and they hold nothing else.
(function() {
  var bar = document.getElementById("cmpFilter");
  if (!bar) return;
  var status = document.getElementById("cmpFilterStatus");
  var allBtn = bar.querySelector(".cmp-all");
  var chips  = Array.prototype.slice.call(
                 bar.querySelectorAll(".cmp-chip[data-cmp]"))
               .filter(function(c){ return c.getAttribute("data-cmp") !== "__all__"; });
  var items  = Array.prototype.slice.call(
                 document.querySelectorAll("[data-compounds]"));

  // -- Preview + search over the chip roster --------------------------------
  // With hundreds of compounds the full chip list is unusable, so show only
  // the first PREVIEW chips (plus any that are currently selected) and let the
  // user type to search or click "Show all" to expand. This only governs which
  // chips are visible; it never changes their selected state.
  var PREVIEW  = 10;
  var expanded = false;
  var search   = document.getElementById("cmpSearch");
  var moreBtn  = document.getElementById("cmpMore");
  var note     = document.getElementById("cmpSearchNote");
  var allChips = Array.prototype.slice.call(bar.querySelectorAll(".cmp-chip"))
                 .filter(function(c){ return !c.classList.contains("cmp-all"); });

  function refreshChips() {
    var q = (search ? search.value : "").trim().toLowerCase();
    var searching = q.length > 0;
    var shown = 0, i = 0;
    allChips.forEach(function(c) {
      var matches = !searching ||
                    c.textContent.toLowerCase().indexOf(q) !== -1;
      var show;
      if (searching)     show = matches;
      else if (expanded) show = true;
      else               show = (i < PREVIEW) || c.classList.contains("is-active");
      c.classList.toggle("cmp-chip-collapsed", !show);
      if (show) shown++;
      i++;
    });
    if (moreBtn) {
      if (searching || allChips.length <= PREVIEW) {
        moreBtn.hidden = true;
      } else {
        moreBtn.hidden = false;
        moreBtn.textContent = expanded
          ? "Show fewer"
          : ("Show all (" + allChips.length + ")");
      }
    }
    if (note) {
      var empty = searching && shown === 0;
      note.hidden = !empty;
      if (empty) note.textContent = "No compounds match";
    }
  }

  function selected() {
    return chips.filter(function(c){ return c.classList.contains("is-active"); })
                .map(function(c){ return c.getAttribute("data-cmp"); });
  }

  function visible(el){ return !el.classList.contains("cmp-hidden"); }

  function collapseEmpties() {
    var grids = document.querySelectorAll(".grid, .vc-grid");
    Array.prototype.forEach.call(grids, function(grid) {
      var filterable = grid.querySelectorAll("[data-compounds]");
      if (!filterable.length) return;
      var anyVisible = Array.prototype.some.call(filterable, visible);
      var hasOther = Array.prototype.some.call(grid.children, function(ch){
        return ch.nodeType === 1 && !ch.hasAttribute("data-compounds") &&
               (ch.classList.contains("card") || ch.classList.contains("vc-card"));
      });
      var hide = !anyVisible && !hasOther;
      // Hide the .grid-capped wrapper when there is one, so the "Show all"
      // button goes with the grid; the row label sits before the wrapper.
      var host = grid.closest(".grid-capped") || grid;
      host.classList.toggle("cmp-hidden", hide);
      var prev = host.previousElementSibling;
      if (prev && prev.classList.contains("vc-row-label"))
        prev.classList.toggle("cmp-hidden", hide);
    });
    var secs = document.querySelectorAll("main section");
    Array.prototype.forEach.call(secs, function(sec) {
      var filterable = sec.querySelectorAll("[data-compounds]");
      if (!filterable.length) return;
      var anyVisible = Array.prototype.some.call(filterable, visible);
      var hasOther = sec.querySelector(
        ".card:not([data-compounds]), .vc-card:not([data-compounds]), " +
        ".data-link:not([data-compounds]), table");
      sec.classList.toggle("cmp-hidden", !anyVisible && !hasOther);
    });
  }

  function apply() {
    var sel  = selected();
    var none = sel.length === 0;
    allBtn.classList.toggle("is-active", none);
    items.forEach(function(el) {
      if (none) { el.classList.remove("cmp-hidden"); return; }
      var tags  = (el.getAttribute("data-compounds") || "").split(" ");
      var match = tags.some(function(t){ return sel.indexOf(t) !== -1; });
      el.classList.toggle("cmp-hidden", !match);
    });
    if (window.__recapGrids) window.__recapGrids(true);
    collapseEmpties();
    if (status)
      status.textContent = none ? "" :
        ("Showing " + sel.length + " compound" + (sel.length > 1 ? "s" : ""));
  }

  bar.addEventListener("click", function(ev) {
    var btn = ev.target.closest(".cmp-chip");
    if (!btn || btn.getAttribute("aria-disabled") === "true") return;
    if (btn.classList.contains("cmp-all")) {
      chips.forEach(function(c){ c.classList.remove("is-active"); });
    } else {
      btn.classList.toggle("is-active");
    }
    apply();
    refreshChips();
  });

  if (search) search.addEventListener("input", refreshChips);
  if (moreBtn) moreBtn.addEventListener("click", function() {
    expanded = !expanded;
    refreshChips();
  });

  apply();
  refreshChips();
})();

// -- Interactivity: lazy images + expanders --------------------------------
// Wrapped so a failure in one piece cannot silently kill the others, and run
// after DOM parse so every card/button exists. Any failure logs a specific
// message to the console instead of failing quietly.
(function() {
  function initLazy() {
    var imgs = Array.prototype.slice.call(document.querySelectorAll("img[data-src]"));
    if (!imgs.length) return;
    function load(img) {
      if (!img.getAttribute("data-src")) return;
      img.addEventListener("load", function() {
        img.classList.add("lazy-loaded");
        if (img.getAttribute("data-raster") === "1") rasterize(img);
      });
      img.src = img.getAttribute("data-src");
      img.removeAttribute("data-src");
    }
    // Flatten a heavy vector image (volcano SVG with no PNG preview) into a
    // static <canvas> bitmap at its on-screen size, once. After that, scrolling
    // only composites pixels instead of re-rasterising tens of thousands of SVG
    // paths. Drawing a file:// image taints the canvas, which only blocks
    // reading pixels back - displaying it is fine. Any failure keeps the <img>.
    function rasterize(img) {
      var nw = img.naturalWidth, nh = img.naturalHeight;
      var box = img.getBoundingClientRect();
      if (!nw || !nh || !box.width || !box.height) return;
      var scale = Math.min(box.width / nw, box.height / nh);
      var cssW = nw * scale, cssH = nh * scale;
      var dpr  = Math.min(window.devicePixelRatio || 1, 2);
      var cv   = document.createElement("canvas");
      cv.width  = Math.max(1, Math.round(cssW * dpr));
      cv.height = Math.max(1, Math.round(cssH * dpr));
      try {
        var ctx = cv.getContext("2d");
        ctx.fillStyle = "#fff";
        ctx.fillRect(0, 0, cv.width, cv.height);
        ctx.drawImage(img, 0, 0, cv.width, cv.height);
      } catch (e) { return; }
      cv.className = img.className;
      cv.setAttribute("role", "img");
      cv.setAttribute("aria-label", img.alt || "");
      cv.style.width = cssW + "px";
      cv.style.height = "auto";          // keeps aspect if the column narrows
      cv.style.maxWidth = "100%";
      img.parentNode.replaceChild(cv, img);
    }
    if (!("IntersectionObserver" in window)) { imgs.forEach(load); return; }
    var io = new IntersectionObserver(function(entries, obs) {
      entries.forEach(function(e) {
        if (e.isIntersecting) { load(e.target); obs.unobserve(e.target); }
      });
    }, { rootMargin: "300px 0px" });
    imgs.forEach(function(img) { io.observe(img); });
    window.__lazyObserve = function(img) { io.observe(img); };
  }

  // One delegated click handler for both expanders — cannot be broken by load
  // order and covers cards revealed later.
  function initClicks() {
    document.addEventListener("click", function(ev) {
      var t = ev.target;
      if (!t || !t.closest) return;

      var pillBtn = t.closest(".pill-more");
      if (pillBtn) {
        var group = pillBtn.parentNode;
        var hidden = group.querySelector(".pill-hidden");
        if (!hidden) return;
        if (pillBtn.getAttribute("data-open") === "1") {
          hidden.setAttribute("hidden", "");
          pillBtn.textContent = "+" + hidden.querySelectorAll(".de-pill-trt").length + " more";
          pillBtn.setAttribute("data-open", "0");
        } else {
          hidden.removeAttribute("hidden");
          pillBtn.textContent = "show fewer";
          pillBtn.setAttribute("data-open", "1");
        }
        return;
      }

      var gridBtn = t.closest(".grid-more");
      if (gridBtn) {
        if (window.__toggleGrid) window.__toggleGrid(gridBtn);
        return;
      }
    });
  }

  function boot() {
    try { initLazy(); }  catch (e) { console.error("ADEPT report: lazy-load init failed:", e); }
    try { initClicks(); } catch (e) { console.error("ADEPT report: click init failed:", e); }
  }
  if (document.readyState === "loading")
    document.addEventListener("DOMContentLoaded", boot);
  else
    boot();
})();
</script>

</body>
</html>'
  )
  
  # -- 10. Write --------------------------------------------------------------
  
  report_name <- paste0("final_report_", run_id, ".html")
  out_file <- file.path(results_dir, report_name)
  writeLines(html, out_file, useBytes = TRUE)
  out_file
}