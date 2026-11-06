# --------------------------------------------------------------
# Font harmonisation (cross-machine, SVG + HTML + plotly)
# --------------------------------------------------------------
#
# A plot's font is only fixed at *render* time for rasters (PNG), but for
# SVG the font-family string is resolved by whatever machine *opens* the
# file. If that machine lacks the named font and the SVG has no generic
# fallback, the browser drops to its default serif — i.e. Times New Roman.
# That is the "fonts change when I open it on another PC" symptom.
#
# These two stacks are the single source of truth for typography. Both end
# in a generic family so a machine missing Helvetica/Arial falls back to its
# own sans-serif, never to a serif.

# Quoted form for author stylesheets / HTML (multi-word name needs quotes).
DRUGSEQ_FONT_STACK <- "'Helvetica Neue', Helvetica, Arial, sans-serif"

# Quote-free form, safe to splice into either single- or double-quoted inline
# SVG style attributes without breaking the attribute delimiter.
DRUGSEQ_FONT_STACK_INLINE <- "Helvetica, Arial, sans-serif"

# Plotly's `layout(font = list(family = ...))` accepts a plain CSS family list.
DRUGSEQ_PLOTLY_FONT <- "Helvetica Neue, Helvetica, Arial, sans-serif"

#' Force every <text>/<tspan> in an SVG onto one font stack.
#'
#' Works regardless of which device wrote the file (ggsave's svg, svglite, or
#' grDevices::svg) by combining two mechanisms:
#'   (a) an injected author stylesheet using !important — per the CSS cascade
#'       this beats both inline `style="font-family:..."` and the SVG
#'       `font-family="..."` presentation attribute, so it wins everywhere a
#'       browser renders the file;
#'   (b) a literal rewrite of any font-family values already in the markup, so
#'       the raw file is self-describing for tools that ignore <style>.
#' Idempotent and non-fatal: a missing file is silently skipped.
harmonize_svg_fonts <- function(svg_path, font_stack = DRUGSEQ_FONT_STACK) {
  if (length(svg_path) != 1L || is.na(svg_path) || !file.exists(svg_path))
    return(invisible(svg_path))
  
  svg <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  
  # (a) Rewrite literal font-family values: CSS-property form (...:...;) and
  #     SVG presentation-attribute form (...="..."). The inline (quote-free)
  #     stack is used so it never collides with the attribute's own quotes.
  svg <- gsub("font-family\\s*:\\s*[^;\"'}]+",
              paste0("font-family:", DRUGSEQ_FONT_STACK_INLINE),
              svg, perl = TRUE)
  svg <- gsub("font-family\\s*=\\s*\"[^\"]*\"",
              paste0("font-family=\"", DRUGSEQ_FONT_STACK_INLINE, "\""),
              svg, perl = TRUE)
  
  # (b) Inject the overriding stylesheet right after the opening <svg> tag,
  #     unless it is already present (keeps the function idempotent).
  if (!grepl("drugseq-font-fix", svg, fixed = TRUE)) {
    style_block <- sprintf(
      "<style type=\"text/css\">/* drugseq-font-fix */text,tspan{font-family:%s !important}</style>",
      font_stack)
    svg <- sub("(<svg[^>]*>)", paste0("\\1\n", style_block), svg, perl = TRUE)
  }
  
  writeLines(svg, svg_path, useBytes = TRUE)
  invisible(svg_path)
}

# --------------------------------------------------------------
# Experiment object
# --------------------------------------------------------------

make_experiment_list <- function(meta, combined_counts) {
  required_meta_cols <- c("plate_barcode", "batch", "well_index")
  
  missing <- setdiff(required_meta_cols, names(meta))
  if (length(missing) > 0) {
    stop(
      "Metadata missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
  
  # Ensure plate identifiers are present in both datasets
  meta_plate_ids <- unique(meta$plate_barcode)
  counts_plate_ids <- names(combined_counts)
  
  missing_counts <- setdiff(meta_plate_ids, counts_plate_ids)
  if (length(missing_counts) > 0) {
    stop(
      "These plates exist in metadata but not in combined_counts: ",
      paste(missing_counts, collapse = ", ")
    )
  }
  
  # Build structure for each batch
  Exp <- lapply(unique(meta$batch), function(batch_id) {
    # Plate barcodes in this batch
    plates <- unique(meta$plate_barcode[meta$batch == batch_id])
    
    # Build per-plate UMI + annotation list
    res <- lapply(plates, function(plate_id) {
      # 1) Extract UMI counts for this plate
      UMI_decode <- as.data.frame(combined_counts[[plate_id]])
      
      # 2) Extract metadata rows for this plate
      annot <- meta[meta$plate_barcode == plate_id, ]
      annot <- as.data.frame(annot)
      
      # 3) Set annotation rownames == "<plate>_<well_index>"
      rownames(annot) <- apply(
        annot[, c("plate_barcode", "well_index")],
        1,
        paste,
        collapse = "_"
      )
      
      # 4) Return a clean list
      list(
        UMI.counts = UMI_decode,
        Annotation = annot
      )
    })
    
    names(res) <- plates
    return(res)
  })
  
  names(Exp) <- unique(meta$batch)
  
  return(Exp)
}

# --------------------------------------------------------------
# Heatmaps DE
# --------------------------------------------------------------

get_top_genes_drugseq <- function(de_table, n = 30) {
  de_table %>%
    dplyr::arrange(padj, desc(abs(LFC) * -log10(padj))) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(ENSEMBLID)
}

#' Static PNG thumbnail for a compound DE heatmap, rendered without a browser.
#'
#' The interactive heatmap is a heatmaply widget; the old thumbnail was a
#' webshot2 screenshot of the iframe wrapper, which renders inconsistently in
#' headless/network-restricted environments (some compounds got a thumbnail,
#' others silently fell back to the generic HTML icon). This redraws the same
#' z-scored matrix with base R's stats::heatmap (row clustering, columns left
#' in dose order to mirror dendrogram = "row"), so every compound gets a
#' deterministic thumbnail offline. Non-fatal: any failure leaves no PNG and
#' the report falls back to its icon.
render_heatmap_thumbnail <- function(m, title, png_file,
                                     palette = grDevices::colorRampPalette(
                                       rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)) {
  grDevices::png(filename = png_file, width = 1000, height = 800, res = 110)
  on.exit(grDevices::dev.off(), add = TRUE)
  ok <- FALSE
  tryCatch({
    if (nrow(m) >= 2 && ncol(m) >= 2) {
      stats::heatmap(
        m,
        Rowv    = TRUE,       # cluster genes (rows)
        Colv    = NA,         # keep columns in dose order
        scale   = "none",     # matrix is already z-scored
        col     = palette,
        margins = c(9, 5),
        labRow  = if (nrow(m) > 50) NA else rownames(m),
        cexRow  = 0.5,
        cexCol  = 0.7,
        main    = title
      )
      ok <- TRUE
    }
  }, error = function(e) {})
  if (!ok) {
    # Fallback for tiny matrices or clustering errors: plain image, no tree.
    op <- graphics::par(mar = c(7, 5, 4, 2))
    graphics::image(t(m[nrow(m):1, , drop = FALSE]),
                    col = palette, axes = FALSE, main = title)
    graphics::par(op)
  }
  invisible(png_file)
}

plot_drugseq_compound_heatmap <- function(
    compound_name,
    de_results,
    experiment_object,
    top_n = 30,
    out_dir,
    aggregate_by = c("none", "concentration", "replicate", "timepoint"),
    meta_experiment = NULL
) {
  
  aggregate_by <- match.arg(aggregate_by)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  header_html <- if (!is.null(meta_experiment)) {
    build_simple_header(meta_experiment)
  } else {
    ""
  }
  
  message("Processing compound: ", compound_name)
  
  # ---- Filter DE results ----
  compound_de <- de_results %>%
    dplyr::filter(treatment == compound_name)
  
  if (nrow(compound_de) == 0) {
    message("No DE results for compound: ", compound_name)
    return(NULL)
  }
  
  top_genes <- unique(get_top_genes_drugseq(compound_de, n = top_n))
  
  expr_list <- list()
  
  # ---- Collect expression from all plates ----
  for (batch in names(experiment_object)) {
    
    for (plate in names(experiment_object[[batch]])) {
      
      umi <- experiment_object[[batch]][[plate]]$UMI.counts
      meta <- experiment_object[[batch]][[plate]]$Annotation
      
      expr_long <- umi %>%
        dplyr::filter(gene %in% top_genes) %>%
        tidyr::pivot_longer(
          cols = -c(gene, sum),
          names_to = "well",
          values_to = "expr"
        ) %>%
        dplyr::mutate(
          well_barcode = sub("^[^_]+_", "", well)
        )
      
      expr_annot <- expr_long %>%
        dplyr::left_join(meta, by = c("well_barcode" = "well_index")) %>%
        dplyr::filter(treatment == compound_name) %>%
        dplyr::mutate(
          plate_barcode = plate,
          batch_id = batch
        )
      
      if (nrow(expr_annot) > 0) {
        expr_list[[length(expr_list) + 1]] <- expr_annot
      }
      
    }
  }
  
  if (length(expr_list) == 0) {
    message("No expression data found for compound: ", compound_name)
    return(NULL)
  }
  
  expr_all <- dplyr::bind_rows(expr_list)
  
  # ---- Aggregation modes ----
  if (aggregate_by == "none") {
    
    expr_summary <- expr_all %>%
      dplyr::mutate(
        colname = paste0(
          treatment_dose,
          treatment_dose_unit,
          "_",
          plate_barcode,
          "_",
          well_barcode
        )
      ) %>%
      dplyr::select(gene, colname, expr)
    
  } else if (aggregate_by == "replicate") {
    
    expr_summary <- expr_all %>%
      dplyr::group_by(
        gene,
        treatment_dose,
        treatment_dose_unit,
        plate_barcode,
        replicate
      ) %>%
      dplyr::summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(
        colname = paste0(
          treatment_dose,
          treatment_dose_unit,
          "_",
          replicate,
          "_",
          plate_barcode
        )
      )
    
  } else if (aggregate_by == "concentration") {
    
    expr_summary <- expr_all %>%
      dplyr::group_by(
        gene,
        treatment_dose,
        treatment_dose_unit
      ) %>%
      dplyr::summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(
        colname = paste0(
          treatment_dose,
          treatment_dose_unit
        )
      )
    
  } else if (aggregate_by == "timepoint") {
    
    expr_summary <- expr_all %>%
      dplyr::group_by(
        gene,
        treatment_dose,
        treatment_dose_unit,
        timepoint
      ) %>%
      dplyr::summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(
        colname = paste0(
          timepoint,
          "_",
          treatment_dose,
          treatment_dose_unit
        )
      )
    
  }
  
  # ---- Pivot to matrix ----
  expr_wide <- expr_summary %>%
    dplyr::select(gene, colname, expr) %>%
    tidyr::pivot_wider(
      names_from = colname,
      values_from = expr,
      values_fn = mean
    )
  
  combined <- as.data.frame(expr_wide)
  
  rownames(combined) <- combined$gene
  combined$gene <- NULL
  
  # ---- Order columns by concentration ----
  col_order <- tibble::tibble(
    colname = colnames(combined)
  ) %>%
    dplyr::mutate(
      dose = readr::parse_number(colname)
    ) %>%
    dplyr::arrange(dose) %>%
    dplyr::pull(colname)
  
  combined <- combined[, col_order, drop = FALSE]
  
  # ---- Map gene symbols ----
  symbol_map <- compound_de %>%
    dplyr::select(ENSEMBLID, SYMBOL) %>%
    dplyr::distinct()
  
  symbol_vec <- symbol_map$SYMBOL[
    match(rownames(combined), symbol_map$ENSEMBLID)
  ]
  
  symbol_vec[is.na(symbol_vec)] <- rownames(combined)[is.na(symbol_vec)]
  
  rownames(combined) <- make.unique(symbol_vec)
  
  # ---- Scale per gene ----
  m_mat <- log2(as.matrix(combined) + 1)
  
  if (ncol(m_mat) <= 2) {
    
    m_scaled <- m_mat
    value_label <- "Log2(UMI + 1)"
    
  } else {
    
    m_scaled <- t(scale(t(m_mat)))
    m_scaled[!is.finite(m_scaled)] <- 0
    value_label <- "Z-score (log2 expression)"
    
  }
  # ---- Title ----
  heatmap_title <- if (!is.null(meta_experiment)) {
    paste0(
      meta_experiment$title,
      " – Compound: ",
      compound_name,
      " (mode: ",
      aggregate_by,
      ")"
    )
  } else {
    paste0(
      "DE Heatmap: ",
      compound_name,
      " (mode: ",
      aggregate_by,
      ")"
    )
  }
  
  # ---- Heatmap ----
  ht <- heatmaply::heatmaply(
    m_scaled,
    column_text_angle = 90,
    main = heatmap_title,
    label_names = c("Gene", "Condition", value_label),
    colors = colorRampPalette(
      rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
    )(100),
    dendrogram = "row"
  )
  
  html_file <- file.path(
    out_dir,
    paste0("heatmap_", compound_name, "_", aggregate_by, ".html")
  )
  
  png_file <- file.path(
    out_dir,
    paste0("heatmap_", compound_name, "_", aggregate_by, ".png")
  )
  
  widget_file <- file.path(
    out_dir,
    paste0("widget_", compound_name, "_", aggregate_by, ".html")
  )
  
  htmlwidgets::saveWidget(
    plotly::as_widget(ht),
    widget_file,
    selfcontained = TRUE
  )
  
  full_html <- glue::glue('
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"/>
<title>{heatmap_title}</title>
</head>

<body style="font-family: \'Helvetica Neue\', Helvetica, Arial, sans-serif; margin:20px;">

{header_html}

<h4>{heatmap_title}</h4>

<iframe src="{basename(widget_file)}"
        style="width:100%; height:800px; border:none;">
</iframe>

</body>
</html>
')
  
  writeLines(full_html, html_file)
  
  # ---- PNG thumbnail (offline, no browser) ----
  # Rendered directly from the z-scored matrix so it always succeeds, instead
  # of screenshotting the interactive iframe (which failed unpredictably).
  output_files <- html_file
  tryCatch({
    render_heatmap_thumbnail(m_scaled, heatmap_title, png_file)
    if (file.exists(png_file)) {
      message("PNG saved to: ", png_file)
      output_files <- c(html_file, png_file)
    }
  }, error = function(e) {
    message("\u26a0 heatmap thumbnail failed (", conditionMessage(e),
            ") \u2014 skipping PNG")
  })
  
  message("Saved to: ", html_file)
  invisible(output_files)
  
}

# --------------------------------------------------------------
# Volcano plots
# --------------------------------------------------------------

#' PNG preview next to a volcano SVG (<stem>.png), ~1100 px on the long side.
#' The final report picks it up for the card thumbnail; the SVG stays the
#' full-size / download version. Non-fatal on any device error.
write_volcano_preview <- function(plot, svg_path, width, height, px = 1100) {
  png_path <- sub("\\.svg$", ".png", svg_path, ignore.case = TRUE)
  if (identical(png_path, svg_path)) return(invisible(NULL))
  tryCatch(
    ggplot2::ggsave(png_path, plot, device = "png", width = width, height = height,
                    dpi = px / max(width, height), bg = "white"),
    error = function(e) {
      message("  volcano PNG preview failed (", conditionMessage(e), ") - ",
              basename(svg_path), " will be rasterised in the browser")
      if (file.exists(png_path)) unlink(png_path)
    }
  )
  invisible(png_path)
}

plot_volcano <- function(
    de_result,
    contrast,
    output_dir,
    PADJ_THRESHOLD = 0.05,
    LFC_THRESHOLD = 1,
    meta_experiment = NULL
) {
  # -------------------------------------------------------
  # 1. Filter DE results for the selected contrast
  # -------------------------------------------------------
  
  safe_name <- function(x) {
    gsub("[^A-Za-z0-9_-]", "-", x)
  }
  
  contrast_safe <- safe_name(contrast)
  
  df_plot <- dplyr::filter(de_result, Sample == contrast)
  
  meta <- df_plot |>
    dplyr::select(
      batch,
      treatment,
      treatment_dose,
      treatment_dose_unit,
      timepoint
    ) |>
    dplyr::distinct()
  
  if (nrow(meta) > 1) {
    meta <- meta[1, ]
  }
  
  title_txt <- if (!is.null(meta_experiment)) {
    glue::glue(
      "{meta_experiment$title} – {meta$treatment} {meta$treatment_dose}{meta$treatment_dose_unit} – {meta$timepoint}"
    )
  } else {
    glue::glue(
      "{meta$treatment} {meta$treatment_dose}{meta$treatment_dose_unit} – {meta$timepoint}"
    )
  }
  
  subtitle_txt <- glue::glue(
    "Batch {meta$batch}"
  )
  
  if (nrow(df_plot) == 0) {
    cli::cli_alert_warning("No DE results found for contrast {contrast}. Volcano skipped.")
    return(character(0))
  }
  
  n_up <- sum(df_plot$padj < PADJ_THRESHOLD & df_plot$LFC > LFC_THRESHOLD, na.rm = TRUE)
  n_down <- sum(df_plot$padj < PADJ_THRESHOLD & df_plot$LFC < -LFC_THRESHOLD, na.rm = TRUE)
  
  caption_txt <- glue::glue(
    "Up: {n_up} | Down: {n_down} | padj < {PADJ_THRESHOLD}, |LFC| > {LFC_THRESHOLD}"
  )
  
  # -------------------------------------------------------
  # 2. Prepare data for EnhancedVolcano
  # -------------------------------------------------------
  df_plot$minus_log10_adjP <- -log10(df_plot$padj)
  
  # To avoid label problems, pick the 50 top genes by adj.P.Val
  top50_genes <- df_plot |>
    dplyr::arrange(padj) |>
    dplyr::slice(1:50) |>
    dplyr::pull(SYMBOL)
  
  # -------------------------------------------------------
  # 3. Create volcano plot (ggplot object)
  # -------------------------------------------------------
  v_plot <- EnhancedVolcano::EnhancedVolcano(
    toptable = df_plot,
    lab = df_plot$SYMBOL,
    x = "LFC",
    y = "padj",
    selectLab = top50_genes,
    pCutoff = PADJ_THRESHOLD,
    FCcutoff = LFC_THRESHOLD,
    # EnhancedVolcano renders the title through ggplot2, which draws text
    # literally — HTML like <br>/<span> would print as raw tags (and ggplotly
    # carries the same string into the interactive title). Keep the title
    # plain; the up/down counts live in the subtitle below.
    title = glue::glue("Volcano plot {contrast} vs DMSO"),
    subtitle = caption_txt,
    drawConnectors = TRUE,
    labSize = 4,
    maxoverlapsConnectors = 50
  ) +
    {
      lfc_max <- max(abs(df_plot$LFC), na.rm = TRUE)
      x_lim   <- ceiling(min(lfc_max + 1, 15))
      list(
        ggplot2::coord_cartesian(xlim = c(-x_lim, x_lim)),
        ggplot2::scale_x_continuous(breaks = seq(-x_lim, x_lim, by = if (x_lim > 10) 2 else 1))
      )
    } +
    ggplot2::labs(
      x = "log2 fold change",
      y = "-log10(padj)"
    )
  
  # -------------------------------------------------------
  # 4. Save SVG version
  # -------------------------------------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  header_html <- if (!is.null(meta_experiment)) {
    build_simple_header(meta_experiment)
  } else {
    ""
  }
  
  svg_path <- file.path(output_dir, paste0(contrast_safe, ".svg"))
  ggplot2::ggsave(svg_path, v_plot, device = "svg", width = 15, height = 15)
  harmonize_svg_fonts(svg_path)
  # Lightweight raster twin shown on the final-report card (the SVG, one path
  # per gene, made the report stutter on scroll). Non-fatal: without it the
  # report falls back to rasterising the SVG in the browser.
  write_volcano_preview(v_plot, svg_path, width = 15, height = 15)
  
  # -------------------------------------------------------
  # 5. Save interactive HTML version (self-contained plotly.js)
  # -------------------------------------------------------
  # Embeds this contrast's gene table as JSON and draws with plotly.js. Two
  # sliders (padj, |LFC|) recolour points live and update the up/down counts,
  # so the reader can re-threshold in-browser without re-running the pipeline.
  # The pipeline's thresholds are the slider start values; only which points are
  # *called* significant changes. Built with paste0 (not glue) because the
  # embedded JavaScript uses { } and [ ] heavily and would collide with any glue
  # delimiter choice. The static SVG above remains the fixed-threshold view.
  gene_df <- df_plot %>%
    dplyr::transmute(
      symbol  = dplyr::coalesce(as.character(SYMBOL), ENSEMBLID),
      ensembl = ENSEMBLID,
      lfc     = round(LFC, 3),
      padj    = signif(padj, 4),
      mlog    = round(-log10(padj), 3),
      wells   = if ("plate_well" %in% names(df_plot)) as.character(plate_well) else NA_character_
    ) %>%
    dplyr::filter(is.finite(lfc), is.finite(mlog))
  
  gene_json <- jsonlite::toJSON(gene_df, dataframe = "columns", na = "null", auto_unbox = TRUE)
  
  html_path <- file.path(output_dir, paste0(contrast_safe, ".html"))
  
  full_html <- paste0(
    '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>', title_txt,
    '</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body { font-family:"Helvetica Neue",Helvetica,Arial,sans-serif; margin:20px; color:#111827; }
  h4 { margin:0 0 4px; }
  .sub { color:#6b7280; margin:0 0 16px; }
  .controls { display:flex; gap:28px; flex-wrap:wrap; align-items:center;
              background:#f9fafb; border:1px solid #e5e7eb; border-radius:10px;
              padding:14px 18px; margin-bottom:14px; }
  .ctrl { display:flex; flex-direction:column; min-width:230px; }
  .ctrl label { font-size:12px; font-weight:600; color:#374151; margin-bottom:6px; }
  .ctrl .val { font-family:ui-monospace,monospace; color:#0077c8; }
  .counts { font-size:13px; }
  .counts b.up { color:#b91c1c; } .counts b.down { color:#1d4ed8; }
  input[type=range] { width:230px; accent-color:#0077c8; }
</style>
</head>
<body>
', header_html,
    '
<h4>', title_txt,
    '</h4>
<p class="sub">', subtitle_txt,
    '</p>

<div class="controls">
  <div class="ctrl">
    <label>padj &le; <span class="val" id="padjVal"></span></label>
    <input type="range" id="padjSlider" min="0" max="1" step="0.005"/>
  </div>
  <div class="ctrl">
    <label>|log2FC| &ge; <span class="val" id="lfcVal"></span></label>
    <input type="range" id="lfcSlider" min="0" max="6" step="0.05"/>
  </div>
  <div class="counts">
    <b class="up" id="nUp">0</b> up &middot; <b class="down" id="nDown">0</b> down
    &middot; <span id="nTot">0</span> tested
  </div>
</div>

<div id="volcano" style="width:100%;height:720px;"></div>

<script>
const G = ', gene_json,
    ';
let padjCut = ', PADJ_THRESHOLD,
    ';
let lfcCut  = ', LFC_THRESHOLD,
    ';
const n = G.lfc.length;
const idx = Array.from({length:n}, (_,i)=>i);
function classify(i) {
  const sig = (G.padj[i] !== null) && (G.padj[i] < padjCut);
  if (sig && G.lfc[i] >=  lfcCut) return "up";
  if (sig && G.lfc[i] <= -lfcCut) return "down";
  if (sig) return "p";
  return "ns";
}
const COLORS = { up:"#b91c1c", down:"#1d4ed8", p:"#059669", ns:"rgba(120,120,120,0.45)" };
function hover(i) {
  return "<b>"+G.symbol[i]+"</b><br>"+G.ensembl[i]+
         "<br>log2FC: "+G.lfc[i]+"<br>padj: "+G.padj[i]+
         (G.wells[i] ? "<br>Wells: "+G.wells[i] : "")+"<extra></extra>";
}
function build() {
  const cls = idx.map(classify);
  const trace = {
    x: G.lfc, y: G.mlog, mode:"markers", type:"scattergl",
    marker:{ size:6, color: cls.map(c=>COLORS[c]) },
    text: idx.map(hover), hovertemplate:"%{text}"
  };
  const shapes = [
    {type:"line", x0:lfcCut,  x1:lfcCut,  yref:"paper", y0:0, y1:1, line:{dash:"dash", width:1, color:"#9ca3af"}},
    {type:"line", x0:-lfcCut, x1:-lfcCut, yref:"paper", y0:0, y1:1, line:{dash:"dash", width:1, color:"#9ca3af"}},
    {type:"line", yref:"y", y0:-Math.log10(padjCut), y1:-Math.log10(padjCut), xref:"paper", x0:0, x1:1, line:{dash:"dash", width:1, color:"#9ca3af"}}
  ];
  const layout = {
    font:{ family:"Helvetica Neue, Helvetica, Arial, sans-serif" },
    xaxis:{ title:"log2 fold change", zeroline:false },
    yaxis:{ title:"-log10(padj)" },
    shapes: shapes, hovermode:"closest", showlegend:false, margin:{ t:20 }
  };
  Plotly.react("volcano", [trace], layout, {responsive:true});
  let up=0, down=0;
  for (const c of cls) { if (c==="up") up++; else if (c==="down") down++; }
  document.getElementById("nUp").textContent = up;
  document.getElementById("nDown").textContent = down;
  document.getElementById("nTot").textContent = n;
}
function syncLabels() {
  document.getElementById("padjVal").textContent = padjCut.toPrecision(2);
  document.getElementById("lfcVal").textContent  = lfcCut.toFixed(2);
}
const padjS = document.getElementById("padjSlider");
const lfcS  = document.getElementById("lfcSlider");
padjS.value = padjCut; lfcS.value = lfcCut;
syncLabels();
padjS.addEventListener("input", e => { padjCut = +e.target.value; syncLabels(); build(); });
lfcS.addEventListener("input",  e => { lfcCut  = +e.target.value; syncLabels(); build(); });
build();
</script>
</body>
</html>
'
  )
  
  writeLines(full_html, html_path)
  
  return(c(svg_path, html_path))
}


# --------------------------------------------------------------
# Volcano plot static
# --------------------------------------------------------------

volcano_plot_static <- function(
    this_contrast_toptabel,
    batch_select,
    this_arms,
    p_val,
    lfc_val,
    save_dir = NULL,
    file_prefix = "volcano"
) {
  # -------------------------
  # SETUP
  # -------------------------
  pcut <- p_val
  
  color_custom <- c(
    "not significant"     = rgb(0, 0, 0, 0.5),
    "upregulated genes"   = rgb(1, 0, 0, 1),
    "downregulated genes" = rgb(0, 0, 1, 1)
  )
  
  # -------------------------
  # FILTER + PLOT + SAVE PER BATCH
  # -------------------------
  this_contrast_toptabel <- this_contrast_toptabel %>%
    dplyr::filter(treatment %in% this_arms, batch %in% batch_select)
  
  if (nrow(this_contrast_toptabel) == 0) {
    warning("No data for treatment ", this_arms, " batches ",
            paste(batch_select, collapse = ", "))
    return(character(0))
  }
  
  safe_name <- function(x) gsub("[^A-Za-z0-9_-]", "-", x)
  out_paths <- character(0)
  
  for (b in batch_select) {
    sub_tab <- this_contrast_toptabel %>% dplyr::filter(batch == b)
    if (nrow(sub_tab) == 0) next
    
    sub_tab <- sub_tab %>%
      dplyr::mutate(
        label = dplyr::case_when(
          padj < pcut & LFC >  lfc_val ~ "upregulated genes",
          padj < pcut & LFC < -lfc_val ~ "downregulated genes",
          TRUE                         ~ "not significant"
        ),
        label = factor(label, levels = c("not significant",
                                         "upregulated genes",
                                         "downregulated genes"))
      )
    
    vp <- sub_tab %>%
      dplyr::arrange(label) %>%
      ggplot2::ggplot(ggplot2::aes(x = LFC, y = -log10(padj), color = label)) +
      ggplot2::geom_point(fill = "white", shape = 21, size = 1, stroke = 0.5) +
      ggplot2::geom_hline(yintercept = -log10(pcut), linewidth = 0.3, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = c(-lfc_val, lfc_val), linewidth = 0.3, linetype = "dashed") +
      ggplot2::xlim(-8, 8) +
      ggplot2::scale_color_manual(values = color_custom) +
      ggplot2::facet_grid(treatment ~ treatment_dose) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title    = this_arms,
        subtitle = glue::glue("Batch {b} | padj < {pcut}, |LFC| > {lfc_val}"),
        x = "log2 fold change", y = "-log10(padj)"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(face = "bold", size = 14),
        axis.text.y = ggplot2::element_text(face = "bold", size = 10),
        axis.title  = ggplot2::element_text(size = 14, face = "bold"),
        legend.text = ggplot2::element_text(size = 12),
        strip.text  = ggplot2::element_text(size = 12, face = "bold")
      )
    
    if (!is.null(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      file_svg <- file.path(
        save_dir,
        paste0(file_prefix, "_", safe_name(this_arms), "_batch", b, ".svg")
      )
      message("Saving volcano plot → ", file_svg)
      ggplot2::ggsave(file_svg, vp, device = "svg", width = 12, height = 12)
      harmonize_svg_fonts(file_svg)
      write_volcano_preview(vp, file_svg, width = 12, height = 12)
      out_paths <- c(out_paths, file_svg)
    }
  }
  
  if (!is.null(save_dir)) return(out_paths)
  vp   # last plot object when not saving
}


# --------------------------------------------------------------
# Upset plot
# --------------------------------------------------------------

make_lt <- function(
    df,
    treatments = NULL,
    doses = NULL,
    logFC_cutoff,
    padj_cutoff
) {
  df_filt <- df %>%
    dplyr::filter(
      abs(LFC) > logFC_cutoff,
      padj < padj_cutoff
    )
  
  if (!is.null(treatments)) {
    df_filt <- df_filt %>%
      dplyr::filter(treatment %in% treatments)
  }
  
  if (!is.null(doses)) {
    df_filt <- df_filt %>%
      dplyr::filter(treatment_dose %in% doses)
  }
  
  lt <- df_filt %>%
    dplyr::group_by(treatment, treatment_dose) %>%
    dplyr::summarise(
      genes = list(unique(ENSEMBLID)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      name = paste(treatment, treatment_dose, sep = "_")
    ) %>%
    dplyr::select(name, genes) %>%
    tibble::deframe()
  
  if (length(lt) == 0) {
    warning("No sets matched the filtering criteria.")
  }
  
  lt
}

build_upset_title <- function(
    treatments,
    padj_cutoff,
    lfc_cutoff
) {
  # Listing every treatment produces an unreadable, page-wide title on large
  # screens (hundreds of conditions). Summarise instead: name them only when
  # there are a few, otherwise give a count.
  n <- length(treatments)
  lead <- if (n <= 3) {
    paste(treatments, collapse = " vs ")
  } else {
    sprintf("%d treatment sets", n)
  }
  sprintf("DE-gene overlap: %s  |  padj < %s, |LFC| > %s",
          lead, padj_cutoff, lfc_cutoff)
}

save_upset_plot <- function(
    lt_list,
    output_dir,
    filename,
    title = NULL,
    width = 12,
    height = 8,
    meta_experiment = NULL,
    max_sets = 40
) {
  if (length(lt_list) == 0) {
    warning("UpSet plot skipped: empty list.")
    return(character(0))
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  
  # ---------------------------------------------------------------
  # Guard against the filesystem's per-filename byte limit (255 on
  # ext4 and most others). With a select-all config the caller joins
  # every treatment ("A_vs_B_vs_...") into `filename`, which can exceed
  # the limit and fail with "cannot open stream". Shorten it here — before
  # file_path, html_file, and the <img src> are derived — so all three stay
  # mutually consistent. upset_plot is a single (non-branched) target, so
  # there is one file per run and no collision risk.
  # ---------------------------------------------------------------
  if (nchar(filename, type = "bytes") > 200L) {
    ext      <- tools::file_ext(filename)
    stem     <- tools::file_path_sans_ext(filename)
    short    <- paste0(substr(stem, 1, 120), "_plus_more_",
                       length(lt_list), "sets")
    filename <- if (nzchar(ext)) paste0(short, ".", ext) else short
    message("  [upset] filename exceeded filesystem limit - shortened to ",
            filename)
  }
  
  # UpSetR is unreadable (and very slow) beyond a few dozen sets. On a large
  # screen with no `selection$treatments` restriction, lt_list is one set per
  # treatment x dose and can run to hundreds. Write a placeholder instead of
  # hanging the pipeline; restrict `selection$treatments` to get a real plot.
  if (length(lt_list) > max_sets) {
    file_path <- file.path(output_dir, filename)
    msg <- sprintf(
      "UpSet plot skipped: %d treatment x dose sets exceed max_sets = %d.\nRestrict selection$treatments in the config to plot a subset.",
      length(lt_list), max_sets)
    message("  [upset] ", gsub("\n", " ", msg))
    ph <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 5, lineheight = 1.1) +
      ggplot2::theme_void()
    ggplot2::ggsave(file_path, ph, device = "svg", width = 10, height = 4)
    harmonize_svg_fonts(file_path)
    return(invisible(c(svg = file_path)))
  }
  
  header_html <- if (!is.null(meta_experiment)) {
    build_simple_header(meta_experiment)
  } else {
    ""
  }
  
  file_path <- file.path(output_dir, filename)
  
  # Use svglite (pure C++, no cairo/X11) rather than grDevices::svg(), which
  # needs a cairo-backed display and fails headless with "unable to start
  # device 'svg'". svglite captures UpSetR's grid graphics and the grid.text
  # title overlay just as well, and matches the device used everywhere else
  # in the pipeline (ggsave device = "svg").
  svglite::svglite(file = file_path, width = width, height = height)
  # Close the device and normalise fonts once the file is fully written.
  on.exit({
    grDevices::dev.off()
    harmonize_svg_fonts(file_path)
  }, add = TRUE)
  
  print(
    UpSetR::upset(
      UpSetR::fromList(lt_list),
      nsets = length(lt_list)
    )
  )
  
  if (!is.null(title)) {
    grid::grid.text(
      title,
      x = grid::unit(0.60, "npc"),
      y = grid::unit(0.95, "npc"),
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  }
  
  html_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(filename), ".html")
  )
  
  full_html <- glue::glue('
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>UpSet Plot</title>
</head>

<body style="font-family: \'Helvetica Neue\', Helvetica, Arial, sans-serif; margin:20px;">

{header_html}

<h4>{title}</h4>

<div style="margin-top:20px;">
  <img src="{basename(filename)}" style="max-width:100%;">
</div>

</body>
</html>
')
  
  writeLines(full_html, html_file)
  
  invisible(c(svg = file_path, html = html_file))
}


# --------------------------------------------------------------
# Correlation plot
# --------------------------------------------------------------

ggpairs_custom_tidy <- function(
    toptable_meta,
    this_arms,
    padj_threshold,
    logFC_threshold,
    corr_size = 3,
    batch_s,
    strip_size = 8,
    label_size = 10,
    exclusive_gene = NULL,
    save_dir = NULL,
    file_prefix = "pairs",
    meta_experiment = NULL,
    treatment_name = NULL
) {
  #---------------------------------------------------------
  # 1. Construct contrast column
  #---------------------------------------------------------
  toptable <- toptable_meta %>%
    dplyr::mutate(
      contrast = paste(treatment, treatment_dose, timepoint, sep = "_")
    )
  
  contrasts <- toptable %>%
    dplyr::filter(
      Sample %in% this_arms,
      batch %in% batch_s
    ) %>%
    dplyr::arrange(treatment_dose) %>%
    dplyr::pull(contrast) %>%
    unique()
  
  #---------------------------------------------------------
  # 2. Identify valid genes
  #---------------------------------------------------------
  valid_genes <- toptable %>%
    dplyr::filter(
      Sample %in% this_arms,
      batch %in% batch_s,
      padj < padj_threshold
    ) %>%
    dplyr::pull(ENSEMBLID) %>%
    unique()
  
  if (!is.null(exclusive_gene)) {
    valid_genes <- setdiff(valid_genes, exclusive_gene)
  }
  
  #---------------------------------------------------------
  # 2b. Early exits — skip gracefully instead of crashing
  #
  # With a "select everything" config, ggpairs_plots_multiple maps over
  # every treatment. Two situations must not abort the whole pipeline:
  #   (a) a treatment with < 2 contrasts  -> a pairwise plot is undefined
  #       (this is the `i may only be in the range from 1:1` error, from
  #        the 2:pm$nrow loop when pm has a single column), and
  #   (b) a treatment with no significant genes -> valid_genes is empty,
  #       so the wide tables carry no contrast columns and the later
  #       all_of(contrasts) select fails with "Elements ... don't exist".
  # In both cases we write a labelled placeholder SVG (so `format = "file"`
  # branches still return a valid path and the final report has something
  # to link) and return early.
  #---------------------------------------------------------
  treat_label <- if (!is.null(treatment_name)) treatment_name else "allTreatments"
  
  save_placeholder <- function(reason) {
    if (is.null(save_dir)) {
      message("  [ggpairs] ", treat_label, ": ", reason, " - skipping.")
      return(invisible(character(0)))
    }
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    safe_name <- function(x) gsub("[^A-Za-z0-9_-]", "-", x)
    file_svg <- file.path(
      save_dir,
      paste0(
        file_prefix, "_",
        safe_name(treat_label), "_",
        "batch", paste(sort(batch_s), collapse = "-"),
        ".svg"
      )
    )
    ph <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 0,
        label = paste0(treat_label, "\n(", reason, ")"),
        size = 5, lineheight = 1.1
      ) +
      ggplot2::theme_void()
    ggplot2::ggsave(file_svg, ph, device = "svg", width = 10, height = 6)
    message("  [ggpairs] ", treat_label, ": ", reason,
            " - wrote placeholder ", file_svg)
    invisible(file_svg)
  }
  
  # (a) a pairwise correlation matrix needs at least two contrasts to compare
  if (length(contrasts) < 2) {
    return(save_placeholder(sprintf(
      "only %d contrast(s) available - need >= 2 for a pairs plot",
      length(contrasts)
    )))
  }
  
  # (b) nothing to plot if no gene is significant anywhere in this treatment
  if (length(valid_genes) == 0) {
    return(save_placeholder(sprintf(
      "no genes pass padj < %g", padj_threshold
    )))
  }
  
  #---------------------------------------------------------
  # 3. Helper for reshaping
  #---------------------------------------------------------
  logfc_wide <- function(df, genes = NULL) {
    df %>%
      dplyr::filter(
        Sample %in% this_arms,
        batch %in% batch_s,
        if (!is.null(genes)) ENSEMBLID %in% genes else TRUE
      ) %>%
      dplyr::select(ENSEMBLID, contrast, LFC) %>%
      dplyr::group_by(ENSEMBLID, contrast) %>%
      dplyr::summarise(LFC = mean(LFC, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = contrast, values_from = LFC)
  }
  
  table_all <- logfc_wide(toptable, valid_genes)
  
  table_sig <- toptable %>%
    dplyr::filter(
      Sample %in% this_arms,
      batch %in% batch_s,
      padj < padj_threshold
    ) %>%
    logfc_wide()
  
  #---------------------------------------------------------
  # 4. Pad missing columns
  #
  # A contrast can still be absent from a wide table when none of its genes
  # survived that contrast's per-DE filter (table_all) or when no gene was
  # significant in it (table_sig). Pad every expected contrast column with
  # NA in BOTH tables so the all_of(contrasts) select below can never fail
  # on a missing element.
  #---------------------------------------------------------
  pad_to_contrasts <- function(tbl) {
    miss <- setdiff(contrasts, colnames(tbl))
    for (m in miss) tbl[[m]] <- NA_real_
    tbl
  }
  
  table_all <- pad_to_contrasts(table_all) %>%
    dplyr::select(ENSEMBLID, dplyr::all_of(contrasts))
  
  table_sig <- pad_to_contrasts(table_sig) %>%
    dplyr::select(ENSEMBLID, dplyr::all_of(contrasts))
  
  #---------------------------------------------------------
  # 5. Merge matrices
  #---------------------------------------------------------
  data_plot <- table_all %>%
    dplyr::inner_join(
      table_sig,
      by = "ENSEMBLID",
      suffix = c(".x", ".y")
    )
  
  df <- data_plot %>% dplyr::select(-ENSEMBLID)
  
  axis_lim <- df %>%
    unlist() %>%
    abs() %>%
    max(na.rm = TRUE)
  
  #---------------------------------------------------------
  # 6. Panel functions
  #---------------------------------------------------------
  diagfun <- function(data, mapping) {
    ggplot2::ggplot(data, mapping) +
      ggplot2::geom_density(color = "darkblue", fill = "lightblue") +
      ggplot2::scale_x_continuous(limits = c(-axis_lim, axis_lim))
  }
  
  corr_heatmap <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    corr <- cor(x, y, use = "pairwise.complete.obs")
    
    colFn <- colorRampPalette(c("blue", "white", "red"))
    fill <- colFn(100)[findInterval(corr, seq(-1, 1, length.out = 100))]
    
    GGally::ggally_cor(data = data, mapping = mapping, ...) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = fill, color = fill)
      )
  }
  
  #---------------------------------------------------------
  # 7. Base ggpairs object
  #---------------------------------------------------------
  pm <- GGally::ggpairs(
    df,
    columns = seq_along(contrasts),
    columnLabels = contrasts,
    diag = list(continuous = GGally::wrap(diagfun)),
    upper = list(continuous = corr_heatmap),
    progress = FALSE
  ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = label_size),
      text = ggplot2::element_text(size = strip_size)
    )
  
  #---------------------------------------------------------
  # 8. Custom lower panels
  #---------------------------------------------------------
  color_vec <- c(
    "not significant" = rgb(0, 0, 0, 0.5),
    "x axis DE genes" = rgb(0, 1, 0, 0.7),
    "y axis DE genes" = rgb(0, 0, 1, 0.7),
    "shared DE genes" = rgb(1, 0, 0, 1)
  )
  
  for (i in 2:pm$nrow) {
    for (j in 1:(i - 1)) {
      panel <- pm[i, j]
      
      xname <- rlang::as_name(panel$mapping$x)
      yname <- rlang::as_name(panel$mapping$y)
      
      base_x <- stringr::str_remove(xname, "\\.[^.]+$")
      base_y <- stringr::str_remove(yname, "\\.[^.]+$")
      
      dt <- data_plot %>%
        dplyr::transmute(
          x = .data[[paste0(base_x, ".x")]],
          y = .data[[paste0(base_y, ".x")]],
          dx = .data[[paste0(base_x, ".y")]],
          dy = .data[[paste0(base_y, ".y")]]
        ) %>%
        dplyr::mutate(
          label = dplyr::case_when(
            abs(dx) > logFC_threshold & abs(dy) > logFC_threshold ~ "shared DE genes",
            abs(dx) > logFC_threshold ~ "x axis DE genes",
            abs(dy) > logFC_threshold ~ "y axis DE genes",
            TRUE ~ "not significant"
          ),
          label = factor(label, levels = names(color_vec))
        )
      
      corr <- cor(dt$x, dt$y, use = "pairwise.complete.obs")
      
      p <- ggplot2::ggplot(dt, ggplot2::aes(x = x, y = y, color = label)) +
        ggplot2::geom_point(size = 2, shape = 21, stroke = 0.5, fill = "white") +
        ggplot2::scale_color_manual(values = color_vec) +
        ggplot2::scale_x_continuous(limits = c(-axis_lim, axis_lim)) +
        ggplot2::scale_y_continuous(limits = c(-axis_lim, axis_lim)) +
        ggplot2::geom_hline(
          yintercept = c(-logFC_threshold, logFC_threshold),
          linetype = "dashed", linewidth = 0.3
        ) +
        ggplot2::geom_vline(
          xintercept = c(-logFC_threshold, logFC_threshold),
          linetype = "dashed", linewidth = 0.3
        ) +
        ggplot2::theme_bw() +
        ggplot2::geom_label(
          data = tibble::tibble(x = -axis_lim, y = axis_lim, lab = round(corr, 3)),
          ggplot2::aes(x = x, y = y, label = lab),
          inherit.aes = FALSE,
          hjust = 0, vjust = 1,
          size = corr_size, fontface = "bold"
        )
      
      pm[i, j] <- p
    }
  }
  
  #---------------------------------------------------------
  # 9. Saving
  #---------------------------------------------------------
  if (!is.null(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Build a small standalone legend that explains the lower-triangle
    # colour scheme. ggpairs/GGally strips per-panel legends, so we render
    # a tiny dummy ggplot purely to extract its legend grob, then compose
    # it with the ggpairs matrix using patchwork.
    legend_df <- data.frame(
      x = seq_along(color_vec),
      y = 0,
      label = factor(names(color_vec), levels = names(color_vec))
    )
    legend_plot <- ggplot2::ggplot(
      legend_df,
      ggplot2::aes(x = x, y = y, color = label)
    ) +
      ggplot2::geom_point(size = 3, shape = 21, stroke = 0.5, fill = "white") +
      ggplot2::scale_color_manual(
        name   = NULL,
        values = color_vec,
        labels = c(
          "not significant"  = "Not significant in either contrast",
          "x axis DE genes"  = "DE only on x-axis contrast",
          "y axis DE genes"  = "DE only on y-axis contrast",
          "shared DE genes"  = "DE in both contrasts (shared)"
        )[names(color_vec)]
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(
        nrow = 1, override.aes = list(size = 4)
      )) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text     = ggplot2::element_text(size = label_size),
        legend.key      = ggplot2::element_blank()
      )
    
    # Compose pm + legend using patchwork if available, else save bare pm.
    final_plot <- if (requireNamespace("patchwork", quietly = TRUE)) {
      # Try to extract just the legend grob via cowplot for a tighter strip;
      # if cowplot is unavailable we embed the whole tiny legend_plot, which
      # still renders correctly (just with a sliver of empty plot area).
      legend_only <- tryCatch(
        cowplot::get_legend(legend_plot),
        error = function(e) NULL
      )
      if (!is.null(legend_only)) {
        patchwork::wrap_elements(GGally::ggmatrix_gtable(pm)) /
          patchwork::wrap_elements(legend_only) +
          patchwork::plot_layout(heights = c(12, 1))
      } else {
        patchwork::wrap_elements(GGally::ggmatrix_gtable(pm)) /
          legend_plot +
          patchwork::plot_layout(heights = c(12, 1))
      }
    } else {
      message("  patchwork not available - saving without legend strip")
      pm
    }
    
    safe_name <- function(x) {
      gsub("[^A-Za-z0-9_-]", "-", x)
    }
    
    treat_part <- if (!is.null(treatment_name)) {
      safe_name(treatment_name)
    } else {
      "allTreatments"
    }
    
    batch_part <- paste(sort(batch_s), collapse = "-")
    
    file_svg <- file.path(
      save_dir,
      paste0(
        file_prefix, "_",
        treat_part, "_",
        "batch", batch_part,
        ".svg"
      )
    )
    
    message("Saving ggpairs plot - ", file_svg)
    
    ggplot2::ggsave(
      file_svg,
      final_plot,
      device = "svg",
      width  = 14,
      height = 14.8     # +0.8 in to accommodate the legend strip
    )
    harmonize_svg_fonts(file_svg)
    
    return(invisible(file_svg))
  }
  
  pm
}


# =============================================================================
# UMAP + clustering utility functions
# =============================================================================


# -----------------------------------------------------------------------------
# Build global UMI matrix (all batches × plates → genes × wells)
# -----------------------------------------------------------------------------
build_global_umi_matrix <- function(Exp) {
  umi_list <- list()
  
  for (b in names(Exp)) {
    for (plate in names(Exp[[b]])) {
      mat <- Exp[[b]][[plate]]$UMI.counts
      mat <- as.data.frame(mat)
      
      rownames(mat) <- mat$gene
      mat <- mat[, !colnames(mat) %in% c("gene", "sum"), drop = FALSE]
      
      # 🔑 Extract well IDs directly from column names
      wells <- sub("_.*$", "", colnames(mat))
      
      # Collapse cell barcodes → wells
      mat_well <- rowsum(t(mat), group = wells)
      mat_well <- t(mat_well)
      
      # Prefix plate ID
      colnames(mat_well) <- paste0(plate, "_", colnames(mat_well))
      
      umi_list[[paste0(b, "_", plate)]] <- mat_well
    }
  }
  
  do.call(cbind, umi_list)
}


# -----------------------------------------------------------------------------
# Align metadata rows to UMI matrix column order
# -----------------------------------------------------------------------------
# NOTE: this OVERWRITES the `plate_id` column that process_plate_metadata()
# created as <plate_barcode>_<plate_well>_<well_index>. Inside the UMAP path
# `plate_id` therefore means <plate_barcode>_<plate_well>, matching the column
# names build_global_umi_matrix() produces. Do not join UMAP tables to the DE
# or DMSO tables on `plate_id`; the two definitions differ.
align_meta_to_umi <- function(UMI_all, meta) {
  meta$plate_id <- paste0(meta$plate_barcode, "_", meta$plate_well)
  meta[match(colnames(UMI_all), meta$plate_id), ]
}


# -----------------------------------------------------------------------------
# Normalise counts → logCPM
#
# Uses edgeR::cpm(log = TRUE, prior.count = 1) for consistency with the
# DE functions (both use log2 + prior count stabilisation).
# TMM normalisation is applied first.
# -----------------------------------------------------------------------------
normalize_logCPM <- function(counts) {
  d <- edgeR::DGEList(counts = counts)
  d <- edgeR::calcNormFactors(d, method = "TMM")
  edgeR::cpm(d, log = TRUE, prior.count = 1)
}


# -----------------------------------------------------------------------------
# Remove plate/batch effect from logCPM matrix before UMAP
#
# @param logCPM       genes × wells logCPM matrix
# @param meta         metadata data.frame aligned to logCPM columns
# @param CTRL         control treatment label (used to set the reference level)
# @param batch_var    column in meta to use as batch (default "plate_barcode")
# @param remove_batch logical — skip correction if FALSE
# -----------------------------------------------------------------------------
remove_umap_batch_effect <- function(
    logCPM,
    meta,
    CTRL         = "DMSO",
    batch_var    = "plate_barcode",
    remove_batch = TRUE
) {
  if (!remove_batch) return(logCPM)
  
  # Use 'treatment' column — consistent with the rest of the pipeline
  if (!"treatment" %in% colnames(meta)) {
    stop(
      "remove_umap_batch_effect: 'treatment' column not found in meta. ",
      "Available columns: ", paste(colnames(meta), collapse = ", ")
    )
  }
  
  treatment_levels <- unique(meta$treatment)
  if (!CTRL %in% treatment_levels) {
    warning(
      "remove_umap_batch_effect: control '", CTRL,
      "' not found in meta$treatment — batch correction applied without reference level"
    )
    groups <- factor(meta$treatment)
  } else {
    groups <- relevel(factor(meta$treatment), ref = CTRL)
  }
  
  batch_factor <- factor(meta[[batch_var]])
  design       <- model.matrix(~ 0 + groups)
  
  limma::removeBatchEffect(logCPM, batch = batch_factor, design = design)
}


# -----------------------------------------------------------------------------
# PCA on transposed logCPM (wells × genes)
# -----------------------------------------------------------------------------
run_umap_pca <- function(logCPM, scale = TRUE) {
  prcomp(t(logCPM), scale. = scale)
}


# -----------------------------------------------------------------------------
# Leiden clustering from PCA embedding
# -----------------------------------------------------------------------------
run_leiden_from_pca <- function(pca_matrix, n_pc = 20, k = 35, resolution = 1.0) {
  set.seed(123)
  
  coords <- pca_matrix[, seq_len(min(n_pc, ncol(pca_matrix))), drop = FALSE]
  n      <- nrow(coords)
  knn    <- RANN::nn2(coords, k = k)$nn.idx
  
  edges <- cbind(
    rep(seq_len(n), each = k),
    as.vector(knn)
  )
  
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  g <- igraph::simplify(g)   # remove duplicated A→B / B→A edges
  
  igraph::cluster_leiden(g, resolution_parameter = resolution)$membership
}


# -----------------------------------------------------------------------------
# UMAP layout from PCA embedding
# -----------------------------------------------------------------------------
run_umap_from_pca <- function(pca_matrix, n_pca = 20, n_neighbors = 20) {
  set.seed(123)
  
  n_pca <- min(n_pca, ncol(pca_matrix))
  
  um <- umap::umap(
    pca_matrix[, seq_len(n_pca), drop = FALSE],
    method      = "umap-learn",
    metric      = "euclidean",
    n_neighbors = n_neighbors
  )
  
  df          <- as.data.frame(um$layout)
  colnames(df) <- c("UMAP1", "UMAP2")
  df$plate_id  <- rownames(pca_matrix)
  
  df
}


# -----------------------------------------------------------------------------
# Annotate wells with spatial position category (center / middle / edge)
#
# Falls back to well_category = "unknown" for unrecognised plate formats
# instead of hard-stopping, so UMAP targets continue even if annotation fails.
# -----------------------------------------------------------------------------
annotate_plate_position <- function(df) {
  df <- df %>%
    dplyr::mutate(
      row = stringr::str_extract(plate_well, "^[A-Z]+"),
      col = as.integer(stringr::str_extract(plate_well, "[0-9]+"))
    )
  
  row_levels <- sort(unique(df$row))
  n_rows     <- length(row_levels)
  n_cols     <- max(df$col, na.rm = TRUE)
  
  plate_type <- dplyr::case_when(
    n_rows == 8  && n_cols <= 12 ~ "96",
    n_rows == 16 && n_cols <= 24 ~ "384",
    n_rows == 32 && n_cols <= 48 ~ "1536",
    TRUE                          ~ "unknown"
  )
  
  if (plate_type == "unknown") {
    warning(
      "annotate_plate_position: unrecognised plate format ",
      "(", n_rows, " rows × ", n_cols, " cols) — ",
      "well_category set to 'unknown' for all wells"
    )
    df$well_category <- "unknown"
    return(df)
  }
  
  df <- df %>% dplyr::mutate(row_index = match(row, row_levels))
  
  df <- switch(
    plate_type,
    
    "96" = df %>% dplyr::mutate(
      well_category = dplyr::case_when(
        row %in% c("C","D","E","F") & col %in% 4:9   ~ "center",
        row %in% c("B","G")         & col %in% 3:10  ~ "middle",
        row %in% c("C","D","E","F") & col %in% c(3,10) ~ "middle",
        TRUE ~ "edge"
      )
    ),
    
    "384" = df %>% dplyr::mutate(
      well_category = dplyr::case_when(
        row %in% c("G","H","I","J") & col %in% 9:16            ~ "center",
        row %in% c("D","E","F")     & col %in% 5:20            ~ "middle",
        row %in% c("G","H","I","J") & col %in% c(5:8, 17:20)   ~ "middle",
        row %in% c("K","L","M")     & col %in% 5:20            ~ "middle",
        TRUE ~ "edge"
      )
    ),
    
    "1536" = df %>% dplyr::mutate(
      well_category = dplyr::case_when(
        row_index %in% 13:20 & col %in% 17:32 ~ "center",
        row_index %in% 9:24  & col %in% 13:36 ~ "middle",
        TRUE ~ "edge"
      )
    )
  )
  
  df
}


# -----------------------------------------------------------------------------
# Assemble final UMAP data.frame (layout + metadata + Leiden clusters)
# -----------------------------------------------------------------------------
assemble_umap_table <- function(umap_layout, meta, leiden_clusters) {
  df                <- dplyr::left_join(umap_layout, meta, by = "plate_id")
  df$cluster_leiden <- factor(leiden_clusters)
  df
}


# -----------------------------------------------------------------------------
# Build a plot-ready table from a prcomp result (wells x PCs).
#
# Used to make PCA scatterplots that mirror the UMAP plots: same metadata
# join, same downstream plotly call. Variance-explained per PC is returned
# as an attribute so axis labels can show e.g. "PC1 (24.3%)".
# -----------------------------------------------------------------------------
assemble_pca_table <- function(pca_obj, meta, leiden_clusters = NULL) {
  scores <- as.data.frame(pca_obj$x[, 1:2, drop = FALSE])
  colnames(scores) <- c("PC1", "PC2")
  scores$plate_id  <- rownames(pca_obj$x)
  
  df <- dplyr::left_join(scores, meta, by = "plate_id")
  if (!is.null(leiden_clusters)) {
    df$cluster_leiden <- factor(leiden_clusters)
  }
  
  # Variance explained for the first two PCs (percent of total variance)
  total_var <- sum(pca_obj$sdev^2)
  var_pct   <- 100 * (pca_obj$sdev^2 / total_var)
  attr(df, "pc1_var_pct") <- round(var_pct[1], 1)
  attr(df, "pc2_var_pct") <- round(var_pct[2], 1)
  
  df
}


# -----------------------------------------------------------------------------
# Save a plotly widget as a self-contained HTML page with optional header
#
# Returns only the final user-facing HTML path (not the inner widget file)
# so targets format = "file" tracks a single stable output.
# -----------------------------------------------------------------------------
save_plotly_html <- function(plot, filename, out_dir, meta_experiment = NULL,
                             write_thumbnail = TRUE,
                             thumbnail_width  = 900,
                             thumbnail_height = 600) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Keep interactive widgets on the same font stack as the static plots.
  plot <- tryCatch(
    plotly::layout(plot, font = list(family = DRUGSEQ_PLOTLY_FONT)),
    error = function(e) plot
  )
  
  widget_file <- file.path(out_dir, paste0("widget_", filename))
  final_html  <- file.path(out_dir, filename)
  
  htmlwidgets::saveWidget(plot, file = widget_file, selfcontained = TRUE)
  
  header_html <- if (!is.null(meta_experiment)) build_simple_header(meta_experiment) else ""
  title_txt   <- tools::file_path_sans_ext(filename)
  
  writeLines(
    glue::glue('
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>{title_txt}</title>
</head>
<body style="font-family: \'Helvetica Neue\', Helvetica, Arial, sans-serif; margin:20px;">

{header_html}

<h4>{title_txt}</h4>

<iframe src="{basename(widget_file)}"
        style="width:100%; height:800px; border:none;">
</iframe>

</body>
</html>
'),
    final_html
  )
  
  # ---- PNG thumbnail (concern: Plotly thumbnails for the report) -----------
  # The final_report.R companion-image lookup looks for a PNG/SVG with the
  # same stem as the HTML, so the screenshot must land at <out_dir>/<stem>.png.
  # We screenshot the widget file (chart only) rather than final_html (chart
  # in an iframe with header), which produces tighter, prettier thumbnails.
  # Failure is non-fatal — the report will fall back to its generic HTML icon.
  if (isTRUE(write_thumbnail)) {
    png_file <- file.path(out_dir,
                          paste0(tools::file_path_sans_ext(filename), ".png"))
    tryCatch({
      webshot2::webshot(
        url     = widget_file,
        file    = png_file,
        vwidth  = thumbnail_width,
        vheight = thumbnail_height,
        delay   = 1.5    # plotly needs time to render the SVG before capture
      )
      message("  Plotly thumbnail saved: ", basename(png_file))
    }, error = function(e) {
      message("  webshot2 thumbnail failed (",
              conditionMessage(e), ") - skipping PNG for ", filename)
    })
  }
  
  # Return only the user-facing file - targets format = "file" tracks this one
  invisible(final_html)
}


# =============================================================================
# DE gene UMAP helpers
# =============================================================================

define_de_gene_universe <- function(de_df, lfc_cutoff, padj_cutoff) {
  unique(de_df$ENSEMBLID[abs(de_df$LFC) > lfc_cutoff & de_df$padj < padj_cutoff])
}


build_gene_contrast_matrix <- function(
    de_df,
    genes,
    value_col = "LFC",
    agg_fun   = mean
) {
  df <- de_df[de_df$ENSEMBLID %in% genes, ]
  
  df_collapsed <- df %>%
    dplyr::group_by(ENSEMBLID, Sample) %>%
    dplyr::summarise(
      value = agg_fun(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  mat           <- tidyr::pivot_wider(
    df_collapsed,
    id_cols     = ENSEMBLID,
    names_from  = Sample,
    values_from = value
  )
  mat           <- as.data.frame(mat)
  rownames(mat) <- mat$ENSEMBLID
  mat$ENSEMBLID <- NULL
  mat[is.na(mat)] <- 0
  
  mat
}


scale_gene_matrix <- function(mat) {
  t(scale(t(mat)))
}


format_threshold_label <- function(lfc, padj) {
  glue::glue(
    "lfc{gsub('\\.', '', format(lfc, trim = TRUE))}_",
    sprintf("%03d", padj * 100)
  )
}

# =============================================================================
# PLATE-HEATMAP DOSE HIGHLIGHTING (shared by QC + DE-count plate heatmaps)
# -----------------------------------------------------------------------------
# Clicking a well selects its compound; every well of that compound gets a
# frame whose colour encodes its concentration (colours assigned in ascending
# dose order), a small dose tag inside the well, and a clickable legend:
#   click a legend chip        -> hide / show that dose
#   shift-click a legend chip  -> show only that dose
# Each frame is drawn twice (white halo underneath, colour on top) so it stays
# visible on any cell colour of either the red QC or the green DE colourscale.
# Used by make_plate_heatmap() (02_input_functions.R) and
# make_all_de_count_heatmaps() (below). The JS expects these globals from the
# host page: compoundMatrix, doseLabelMatrix, doseKeyMatrix, rows, cols, and a
# #clearBtn button. Returned as plain strings (not glue templates), so the JS
# uses single braces; glue does not re-interpolate inserted values.
# =============================================================================

#' Numeric sort key for a dose, normalised to uM when the unit is recognised.
#' Unknown units (%, ug/mL, ...) keep their raw value, which still orders
#' correctly within a compound as long as that compound uses one unit.
plate_dose_sort_key <- function(dose, unit) {
  d <- suppressWarnings(as.numeric(as.character(dose)))
  u <- tolower(trimws(as.character(unit)))
  u <- gsub("\u00b5|\u03bc", "u", u)
  f <- c(pm = 1e-6, nm = 1e-3, um = 1, mm = 1e3, m = 1e6)[u]
  f[is.na(f)] <- 1
  unname(d * f)
}

#' Human-readable dose label ("10 uM"); "" when the dose is missing.
plate_dose_label <- function(dose, unit) {
  d   <- as.character(dose)
  u   <- as.character(unit)
  u[is.na(u)] <- ""
  lab <- trimws(paste(d, u))
  lab[is.na(d) | !nzchar(trimws(d))] <- ""
  lab
}

#' Legend container + label toggle, dropped into the page's #controls div.
plate_dose_controls_html <- function() {
  '<label style="font-size:13px;color:#374151;margin-left:10px;">
    <input type="checkbox" id="doseLabelsChk"> dose tags in wells
  </label>
  <div id="doseLegend" style="display:flex;flex-wrap:wrap;align-items:center;gap:6px;margin-top:8px;min-height:26px;font-size:12px;font-family:Helvetica,Arial,sans-serif;"></div>'
}

#' Click / legend logic. show_labels sets the initial state of the tag toggle
#' (off by default on 1536-well plates, where the wells are too small).
plate_dose_highlight_js <- function(div_id = "heatmap", show_labels = TRUE) {
  js <- r"---(
(function() {
  // Ordered, colour-blind-aware palette with no pure red/green, so frames do
  // not blend into either heatmap colourscale.
  const DOSE_PALETTE = ["#0072B2", "#E69F00", "#7B2CBF", "#00A6C8", "#D6007A",
                        "#111827", "#8C6D1F", "#4B6CFF", "#6B4226", "#7F7F7F"];
  const hm     = document.getElementById("__DIV_ID__");
  const legend = document.getElementById("doseLegend");
  const lblChk = document.getElementById("doseLabelsChk");
  if (lblChk) lblChk.checked = __SHOW_LABELS__;

  let selCompound = null;
  let doseOrder   = [];            // [{label, key, n, color}] ascending dose
  const hidden    = new Set();     // dose labels currently hidden

  function rect(r, c, color, width) {
    return { type: "rect", xref: "x", yref: "y",
             x0: c + 0.5, x1: c + 1.5, y0: r + 0.5, y1: r + 1.5,
             line: { color: color, width: width },
             fillcolor: "rgba(0,0,0,0)", layer: "above" };
  }

  function wellLabel(r, c) { return doseLabelMatrix[r][c] || "n/a"; }

  function buildDoseOrder(cmp) {
    const m = new Map();
    for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
      if (compoundMatrix[r][c] !== cmp) continue;
      const lab = wellLabel(r, c);
      const key = doseKeyMatrix[r][c];
      if (!m.has(lab)) m.set(lab, { label: lab, key: (key === null ? Infinity : key), n: 0 });
      m.get(lab).n++;
    }
    doseOrder = Array.from(m.values()).sort(function(a, b) {
      return (a.key - b.key) || a.label.localeCompare(b.label);
    });
    doseOrder.forEach(function(d, i) { d.color = DOSE_PALETTE[i % DOSE_PALETTE.length]; });
  }

  function renderLegend() {
    if (!legend) return;
    legend.innerHTML = "";
    if (!selCompound) return;
    const head = document.createElement("span");
    head.textContent = selCompound + " \u2014 low \u2192 high:";
    head.style.cssText = "font-weight:600;color:#374151;margin-right:4px;";
    legend.appendChild(head);
    doseOrder.forEach(function(d) {
      const off  = hidden.has(d.label);
      const chip = document.createElement("button");
      chip.type  = "button";
      chip.title = "Click: hide/show this dose \u00b7 Shift-click: show only this dose";
      chip.style.cssText =
        "display:inline-flex;align-items:center;gap:5px;padding:2px 8px;cursor:pointer;" +
        "border-radius:12px;background:#fff;font-size:12px;border:2px solid " + d.color + ";" +
        (off ? "opacity:.35;text-decoration:line-through;" : "");
      const sw = document.createElement("span");
      sw.style.cssText = "width:10px;height:10px;border-radius:2px;background:" + d.color + ";";
      chip.appendChild(sw);
      chip.appendChild(document.createTextNode(d.label + " (" + d.n + ")"));
      chip.addEventListener("click", function(ev) {
        if (ev.shiftKey) {
          hidden.clear();
          doseOrder.forEach(function(o) { if (o.label !== d.label) hidden.add(o.label); });
        } else if (hidden.has(d.label)) {
          hidden.delete(d.label);
        } else {
          hidden.add(d.label);
        }
        draw();
      });
      legend.appendChild(chip);
    });
    if (hidden.size) {
      const all = document.createElement("button");
      all.type = "button";
      all.textContent = "show all doses";
      all.style.cssText = "font-size:12px;padding:2px 8px;cursor:pointer;";
      all.addEventListener("click", function() { hidden.clear(); draw(); });
      legend.appendChild(all);
    }
  }

  function draw() {
    const shapes = [], ann = [];
    if (selCompound) {
      const byLab   = new Map(doseOrder.map(function(d) { return [d.label, d]; }));
      const showLbl = lblChk ? lblChk.checked : false;
      for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
        if (compoundMatrix[r][c] !== selCompound) continue;
        const lab = wellLabel(r, c);
        if (hidden.has(lab)) continue;
        const d = byLab.get(lab);
        shapes.push(rect(r, c, "#ffffff", 6));   // halo
        shapes.push(rect(r, c, d.color, 3));     // dose colour
        if (showLbl) ann.push({
          x: c + 1, y: r + 1, xref: "x", yref: "y", showarrow: false,
          text: lab.split(" ")[0], bgcolor: d.color, borderpad: 1,
          font: { size: 9, color: "#ffffff" }
        });
      }
    }
    Plotly.relayout(hm, { shapes: shapes, annotations: ann });
    renderLegend();
  }

  hm.on("plotly_click", function(ev) {
    const pt  = ev.points[0];
    const cmp = compoundMatrix[pt.y - 1][pt.x - 1];
    if (!cmp) return;
    if (cmp !== selCompound) { selCompound = cmp; hidden.clear(); buildDoseOrder(cmp); }
    draw();
  });
  document.getElementById("clearBtn").addEventListener("click", function() {
    selCompound = null; hidden.clear(); draw();
  });
  if (lblChk) lblChk.addEventListener("change", draw);
})();
)---"
  js <- gsub("__DIV_ID__", div_id, js, fixed = TRUE)
  gsub("__SHOW_LABELS__", if (isTRUE(show_labels)) "true" else "false", js, fixed = TRUE)
}

# =============================================================================
# DE-GENE-COUNT PLATE HEATMAP
# -----------------------------------------------------------------------------
# Mirrors the QC plate heatmap in 02_input_functions.R, but instead of colouring
# each well by its UMI count it colours each well by the number of DE genes that
# the well's compound x dose contrast produced (padj < padj_cutoff & |LFC| >
# lfc_cutoff). DE is defined per CONTRAST, not per well, so every well of a
# given compound x dose shows that contrast's count; DMSO wells have no contrast
# and are shown blank (grey). One interactive HTML + one PNG twin per plate,
# so the final report picks them up with the same companion-image lookup used
# for the QC heatmaps.
# =============================================================================

#' @param de_results  de_results_collapse (must contain Sample, treatment,
#'                     treatment_dose, timepoint, batch, LFC, padj)
#' @param plate_layout process_plate_metadata() output (well -> plate mapping)
#' @param meta_experiment  experiment metadata bundle (title etc.), optional
#' @param outdir       directory for plate_de_<barcode>.html/.png
#' @param plate_format 96 / 384 / 1536
#' @param lfc_cutoff,padj_cutoff  DE thresholds used to count genes per contrast
#' @param CTRL         control compound name (blanked on the plate)
make_all_de_count_heatmaps <- function(de_results,
                                       plate_layout,
                                       meta_experiment = NULL,
                                       outdir,
                                       plate_format = 384,
                                       lfc_cutoff   = 1,
                                       padj_cutoff  = 0.1,
                                       CTRL         = "DMSO") {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- 1. DE-gene count per contrast (Sample) ------------------------------
  # Key on (batch, Sample), not Sample alone: reference/QC compounds (e.g.
  # staurosporine, brefeldina) are plated on EVERY plate, so the same Sample
  # string (compound_dose_tp, no plate in it) legitimately appears once per
  # batch in de_results_collapse. Counting by Sample alone pools all 14 plates
  # into one well ("45k DE genes"); keying on batch gives each well its own
  # plate's count. distinct() also guards against any per-gene duplication.
  # Falls back to Sample-only if no batch column is present.
  has_batch <- "batch" %in% names(de_results) && "batch" %in% names(plate_layout)
  
  de_flt <- de_results %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff, abs(LFC) > lfc_cutoff)
  
  if (has_batch) {
    de_counts <- de_flt %>%
      dplyr::mutate(batch = as.character(batch)) %>%
      dplyr::distinct(batch, Sample, ENSEMBLID) %>%
      dplyr::count(batch, Sample, name = "n_de")
  } else {
    de_counts <- de_flt %>%
      dplyr::distinct(Sample, ENSEMBLID) %>%
      dplyr::count(Sample, name = "n_de")
  }
  
  # sample_clean is how a well maps to a contrast; rebuild it on the layout the
  # same way the DE functions do (treatment_dose + unit + timepoint), because
  # de_results carries the human-readable Sample but wells carry the parts.
  lay <- plate_layout %>%
    dplyr::mutate(
      tp_clean     = paste0(stringr::str_extract(as.character(timepoint), "\\d+"), "hr"),
      dose_clean   = paste0(treatment_dose, treatment_dose_unit),
      Sample       = paste0(treatment, "_", dose_clean, "_", tp_clean)
    )
  if (has_batch) lay <- dplyr::mutate(lay, batch = as.character(batch))
  lay <- lay %>%
    dplyr::left_join(de_counts, by = if (has_batch) c("batch", "Sample") else "Sample") %>%
    dplyr::mutate(
      n_de = dplyr::if_else(treatment == CTRL, NA_real_, as.numeric(n_de)),
      # a treated contrast that ran but called nothing is a real 0, not missing
      n_de = dplyr::if_else(treatment != CTRL & is.na(n_de), 0, n_de)
    )
  
  plate_dims <- switch(
    as.character(plate_format),
    "96"   = list(rows = LETTERS[1:8],  cols = 1:12),
    "384"  = list(rows = LETTERS[1:16], cols = 1:24),
    "1536" = list(rows = LETTERS[1:32], cols = 1:48),
    stop("Unsupported plate format: ", plate_format)
  )
  row_levels <- plate_dims$rows
  col_levels <- plate_dims$cols
  n_rows <- length(row_levels); n_cols <- length(col_levels)
  
  header_html <- if (!is.null(meta_experiment)) build_simple_header(meta_experiment) else ""
  exp_title   <- if (!is.null(meta_experiment)) meta_experiment$title else "DRUG-seq"
  
  out_files <- character(0)
  
  for (plate_name in unique(lay$plate_barcode)) {
    wi <- lay %>% dplyr::filter(plate_barcode == plate_name) %>%
      dplyr::mutate(
        row_id    = stringr::str_extract(plate_well, "^[A-Z]+"),
        column_id = as.numeric(stringr::str_extract(plate_well, "\\d+"))
      )
    
    de_matrix    <- matrix(NA_real_, n_rows, n_cols, dimnames = list(row_levels, col_levels))
    cmp_matrix   <- matrix("",       n_rows, n_cols, dimnames = list(row_levels, col_levels))
    hover_matrix <- matrix("",       n_rows, n_cols, dimnames = list(row_levels, col_levels))
    dose_lab_matrix <- matrix("",       n_rows, n_cols, dimnames = list(row_levels, col_levels))
    dose_key_matrix <- matrix(NA_real_, n_rows, n_cols, dimnames = list(row_levels, col_levels))
    wi_dose_lab <- plate_dose_label(wi$treatment_dose, wi$treatment_dose_unit)
    wi_dose_key <- plate_dose_sort_key(wi$treatment_dose, wi$treatment_dose_unit)
    
    for (i in seq_len(nrow(wi))) {
      r <- match(wi$row_id[i], row_levels)
      c <- match(wi$column_id[i], col_levels)
      if (is.na(r) || is.na(c)) next
      de_matrix[r, c]  <- wi$n_de[i]
      cmp_matrix[r, c] <- wi$treatment[i]
      dose_lab_matrix[r, c] <- wi_dose_lab[i]
      dose_key_matrix[r, c] <- wi_dose_key[i]
      hover_matrix[r, c] <- paste0(
        "Well: ", wi$plate_well[i], "<br>",
        "Compound: ", wi$treatment[i], "<br>",
        "Conc: ", wi$treatment_dose[i], " ", wi$treatment_dose_unit[i], "<br>",
        if (identical(wi$treatment[i], CTRL))
          "Control well (no contrast)"
        else
          paste0("DE genes: ", ifelse(is.na(wi$n_de[i]), "n/a", wi$n_de[i]))
      )
    }
    
    cmp_json    <- jsonlite::toJSON(cmp_matrix,    auto_unbox = TRUE)
    de_json     <- jsonlite::toJSON(de_matrix,     auto_unbox = TRUE, na = "null")
    hover_json  <- jsonlite::toJSON(hover_matrix,  auto_unbox = TRUE)
    dose_lab_json <- jsonlite::toJSON(dose_lab_matrix, auto_unbox = TRUE)
    dose_key_json <- jsonlite::toJSON(dose_key_matrix, auto_unbox = TRUE, na = "null", digits = NA)
    dose_controls <- plate_dose_controls_html()
    dose_js       <- plate_dose_highlight_js("heatmap", show_labels = n_rows <= 16)
    
    html <- glue::glue('
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>DE-gene count \u2014 Plate {plate_name}</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
</head>
<body>
<div style="margin:20px;">
  {header_html}
  <div id="controls" style="margin-bottom:10px;">
    <button id="clearBtn">Clear selection</button>
    <span style="color:#6b7280;font-size:13px;margin-left:8px;">
      Colour = DE genes for that well\'s compound \u00d7 dose (padj &lt; {padj_cutoff}, |log2FC| &gt; {lfc_cutoff}). Grey = control / no contrast.
    </span>
    {dose_controls}
  </div>
  <div id="heatmap" style="width:100%;height:80vh;min-height:600px;"></div>
</div>
<script>
const compoundMatrix = {cmp_json};
const deMatrix       = {de_json};
const hoverMatrix    = {hover_json};
const doseLabelMatrix = {dose_lab_json};
const doseKeyMatrix   = {dose_key_json};
const rows = compoundMatrix.length;
const cols = compoundMatrix[0].length;
const xVals = Array.from({{length: cols}}, (_, i) => i + 1);
const yVals = Array.from({{length: rows}}, (_, i) => i + 1);
const yTickText = Array.from({{length: rows}}, (_,i) => String.fromCharCode(65+i));
const baseTrace = {{
  z: deMatrix, x: xVals, y: yVals, text: hoverMatrix,
  customdata: compoundMatrix,
  hovertemplate: "%{{text}}<extra></extra>",
  type: "heatmap",
  colorscale: [
   [0.00, "rgb(237,248,233)"],
   [0.33, "rgb(161,217,155)"],
   [0.66, "rgb(65,171,93)"],
   [1.00, "rgb(0,90,50)"]
  ],
  zauto: true, showscale: true, colorbar: {{ title: "DE genes" }}
}};
const layout = {{
  title: "{exp_title} \u2013 DE-gene count \u2013 Plate {plate_name}",
  font: {{ family: "Helvetica Neue, Helvetica, Arial, sans-serif" }},
  xaxis: {{ tickmode:"array", tickvals:xVals, ticktext:xVals.map(String), side:"top" }},
  yaxis: {{ tickmode:"array", tickvals:yVals, ticktext:yTickText, autorange:"reversed" }},
  shapes: []
}};
Plotly.newPlot("heatmap", [baseTrace], layout);
{dose_js}
</script>
</body>
</html>
')
    
    html_path <- file.path(outdir, glue::glue("plate_de_{plate_name}.html"))
    writeLines(html, html_path)
    out_files <- c(out_files, html_path)
    
    png_path <- file.path(outdir, glue::glue("plate_de_{plate_name}.png"))
    tryCatch(
      render_plate_thumbnail(
        counts_matrix = de_matrix,
        row_levels    = row_levels,
        col_levels    = col_levels,
        title         = glue::glue("{exp_title} \u2013 DE genes \u2013 Plate {plate_name}"),
        out_png       = png_path,
        palette       = c("#edf8e9", "#a1d99b", "#41ab5d", "#005a32"),
        legend_name   = "DE genes"
      ),
      error = function(e)
        message("  DE-count thumbnail failed (", conditionMessage(e),
                ") - skipping PNG for ", plate_name)
    )
  }
  
  invisible(out_files)
}