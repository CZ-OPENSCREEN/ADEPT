read_config <- function(path = "config.yml") {
  if (!file.exists(path)) {
    stop("Config file not found: ", path)
  }
  
  config <- yaml::read_yaml(path)
  
  # Top-level blocks the pipeline actually reads. `qc`, `selection` and `umap`
  # are optional (every consumer tolerates NULL). An `outputs:` block, present
  # in older configs, is ignored.
  required_fields <- c("experiment", "analysis", "input_files", "filters")
  
  if (is.null(config$paths$data_dir))
    stop("config$paths$data_dir is missing")
  
  if (is.null(config$paths$output_dir))
    stop("config$paths$output_dir is missing")
  
  if (is.null(config$paths$cache_dir))
    stop("config$paths$cache_dir is missing")
  
  missing <- setdiff(required_fields, names(config))
  
  if (length(missing) > 0) {
    stop(
      "Config file is missing required fields: ",
      paste(missing, collapse = ", ")
    )
  }
  
  # Defaults for fields older configs may lack, and a check on their values.
  # top_n: "all" keeps every non-flagged DMSO well (with a floor of min_keep
  # per plate); a number keeps that many best-ranked wells per plate.
  config$analysis$dmso$top_n    <- config$analysis$dmso$top_n %||% "all"
  config$analysis$dmso$min_keep <- config$analysis$dmso$min_keep %||% 6
  config$analysis$de$method     <- config$analysis$de$method %||% "auto"

  top_n <- config$analysis$dmso$top_n
  if (is.character(top_n)) top_n <- tolower(top_n)
  if (!(identical(top_n, "all") || (is.numeric(top_n) && length(top_n) == 1 && top_n >= 1)))
    stop("analysis$dmso$top_n must be \"all\" or a positive number, not: ", format(top_n))
  config$analysis$dmso$top_n <- top_n

  if (!(is.numeric(config$analysis$dmso$min_keep) && config$analysis$dmso$min_keep >= 1))
    stop("analysis$dmso$min_keep must be a positive number")
  if (!config$analysis$de$method %in% c("auto", "trend", "voom"))
    stop("analysis$de$method must be one of auto, trend, voom; not: ", config$analysis$de$method)

  message("Loaded config from: ", path)
  config
}


load_experiment_meta <- function(config) {
  list(
    title = config$experiment$title,
    description = config$experiment$description,
    cell_line = config$experiment$cell_line,
    citation = config$experiment$citation
  )
}


create_output_dirs <- function(dirs) {
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  message("Ensured output directories exist.")
  dirs
}


build_simple_header <- function(meta_experiment) {
  glue::glue("
<div style='margin-bottom:20px; font-family: \"Helvetica Neue\", Helvetica, Arial, sans-serif;'>

<h3>{meta_experiment$title}</h3>

<p>
<b>Cell line:</b> {meta_experiment$cell_line}
</p>

<p style='max-width:900px;'>
{meta_experiment$description}
</p>

<p style='font-size:0.9em; color:gray;'>
{meta_experiment$citation}
</p>

<hr>

</div>
")
}