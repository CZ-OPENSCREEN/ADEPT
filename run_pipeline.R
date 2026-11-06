args <- commandArgs(trailingOnly = TRUE)

# Default config location
config_path <- Sys.getenv("DRUGSEQ_CONFIG", "/data/config.yml")
Sys.setenv(DRUGSEQ_CONFIG = config_path)

if ("--configure" %in% args) {
  source("/app/config_builder.R")
  quit(save = "no")
}

library(targets)

# The targets store lives where config$paths$cache_dir says (the configurator
# writes /data/_targets). Read it here so _targets.R does not need to call
# tar_config_set(), which would write a stray _targets.yaml into /app.
cfg   <- yaml::read_yaml(config_path)
store <- cfg$paths$cache_dir
if (is.null(store) || !nzchar(store)) store <- "/data/_targets"

tar_make(
  script = "/app/_targets.R",
  store  = store
)