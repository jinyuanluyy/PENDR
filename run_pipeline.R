#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
stage <- if (length(arguments) == 0L) "help" else tolower(arguments[[1]])

run_script <- function(path, script_arguments = character()) {
  status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = c(path, script_arguments)
  )
  if (!identical(status, 0L)) {
    stop("Pipeline stage failed: ", path, call. = FALSE)
  }
}

score_all_cohorts <- function() {
  source("config.R")
  for (cohort_id in pendr_config$cohorts$id) {
    run_script(
      "scripts/01_survival_associated_edge_scores.R",
      cohort_id
    )
  }
}

if (stage == "validate") {
  run_script("scripts/validate_inputs.R")
} else if (stage == "score") {
  if (length(arguments) != 2L) {
    stop(
      paste(
        "Usage: Rscript run_pipeline.R score",
        "<CPTAC-PDAC|TCGA-PDAC|ICGC-PDAC-CA>"
      ),
      call. = FALSE
    )
  }
  run_script(
    "scripts/01_survival_associated_edge_scores.R",
    arguments[[2]]
  )
} else if (stage == "meta-analysis") {
  run_script("scripts/02_rank_product_meta_analysis.R")
} else if (stage == "optimize") {
  run_script("scripts/03_exhaustive_combinatorial_optimization.R")
} else if (stage == "all") {
  run_script("scripts/validate_inputs.R")
  score_all_cohorts()
  run_script("scripts/02_rank_product_meta_analysis.R")
  run_script("scripts/03_exhaustive_combinatorial_optimization.R")
} else if (stage == "help") {
  cat(
    paste(
      "PENDR edge-perturbation analysis",
      "",
      "Commands:",
      "  Rscript run_pipeline.R validate",
      "  Rscript run_pipeline.R score CPTAC-PDAC",
      "  Rscript run_pipeline.R score TCGA-PDAC",
      "  Rscript run_pipeline.R score ICGC-PDAC-CA",
      "  Rscript run_pipeline.R meta-analysis",
      "  Rscript run_pipeline.R optimize",
      "  Rscript run_pipeline.R all",
      "",
      sep = "\n"
    )
  )
} else {
  stop("Unknown pipeline stage: ", stage, call. = FALSE)
}
