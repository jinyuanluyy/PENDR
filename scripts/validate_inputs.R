#!/usr/bin/env Rscript

source("config.R")
source("R/pipeline_utils.R")
source("R/simms_compat.R")

validate_pathway_guided_submodules <- function() {
  submodules <- get.adjacency.matrix(
    pendr_config$pathway_guided_submodules_file
  )
  empty_submodules <- names(submodules)[
    vapply(submodules, nrow, integer(1)) < 2L
  ]
  if (length(empty_submodules) > 0L) {
    stop(
      "Pathway-guided submodules with fewer than two nodes: ",
      paste(empty_submodules, collapse = ", "),
      call. = FALSE
    )
  }

  expected <- pendr_config$expected_pathway_guided_submodules
  if (length(submodules) != expected) {
    warning(
      "The manuscript reports ",
      expected,
      " pathway-guided submodules, but the input contains ",
      length(submodules),
      "."
    )
  }
  submodules
}

validate_cohort_manifest <- function() {
  manifest_path <- file.path(
    pendr_config$cohort_data_directory,
    pendr_config$datasets_file
  )
  manifest <- utils::read.delim(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c(
    "dataset",
    "annotation",
    "survstat",
    "survtime",
    "survtime.unit",
    pendr_config$data_type
  )
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns) > 0L) {
    stop(
      "Cohort manifest is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  missing_cohorts <- setdiff(
    pendr_config$cohorts$dataset_name,
    manifest$dataset
  )
  if (length(missing_cohorts) > 0L) {
    stop(
      "Cohort manifest is missing dataset entries: ",
      paste(missing_cohorts, collapse = ", "),
      call. = FALSE
    )
  }

  for (cohort_id in pendr_config$cohorts$dataset_name) {
    cohort <- manifest[manifest$dataset == cohort_id, , drop = FALSE]
    cohort_directory <- file.path(
      pendr_config$cohort_data_directory,
      cohort_id
    )
    assert_files_exist(
      c(
        file.path(cohort_directory, cohort$annotation),
        file.path(cohort_directory, cohort[[pendr_config$data_type]])
      ),
      paste0(cohort_id, " input file")
    )
  }
  invisible(manifest)
}

validate_optimization_data <- function() {
  expression <- utils::read.delim(
    pendr_config$optimization_expression_file,
    row.names = 1L,
    check.names = FALSE
  )
  samples <- utils::read.delim(
    pendr_config$optimization_sample_metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c("sample_id", "class")
  if (!all(required_columns %in% names(samples))) {
    stop(
      "Optimization sample metadata must contain sample_id and class.",
      call. = FALSE
    )
  }
  if (anyDuplicated(samples$sample_id)) {
    stop("Optimization sample IDs must be unique.", call. = FALSE)
  }
  if (!all(colnames(expression) %in% samples$sample_id)) {
    stop(
      "Every optimization expression column must have sample metadata.",
      call. = FALSE
    )
  }

  expected_classes <- c(
    pendr_config$tumor_class_label,
    pendr_config$normal_class_label
  )
  observed_classes <- unique(samples$class)
  if (!all(expected_classes %in% observed_classes)) {
    stop(
      "Optimization metadata must include Tumor and Normal samples.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

main <- function() {
  assert_files_exist(
    c(
      pendr_config$pathway_guided_submodules_file,
      file.path(
        pendr_config$cohort_data_directory,
        pendr_config$datasets_file
      ),
      pendr_config$optimization_expression_file,
      pendr_config$optimization_sample_metadata_file
    )
  )
  validate_pathway_guided_submodules()
  validate_cohort_manifest()
  validate_optimization_data()
  message("All PENDR analysis inputs passed validation.")
}

main()
