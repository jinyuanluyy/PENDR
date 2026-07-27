#!/usr/bin/env Rscript

source("config.R")
source("R/pipeline_utils.R")

read_submodule_score_matrix <- function(score_files, cohort_labels) {
  cohort_scores <- lapply(score_files, readRDS)
  reference_submodules <-
    cohort_scores[[1]]$pathway_guided_submodule

  aligned_scores <- lapply(seq_along(cohort_scores), function(index) {
    scores <- cohort_scores[[index]]
    if (anyDuplicated(scores$pathway_guided_submodule)) {
      stop(
        "Duplicated pathway-guided submodule names in ",
        score_files[[index]],
        call. = FALSE
      )
    }

    scores <- scores[
      match(
        reference_submodules,
        scores$pathway_guided_submodule
      ),
      ,
      drop = FALSE
    ]
    if (anyNA(scores$pathway_guided_submodule)) {
      stop(
        "Cohort score tables contain different submodules.",
        call. = FALSE
      )
    }
    as.numeric(scores$edge_based_submodule_score)
  })

  score_matrix <- do.call(cbind, aligned_scores)
  rownames(score_matrix) <- reference_submodules
  colnames(score_matrix) <- cohort_labels
  if (anyNA(score_matrix) || any(!is.finite(score_matrix))) {
    stop(
      "The submodule-by-cohort score matrix contains non-finite values.",
      call. = FALSE
    )
  }
  score_matrix
}

calculate_rank_product <- function(score_matrix) {
  cohort_ranks <- apply(
    score_matrix,
    2L,
    function(scores) rank(-scores, ties.method = "first")
  )
  apply(cohort_ranks, 1L, prod)
}

estimate_empirical_significance <- function(
  score_matrix,
  observed_rank_product,
  iterations
) {
  extreme_counts <- integer(nrow(score_matrix))

  for (iteration in seq_len(iterations)) {
    permuted_scores <- apply(score_matrix, 2L, sample)
    permuted_rank_product <- calculate_rank_product(permuted_scores)
    extreme_counts <- extreme_counts +
      as.integer(permuted_rank_product <= observed_rank_product)
  }
  extreme_counts / iterations
}

main <- function() {
  set.seed(pendr_config$random_seed)
  score_files <- file.path(
    "results",
    "cohort_submodule_scores",
    paste0(
      pendr_config$cohorts$id,
      "_edge_based_submodule_scores.rds"
    )
  )
  assert_files_exist(score_files, "cohort submodule score file")

  score_matrix <- read_submodule_score_matrix(
    score_files,
    pendr_config$cohorts$label
  )
  observed_rank_product <- calculate_rank_product(score_matrix)
  empirical_significance <- estimate_empirical_significance(
    score_matrix,
    observed_rank_product,
    iterations = pendr_config$rank_product_permutations
  )

  meta_analysis_results <- data.frame(
    pathway_guided_submodule = rownames(score_matrix),
    observed_rank_product = observed_rank_product,
    empirical_significance = empirical_significance,
    adjusted_q_value = stats::p.adjust(
      empirical_significance,
      method = "BH"
    ),
    stringsAsFactors = FALSE
  )
  meta_analysis_results <- meta_analysis_results[
    order(
      meta_analysis_results$adjusted_q_value,
      meta_analysis_results$observed_rank_product
    ),
    ,
    drop = FALSE
  ]
  meta_analysis_results$consensus_rank <-
    seq_len(nrow(meta_analysis_results))

  score_table <- data.frame(
    pathway_guided_submodule = rownames(score_matrix),
    score_matrix,
    check.names = FALSE
  )
  output_directory <- file.path("results", "rank_product_meta_analysis")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    score_table,
    file.path(output_directory, "submodule_by_cohort_score_matrix.tsv")
  )
  write_tsv(
    meta_analysis_results,
    file.path(output_directory, "rank_product_meta_analysis_results.tsv")
  )
  saveRDS(
    meta_analysis_results,
    file.path(output_directory, "rank_product_meta_analysis_results.rds")
  )

  message("Saved rank product meta-analysis results: ", output_directory)
}

main()
