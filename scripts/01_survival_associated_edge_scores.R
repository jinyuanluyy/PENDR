#!/usr/bin/env Rscript

source("config.R")
source("R/pipeline_utils.R")
source("R/simms_compat.R")

resolve_cohort <- function(cohort_id) {
  cohort <- pendr_config$cohorts[
    pendr_config$cohorts$id == cohort_id,
    ,
    drop = FALSE
  ]
  if (nrow(cohort) != 1L) {
    stop("Unknown cohort ID: ", cohort_id, call. = FALSE)
  }
  cohort
}

extract_unique_edges <- function(adjacency_matrices) {
  edge_tables <- lapply(adjacency_matrices, function(adjacency) {
    edge_indices <- which(
      upper.tri(adjacency) & adjacency == 1,
      arr.ind = TRUE
    )
    if (nrow(edge_indices) == 0L) {
      return(NULL)
    }
    data.frame(
      GeneID1 = rownames(adjacency)[edge_indices[, 1]],
      GeneID2 = colnames(adjacency)[edge_indices[, 2]],
      stringsAsFactors = FALSE
    )
  })
  unique(do.call(rbind, Filter(Negate(is.null), edge_tables)))
}

fit_edge_models <- function(gene_pairs, cancer_data) {
  genes <- unique(c(gene_pairs$GeneID1, gene_pairs$GeneID2))
  coefficients <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(genes),
    dimnames = list(genes, genes)
  )
  p_values <- coefficients

  for (edge_index in seq_len(nrow(gene_pairs))) {
    gene_1 <- gene_pairs$GeneID1[[edge_index]]
    gene_2 <- gene_pairs$GeneID2[[edge_index]]
    feature_1 <- paste0(gene_1, pendr_config$expression_feature_suffix)
    feature_2 <- paste0(gene_2, pendr_config$expression_feature_suffix)
    fit <- fit.interaction.model(
      feature1 = feature_1,
      feature2 = feature_2,
      expression.data = cancer_data[["all.data"]][[pendr_config$data_type]],
      survival.data = cancer_data$all.survobj,
      data.type.ordinal = FALSE,
      centre.data = "median"
    )

    interaction <- fit[["cox.int"]]
    if (is.null(interaction) || length(interaction) < 2L) {
      next
    }

    hazard_ratio <- suppressWarnings(as.numeric(interaction[[1]]))
    p_value <- suppressWarnings(as.numeric(interaction[[2]]))
    if (
      !is.finite(hazard_ratio) ||
      hazard_ratio <= 0 ||
      !is.finite(p_value)
    ) {
      next
    }

    coefficient <- log2(hazard_ratio)
    coefficients[gene_1, gene_2] <- coefficient
    coefficients[gene_2, gene_1] <- coefficient
    p_values[gene_1, gene_2] <- p_value
    p_values[gene_2, gene_1] <- p_value
  }

  list(coefficients = coefficients, p_values = p_values)
}

calculate_edge_based_submodule_score <- function(
  adjacency,
  edge_models
) {
  nodes <- intersect(
    rownames(adjacency),
    rownames(edge_models$coefficients)
  )
  if (length(nodes) < 2L) {
    return(0)
  }

  adjacency <- adjacency[nodes, nodes, drop = FALSE]
  edge_indices <- which(
    upper.tri(adjacency) & adjacency == 1,
    arr.ind = TRUE
  )
  if (nrow(edge_indices) == 0L) {
    return(0)
  }

  edge_contributions <- vapply(seq_len(nrow(edge_indices)), function(index) {
    gene_1 <- nodes[[edge_indices[index, 1]]]
    gene_2 <- nodes[[edge_indices[index, 2]]]
    p_value <- edge_models$p_values[gene_1, gene_2]
    coefficient <- edge_models$coefficients[gene_1, gene_2]

    if (
      is.finite(p_value) &&
      p_value < pendr_config$interaction_p_value_threshold &&
      is.finite(coefficient)
    ) {
      abs(coefficient)
    } else {
      0
    }
  }, numeric(1))

  sum(edge_contributions) / length(nodes)
}

main <- function() {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 1L) {
    stop(
      paste(
        "Usage: Rscript scripts/01_survival_associated_edge_scores.R",
        "<CPTAC-PDAC|TCGA-PDAC|ICGC-PDAC-CA>"
      ),
      call. = FALSE
    )
  }

  cohort_id <- arguments[[1]]
  cohort <- resolve_cohort(cohort_id)
  assert_files_exist(pendr_config$pathway_guided_submodules_file)

  adjacency_matrices <- get.adjacency.matrix(
    pendr_config$pathway_guided_submodules_file
  )
  gene_pairs <- extract_unique_edges(adjacency_matrices)
  cancer_data <- load.cancer.datasets(
    tumour.only = TRUE,
    with.survival.only = TRUE,
    truncate.survival = pendr_config$truncate_survival_years,
    datasets.to.load = cohort$dataset_name,
    data.types = pendr_config$data_type,
    datasets.file = pendr_config$datasets_file,
    data.directory = pendr_config$cohort_data_directory,
    verbose = FALSE,
    subset = NULL
  )
  edge_models <- fit_edge_models(gene_pairs, cancer_data)
  submodule_scores <- vapply(
    adjacency_matrices,
    calculate_edge_based_submodule_score,
    numeric(1),
    edge_models = edge_models
  )

  score_table <- data.frame(
    pathway_guided_submodule = normalize_subnetwork_name(
      names(submodule_scores)
    ),
    edge_based_submodule_score = unname(submodule_scores),
    stringsAsFactors = FALSE
  )
  output_file <- file.path(
    "results",
    "cohort_submodule_scores",
    paste0(cohort_id, "_edge_based_submodule_scores.rds")
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(score_table, output_file)
  message("Saved survival-associated edge scores: ", output_file)
}

main()
