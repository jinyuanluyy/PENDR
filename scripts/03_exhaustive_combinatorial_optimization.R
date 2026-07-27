#!/usr/bin/env Rscript

source("config.R")
source("R/pipeline_utils.R")
source("R/simms_compat.R")

required_packages <- c("AnnotationDbi", "org.Hs.eg.db", "pROC")

map_submodule_genes <- function(submodule_names, adjacency_matrices) {
  gene_lists <- lapply(submodule_names, function(submodule_name) {
    adjacency <- adjacency_matrices[[submodule_name]]
    if (is.null(adjacency)) {
      return(character())
    }

    entrez_ids <- rownames(adjacency)
    symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = entrez_ids,
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )
    unique(stats::na.omit(unname(symbols)))
  })
  names(gene_lists) <- submodule_names
  gene_lists[lengths(gene_lists) > 0L]
}

select_candidate_submodules <- function(meta_analysis_results) {
  significant <- meta_analysis_results[
    meta_analysis_results$adjusted_q_value <
      pendr_config$submodule_fdr_threshold,
    ,
    drop = FALSE
  ]
  head(
    significant$pathway_guided_submodule,
    pendr_config$max_candidate_submodules
  )
}

read_discrimination_data <- function() {
  expression <- as.matrix(utils::read.delim(
    pendr_config$optimization_expression_file,
    row.names = 1L,
    check.names = FALSE
  ))
  storage.mode(expression) <- "numeric"
  metadata <- utils::read.delim(
    pendr_config$optimization_sample_metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  metadata <- metadata[
    match(colnames(expression), metadata$sample_id),
    ,
    drop = FALSE
  ]
  if (anyNA(metadata$sample_id)) {
    stop(
      "Optimization expression samples are missing from metadata.",
      call. = FALSE
    )
  }

  classes <- ifelse(
    metadata$class == pendr_config$tumor_class_label,
    2L,
    ifelse(
      metadata$class == pendr_config$normal_class_label,
      1L,
      NA_integer_
    )
  )
  if (anyNA(classes) || length(unique(classes)) != 2L) {
    stop(
      "Optimization classes must be exactly Tumor and Normal.",
      call. = FALSE
    )
  }
  names(classes) <- metadata$sample_id
  list(expression = expression, classes = classes)
}

safe_ratio <- function(numerator, denominator) {
  if (denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

evaluate_submodule_combination <- function(
  submodule_names,
  gene_lists,
  discrimination_data
) {
  genes <- unique(unlist(gene_lists[submodule_names], use.names = FALSE))
  selected_expression <- discrimination_data$expression[
    rownames(discrimination_data$expression) %in% genes,
    ,
    drop = FALSE
  ]
  if (nrow(selected_expression) < 2L) {
    return(NULL)
  }

  principal_components <- stats::prcomp(
    t(selected_expression),
    center = TRUE,
    scale. = FALSE
  )
  prediction <- principal_components$x[, 1]
  labels <- discrimination_data$classes[names(prediction)]
  roc_auc <- as.numeric(pROC::auc(
    pROC::roc(labels, prediction, quiet = TRUE)
  ))

  predicted_classes <- ifelse(prediction > 0, 2L, 1L)
  true_positives <- sum(predicted_classes == 2L & labels == 2L)
  false_positives <- sum(predicted_classes == 2L & labels == 1L)
  false_negatives <- sum(predicted_classes == 1L & labels == 2L)
  precision <- safe_ratio(
    true_positives,
    true_positives + false_positives
  )
  recall <- safe_ratio(
    true_positives,
    true_positives + false_negatives
  )
  f1_score <- if (
    is.finite(precision) &&
    is.finite(recall) &&
    precision + recall > 0
  ) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }

  data.frame(
    submodule_combination = paste(submodule_names, collapse = " | "),
    submodule_count = length(submodule_names),
    gene_count = nrow(selected_expression),
    AUROC = roc_auc,
    Accuracy = mean(predicted_classes == labels),
    F1_score = f1_score,
    stringsAsFactors = FALSE
  )
}

evaluate_all_combinations <- function(
  candidate_names,
  gene_lists,
  discrimination_data
) {
  results <- list()
  result_index <- 1L

  for (combination_size in seq_along(candidate_names)) {
    combinations <- utils::combn(
      candidate_names,
      combination_size,
      simplify = FALSE
    )
    for (submodule_names in combinations) {
      result <- evaluate_submodule_combination(
        submodule_names,
        gene_lists,
        discrimination_data
      )
      if (!is.null(result)) {
        results[[result_index]] <- result
        result_index <- result_index + 1L
      }
    }
  }

  if (length(results) == 0L) {
    stop("No submodule combination could be evaluated.", call. = FALSE)
  }
  do.call(rbind, results)
}

calculate_composite_scores <- function(results) {
  metric_columns <- c("AUROC", "Accuracy", "F1_score")
  results <- results[
    stats::complete.cases(results[metric_columns]),
    ,
    drop = FALSE
  ]
  if (nrow(results) == 0L) {
    stop("All optimization metrics contain missing values.", call. = FALSE)
  }

  metric_ranks <- apply(
    results[metric_columns],
    2L,
    function(values) rank(-values, ties.method = "min")
  )
  results$composite_score <- apply(
    metric_ranks,
    1L,
    prod
  )^(1 / length(metric_columns))
  results <- results[
    order(results$composite_score),
    ,
    drop = FALSE
  ]
  results$optimization_rank <- seq_len(nrow(results))
  results
}

write_optimization_outputs <- function(results, gene_lists) {
  output_directory <- file.path(
    "results",
    "exhaustive_combinatorial_optimization"
  )
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  write_tsv(
    results,
    file.path(output_directory, "submodule_combination_metrics.tsv")
  )
  saveRDS(
    results,
    file.path(output_directory, "submodule_combination_metrics.rds")
  )

  optimized_submodules <- strsplit(
    results$submodule_combination[[1]],
    " | ",
    fixed = TRUE
  )[[1]]
  optimized_gene_lists <- gene_lists[optimized_submodules]
  optimized_genes <- sort(unique(unlist(
    optimized_gene_lists,
    use.names = FALSE
  )))
  writeLines(
    optimized_submodules,
    file.path(output_directory, "optimized_pdac_driver_submodules.txt")
  )
  writeLines(
    optimized_genes,
    file.path(output_directory, "optimized_pdac_driver_genes.txt")
  )
}

main <- function() {
  assert_packages_installed(required_packages)
  meta_analysis_file <- file.path(
    "results",
    "rank_product_meta_analysis",
    "rank_product_meta_analysis_results.rds"
  )
  assert_files_exist(
    c(
      pendr_config$pathway_guided_submodules_file,
      meta_analysis_file,
      pendr_config$optimization_expression_file,
      pendr_config$optimization_sample_metadata_file
    )
  )

  meta_analysis_results <- readRDS(meta_analysis_file)
  candidate_names <- select_candidate_submodules(
    meta_analysis_results
  )
  if (length(candidate_names) == 0L) {
    stop(
      "No submodule passed the configured Q-value threshold.",
      call. = FALSE
    )
  }

  adjacency_matrices <- get.adjacency.matrix(
    pendr_config$pathway_guided_submodules_file
  )
  adjacency_matrices <- stats::setNames(
    adjacency_matrices,
    normalize_subnetwork_name(names(adjacency_matrices))
  )
  gene_lists <- map_submodule_genes(
    candidate_names,
    adjacency_matrices
  )
  candidate_names <- intersect(candidate_names, names(gene_lists))
  if (length(candidate_names) == 0L) {
    stop(
      "Candidate submodules could not be mapped to gene symbols.",
      call. = FALSE
    )
  }

  discrimination_data <- read_discrimination_data()
  optimization_results <- evaluate_all_combinations(
    candidate_names,
    gene_lists,
    discrimination_data
  )
  optimization_results <- calculate_composite_scores(
    optimization_results
  )
  write_optimization_outputs(optimization_results, gene_lists)
  message(
    "Saved exhaustive combinatorial optimization results under results/."
  )
}

main()
