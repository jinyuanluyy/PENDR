#!/usr/bin/env Rscript

source("R/simms_compat.R")

temporary_directory <- tempfile("pendr_test_")
dir.create(temporary_directory, recursive = TRUE)
on.exit(unlink(temporary_directory, recursive = TRUE), add = TRUE)

# Test adjacency parsing.
network_file <- file.path(temporary_directory, "subnetworks.txt")
writeLines(
  c(
    "#Example network",
    "1\tgeneA\tgeneB",
    "2\tgeneB\tgeneC"
  ),
  network_file
)
adjacency <- get.adjacency.matrix(network_file)
stopifnot(
  length(adjacency) == 1L,
  adjacency[[1]]["geneA", "geneB"] == 1,
  adjacency[[1]]["geneB", "geneC"] == 1,
  adjacency[[1]]["geneA", "geneC"] == 0
)

# Test the XNOR interaction model with one in-memory cohort.
set.seed(1)
sample_count <- 80L
gene1_group <- rep(c(0, 0, 1, 1), each = sample_count / 4)
gene2_group <- rep(c(0, 1, 0, 1), each = sample_count / 4)
concordant <- as.numeric(gene1_group == gene2_group)
expression_matrix <- rbind(
  gene1 = gene1_group + stats::rnorm(sample_count, sd = 0.05),
  gene2 = gene2_group + stats::rnorm(sample_count, sd = 0.05)
)
colnames(expression_matrix) <- paste0("sample", seq_len(sample_count))
survival_time <- stats::rexp(
  sample_count,
  rate = exp(0.8 * concordant)
)
survival_status <- rep(1, sample_count)

interaction_result <- fit.interaction.model(
  feature1 = "gene1",
  feature2 = "gene2",
  expression.data = list(test = expression_matrix),
  survival.data = list(
    test = survival::Surv(survival_time, survival_status)
  )
)
stopifnot(
  is.list(interaction_result),
  all(c("cox.uv.1", "cox.uv.2", "cox.int") %in%
    names(interaction_result)),
  length(interaction_result$cox.int) == 2L,
  is.finite(interaction_result$cox.int[["HR"]]),
  is.finite(interaction_result$cox.int[["P"]])
)

message("All local SIMMS-compatibility tests passed.")
