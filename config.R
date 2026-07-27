pendr_config <- list(
  data_type = "mRNA",
  interaction_p_value_threshold = 0.05,
  truncate_survival_years = 10,
  rank_product_permutations = 10000L,
  random_seed = 2026L,
  submodule_fdr_threshold = 0.01,
  max_candidate_submodules = 15L,
  expected_pathway_guided_submodules = 2732L,
  pathway_guided_submodules_file = "data/pathway_guided_submodules.txt",
  expression_feature_suffix = "_at",
  cohort_data_directory = "testdata",
  datasets_file = "datasets.txt",
  optimization_expression_file =
    "data/optimization/pdac_discrimination_expression.tsv",
  optimization_sample_metadata_file =
    "data/optimization/pdac_discrimination_samples.tsv",
  tumor_class_label = "Tumor",
  normal_class_label = "Normal",
  cohorts = data.frame(
    id = c("CPTAC-PDAC", "TCGA-PDAC", "ICGC-PDAC-CA"),
    label = c("CPTAC-PDAC", "TCGA-PDAC", "ICGC-PDAC-CA"),
    dataset_name = c(
      "Pancreatic ductal adenocarcinoma1",
      "Pancreatic ductal adenocarcinoma2",
      "Pancreatic ductal adenocarcinoma3"
    ),
    stringsAsFactors = FALSE
  )
)
