# Local compatibility functions for the subset of SIMMS used by PENDR.
#
# These functions reproduce the documented behaviour required by the
# edge-only pipeline without requiring installation of the SIMMS package.
# The interaction term is concordance/XNOR:
#   1 when both genes are in the same expression group, otherwise 0.

.require_survival <- function() {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The CRAN package 'survival' is required.")
  }
}

.empty_cox_result <- function() {
  list(
    cox.stats = setNames(
      rep(NA_real_, 5L),
      c("HR", "CI95L", "CI95U", "P", "n")
    ),
    cox.obj = NULL
  )
}

.create_survobj <- function(annotation, truncate.survival = 100) {
  .require_survival()

  required_columns <- c("survtime", "survstat", "survtime.unit")
  missing_columns <- setdiff(required_columns, colnames(annotation))
  if (length(missing_columns) > 0L) {
    stop(
      "Patient annotation is missing survival column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  survival_time <- suppressWarnings(as.numeric(annotation$survtime))
  survival_status <- suppressWarnings(as.numeric(annotation$survstat))
  survival_unit <- tolower(trimws(as.character(annotation$survtime.unit)))

  supported_units <- c("days", "weeks", "months", "years")
  unsupported_units <- unique(
    survival_unit[!is.na(survival_unit) & !(survival_unit %in% supported_units)]
  )
  if (length(unsupported_units) > 0L) {
    stop(
      "Unsupported survival-time unit(s): ",
      paste(unsupported_units, collapse = ", ")
    )
  }

  survival_time[survival_unit == "days"] <-
    survival_time[survival_unit == "days"] / 365.25
  survival_time[survival_unit == "weeks"] <-
    survival_time[survival_unit == "weeks"] / 52.18
  survival_time[survival_unit == "months"] <-
    survival_time[survival_unit == "months"] / 12

  survival_time[survival_time <= 0] <- 1e-05
  over_limit <- is.finite(survival_time) &
    survival_time > truncate.survival
  survival_status[over_limit] <- 0
  survival_time[over_limit] <- truncate.survival

  survival::Surv(survival_time, survival_status)
}

.resolve_annotation_column <- function(annotation, dataset_row, canonical_name) {
  if (canonical_name %in% colnames(annotation)) {
    return(annotation[[canonical_name]])
  }

  source_name <- as.character(dataset_row[[canonical_name]])
  if (
    length(source_name) != 1L ||
    is.na(source_name) ||
    !nzchar(source_name) ||
    !(source_name %in% colnames(annotation))
  ) {
    stop(
      "Cannot resolve annotation column '",
      canonical_name,
      "' for dataset '",
      dataset_row$dataset,
      "'."
    )
  }

  annotation[[source_name]]
}

load.cancer.datasets <- function(
  tumour.only = TRUE,
  with.survival.only = TRUE,
  truncate.survival = 100,
  datasets.to.load = "all",
  data.types = c("mRNA"),
  datasets.file = "datasets.txt",
  data.directory = ".",
  verbose = FALSE,
  subset = NULL
) {
  .require_survival()

  datasets_path <- file.path(data.directory, datasets.file)
  if (!file.exists(datasets_path)) {
    stop("Dataset manifest not found: ", datasets_path)
  }

  datasets <- read.delim(
    datasets_path,
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_manifest_columns <- c(
    "dataset",
    "annotation",
    "survstat",
    "survtime",
    "survtime.unit",
    data.types
  )
  missing_manifest_columns <- setdiff(
    required_manifest_columns,
    colnames(datasets)
  )
  if (length(missing_manifest_columns) > 0L) {
    stop(
      "datasets.txt is missing column(s): ",
      paste(missing_manifest_columns, collapse = ", ")
    )
  }

  if ("all" %in% datasets.to.load) {
    datasets.to.load <- datasets$dataset
  }
  unknown_datasets <- setdiff(datasets.to.load, datasets$dataset)
  if (length(unknown_datasets) > 0L) {
    stop(
      "Unknown dataset(s) requested: ",
      paste(unknown_datasets, collapse = ", ")
    )
  }

  all.data <- setNames(
    replicate(length(data.types), list(), simplify = FALSE),
    data.types
  )
  all.survobj <- list()
  all.probesets <- character()

  for (dataset_name in datasets.to.load) {
    if (verbose) {
      message("Reading dataset: ", dataset_name)
    }

    dataset_row <- datasets[
      match(dataset_name, datasets$dataset),
      ,
      drop = FALSE
    ]
    dataset_directory <- file.path(data.directory, dataset_name)
    annotation_path <- file.path(
      dataset_directory,
      as.character(dataset_row$annotation)
    )
    if (!file.exists(annotation_path)) {
      stop("Patient annotation not found: ", annotation_path)
    }

    annotation <- read.delim(
      annotation_path,
      header = TRUE,
      row.names = 1,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (tumour.only) {
      if (!("Tumour" %in% colnames(annotation))) {
        stop(
          "tumour.only=TRUE but annotation has no 'Tumour' column: ",
          annotation_path
        )
      }
      annotation <- annotation[
        !is.na(annotation$Tumour) & annotation$Tumour == "Yes",
        ,
        drop = FALSE
      ]
    }

    if (!is.null(subset)) {
      if (!all(c("Field", "Entry") %in% names(subset))) {
        stop("subset must contain named elements 'Field' and 'Entry'.")
      }
      subset_field <- as.character(subset$Field)
      if (!(subset_field %in% colnames(annotation))) {
        stop("Subset field not found in annotation: ", subset_field)
      }
      annotation <- annotation[
        annotation[[subset_field]] %in% subset$Entry,
        ,
        drop = FALSE
      ]
    }

    annotation$survstat <- .resolve_annotation_column(
      annotation,
      dataset_row,
      "survstat"
    )
    annotation$survtime <- .resolve_annotation_column(
      annotation,
      dataset_row,
      "survtime"
    )
    annotation$survtime.unit <- .resolve_annotation_column(
      annotation,
      dataset_row,
      "survtime.unit"
    )

    annotation$survtime <- suppressWarnings(
      as.numeric(annotation$survtime)
    )
    annotation$survstat <- suppressWarnings(
      as.numeric(annotation$survstat)
    )
    annotation$survtime[annotation$survtime <= 0] <- 1e-05

    if (with.survival.only) {
      annotation <- annotation[
        !is.na(annotation$survtime) & !is.na(annotation$survstat),
        ,
        drop = FALSE
      ]
    }

    rownames(annotation) <- make.names(tolower(rownames(annotation)))
    common_samples <- rownames(annotation)

    for (data_type in data.types) {
      data_file <- as.character(dataset_row[[data_type]])
      data_path <- file.path(dataset_directory, data_file)
      if (!file.exists(data_path)) {
        stop(
          "Data file for ",
          data_type,
          " not found: ",
          data_path
        )
      }

      data_matrix <- read.delim(
        data_path,
        header = TRUE,
        row.names = 1,
        check.names = TRUE
      )
      data_matrix <- as.matrix(data_matrix)
      storage.mode(data_matrix) <- "numeric"
      colnames(data_matrix) <- make.names(tolower(colnames(data_matrix)))

      nonempty_samples <- apply(
        data_matrix,
        2,
        function(x) any(!is.na(x))
      )
      data_matrix <- data_matrix[
        ,
        nonempty_samples,
        drop = FALSE
      ]

      all.data[[data_type]][[dataset_name]] <- data_matrix
      common_samples <- intersect(common_samples, colnames(data_matrix))
    }

    if (length(common_samples) == 0L) {
      warning(
        "Dataset '",
        dataset_name,
        "' has no samples shared by annotation and molecular data."
      )
      for (data_type in data.types) {
        all.data[[data_type]][[dataset_name]] <- NULL
      }
      next
    }

    annotation <- annotation[
      common_samples,
      ,
      drop = FALSE
    ]
    for (data_type in data.types) {
      all.data[[data_type]][[dataset_name]] <-
        all.data[[data_type]][[dataset_name]][
          ,
          common_samples,
          drop = FALSE
        ]
      all.probesets <- union(
        all.probesets,
        rownames(all.data[[data_type]][[dataset_name]])
      )
    }

    all.survobj[[dataset_name]] <- .create_survobj(
      annotation,
      truncate.survival = truncate.survival
    )
  }

  list(
    all.data = all.data,
    all.survobj = all.survobj,
    all.probesets = all.probesets
  )
}

.dichotomize_dataset <- function(x, split.at = "median") {
  if (length(split.at) != 1L || is.na(split.at)) {
    stop("split.at must be 'median', 'mean', or one numeric threshold.")
  }

  if (is.character(split.at) && split.at == "median") {
    threshold <- stats::median(x, na.rm = TRUE)
  } else if (is.character(split.at) && split.at == "mean") {
    threshold <- mean(x, na.rm = TRUE)
  } else {
    threshold <- suppressWarnings(as.numeric(split.at))
    if (!is.finite(threshold)) {
      stop("Invalid split.at value.")
    }
  }

  as.numeric(x > threshold)
}

.dichotomize_meta_dataset <- function(
  feature.name,
  expression.data,
  survival.data,
  data.type.ordinal = FALSE,
  centre.data = "median"
) {
  if (!is.list(expression.data) || !is.list(survival.data)) {
    stop("expression.data and survival.data must both be lists.")
  }
  if (length(expression.data) != length(survival.data)) {
    stop("Expression and survival lists must have the same length.")
  }

  groups <- numeric()
  survival_time <- numeric()
  survival_status <- numeric()

  for (dataset_index in seq_along(expression.data)) {
    expression_matrix <- expression.data[[dataset_index]]
    survival_object <- survival.data[[dataset_index]]

    if (ncol(expression_matrix) != nrow(survival_object)) {
      stop(
        "Expression and survival sample counts differ in dataset ",
        dataset_index,
        "."
      )
    }

    if (!(feature.name %in% rownames(expression_matrix))) {
      expression_values <- rep(NA_real_, ncol(expression_matrix))
    } else {
      expression_values <- as.numeric(
        expression_matrix[feature.name, , drop = TRUE]
      )
    }

    if (data.type.ordinal) {
      feature_groups <- expression_values
    } else {
      feature_groups <- .dichotomize_dataset(
        expression_values,
        split.at = centre.data
      )
    }

    groups <- c(groups, feature_groups)
    survival_time <- c(survival_time, survival_object[, 1])
    survival_status <- c(survival_status, survival_object[, 2])
  }

  list(
    groups = groups,
    survtime = survival_time,
    survstat = survival_status
  )
}

.fit_coxmodel <- function(
  groups,
  survival_time,
  survival_status,
  rounding = 3,
  data.type.ordinal = FALSE
) {
  .require_survival()

  valid_groups <- unique(groups[!is.na(groups)])
  if (
    length(groups) != length(survival_time) ||
    length(groups) != length(survival_status) ||
    length(valid_groups) < 2L
  ) {
    return(.empty_cox_result())
  }

  if (data.type.ordinal) {
    if (0 %in% valid_groups) {
      groups <- factor(groups, levels = c(0, setdiff(valid_groups, 0)))
    } else {
      groups <- factor(groups, levels = sort(valid_groups))
    }
  }

  model_data <- data.frame(
    survival_time = survival_time,
    survival_status = survival_status,
    groups = groups
  )

  cox_fit <- tryCatch(
    survival::coxph(
      survival::Surv(survival_time, survival_status) ~ groups,
      data = model_data
    ),
    error = function(error) NULL
  )
  if (is.null(cox_fit)) {
    return(.empty_cox_result())
  }

  cox_summary <- summary(cox_fit)
  if (
    is.null(cox_summary$conf.int) ||
    nrow(cox_summary$conf.int) < 1L ||
    is.null(cox_summary$coefficients)
  ) {
    return(.empty_cox_result())
  }

  if (data.type.ordinal) {
    result <- cbind(
      HR = cox_summary$conf.int[, 1],
      CI95L = cox_summary$conf.int[, 3],
      CI95U = cox_summary$conf.int[, 4],
      P = cox_summary$coefficients[, 5],
      n = rep(sum(!is.na(groups)), nrow(cox_summary$conf.int))
    )
    rownames(result) <- sub("^groups", "", rownames(result))
    return(list(cox.stats = result, cox.obj = cox_fit))
  }

  result <- c(
    HR = cox_summary$conf.int[1, 1],
    CI95L = cox_summary$conf.int[1, 3],
    CI95U = cox_summary$conf.int[1, 4],
    P = cox_summary$coefficients[1, 5],
    n = sum(!is.na(groups))
  )
  if (!isFALSE(rounding) && is.finite(rounding)) {
    result[c("HR", "CI95L", "CI95U")] <- round(
      result[c("HR", "CI95L", "CI95U")],
      digits = rounding
    )
  }

  list(cox.stats = result, cox.obj = cox_fit)
}

.calculate_meta_survival <- function(
  feature.name,
  expression.data,
  survival.data,
  rounding = 3,
  data.type.ordinal = FALSE,
  centre.data = "median"
) {
  if (
    !is.list(expression.data) ||
    !is.list(survival.data) ||
    length(expression.data) != length(survival.data)
  ) {
    return(.empty_cox_result())
  }

  dichotomized <- .dichotomize_meta_dataset(
    feature.name = feature.name,
    expression.data = expression.data,
    survival.data = survival.data,
    data.type.ordinal = data.type.ordinal,
    centre.data = centre.data
  )

  if (
    length(dichotomized$groups) == 0L ||
    all(is.na(dichotomized$groups))
  ) {
    return(.empty_cox_result())
  }

  .fit_coxmodel(
    groups = dichotomized$groups,
    survival_time = dichotomized$survtime,
    survival_status = dichotomized$survstat,
    rounding = rounding,
    data.type.ordinal = data.type.ordinal
  )
}

fit.interaction.model <- function(
  feature1,
  feature2,
  expression.data,
  survival.data,
  data.type.ordinal = FALSE,
  centre.data = "median"
) {
  .require_survival()

  groups1 <- .dichotomize_meta_dataset(
    feature.name = feature1,
    expression.data = expression.data,
    survival.data = survival.data,
    data.type.ordinal = data.type.ordinal,
    centre.data = centre.data
  )
  groups2 <- .dichotomize_meta_dataset(
    feature.name = feature2,
    expression.data = expression.data,
    survival.data = survival.data,
    data.type.ordinal = data.type.ordinal,
    centre.data = centre.data
  )

  interaction_hr <- NA_real_
  interaction_p <- NA_real_

  valid_input <-
    length(groups1$groups) == length(groups2$groups) &&
    length(groups1$groups) == length(groups1$survtime) &&
    length(unique(groups1$groups[!is.na(groups1$groups)])) >= 2L &&
    length(unique(groups2$groups[!is.na(groups2$groups)])) >= 2L

  if (valid_input) {
    model_data <- data.frame(
      survival_time = groups1$survtime,
      survival_status = groups1$survstat,
      group1 = factor(groups1$groups),
      group2 = factor(groups2$groups)
    )
    model_data$concordant <- as.numeric(
      model_data$group1 == model_data$group2
    )

    interaction_fit <- tryCatch(
      survival::coxph(
        survival::Surv(survival_time, survival_status) ~
          group1 + group2 + concordant,
        data = model_data
      ),
      error = function(error) NULL
    )

    if (!is.null(interaction_fit)) {
      interaction_summary <- summary(interaction_fit)
      coefficients <- interaction_summary$coefficients
      if (!is.null(coefficients) && "concordant" %in% rownames(coefficients)) {
        interaction_hr <- as.numeric(coefficients["concordant", 2])
        interaction_p <- as.numeric(coefficients["concordant", 5])
      }
    }
  }

  list(
    cox.uv.1 = .calculate_meta_survival(
      feature.name = feature1,
      expression.data = expression.data,
      survival.data = survival.data,
      data.type.ordinal = data.type.ordinal,
      centre.data = centre.data
    ),
    cox.uv.2 = .calculate_meta_survival(
      feature.name = feature2,
      expression.data = expression.data,
      survival.data = survival.data,
      data.type.ordinal = data.type.ordinal,
      centre.data = centre.data
    ),
    cox.int = c(HR = interaction_hr, P = interaction_p)
  )
}

get.adjacency.matrix <- function(subnets.file = NULL) {
  if (is.null(subnets.file) || !file.exists(subnets.file)) {
    stop("Subnetwork file not found: ", subnets.file)
  }

  lines <- readLines(subnets.file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  header_indices <- which(grepl("^#", lines))
  if (length(header_indices) == 0L) {
    stop("No subnetwork headers beginning with '#' were found.")
  }

  adjacency_matrices <- list()
  for (header_position in seq_along(header_indices)) {
    start_index <- header_indices[[header_position]]
    end_index <- if (header_position < length(header_indices)) {
      header_indices[[header_position + 1L]] - 1L
    } else {
      length(lines)
    }

    raw_name <- sub("\t$", "", lines[[start_index]])
    graph_name <- make.names(raw_name)
    edge_lines <- if (end_index > start_index) {
      lines[(start_index + 1L):end_index]
    } else {
      character()
    }

    parsed_edges <- lapply(edge_lines, function(edge_line) {
      fields <- strsplit(edge_line, "\t", fixed = TRUE)[[1]]
      fields <- trimws(fields)
      fields <- fields[nzchar(fields)]
      if (length(fields) >= 3L) {
        genes <- fields[2:3]
      } else if (length(fields) == 2L) {
        genes <- fields[1:2]
      } else {
        return(NULL)
      }
      gsub("[()]", "-", genes)
    })
    parsed_edges <- Filter(Negate(is.null), parsed_edges)

    if (length(parsed_edges) == 0L) {
      adjacency_matrices[[graph_name]] <- matrix(
        numeric(),
        nrow = 0L,
        ncol = 0L
      )
      next
    }

    edge_matrix <- do.call(rbind, parsed_edges)
    vertices <- unique(c(edge_matrix[, 1], edge_matrix[, 2]))
    adjacency <- matrix(
      0,
      nrow = length(vertices),
      ncol = length(vertices),
      dimnames = list(vertices, vertices)
    )
    for (edge_index in seq_len(nrow(edge_matrix))) {
      gene1 <- edge_matrix[edge_index, 1]
      gene2 <- edge_matrix[edge_index, 2]
      adjacency[gene1, gene2] <- 1
      adjacency[gene2, gene1] <- 1
    }

    adjacency_matrices[[graph_name]] <- adjacency
  }

  adjacency_matrices
}
