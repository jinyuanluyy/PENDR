assert_files_exist <- function(paths, label = "required file") {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "Missing ",
      label,
      if (length(missing_paths) > 1L) "s" else "",
      ": ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(paths)
}

assert_packages_installed <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      "Missing R package",
      if (length(missing_packages) > 1L) "s" else "",
      ": ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(packages)
}

write_tsv <- function(data, path, row.names = FALSE) {
  output_directory <- dirname(path)
  if (!identical(output_directory, ".")) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  }

  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = row.names,
    col.names = TRUE
  )
  invisible(path)
}

normalize_subnetwork_name <- function(x) {
  normalized <- make.names(trimws(x))
  sub("^X\\.", "", normalized)
}
