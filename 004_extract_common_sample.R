extract_common_columns_keep_first <- function(
  input_xlsx,
  sheets = c("Cecal", "Serum", "Species"),
  output_dir,
  output_xlsx = "common_columns_with_ids.xlsx"
) {

  # -------------------------------
  # Load libraries
  # -------------------------------
  suppressPackageStartupMessages({
    library(openxlsx)
    library(dplyr)
    library(purrr)
  })

  # -------------------------------
  # Create output directory
  # -------------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  out_file <- file.path(output_dir, output_xlsx)

  # -------------------------------
  # Read sheets
  # -------------------------------
  message("Reading sheets...")
  sheet_data <- map(sheets, ~ read.xlsx(input_xlsx, sheet = .x))
  names(sheet_data) <- sheets

  # -------------------------------
  # Identify common NON-ID columns
  # -------------------------------
  message("Identifying common columns (excluding first column)...")

  non_id_cols <- map(
    sheet_data,
    ~ colnames(.x)[-1]
  )

  common_cols <- reduce(non_id_cols, intersect)

  if (length(common_cols) == 0) {
    stop("No common data columns found across sheets.")
  }

  message("Common data columns: ", paste(common_cols, collapse = ", "))

  # -------------------------------
  # Subset each sheet:
  #   first column + common columns
  # -------------------------------
  filtered_data <- imap(
    sheet_data,
    ~ {
      id_col <- colnames(.x)[1]
      dplyr::select(.x, all_of(c(id_col, common_cols)))
    }
  )

  # -------------------------------
  # Write output workbook
  # -------------------------------
  message("Writing output workbook...")

  wb <- createWorkbook()

  walk(names(filtered_data), function(sheet_nm) {
    addWorksheet(wb, sheet_nm)
    writeData(wb, sheet_nm, filtered_data[[sheet_nm]])
  })

  saveWorkbook(wb, out_file, overwrite = TRUE)

  message("Output saved at: ", out_file)
  message("Extraction completed successfully ✔")

  invisible(filtered_data)
}

extract_common_columns_keep_first(
  input_xlsx = "data/processed_data/004_integrated/control_cecal_serum_species_merged.xlsx",
  sheets     = c("Cecal", "Serum", "Species"),
  output_dir = "data/processed_data/005_common_samples",
  output_xlsx = "control_common_sample.xlsx"
)

extract_common_columns_keep_first(
  input_xlsx = "data/processed_data/004_integrated/treatment_cecal_serum_species_merged.xlsx",
  sheets     = c("Cecal", "Serum", "Species"),
  output_dir = "data/processed_data/005_common_samples",
  output_xlsx = "treatment_common_sample.xlsx"
)
