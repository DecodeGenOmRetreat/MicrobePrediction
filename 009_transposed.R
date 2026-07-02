transpose_sheets_pipeline <- function(
  input_file,
  output_dir,
  output_xlsx,
  sheets = c("Cecal", "Serum", "Species")
) {
  # -------------------------------
  # Step 1: Load libraries
  # -------------------------------
  suppressPackageStartupMessages({
    library(openxlsx)
    library(dplyr)
    library(readxl)
  })

  # -------------------------------
  # Step 2: Create output directory
  # -------------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  } else {
    message("Using existing output directory: ", output_dir)
  }

  out_file <- file.path(output_dir, output_xlsx)

  # -------------------------------
  # Step 3: Read and transpose sheets
  # -------------------------------
  transpose_sheet <- function(sheet_name) {
    message("Processing sheet: ", sheet_name)
    df <- read_excel(input_file, sheet = sheet_name)

    # Keep first column as rownames
    row_names <- df[[1]]

    df_data <- df[, -1, drop = FALSE]

    df_t <- as.data.frame(t(df_data))
    colnames(df_t) <- row_names
    df_t <- tibble::rownames_to_column(df_t, var = "Sample")
    return(df_t)
  }

  # Apply to all sheets
  transposed_list <- lapply(sheets, transpose_sheet)
  names(transposed_list) <- sheets

  # -------------------------------
  # Step 4: Write to Excel
  # -------------------------------
  wb <- createWorkbook()
  for (sheet_name in sheets) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, transposed_list[[sheet_name]])
  }

  saveWorkbook(wb, out_file, overwrite = TRUE)
  message("Output saved to: ", out_file)

  invisible(transposed_list)
}

# -------------------------------
# Example usage
# -------------------------------
transpose_sheets_pipeline(
  input_file  = "data/processed_data/009_prevalence_filtered/control_kegg_minsample_filtered.xlsx",
  output_dir  = "data/processed_data/010_transposed",
  output_xlsx = "control_transposed.xlsx"
)

transpose_sheets_pipeline(
  input_file  = "data/processed_data/009_prevalence_filtered/treatment_kegg_minsample_filtered.xlsx",
  output_dir  = "data/processed_data/010_transposed",
  output_xlsx = "treatment_transposed.xlsx"
)

