fill_blank_cells_with_zero <- function(
  input_xlsx,
  output_dir,
  output_xlsx,
  sheets = NULL
) {

  suppressPackageStartupMessages({
    library(openxlsx)
    library(dplyr)
    library(purrr)
  })

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  out_file <- file.path(output_dir, output_xlsx)

  # Detect sheets if not provided
  if (is.null(sheets)) {
    sheets <- getSheetNames(input_xlsx)
  }

  message("Processing sheets: ", paste(sheets, collapse = ", "))

  # Read → replace blanks → store
  cleaned_data <- map(
    sheets,
    function(sheet) {
      df <- read.xlsx(input_xlsx, sheet = sheet)

      # Replace NA or empty strings with 0
      df[] <- lapply(df, function(x) {
        ifelse(is.na(x) | x == "", 0, x)
      })

      df
    }
  )
  names(cleaned_data) <- sheets

  # Write to Excel
  wb <- createWorkbook()

  walk(sheets, function(sheet) {
    addWorksheet(wb, sheet)
    writeData(wb, sheet, cleaned_data[[sheet]])
  })

  saveWorkbook(wb, out_file, overwrite = TRUE)

  message("Output saved at: ", out_file)
  message("Blank cells replaced with 0 successfully ✔")

  invisible(cleaned_data)
}

fill_blank_cells_with_zero(
  input_xlsx = "data/processed_data/006_kegg_mapped/control_kegg.xlsx",
  output_dir = "data/processed_data/007_zero_filled",
  output_xlsx = "control_zerofilled.xlsx"
)

fill_blank_cells_with_zero(
  input_xlsx = "data/processed_data/006_kegg_mapped/treatment_kegg.xlsx",
  output_dir = "data/processed_data/007_zero_filled",
  output_xlsx = "treatment_zerofilled.xlsx"
)