impute_sheets_pipeline <- function(
  input_file,
  species_sheet = "Species",
  cecal_sheet = "Cecal",
  serum_sheet = "Serum",
  output_dir = "data/processed_data/imputed",
  output_xlsx = "imputed_cecal_serum_species.xlsx",
  impute_method = c("none", "knn", "median", "mean", "missforest")
) {
  # -------------------------------
  # Load libraries
  # -------------------------------
  suppressPackageStartupMessages({
    library(openxlsx)
    library(dplyr)
    library(missForest)
  })

  impute_method <- match.arg(impute_method)

  # -------------------------------
  # Create output directory
  # -------------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  out_file <- file.path(output_dir, output_xlsx)

  # -------------------------------
  # Helper: imputation (preserve row names)
  # -------------------------------
  impute_data <- function(df, method) {
    # Preserve row names as a column
    df <- df %>% tibble::rownames_to_column(var = "Sample")
    numeric_cols <- names(df)[sapply(df, is.numeric)]

    if (method == "none") {
      return(df)
    }

    if (method == "knn") {
      message("Applying KNN imputation...")
      if (!requireNamespace("impute", quietly = TRUE)) {
        stop("Package 'impute' needed for KNN. Install it first.")
      }
      df[numeric_cols] <- as.data.frame(
        impute::impute.knn(as.matrix(df[numeric_cols]))$data
      )
      return(df)
    }

    if (method == "median") {
      message("Applying median imputation...")
      df[numeric_cols] <- lapply(df[numeric_cols], function(col) {
        col[is.na(col)] <- median(col, na.rm = TRUE)
        col
      })
      return(df)
    }

    if (method == "mean") {
      message("Applying mean imputation...")
      df[numeric_cols] <- lapply(df[numeric_cols], function(col) {
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
      })
      return(df)
    }

    if (method == "missforest") {
      message("Applying Random Forest imputation (missForest)...")
      if (!requireNamespace("missForest", quietly = TRUE)) {
        stop("Package 'missForest' needed. Install it first.")
      }
      df[numeric_cols] <- as.data.frame(
        missForest::missForest(as.matrix(df[numeric_cols]))$ximp
      )
      return(df)
    }
  }

  # -------------------------------
  # Read sheets
  # -------------------------------
  message("Reading Excel sheets...")
  species <- read.xlsx(input_file, sheet = species_sheet, rowNames = TRUE,sep.names = " ")
  cecal <- read.xlsx(input_file, sheet = cecal_sheet, rowNames = TRUE)
  serum <- read.xlsx(input_file, sheet = serum_sheet, rowNames = TRUE)

  # -------------------------------
  # Apply imputation
  # -------------------------------
  message("Applying imputation using method: ", impute_method)
  species_imp <- impute_data(species, impute_method)
  cecal_imp <- impute_data(cecal, impute_method)
  serum_imp <- impute_data(serum, impute_method)

  # -------------------------------
  # Write imputed sheets to Excel
  # -------------------------------
  wb <- createWorkbook()

  addWorksheet(wb, cecal_sheet)
  writeData(wb, cecal_sheet, cecal_imp)

  addWorksheet(wb, serum_sheet)
  writeData(wb, serum_sheet, serum_imp)

  addWorksheet(wb, species_sheet)
  writeData(wb, species_sheet, species_imp)

  saveWorkbook(wb, out_file, overwrite = TRUE)
  message("Imputed data saved to: ", out_file)

  invisible(list(
    species = species_imp,
    cecal = cecal_imp,
    serum = serum_imp
  ))
}

# -------------------------------
# Example usage
# -------------------------------
impute_sheets_pipeline(
  input_file = "data/processed_data/010_transposed/control_transposed.xlsx",
  output_dir = "data/processed_data/011_imputed",
  output_xlsx = "control_imputed.xlsx",
  impute_method = "missforest"
)

impute_sheets_pipeline(
  input_file = "data/processed_data/010_transposed/treatment_transposed.xlsx",
  output_dir = "data/processed_data/011_imputed",
  output_xlsx = "treatment_imputed.xlsx",
  impute_method = "missforest"
)
