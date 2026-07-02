filter_kegg_by_sample_presence <- function(
  input_xlsx,
  output_dir,
  output_xlsx,
  min_samples = 10,
  cecal_sheet = "Cecal",
  serum_sheet = "Serum",
  species_sheet = "Species"
) {

  suppressPackageStartupMessages({
    library(openxlsx)
    library(dplyr)
  })

  # -----------------------------
  # Create output directory
  # -----------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  out_file <- file.path(output_dir, output_xlsx)

  # -----------------------------
  # Read sheets (keep rownames)
  # -----------------------------
  message("Reading Excel sheets...")

  cecal_df   <- read.xlsx(input_xlsx, sheet = cecal_sheet, rowNames = TRUE)
  serum_df   <- read.xlsx(input_xlsx, sheet = serum_sheet, rowNames = TRUE)
  species_df <- read.xlsx(input_xlsx, sheet = species_sheet, rowNames = TRUE)

  # -----------------------------
  # Helper: prevalence filter
  # -----------------------------
  filter_rows <- function(df, min_samples) {
    keep <- apply(
      df,
      1,
      function(x) sum(!is.na(x) & x != 0) >= min_samples
    )
    df[keep, , drop = FALSE]
  }

  # -----------------------------
  # Apply filtering
  # -----------------------------
  message("Filtering Cecal KEGG IDs...")
  cecal_filt <- filter_rows(cecal_df, min_samples)

  message("Filtering Serum KEGG IDs...")
  serum_filt <- filter_rows(serum_df, min_samples)

  message("Filtering Species features...")
  species_filt <- filter_rows(species_df, min_samples)

  # -----------------------------
  # Summary
  # -----------------------------
  message("Rows retained after filtering:")
  message("  Cecal   : ", nrow(cecal_filt))
  message("  Serum   : ", nrow(serum_filt))
  message("  Species : ", nrow(species_filt))

  # -----------------------------
  # Write output workbook
  # -----------------------------
  message("Writing filtered workbook...")

  wb <- createWorkbook()

  addWorksheet(wb, cecal_sheet)
  writeData(wb, cecal_sheet, cecal_filt, rowNames = TRUE)

  addWorksheet(wb, serum_sheet)
  writeData(wb, serum_sheet, serum_filt, rowNames = TRUE)

  addWorksheet(wb, species_sheet)
  writeData(wb, species_sheet, species_filt, rowNames = TRUE)

  saveWorkbook(wb, out_file, overwrite = TRUE)

  message("Output saved at: ", out_file)
  message("Prevalence filtering completed successfully ✔")

  invisible(
    list(
      cecal   = cecal_filt,
      serum   = serum_filt,
      species = species_filt
    )
  )
}

filter_kegg_by_sample_presence(
  input_xlsx = "data/processed_data/008_kegg_summed/control_kegg_summed.xlsx",
  output_dir = "data/processed_data/009_prevalence_filtered",
  output_xlsx = "control_kegg_minsample_filtered.xlsx",
  min_samples = 10
)

filter_kegg_by_sample_presence(
  input_xlsx = "data/processed_data/008_kegg_summed/treatment_kegg_summed.xlsx",
  output_dir = "data/processed_data/009_prevalence_filtered",
  output_xlsx = "treatment_kegg_minsample_filtered.xlsx",
  min_samples = 10
)
