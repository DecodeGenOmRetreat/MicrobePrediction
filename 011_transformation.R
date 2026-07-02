full_transform_pipeline <- function(
  input_file,
  output_dir,
  output_xlsx,
  species_sheet = "Species",
  cecal_sheet = "Cecal",
  serum_sheet = "Serum"
) {
  suppressPackageStartupMessages({
    library(openxlsx)
    library(microbiome)
    library(tibble)
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
  # Read sheets (row names = samples)
  # -----------------------------
  message("Reading input sheets...")

  species <- read.xlsx(input_file, sheet = species_sheet, rowNames = TRUE,sep.names = " ")
  cecal <- read.xlsx(input_file, sheet = cecal_sheet, rowNames = TRUE)
  serum <- read.xlsx(input_file, sheet = serum_sheet, rowNames = TRUE)

  # -----------------------------
  # Sanity checks
  # -----------------------------
  if (
    !all(rownames(species) == rownames(cecal)) ||
      !all(rownames(species) == rownames(serum))
  ) {
    stop("Sample names (rownames) do not match across sheets")
  }

  # -----------------------------
  # Transformations
  # -----------------------------
  message("Applying transformations...")

  species_clr <- microbiome::transform(species, "clr")
  cecal_log <- microbiome::transform(cecal, "log10p")
  serum_log <- microbiome::transform(serum, "log10p")

  # -----------------------------
  # Add Sample column explicitly
  # -----------------------------
  species_out <- rownames_to_column(
    as.data.frame(species_clr),
    var = "Sample"
  )

  cecal_out <- rownames_to_column(
    as.data.frame(cecal_log),
    var = "Sample"
  )

  serum_out <- rownames_to_column(
    as.data.frame(serum_log),
    var = "Sample"
  )

  # -----------------------------
  # Write Excel output
  # -----------------------------
  message("Writing transformed data to Excel...")

  wb <- createWorkbook()

  addWorksheet(wb, "Cecal_log10p")
  writeData(wb, "Cecal_log10p", cecal_out)

  addWorksheet(wb, "Serum_log10p")
  writeData(wb, "Serum_log10p", serum_out)

  addWorksheet(wb, "Species_CLR")
  writeData(wb, "Species_CLR", species_out)

  saveWorkbook(wb, out_file, overwrite = TRUE)

  message("Pipeline completed successfully ✔")
  message("Output saved at: ", out_file)

  invisible(list(
    species = species_out,
    cecal = cecal_out,
    serum = serum_out
  ))
}

# -------------------------------
# Example usage
# -------------------------------
full_transform_pipeline(
  input_file = "data/processed_data/011_imputed/control_imputed.xlsx",
  output_dir = "data/processed_data/012_transformed",
  output_xlsx = "control_transformed.xlsx"
)

full_transform_pipeline(
  input_file = "data/processed_data/011_imputed/treatment_imputed.xlsx",
  output_dir = "data/processed_data/012_transformed",
  output_xlsx = "treatment_transformed.xlsx"
)
