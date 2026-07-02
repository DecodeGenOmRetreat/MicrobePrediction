map_kegg_cecal_serum_pipeline <- function(
  input_xlsx,
  output_dir,
  output_xlsx,
  metabolite_col = "METABOLITE",
  cecal_sheet = "Cecal",
  serum_sheet = "Serum",
  species_sheet = "Species"
) {
  # -------------------------------
  # Step 1: Load libraries
  # -------------------------------
  suppressPackageStartupMessages({
    library(openxlsx)
    library(KEGGREST)
    library(dplyr)
    library(tidyr)
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
  # Step 3: Read Excel sheets
  # -------------------------------
  message("Reading Excel sheets...")

  cecal_df <- read.xlsx(input_xlsx, sheet = cecal_sheet)
  serum_df <- read.xlsx(input_xlsx, sheet = serum_sheet)
  species_df <- read.xlsx(input_xlsx, sheet = species_sheet)

  # -------------------------------
  # Step 4: Normalize metabolite names
  # -------------------------------
  cecal_df[[metabolite_col]] <- tolower(trimws(cecal_df[[metabolite_col]]))
  serum_df[[metabolite_col]] <- tolower(trimws(serum_df[[metabolite_col]]))

  # -------------------------------
  # Step 5: Download KEGG database
  # -------------------------------
  message("Downloading KEGG compound database...")
  kegg_compounds <- keggList("compound")

  kegg_db <- data.frame(
    KEGG_ID = names(kegg_compounds),
    Name = as.character(kegg_compounds),
    stringsAsFactors = FALSE
  )

  kegg_long <- kegg_db %>%
    separate_rows(Name, sep = ";") %>%
    mutate(Name = tolower(trimws(Name)))

  # -------------------------------
  # Step 6: Map Cecal metabolites
  # -------------------------------
  message("Mapping Cecal metabolites to KEGG IDs...")
  cecal_kegg <- cecal_df %>%
    inner_join(
      kegg_long,
      by = setNames("Name", metabolite_col),
      relationship = "many-to-many"
    ) %>%
    relocate(KEGG_ID, .after = all_of(metabolite_col))

  # -------------------------------
  # Step 7: Map Serum metabolites
  # -------------------------------
  message("Mapping Serum metabolites to KEGG IDs...")
  serum_kegg <- serum_df %>%
    inner_join(
      kegg_long,
      by = setNames("Name", metabolite_col),
      relationship = "many-to-many"
    ) %>%
    relocate(KEGG_ID, .after = all_of(metabolite_col))

  # -------------------------------
  # Step 8: Write Excel workbook
  # -------------------------------
  message("Writing output Excel file...")

  wb <- createWorkbook()

  addWorksheet(wb, "Cecal")
  writeData(wb, "Cecal", cecal_kegg)

  addWorksheet(wb, "Serum")
  writeData(wb, "Serum", serum_kegg)

  addWorksheet(wb, "Species")
  writeData(wb, "Species", species_df) # unchanged

  saveWorkbook(wb, out_file, overwrite = TRUE)

  message("Output saved at: ", out_file)
  message("KEGG mapping pipeline completed successfully ✔")

  invisible(
    list(
      cecal = cecal_kegg,
      serum = serum_kegg,
      species = species_df
    )
  )
}

map_kegg_cecal_serum_pipeline(
  input_xlsx = "data/processed_data/005_common_samples/control_common_sample.xlsx",
  output_dir = "data/processed_data/006_kegg_mapped",
  output_xlsx = "control_kegg.xlsx"
)

map_kegg_cecal_serum_pipeline(
  input_xlsx = "data/processed_data/005_common_samples/treatment_common_sample.xlsx",
  output_dir = "data/processed_data/006_kegg_mapped",
  output_xlsx = "treatment_kegg.xlsx"
)