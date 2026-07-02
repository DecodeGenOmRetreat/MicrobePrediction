merge_treatment_cecal_serum_species_files <- function(
  input_dir,
  output_dir,
  metabolite_col,
  species_col,
  output_xlsx = "treatment_cecal_serum_species_merged.xlsx"
) {

  # Load libraries
  library(data.table)
  library(dplyr)
  library(purrr)
  library(openxlsx)

  message("Step 1: Listing files...")

  cecal_files <- list.files(
    input_dir,
    pattern = "^treatment_cecal_.*_merged.tsv",
    full.names = TRUE,
    recursive = TRUE
  )

  serum_files <- list.files(
    input_dir,
    pattern = "^treatment_serum_.*_merged.tsv",
    full.names = TRUE,
    recursive = TRUE
  )

  species_files <- list.files(
    input_dir,
    pattern = "^treatment_species_.*_merged.tsv",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(cecal_files) < 2)   stop("Need at least two treatment cecal files.")
  if (length(serum_files) < 2)   stop("Need at least two treatment serum files.")
  if (length(species_files) < 2) stop("Need at least two treatment species files.")

  message("  Cecal files   : ", length(cecal_files))
  message("  Serum files   : ", length(serum_files))
  message("  Species files : ", length(species_files))

  # Helper: read TSV files
  read_tsv <- function(files) {
    map(files, ~ fread(.x, sep = "\t"))
  }

  message("Step 2: Reading files...")
  cecal_list   <- read_tsv(cecal_files)
  serum_list   <- read_tsv(serum_files)
  species_list <- read_tsv(species_files)

  # Helper: inner merge
 inner_merge <- function(df_list, by_col) {
  df_list <- map(df_list, ~ distinct(.x, .data[[by_col]], .keep_all = TRUE))
  reduce(df_list, inner_join, by = by_col)
}

  message("Step 3: Merging cecal files (by metabolite)...")
  cecal_merged <- inner_merge(cecal_list, metabolite_col)

  message("Step 4: Merging serum files (by metabolite)...")
  serum_merged <- inner_merge(serum_list, metabolite_col)

  message("Step 5: Merging species files (by species)...")
  species_merged <- inner_merge(species_list, species_col)

  message("Step 6: Merge summaries")
  message("  Cecal   : ", nrow(cecal_merged),   " rows, ", ncol(cecal_merged),   " columns")
  message("  Serum   : ", nrow(serum_merged),   " rows, ", ncol(serum_merged),   " columns")
  message("  Species : ", nrow(species_merged), " rows, ", ncol(species_merged), " columns")

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Write Excel workbook
  message("Step 7: Writing Excel output...")

  wb <- createWorkbook()

  addWorksheet(wb, "Cecal")
  writeData(wb, "Cecal", cecal_merged)

  addWorksheet(wb, "Serum")
  writeData(wb, "Serum", serum_merged)

  addWorksheet(wb, "Species")
  writeData(wb, "Species", species_merged)

  out_file <- file.path(output_dir, output_xlsx)
  saveWorkbook(wb, out_file, overwrite = TRUE)

  message("Output saved as: ", out_file)
  message("Process completed successfully!")

  invisible(
    list(
      cecal   = cecal_merged,
      serum   = serum_merged,
      species = species_merged
    )
  )
}

merge_treatment_cecal_serum_species_files(
  input_dir      = "data/processed_data",
  output_dir     = "data/processed_data/004_integrated",
  metabolite_col = "METABOLITE",
  species_col    = "Species"
)
