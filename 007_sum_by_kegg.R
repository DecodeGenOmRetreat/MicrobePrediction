sum_by_kegg <- function(
  input_file,
  cecum_sheet,
  serum_sheet,
  species_sheet,
  metabolite_col,
  kegg_col,
  output_dir,
  output_file
) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(readxl)
    library(openxlsx)
  })

  # -----------------------------
  # Create output directory
  # -----------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  } else {
    message("Using existing output directory: ", output_dir)
  }

  out_path <- file.path(output_dir, output_file)

  # -----------------------------
  # Read sheets
  # -----------------------------
  cecum   <- read_excel(input_file, sheet = cecum_sheet)
  serum   <- read_excel(input_file, sheet = serum_sheet)
  species <- read_excel(input_file, sheet = species_sheet)

  # -----------------------------
  # Column checks
  # -----------------------------
  if (!metabolite_col %in% colnames(cecum))
    stop("Metabolite column missing in cecum")

  if (!metabolite_col %in% colnames(serum))
    stop("Metabolite column missing in serum")

  if (!kegg_col %in% colnames(cecum))
    stop("KEGG column missing in cecum")

  if (!kegg_col %in% colnames(serum))
    stop("KEGG column missing in serum")

  # -------------------------------------------------------
  # Helper: Remove duplicate metabolites → Sum by KEGG
  # -------------------------------------------------------
  process_sheet <- function(df) {
    df %>%
      distinct(.data[[metabolite_col]], .keep_all = TRUE) %>%
      group_by(.data[[kegg_col]]) %>%
      summarise(
        across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::select(all_of(kegg_col), everything())
  }

  # -----------------------------
  # Apply processing
  # -----------------------------
  message("Summing Cecal by KEGG...")
  cecum_summed <- process_sheet(cecum)

  message("Summing Serum by KEGG...")
  serum_summed <- process_sheet(serum)

  # Species unchanged
  species_final <- species

  # -----------------------------
  # Write Excel output
  # -----------------------------
  wb <- createWorkbook()

  addWorksheet(wb, "Cecal")
  writeData(wb, "Cecal", cecum_summed)

  addWorksheet(wb, "Serum")
  writeData(wb, "Serum", serum_summed)

  addWorksheet(wb, "Species")
  writeData(wb, "Species", species_final)

  saveWorkbook(wb, out_path, overwrite = TRUE)

  message("File saved to: ", out_path)

  invisible(
    list(
      cecal   = cecum_summed,
      serum   = serum_summed,
      species = species_final
    )
  )
}

sum_by_kegg(
  input_file   = "data/processed_data/007_zero_filled/control_zerofilled.xlsx",
  cecum_sheet  = "Cecal",
  serum_sheet  = "Serum",
  species_sheet= "Species",
  metabolite_col = "METABOLITE",
  kegg_col       = "KEGG_ID",
  output_dir     = "data/processed_data/008_kegg_summed",
  output_file    = "control_kegg_summed.xlsx"
)

sum_by_kegg(
  input_file   = "data/processed_data/007_zero_filled/treatment_zerofilled.xlsx",
  cecum_sheet  = "Cecal",
  serum_sheet  = "Serum",
  species_sheet= "Species",
  metabolite_col = "METABOLITE",
  kegg_col       = "KEGG_ID",
  output_dir     = "data/processed_data/008_kegg_summed",
  output_file    = "treatment_kegg_summed.xlsx"
)