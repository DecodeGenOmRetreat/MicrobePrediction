merge_sheets <- function(file_name, column_name,
                         out_dir) {

  # Step 1: Load libraries
  message("Step 1: Loading required libraries...")
  library(readxl)
  library(dplyr)
  library(purrr)
  library(data.table)
  message("   Libraries loaded successfully!")

  # Step 2: Check if file exists
  message("Step 2: Checking if file exists...")
  if (!file.exists(file_name)) {
    stop("   Error: File not found - ", file_name)
  }
  message("   File found: ", file_name)

  # Step 3: Get all sheet names
  message("Step 3: Reading sheet names...")
  sheet_names <- excel_sheets(file_name)
  message("   Found ", length(sheet_names), " sheets: ",
          paste(sheet_names, collapse = ", "))

  # Step 4: Read all sheets
  message("Step 4: Reading data from each sheet...")
  all_sheets <- lapply(sheet_names, function(sheet) {
    data <- read_excel(file_name, sheet = sheet)
    message("   - Sheet '", sheet, "': ",
            nrow(data), " rows, ", ncol(data), " columns")
    data
  })
  message("   All sheets read successfully!")

  # Step 5: Merge sheets
  message("Step 5: Merging sheets by column '", column_name, "'...")
  merged_data <- reduce(all_sheets, full_join, by = column_name)
  message("   Sheets merged successfully!")

  # Step 6: Summary
  message("Step 6: Final Summary")
  message("   Total rows: ", nrow(merged_data))
  message("   Total columns: ", ncol(merged_data))

  # Step 7: Save output (always in processed_data)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message("Created output directory: ", out_dir)
  }

  out_file <- file.path(
    out_dir,
    paste0(
      tools::file_path_sans_ext(basename(file_name)),
      "_merged.tsv"
    )
  )

  fwrite(merged_data, out_file, sep = "\t")
  message("Step 7: Output saved as: ", out_file)

  message("Process completed successfully!")

  return(merged_data)
}

# ankita
# cecal
# Control
file <- "data/raw_data/ankita/001_cecal_ankita/control_cecal_ankita.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/001_cecal")
# Treatment
file <- "data/raw_data/ankita/001_cecal_ankita/treatment_cecal_ankita.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/001_cecal")
# serum
# Control
file <- "data/raw_data/ankita/002_serum_ankita/control_serum_ankita.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/002_serum")
# Treatment
file <- "data/raw_data/ankita/002_serum_ankita/treatment_serum_ankita.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/002_serum")

# Arko
# cecal
# Control
file <- "data/raw_data/arko/001_cecal_arko/control_cecal_arko.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/001_cecal")
# Treatment
file <- "data/raw_data/arko/001_cecal_arko/treatment_cecal_arko.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/001_cecal")
# serum
# Control
file <- "data/raw_data/arko/002_serum_arko/control_serum_arko.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/002_serum")
# Treatment
file <- "data/raw_data/arko/002_serum_arko/treatment_serum_arko.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/002_serum")


# Tanuja
# cecal
# Control
file <- "data/raw_data/tanuja/001_cecal_tanuja/control_cecal_tanuja.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/001_cecal")
# Treatment
file <- "data/raw_data/tanuja/001_cecal_tanuja/treatment_cecal_tanuja.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/001_cecal")
# serum
# Control
file <- "data/raw_data/tanuja/002_serum_tanuja/control_serum_tanuja.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/002_serum")
# Treatment
file <- "data/raw_data/tanuja/002_serum_tanuja/treatment_serum_tanuja.xlsx"
result <- merge_sheets(file, "METABOLITE","data/processed_data/002_serum")


# ankita
# species
# Control
file <- "data/raw_data/ankita/003_species_ankita/control_species_ankita.xlsx"
result <- merge_sheets(file, "Species","data/processed_data/003_species")
# Treatment
file <- "data/raw_data/ankita/003_species_ankita/treatment_species_ankita.xlsx"
result <- merge_sheets(file, "Species","data/processed_data/003_species")

# arko
# species
# Control
file <- "data/raw_data/arko/003_species_arko/control_species_arko.xlsx"
result <- merge_sheets(file, "Species","data/processed_data/003_species")
# Treatment
file <- "data/raw_data/arko/003_species_arko/treatment_species_arko.xlsx"
result <- merge_sheets(file, "Species","data/processed_data/003_species")

# Tanuja
# species
# Control
file <- "data/raw_data/tanuja/003_species_tanuja/control_species_tanuja.xlsx"
result <- merge_sheets(file, "Species","data/processed_data/003_species")
# Treatment
file <- "data/raw_data/tanuja/003_species_tanuja/treatment_species_tanuja.xlsx"
result <- merge_sheets(file, "Species","data/processed_data/003_species")


