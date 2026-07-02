variance_filter_top200 <- function(
  input_file,
  output_xlsx,
  output_dir = "data/processed_data/013_variance_filter",
  species_sheet = "Species_CLR",
  cecal_sheet = "Cecal_log10p",
  serum_sheet = "Serum_log10p",
  top_n = 200
) {
  suppressPackageStartupMessages({
    library(openxlsx)
    library(tibble)
  })

  #----------------------------------------------------#
  # CREATE OUTPUT DIRECTORY
  #----------------------------------------------------#
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  output_file <- file.path(output_dir, output_xlsx)

  #----------------------------------------------------#
  # READ TRANSFORMED DATA
  #----------------------------------------------------#
  message("Step 1: Reading transformed sheets")

  species <- read.xlsx(input_file, sheet = species_sheet, rowNames = TRUE,sep.names = " ")
  cecal <- read.xlsx(input_file, sheet = cecal_sheet, rowNames = TRUE)
  serum <- read.xlsx(input_file, sheet = serum_sheet, rowNames = TRUE)

  stopifnot(
    identical(rownames(species), rownames(cecal)),
    identical(rownames(species), rownames(serum))
  )

  message("✔ Sample IDs aligned")

  #----------------------------------------------------#
  # TOP FEATURES BY VARIANCE
  #----------------------------------------------------#
  top_by_variance <- function(df, n) {
    vars <- apply(df, 2, var, na.rm = TRUE)
    df[,
      order(vars, decreasing = TRUE)[seq_len(min(n, ncol(df)))],
      drop = FALSE
    ]
  }

  message("Step 2: Selecting top ", top_n, " variance features")

  cecal_top200 <- top_by_variance(cecal, top_n)
  serum_top200 <- top_by_variance(serum, top_n)

  #----------------------------------------------------#
  # ADD SAMPLE COLUMN (YOUR EXACT STYLE)
  #----------------------------------------------------#
  species_out <- rownames_to_column(
    as.data.frame(species),
    var = "Sample"
  )

  cecal_out <- rownames_to_column(
    as.data.frame(cecal_top200),
    var = "Sample"
  )

  serum_out <- rownames_to_column(
    as.data.frame(serum_top200),
    var = "Sample"
  )

  #----------------------------------------------------#
  # WRITE OUTPUT WORKBOOK
  #----------------------------------------------------#
  message("Step 3: Writing Excel output")

  wb <- createWorkbook()

  addWorksheet(wb, "Cecal")
  writeData(wb, "Cecal", cecal_out, colNames = TRUE)

  addWorksheet(wb, "Serum")
  writeData(wb, "Serum", serum_out, colNames = TRUE)

  addWorksheet(wb, "Species")
  writeData(wb, "Species", species_out, colNames = TRUE)

  saveWorkbook(wb, output_file, overwrite = TRUE)

  message(" Saved: ", output_file)

  invisible(list(
    cecal_top200 = cecal_out,
    serum_top200 = serum_out,
    species = species_out
  ))
}


variance_filter_top200(
  input_file = "data/processed_data/012_transformed/control_transformed.xlsx",
  output_xlsx = "control_top200.xlsx"
)

variance_filter_top200(
  input_file = "data/processed_data/012_transformed/treatment_transformed.xlsx",
  output_xlsx = "treatment_top200.xlsx"
)
