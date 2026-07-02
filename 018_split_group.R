split_groups <- function(
  input_file,
  control_samples,
  treatment_samples,
  output_dir = "."
) {
  source("code/requirement.R")
   if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  message("Reading input file: ", input_file)
  data <- fread(input_file, header = TRUE)

  # Check if "sample" column exists
  if (!"sample" %in% colnames(data)) {
    stop("The input file must contain a column named 'sample'.")
  }

  message("Filtering data for control and treatment groups...")

  # Filter control and treatment
  control_data <- data[data$sample %in% control_samples, ]
  treatment_data <- data[data$sample %in% treatment_samples, ]

  # Select only desired columns
  control_selected <- control_data[, .(
    Taxon,
    PathwayId,
    PathwayName
  )]
  treatment_selected <- treatment_data[, .(
    Taxon,
    PathwayId,
    PathwayName
  )]

  # Remove duplicates
  message("🧹 Removing duplicates from both control and treatment datasets...")
  control_selected <- unique(control_selected)
  treatment_selected <- unique(treatment_selected)

  # Define output file paths
  control_file <- file.path(output_dir, "control_species_pathway.tsv")
  treatment_file <- file.path(output_dir, "treatment_species_pathway.tsv")

  message("Writing filtered files to: ", output_dir)
  fwrite(control_selected, control_file, sep = "\t")
  fwrite(treatment_selected, treatment_file, sep = "\t")

  message("Control data saved as: ", control_file)
  message("Treatment data saved as: ", treatment_file)
  message("Done!")
  # Optionally return list of data tables
  #return(list(control = control_data, treatment = treatment_data))
}

# Run
input_dir <- "data/raw_data/pathway"
input_file <- file.path(input_dir, "pathway_taxon_abundance.tsv")
output_dir <- "data/processed_data/015_pathway"

control_samples <- c(
  "CST1D3",
  "CST2D3",
  "CST3D3",
  "CST4D3",
  "CST5D3",
  "CST1D5",
  "CST2D5",
  "CST3D5",
  "CST4D5",
  "CST5D5",
  "CST1D7",
  "CST2D7",
  "CST3D7",
  "CST4D7",
  "CST5D7",
  "CST1D10",
  "CST2D10",
  "CST3D10",
  "CST4D10",
  "MC2",
  "2DC4",
  "JAC1",
  "JAC2",
  "DC2",
  "2DC2",
  "2DC3",
  "EC1",
  "MC3",
  "MC4",
  "2NC4",
  "2NC5",
  "MC6",
  "2NC1",
  "2NC2",
  "OC4",
  "S2C3",
  "JC4",
  "MC8",
  "2NC3",
  "OB3",
  "2DC1",
  "OC6",
  "OC8",
  "C1",
  "C2",
  "C3"
)
treatment_samples <- c(
  "MST1D3",
  "MST2D3",
  "MST3D3",
  "MST4D3",
  "MST5D3",
  "MST1D5",
  "MST2D5",
  "MST3D5",
  "MST4D5",
  "MST5D5",
  "MST4D7",
  "MST5D7",
  "MST1D10",
  "MST2D10",
  "MST3D10",
  "MST4D10",
  "MST5D10",
  "NF5",
  "M1F4",
  "JAF3",
  "JAF2",
  "NF1",
  "SEPF1",
  "2NF3",
  "2NF4",
  "S2F1",
  "S3F1",
  "SEPF3",
  "2DF4",
  "2DF1",
  "2NF1",
  "2NF2",
  "M1F2",
  "Z16F",
  "EF3",
  "NF4",
  "S3F3",
  "M1F1",
  "EF7",
  "T1",
  "T2",
  "T3"
)

split_groups(
  input_file = input_file,
  control_samples = control_samples,
  treatment_samples = treatment_samples,
  output_dir = output_dir
)
