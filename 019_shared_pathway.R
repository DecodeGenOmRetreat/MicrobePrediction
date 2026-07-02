map_metabolite_microbe_pathways <- function(
  serum_cecal_file,
  cecal_species_file,
  species_pathway_file,
  output_dir = ".",
  study_name
) {
  # Load dependencies
  source("code/requirement.R") # should load data.table, dplyr, etc.
  source("code/requirement.R")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  message(" Reading input files...")
  message(" Serum-Cecal metabolite interaction file...")
  #  Serum-Cecal metabolite interaction file
  sc <- fread(serum_cecal_file, sep = "\t", header = TRUE)
  sc <- sc[, .(Serum = Serum, SerumId = SerumId, Cecal = Cecal)]
  message(" Cecal metabolite-species interaction file...")
  #  Cecal metabolite-species interaction file
  cs <- fread(cecal_species_file, sep = "\t", header = TRUE)
  cs <- cs[, .(Cecal = Cecal, CecalId = CecalId, Species = Species)]
  message(" Species-pathway mapping file...")
  #  Species-pathway mapping file
  sp <- fread(species_pathway_file, sep = "\t", header = TRUE)

  message(" Downloading KEGG compound → pathway mapping...")
  kegg_map <- read.delim(
    "https://rest.kegg.jp/link/pathway/compound",
    header = FALSE,
    sep = "\t"
  )
  colnames(kegg_map) <- c("CompoundId", "PathwayId")

  # Clean IDs
  kegg_map$CompoundId <- sub("cpd:", "", kegg_map$CompoundId)
  kegg_map$PathwayId <- sub("path:", "", kegg_map$PathwayId)

  message(" Mapping Serum metabolites to KEGG pathways...")
  metabolite_serum_mapped <- sc %>%
    inner_join(
      kegg_map,
      by = c("SerumId" = "CompoundId"),
      relationship = "many-to-many"
    ) %>%
    distinct()

  message(" Mapping Cecal metabolites to KEGG pathways...")
  metabolite_cecal_mapped <- cs %>%
    inner_join(
      kegg_map,
      by = c("CecalId" = "CompoundId"),
      relationship = "many-to-many"
    ) %>%
    distinct()

  message(" Finding shared pathways between Serum and Cecal metabolites...")
  metabolite_merged <- inner_join(
    metabolite_serum_mapped,
    metabolite_cecal_mapped,
    by = c("Cecal", "PathwayId"),
    relationship = "many-to-many"
  ) %>%
    distinct()

  message(
    "Mapping shared metabolite pathways to microbial species pathways..."
  )
  metabolite_species_merged <- metabolite_merged %>%
    inner_join(
      sp,
      by = c("Species", "PathwayId"),
      relationship = "many-to-many"
    ) %>%
    distinct()

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Define output paths
  #out1 <- file.path(output_dir, paste0(study_name, "_metabolite_serum_pathway.tsv"))
  #out2 <- file.path(output_dir, paste0(study_name, "_metabolite_cecal_pathway.tsv"))
  #out3 <- file.path(output_dir, paste0(study_name, "_shared_metabolite_pathway.tsv"))
  out4 <- file.path(
    output_dir,
    paste0(study_name, "_metabolite_species_shared_pathway.tsv")
  )

  message(" Writing output files...")
  #fwrite(metabolite_serum_mapped, out1, sep = "\t")
  #fwrite(metabolite_cecal_mapped, out2, sep = "\t")
  #fwrite(metabolite_merged, out3, sep = "\t")
  fwrite(metabolite_species_merged, out4, sep = "\t")

  message(" Pathway mapping completed successfully!")
  message(
    "Total shared pathways (metabolite ↔ species): ",
    nrow(metabolite_species_merged)
  )

  #return(list(
  #   Serum_Pathways = metabolite_serum_mapped,
  #   Cecal_Pathways = metabolite_cecal_mapped,
  #   Shared_Metabolite_Pathways = metabolite_merged,
  #   Shared_Metabolite_Species_Pathways = metabolite_species_merged
  # ))
}
# Control
input_dir <- "results/001_multi_omics/001_control"
serum_cecal_file <- file.path(
  input_dir,
  "control_edges_Serum_Cecal_id_edited.txt"
)
cecal_species_file <- file.path(
  input_dir,
  "control_edges_Cecal_Species_id_edited.txt"
)
pathway_dir <- "data/processed_data/015_pathway"
species_pathway_file <- file.path(
  pathway_dir,
  "control_species_pathway_edited.tsv"
)
map_metabolite_microbe_pathways(
  serum_cecal_file = serum_cecal_file,
  cecal_species_file = cecal_species_file,
  species_pathway_file = species_pathway_file,
  output_dir = "results/004_pathway",
  study_name = "control"
)
#Treatment
input_dir <- "results/001_multi_omics/002_treatment"
serum_cecal_file <- file.path(input_dir, "treatment_edges_Serum_Cecal_id_edited.txt")
cecal_species_file <- file.path(
  input_dir,
  "treatment_edges_Cecal_Species_id_edited.txt"
)
pathway_dir <- "data/processed_data/015_pathway"
species_pathway_file <- file.path(
  pathway_dir,
  "treatment_species_pathway_edited.tsv"
)
map_metabolite_microbe_pathways(
  serum_cecal_file = serum_cecal_file,
  cecal_species_file = cecal_species_file,
  species_pathway_file = species_pathway_file,
  output_dir = "results/004_pathway",
  study_name = "treatment"
)
