annotate_and_stack_nodes <- function(
  input_file,
  species_annot_file,
  output_dir = "node_tables",
  study
) {
  library(data.table)
  library(dplyr)
  message("📂 Creating output directory...")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  #==================================================
  # STEP 1: Read node file
  #==================================================
  message("📖 Reading node file...")
  df <- fread(input_file)

  #==================================================
  # STEP 2: Assign SourceNode
  #==================================================
  message("🏷️ Assigning SourceNode labels...")
  df[,
    SourceNode := fifelse(
      grepl("^serum_", Node, ignore.case = TRUE),
      2L,
      fifelse(grepl("^cecal_", Node, ignore.case = TRUE), 1L, 3L)
    )
  ]

  #==================================================
  # STEP 3: Split node tables
  #==================================================
  message("✂️ Splitting nodes into serum / cecal / species...")
  serum_df <- df[SourceNode == 2]
  cecal_df <- df[SourceNode == 1]
  species_df <- df[SourceNode == 3]

  #==================================================
  # STEP 4: Create CleanName column
  #==================================================
  message("🧹 Creating CleanName column...")
  serum_df[, CleanName := gsub("^serum_", "", Node, ignore.case = TRUE)]
  cecal_df[, CleanName := gsub("^cecal_", "", Node, ignore.case = TRUE)]
  species_df[, CleanName := Node]

  #==================================================
  # STEP 5: Download KEGG annotations
  #==================================================
  message("🔗 Downloading KEGG Compound → Pathway mapping...")
  kegg_map <- fread(
    "https://rest.kegg.jp/link/pathway/compound",
    header = FALSE
  )
  setnames(kegg_map, c("CompoundId", "PathwayId"))
  kegg_map[, CompoundId := gsub("^cpd:", "", CompoundId)]
  kegg_map[, PathwayId := gsub("^path:", "", PathwayId)]

  message("🔗 Downloading KEGG Pathway names...")
  kegg_pathways <- fread("https://rest.kegg.jp/list/pathway", header = FALSE)
  setnames(kegg_pathways, c("PathwayId", "PathwayName"))
  kegg_pathways[, PathwayId := gsub("^path:", "", PathwayId)]

  kegg_annot <- kegg_map %>%
    left_join(kegg_pathways, by = "PathwayId")

  #==================================================
  # STEP 6: Annotate serum & cecal nodes
  #==================================================
  message("🧬 Annotating serum metabolites...")
  serum_annot <- serum_df %>%
    left_join(kegg_annot, by = c("CleanName" = "CompoundId"))

  message("🧬 Annotating cecal metabolites...")
  cecal_annot <- cecal_df %>%
    left_join(kegg_annot, by = c("CleanName" = "CompoundId"))

  #==================================================
  # STEP 7: Annotate species nodes
  #==================================================
  message("🦠 Annotating species nodes...")
  species_annot <- fread(species_annot_file)

  species_annotated <- species_df %>%
    left_join(species_annot, by = c("CleanName" = "Species"))

  #==================================================
  # STEP 8: Stack all nodes
  #==================================================
  all_nodes_stacked <- rbindlist(
    list(serum_annot, cecal_annot, species_annotated),
    use.names = TRUE,
    fill = TRUE
  )

  #==================================================
  # STEP 9: Save outputs
  #==================================================
  message("💾 Writing output files...")
  fwrite(
    all_nodes_stacked,
    file.path(output_dir, paste0(study, "_pathway.txt")),
    sep = "\t"
  )

  message("✅ Node annotation pipeline completed successfully")

  invisible(all_nodes_stacked)
}


annotate_and_stack_nodes(
  input_file = "results/003_network_comp/treatment_unique_nodes.txt",
  species_annot_file = "data/processed_data/015_pathway/treatment_species_pathway_edited.tsv",
  output_dir = "results/004_pathway",
  study = "treatment_unique_node"
)

annotate_and_stack_nodes(
  input_file = "results/003_network_comp/control_unique_nodes.txt",
  species_annot_file = "data/processed_data/015_pathway/control_species_pathway_edited.tsv",
  output_dir = "results/004_pathway",
  study = "control_unique_node"
)

annotate_and_stack_nodes(
  input_file = "results/003_network_comp/shared_nodes.txt",
  species_annot_file = "data/processed_data/015_pathway/treatment_species_pathway_edited.tsv",
  output_dir = "results/004_pathway",
  study = "shared_nodes"
)

