# =====================================================
# Convert PICRUSt2 KO output to KEGG Pathway abundances
# =====================================================

ko_to_pathway <- function(ko_abundance, resultPath = ".") {
  source("code/requirement.R")
  message("📂 Reading KO abundance file...")

  # 1️⃣ Read KO abundance table
  abundance <- fread(ko_abundance, sep = "\t", header = TRUE)

  # Check that 'function' column exists
  if (!"function" %in% colnames(abundance)) {
    stop("❌ Input file must contain a column named 'function' with KO IDs.")
  }

  message("🔗 Downloading KO → Pathway mapping from KEGG...")
  kegg_map <- fread(
    "https://rest.kegg.jp/link/pathway/ko",
    header = FALSE,
    sep = "\t"
  )
  colnames(kegg_map) <- c("KoId", "PathwayId")

  # Clean IDs
  kegg_map$KoId <- sub("ko:", "", kegg_map$KoId)
  kegg_map$PathwayId <- sub("path:", "", kegg_map$PathwayId)

  message("🧬 Mapping KO IDs to KEGG pathways...")
  ko_mapped <- abundance %>%
    inner_join(kegg_map, by = c("function" = "KoId"),relationship = "many-to-many)

  message("📖 Downloading KEGG pathway names...")
  pathway_names <- fread(
    "https://rest.kegg.jp/list/pathway",
    header = FALSE,
    sep = "\t"
  )
  colnames(pathway_names) <- c("PathwayId", "PathwayName")
  pathway_names$PathwayId <- sub("path:", "", pathway_names$PathwayId)

  # Join pathway names
  pathway_final <- ko_mapped %>%
    inner_join(pathway_names, by = "PathwayId",relationship = "many-to-many)

  # Define output path
  outFile <- file.path(resultPath, "pathway_taxon_abundance.tsv")

  message("💾 Writing combined KEGG pathway table to: ", outFile)
  fwrite(pathway_final, outFile, sep = "\t")

  message("✅ KEGG pathway mapping completed successfully!")
  message(
    "📊 Output contains ",
    nrow(pathway_final),
    " rows and ",
    ncol(pathway_final),
    " columns."
  )

  return(pathway_final)
}


ko_to_pathway(
  ko_abundance = "data/raw_data/picrust/picrust_ko_all_stacked.txt",
  resultPath = "data/raw_data/pathway"
)
