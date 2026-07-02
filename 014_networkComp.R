source("code/requirement.R")
#---------------------------
# Helper: Export node properties (Cytoscape-like)
#---------------------------
export_node_properties <- function(g, resultDir, network_name = "network") {
  message(" Calculating node-level properties for: ", network_name)

  node_data <- data.frame(
    Node = V(g)$name,
    Degree = degree(g),
    InDegree = if (is_directed(g)) degree(g, mode = "in") else NA,
    OutDegree = if (is_directed(g)) degree(g, mode = "out") else NA,
    Betweenness = betweenness(g, normalized = TRUE),
    Closeness = closeness(g, normalized = TRUE),
    Eigenvector = eigen_centrality(g)$vector,
    PageRank = page_rank(g)$vector,
    ClusteringCoeff = transitivity(g, type = "local", isolates = "zero")
  )

  # Include vertex attributes (color, shape, etc.)
  vertex_attrs <- list.vertex.attributes(g)
  for (attr in vertex_attrs) {
    node_data[[attr]] <- get.vertex.attribute(g, attr)
  }

  # Save node-level CSV
  out_path <- file.path(resultDir, paste0(network_name, "_node_properties.txt"))
  fwrite(node_data, out_path, sep = "\t")
  message(" Node properties saved to: ", out_path)

  return(node_data)
}

#---------------------------
# Helper: Extract global network metrics
#---------------------------
get_network_properties <- function(g, label) {
  message(" Calculating global network metrics for: ", label)
  props <- list(
    Network = label,
    Nodes = gorder(g),
    Edges = gsize(g),
    Density = edge_density(g),
    Diameter = diameter(g, directed = FALSE, weights = NA),
    Avg_Degree = mean(degree(g)),
    Avg_Clustering = transitivity(g, type = "average"),
    Avg_Path_Length = mean_distance(g, directed = FALSE),
    Modularity = modularity(cluster_louvain(g)),
    Degree_Centralization = centr_degree(g)$centralization,
    Betweenness_Centralization = centr_betw(g)$centralization
  )
  return(as.data.frame(props))
}

#---------------------------
# Helper: Compare two networks
#---------------------------
compare_networks <- function(g1, g2, label1, label2, resultDir) {
  message(" Comparing networks: ", label1, " vs ", label2)

  # Node overlap
  nodes1 <- V(g1)$name
  nodes2 <- V(g2)$name
  shared_nodes <- intersect(nodes1, nodes2)
  unique_nodes1 <- setdiff(nodes1, nodes2)
  unique_nodes2 <- setdiff(nodes2, nodes1)

  # Edge overlap
  e1 <- igraph::as_data_frame(g1, what = "edges")[, 1:2]
  e2 <- igraph::as_data_frame(g2, what = "edges")[, 1:2]
  e1 <- t(apply(e1, 1, sort))
  e2 <- t(apply(e2, 1, sort))
  e1_df <- as.data.frame(e1)
  e2_df <- as.data.frame(e2)

  shared_edges <- merge(e1_df, e2_df, by = c("V1", "V2"))
  unique_edges1 <- anti_join(e1_df, e2_df, by = c("V1", "V2"))
  unique_edges2 <- anti_join(e2_df, e1_df, by = c("V1", "V2"))

  # Jaccard similarities
  jaccard_edges <- nrow(shared_edges) / length(unique(rbind(e1, e2)))
  jaccard_nodes <- length(shared_nodes) / length(unique(c(nodes1, nodes2)))

  # Centrality correlation on shared nodes
  deg1 <- degree(g1)
  deg2 <- degree(g2)
  common <- intersect(names(deg1), names(deg2))
  deg_corr <- cor(deg1[common], deg2[common])

  betw1 <- betweenness(g1)
  betw2 <- betweenness(g2)
  betw_corr <- cor(betw1[common], betw2[common])

  # Save shared and unique nodes
  fwrite(
    data.frame(Node = shared_nodes),
    file.path(resultDir, "shared_nodes.txt"),
    sep = "\t"
  )

  fwrite(
    data.frame(Node = unique_nodes1),
    file.path(resultDir, paste0(label1, "_unique_nodes.txt")),
    sep = "\t"
  )

  fwrite(
    data.frame(Node = unique_nodes2),
    file.path(resultDir, paste0(label2, "_unique_nodes.txt")),
    sep = "\t"
  )

  # Save shared and unique edges
  fwrite(shared_edges, file.path(resultDir, "shared_edges.txt"), sep = "\t")
  fwrite(
    unique_edges1,
    file.path(resultDir, paste0(label1, "_unique_edges.txt")),
    sep = "\t"
  )
  fwrite(
    unique_edges2,
    file.path(resultDir, paste0(label2, "_unique_edges.txt")),
    sep = "\t"
  )

  # Summary
  summary_df <- data.frame(
    Metric = c(
      "Nodes (Network 1)",
      "Nodes (Network 2)",
      "Shared Nodes",
      "Shared Edges",
      "Jaccard Node Similarity",
      "Jaccard Edge Similarity",
      "Degree Correlation",
      "Betweenness Correlation"
    ),
    Value = c(
      length(nodes1),
      length(nodes2),
      length(shared_nodes),
      nrow(shared_edges),
      round(jaccard_nodes, 3),
      round(jaccard_edges, 3),
      round(deg_corr, 3),
      round(betw_corr, 3)
    )
  )

  fwrite(
    summary_df,
    file.path(resultDir, "network_comparison_summary.txt"),
    sep = "\t"
  )

  message(" Network comparison completed. Results saved in: ", resultDir)
  return(summary_df)
}


#---------------------------
# Helper: Plot shared and unique edges with preserved colors/shapes
#---------------------------
plot_shared_unique_edges <- function(
  g1,
  g2,
  label1,
  label2,
  resultDir,
  layout
) {
  message(" Plotting shared and unique edges...")

  # Convert edges to undirected pairs (sorted)
  e1 <- igraph::as_data_frame(g1, what = "edges")[, 1:2]
  e2 <- igraph::as_data_frame(g2, what = "edges")[, 1:2]
  e1 <- t(apply(e1, 1, sort))
  e2 <- t(apply(e2, 1, sort))
  e1_df <- as.data.frame(e1)
  e2_df <- as.data.frame(e2)

  # Identify shared and unique edges
  shared_edges <- merge(e1_df, e2_df, by = c("V1", "V2"))
  unique_edges1 <- anti_join(e1_df, e2_df, by = c("V1", "V2"))
  unique_edges2 <- anti_join(e2_df, e1_df, by = c("V1", "V2"))

  # Combine for plotting
  combined_edges <- bind_rows(
    shared_edges %>% mutate(Type = "Shared"),
    unique_edges1 %>% mutate(Type = paste0("Unique_", label1)),
    unique_edges2 %>% mutate(Type = paste0("Unique_", label2))
  )
  combined_nodes <- data.frame(
    name = unique(c(combined_edges$V1, combined_edges$V2))
  )

  # Create combined graph
  g_combined <- graph_from_data_frame(
    combined_edges,
    vertices = combined_nodes,
    directed = FALSE
  )

  # Transfer node colors/shapes
  V(g_combined)$color <- coalesce(
    V(g1)$color[match(V(g_combined)$name, V(g1)$name)],
    V(g2)$color[match(V(g_combined)$name, V(g2)$name)],
    "grey"
  )

  V(g_combined)$shape <- coalesce(
    V(g1)$shape[match(V(g_combined)$name, V(g1)$name)],
    V(g2)$shape[match(V(g_combined)$name, V(g2)$name)],
    "circle"
  )

  # Assign edge colors
  E(g_combined)$color <- case_when(
    E(g_combined)$Type == "Shared" ~ "darkgreen",
    E(g_combined)$Type == paste0("Unique_", label1) ~ "red",
    E(g_combined)$Type == paste0("Unique_", label2) ~ "blue"
  )
  shared_unique <- igraph::as_data_frame(g_combined, what = "edges")
  fwrite(
    shared_unique,
    file.path(resultDir, "shared_unique_network_comparison.txt"),
    sep = "\t"
  )
  # Plot combined network
  tiff(
    file.path(
      resultDir,
      paste0("Shared_Unique_", label1, "_vs_", label2, ".tiff")
    ),
    width = 3000,
    height = 3000,
    res = 300
  )
  plot(
    g_combined,
    layout = layout_with_drl,
    vertex.color = V(g_combined)$color,
    vertex.shape = V(g_combined)$shape,
    vertex.label = NA,
    vertex.size = 4,
    edge.width = 1,
    edge.color = E(g_combined)$color,
    main = paste("Shared (green) and Unique Edges\n", label1, "vs", label2)
  )
  legend(
    "topleft",
    legend = c("Shared", paste("Unique -", label1), paste("Unique -", label2)),
    col = c("darkgreen", "red", "blue"),
    lwd = 2,
    bty = "n"
  )
  dev.off()

  message(" Shared/Unique edge plot saved to: ", resultDir)
}

#---------------------------
# Helper: Detect and export clusters + clustered edges
#---------------------------
export_clusters <- function(g, resultDir, label) {
  message(" 🧩 Detecting clusters for: ", label)

  # Perform community detection (Louvain)
  community <- cluster_louvain(g)

  # Node-to-cluster mapping
  membership_df <- data.frame(
    Node = names(membership(community)),
    Cluster = membership(community)
  )

  # Save node-level cluster membership
  node_out_path <- file.path(resultDir, paste0(label, "_clusters.txt"))
  fwrite(membership_df, node_out_path, sep = "\t")

  # Cluster summary (number of nodes per cluster)
  cluster_summary <- data.frame(
    Cluster = as.numeric(names(table(membership_df$Cluster))),
    Size = as.numeric(table(membership_df$Cluster))
  )
  fwrite(
    cluster_summary,
    file.path(resultDir, paste0(label, "_cluster_summary.txt")),
    sep = "\t"
  )

  #---------------------------
  # Annotate edges with cluster info
  #---------------------------
  edges <- igraph::as_data_frame(g, what = "edges")
  edges$Cluster_Source <- membership_df$Cluster[
    match(edges$from, membership_df$Node)
  ]
  edges$Cluster_Target <- membership_df$Cluster[
    match(edges$to, membership_df$Node)
  ]
  edges$Interaction_Type <- ifelse(
    edges$Cluster_Source == edges$Cluster_Target,
    "Intra-cluster",
    "Inter-cluster"
  )

  # Save edge file annotated with cluster info
  edge_out_path <- file.path(resultDir, paste0(label, "_clustered_edges.txt"))
  fwrite(edges, edge_out_path, sep = "\t")

  # Summary message
  message(" ✅ Clusters detected and saved for ", label)
  message(" 📂 Clustered edges saved to: ", edge_out_path)

  return(list(community = community, clusters = membership_df, edges = edges))
}


#---------------------------
# Main Function
#---------------------------
mynetwork <- function(
  treatment,
  control,
  label_treatment,
  label_control,
  direct = FALSE,
  resultDir,
  layout
) {
  if (!dir.exists(resultDir)) {
    dir.create(resultDir, recursive = TRUE)
    message(" Created result directory: ", resultDir)
  }

  # Read edge files
  message(" Reading treatment and control edge files...")
  edges_treatment <- fread(treatment, header = TRUE, sep = "\t")
  edges_control <- fread(control, header = TRUE, sep = "\t")

  #-------------------------
  # 🧬 Treatment network
  #-------------------------
  message(" Building treatment network...")
  g_treatment <- graph_from_data_frame(edges_treatment, directed = direct)

  message(" Setting node colors and shapes for treatment network...")
  V(g_treatment)$source <- edges_treatment$SourceNode[match(
    V(g_treatment)$name,
    edges_treatment$SerumCecal
  )]

  V(g_treatment)$color <- ifelse(
    is.na(V(g_treatment)$source),
    "red",
    ifelse(
      V(g_treatment)$source == 1,
      "green",
      ifelse(V(g_treatment)$source == 2, "yellow", "red")
    )
  )

  V(g_treatment)$shape <- ifelse(
    is.na(V(g_treatment)$source),
    "circle",
    ifelse(
      V(g_treatment)$source == 1,
      "square",
      ifelse(V(g_treatment)$source == 2, "square", "circle")
    )
  )

  # Plot treatment network
  tiff(
    file.path(resultDir, paste0(label_treatment, "_network.tiff")),
    width = 3000,
    height = 3000,
    res = 300
  )
  plot(
    g_treatment,
    layout = layout,
    vertex.label = NA,
    vertex.size = 4,
    edge.width = 0.5,
    vertex.color = V(g_treatment)$color,
    vertex.shape = V(g_treatment)$shape
  )
  dev.off()

  export_node_properties(g_treatment, resultDir, label_treatment)
  treatment_props <- get_network_properties(g_treatment, label_treatment)
  community_treatment <- export_clusters(
    g_treatment,
    resultDir,
    label_treatment
  )
  #-------------------------
  # 🧪 Control network
  #-------------------------
  message(" Building control network...")
  g_control <- graph_from_data_frame(edges_control, directed = direct)

  V(g_control)$source <- edges_control$SourceNode[match(
    V(g_control)$name,
    edges_control$SerumCecal
  )]

  V(g_control)$color <- ifelse(
    is.na(V(g_control)$source),
    "red",
    ifelse(
      V(g_control)$source == 1,
      "green",
      ifelse(V(g_control)$source == 2, "yellow", "red")
    )
  )

  V(g_control)$shape <- ifelse(
    is.na(V(g_control)$source),
    "circle",
    ifelse(
      V(g_control)$source == 1,
      "square",
      ifelse(V(g_control)$source == 2, "square", "circle")
    )
  )

  # Plot control network
  tiff(
    file.path(resultDir, paste0(label_control, "_network.tiff")),
    width = 3000,
    height = 3000,
    res = 300
  )
  plot(
    g_control,
    layout = layout,
    vertex.label = NA,
    vertex.size = 4,
    edge.width = 0.5,
    vertex.color = V(g_control)$color,
    vertex.shape = V(g_control)$shape
  )
  dev.off()

  export_node_properties(g_control, resultDir, label_control)
  control_props <- get_network_properties(g_control, label_control)
  community_control <- export_clusters(g_control, resultDir, label_control)
  #-------------------------
  # 🧾 Combine and save summary
  #-------------------------
  summary_table <- rbind(treatment_props, control_props)
  summary_path <- file.path(resultDir, "network_summary_properties.txt")
  fwrite(summary_table, summary_path, sep = "\t")

  #-------------------------
  # 🔍 Compare treatment vs control
  #-------------------------
  comparison_summary <- compare_networks(
    g_treatment,
    g_control,
    label_treatment,
    label_control,
    resultDir
  )

  # fwrite(
  #   comparison_summary,
  #   file.path(resultDir, "network_prop_comparison.txt"),
  #   sep = "\t"
  # )

  plot_shared_unique_edges(
    g_treatment,
    g_control,
    label_treatment,
    label_control,
    resultDir,
    layout
  )

  message(" Final comparison summary saved!")

  message(" All plots, metrics, and comparisons saved in: ", resultDir)
}

#---------------------------
# Example Run
#---------------------------
inputDir <- "results/001_multi_omics"
treatment <- file.path(inputDir, "002_treatment/treatment_edges_Serum_Cecal_Species.txt"
)
control <- file.path(inputDir, "001_control/control_edges_Serum_Cecal_Species.txt")
resultDir <- "results/003_network_comp"
layout <- layout.sphere
mynetwork(
  treatment = treatment,
  control = control,
  label_treatment = "treatment",
  label_control = "control",
  resultDir = resultDir,
  layout = layout
)
