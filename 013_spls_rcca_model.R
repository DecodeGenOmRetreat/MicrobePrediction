# serum_microbe_interaction + train/test projection
#source("code/requirement.R")

serum_microbe_interaction <- function(
  input_file,
  species_sheet,
  cecal_sheet,
  serum_sheet,
  ncomp = 3,
  cutoff = 0.3,
  study,
  result_path,
  train_frac = 0.7,
  seed = 123
) {
  set.seed(seed)

  library(mixOmics)
  library(openxlsx)
  library(data.table)
  library(dplyr)
  library(igraph)

  if (!dir.exists(result_path)) {
    dir.create(result_path, recursive = TRUE)
  }

  message("Step 1: Reading data")

  species <- read.xlsx(input_file, sheet = species_sheet, rowNames = TRUE,sep.names = " ")
  cecal <- read.xlsx(input_file, sheet = cecal_sheet, rowNames = TRUE)
  colnames(cecal) <- paste0("cecal_", colnames(cecal))
  serum <- read.xlsx(input_file, sheet = serum_sheet, rowNames = TRUE)
  colnames(serum) <- paste0("serum_", colnames(serum))

  #----------------------------------------------------#
  # ALIGN SAMPLES
  #----------------------------------------------------#
  common_samples <- Reduce(
    intersect,
    list(rownames(species), rownames(cecal), rownames(serum))
  )

  species <- species[common_samples, , drop = FALSE]
  cecal <- cecal[common_samples, , drop = FALSE]
  serum <- serum[common_samples, , drop = FALSE]

  #----------------------------------------------------#
  # TRAIN / TEST SPLIT
  #----------------------------------------------------#
  message("Step 2: Train / test split")

  n <- length(common_samples)
  train_idx <- sample(seq_len(n), floor(train_frac * n))

  species_train <- species[train_idx, , drop = FALSE]
  cecal_train <- cecal[train_idx, , drop = FALSE]
  serum_train <- serum[train_idx, , drop = FALSE]

  species_test <- species[-train_idx, , drop = FALSE]
  cecal_test <- cecal[-train_idx, , drop = FALSE]
  serum_test <- serum[-train_idx, , drop = FALSE]

  library(openxlsx)

  # Create workbook
  wb <- createWorkbook()

  # ---- Add sheets ----
  addWorksheet(wb, "species_train")
  addWorksheet(wb, "cecal_train")
  addWorksheet(wb, "serum_train")
  addWorksheet(wb, "species_test")
  addWorksheet(wb, "cecal_test")
  addWorksheet(wb, "serum_test")

  # ---- Write data ----
  writeData(wb, "species_train", species_train, rowNames = TRUE)
  writeData(wb, "cecal_train", cecal_train, rowNames = TRUE)
  writeData(wb, "serum_train", serum_train, rowNames = TRUE)

  writeData(wb, "species_test", species_test, rowNames = TRUE)
  writeData(wb, "cecal_test", cecal_test, rowNames = TRUE)
  writeData(wb, "serum_test", serum_test, rowNames = TRUE)

  # ---- Save workbook ----
  saveWorkbook(
    wb,
    file = file.path(result_path, paste0(study, "_train_test_split_data.xlsx")),
    overwrite = TRUE
  )

  #----------------------------------------------------#
  # sPLS (TRAIN)
  #----------------------------------------------------#
  message("Step 3: sPLS (cecal → species)")

  tune.res <- tune.spls(
    X = cecal_train,
    Y = species_train,
    ncomp = ncomp,
    test.keepX = c(50, 100, ncol(cecal_train)),
    test.keepY = c(20, 50, ncol(species_train)),
    validation = "loo",
    measure = "cor",
    scale = TRUE
  )

  spls_model <- spls(
    X = cecal_train,
    Y = species_train,
    ncomp = ncomp,
    keepX = tune.res$choice.keepX,
    keepY = tune.res$choice.keepY,
    scale = TRUE
  )

  saveRDS(spls_model, file.path(result_path, paste0(study, "_spls_model.rds")))

  perf.spls <- perf(
    spls_model,
    validation = "loo",
    progressBar = TRUE
  )

  avg_mse <- mean(perf.spls$measures$MSEP$values$value)
  message("Average MSE: ", avg_mse)
  #----------------------------------------------------#
  # sPLS TEST PROJECTION
  #----------------------------------------------------#
  message("Step 4: sPLS test projection")

  spls_pred <- predict(spls_model, newdata = cecal_test)
  pred_Y <- spls_pred$predict[,, 1, drop = FALSE]

  test_cor_spls <- sapply(
    seq_len(ncol(species_test)),
    function(j) cor(pred_Y[, j, 1], species_test[, j])
  )
  names(test_cor_spls) <- colnames(species_test)

  fwrite(
    data.frame(Species = names(test_cor_spls), Correlation = test_cor_spls),
    file.path(result_path, paste0(study, "_spls_test_correlation.txt")),
    sep = "\t"
  )

  #----------------------------------------------------#
  # sPLS NETWORK
  #----------------------------------------------------#
  jpeg(
    file.path(result_path, paste0(study, "_cecal_species_network.jpeg")),
    width = 4000,
    height = 4000,
    res = 600
  )

  net_cs <- network(spls_model, comp = 1, cutoff = cutoff)
  dev.off()

  edges_cs <- as.data.frame(igraph::as_data_frame(net_cs$gR, "edges"))
  colnames(edges_cs)[1:3] <- c("cecal", "species", "weight_cs")
  edges_cs$weight_cs_abs <- abs(edges_cs$weight_cs)
  edges_cs$cecalNode <- 1

  fwrite(
    edges_cs,
    file.path(result_path, paste0(study, "_edges_cecal_species.txt")),
    sep = "\t"
  )

  #----------------------------------------------------#
  # rCCA (TRAIN)
  #----------------------------------------------------#

  message("Step 5: rCCA (serum ↔ cecal)")

  # Step 5a: Tune regularization parameters (lambda1 and lambda2)
  message("  Tuning regularization parameters...")

  rcca_tune <- tune.rcc(
    X = serum_train,
    Y = cecal_train,
    grid1 = seq(0.001, 1, length.out = 10), # lambda range for X (serum)
    grid2 = seq(0.001, 1, length.out = 10), # lambda range for Y (cecal)
    validation = "loo" # leave-one-out cross-validation
  )

  # View optimal lambda values
  message("  Optimal lambda1 (serum): ", rcca_tune$opt.lambda1)
  message("  Optimal lambda2 (cecal): ", rcca_tune$opt.lambda2)

  # Step 5b: Run rCCA with tuned parameters
  rcca_model <- rcc(
    X = serum_train,
    Y = cecal_train,
    ncomp = ncomp,
    lambda1 = rcca_tune$opt.lambda1,
    lambda2 = rcca_tune$opt.lambda2
  )
  saveRDS(rcca_model, file.path(result_path, paste0(study, "_rcca_model.rds")))

  #----------------------------------------------------#
  # rCCA TEST PROJECTION
  #----------------------------------------------------#
  message("Step 6: rCCA test projection")

  X_test_scores <- as.matrix(serum_test) %*% rcca_model$loadings$X
  Y_test_scores <- as.matrix(cecal_test) %*% rcca_model$loadings$Y

  test_cor_rcca <- diag(cor(X_test_scores, Y_test_scores))
  names(test_cor_rcca) <- paste0("Comp", seq_along(test_cor_rcca))

  fwrite(
    data.frame(Component = names(test_cor_rcca), Correlation = test_cor_rcca),
    file.path(result_path, paste0(study, "_rcca_test_correlation.txt")),
    sep = "\t"
  )

  #----------------------------------------------------#
  # rCCA NETWORK
  #----------------------------------------------------#
  jpeg(
    file.path(result_path, paste0(study, "_serum_cecal_network.jpeg")),
    width = 4000,
    height = 4000,
    res = 600
  )

  net_sc <- network(rcca_model, comp = 1, cutoff = cutoff)
  dev.off()

  edges_sc <- as.data.frame(igraph::as_data_frame(net_sc$gR, "edges"))
  colnames(edges_sc)[1:3] <- c("serum", "cecal", "weight_sc")
  edges_sc$weight_sc_abs <- abs(edges_sc$weight_sc)
  edges_sc$serumNode <- 2

  fwrite(
    edges_sc,
    file.path(result_path, paste0(study, "_edges_serum_cecal.txt")),
    sep = "\t"
  )

  #----------------------------------------------------#
  # MULTI-LAYER INTEGRATION
  #----------------------------------------------------#
  message("Step 7: Multi-layer integration")

  integrated_edges <- inner_join(
    edges_sc,
    edges_cs,
    by = "cecal",
    relationship = "many-to-many"
  )

  serum_unique <- integrated_edges |>
    dplyr::select(serum, cecal, weight_sc, weight_sc_abs, serumNode) |>
    distinct()

  cecal_unique <- integrated_edges |>
    dplyr::select(cecal, species, weight_cs, weight_cs_abs, cecalNode) |>
    distinct()

  combined_edges <- rbind(as.matrix(serum_unique), as.matrix(cecal_unique))

  colnames(combined_edges) <- c(
    "SerumCecal",
    "CecalSpecies",
    "Weight",
    "AbsWeight",
    "SourceNode"
  )

  fwrite(
    as.data.frame(combined_edges),
    file.path(
      result_path,
      paste0(
        study,
        "_edges_",
        serum_sheet,
        "_",
        cecal_sheet,
        "_",
        species_sheet,
        ".txt"
      )
    ),
    sep = "\t"
  )
  message("✔ Analysis complete")

  invisible(list(
    spls_model = spls_model,
    rcca_model = rcca_model,
    spls_test_cor = test_cor_spls,
    rcca_test_cor = test_cor_rcca,
    integrated_edges = integrated_edges
  ))
}


serum_microbe_interaction(
  input_file = "data/processed_data/013_variance_filter/control_top200.xlsx",
  species_sheet = "Species",
  cecal_sheet = "Cecal",
  serum_sheet = "Serum",
  ncomp = 3,
  cutoff = 0.3,
  study = "control",
  result_path = "results/001_multi_omics/001_control",
  train_frac = 0.8,
  seed = 123
)

serum_microbe_interaction(
  input_file = "data/processed_data/013_variance_filter/treatment_top200.xlsx",
  species_sheet = "Species",
  cecal_sheet = "Cecal",
  serum_sheet = "Serum",
  ncomp = 3,
  cutoff = 0.3,
  study = "treatment",
  result_path = "results/001_multi_omics/002_treatment",
  train_frac = 0.8,
  seed = 123
)
