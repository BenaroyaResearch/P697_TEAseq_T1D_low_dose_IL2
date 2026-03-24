# =============================================================================
# P697_functions.R
# Shared functions for P650 analysis pipeline
# =============================================================================
# This file contains all helper functions used across the P697 analysis scripts.
# Source this file at the beginning of each script after loading libraries.
# =============================================================================

# Helper function to order genes by their peak expression group for dotplots
orderGenesByExpression <- function(seuratObj, genes, groupBy) {
  # Get average expression per gene per group
  avgExp <- AverageExpression(
    seuratObj,
    features = genes,
    group.by = groupBy,
    assays = "RNA",
    slot = "data"
  )$RNA

  # For each gene, find which group has the highest expression
  peakGroup <- apply(avgExp, 1, which.max)
  peakValue <- apply(avgExp, 1, max)

  # Create ordering dataframe
  orderDf <- data.frame(
    gene = rownames(avgExp),
    peakGroup = peakGroup,
    peakValue = peakValue,
    stringsAsFactors = FALSE
  )

  # Order genes by peak group, then by peak value within group (descending)
  orderDf <- orderDf %>%
    dplyr::arrange(peakGroup, dplyr::desc(peakValue))

  return(orderDf$gene)
}

# Helper function to create styled dotplot for a gene list
createGeneDotPlot <- function(seuratObj, genes, groupBy, title, xlab = "Genes",
                              ylab = groupBy, plotDir = NULL, filename = NULL,
                              orderGenes = TRUE, dotScale = 3.5,
                              height = 5.5, width = 12, units = "in",
                              dpi = 600, formats = c("pdf", "png")) {
  # Filter to genes present in the dataset
  genes <- genes[!is.na(genes) & genes != ""]
  genesPresent <- genes[genes %in% rownames(seuratObj)]

  if (length(genesPresent) == 0) {
    warning(paste0("No genes found in dataset for: ", title))
    return(NULL)
  }

  cat(paste0("Creating dotplot: ", title, " (", length(genesPresent), " genes)\n"))

  # Order genes by expression pattern across groups if requested
  if (orderGenes && length(genesPresent) > 1) {
    genesOrdered <- orderGenesByExpression(seuratObj, genesPresent, groupBy)
  } else {
    genesOrdered <- genesPresent
  }

  dotPlot <- DotPlot(
    object = seuratObj,
    features = genesOrdered,
    group.by = groupBy,
    dot.scale = dotScale
  ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
      text = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9)
    ) +
    labs(
      y = ylab,
      x = xlab,
      title = title
    ) +
    guides(
      color = guide_colorbar(title = "Average\nExpression\nz-score"),
      size = guide_legend(title = "Percent\nExpressed")
    ) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red")

  # Save plot if plotDir and filename are provided
  if (!is.null(plotDir) && !is.null(filename)) {
    # Calculate appropriate width based on number of genes if not provided
    if (is.null(width)) {
      width <- max(8, 4 + length(genesPresent) * 0.25)
    }

    savePlot(
      plot = dotPlot,
      plotDir = plotDir,
      filename = filename,
      height = height,
      width = width,
      units = units,
      dpi = dpi,
      formats = formats
    )
  }

  return(dotPlot)
}

# RNA cluster optimization
optimizeRNAseqClustering <- function(seuratObject,
                                     hyperParamGridPCA.df,
                                     min.dist.vector,
                                     n_neighbors.vector,
                                     reduction = "pca",
                                     integrated = FALSE,
                                     useParallel = FALSE,
                                     computeSilhouette = TRUE,
                                     silhouetteSample = NULL,
                                     maxNClusters = 10,
                                     minNClusters = 2) {
  # ------------------------------------------------INPUTS------------------------------------------------
  # seuratObject is... a Seurat object
  # hyperParamGridPCA.df is an expanded grid for a full grid search for PCA clustering. type: dataframe
  # min.dist.vector is a vector of min.dist values to pass to UMAP (via scDEED())
  # n_neighbors.vector is a vector of n_neighbors values to pass to UMAP (via scDEED())
  # maxNClusters is the maximum number of clusters to be considered (will evaluate 2 through maxNClusters)
  #-------------------------------------------------------------------------------------------------------

  #----------------------------------------------OUTPUTS--------------------------------------------------
  # result contains $hyperParamGridOptClusterAndUMAP.df and $hyperParamGridClusterOnly.df
  #   Extract like:
  #     hyperParamGridOptClusterAndUMAP.df <- result$hyperParamGridOptClusterAndUMAP.df (smaller output with opt PCA clustering and UMAP embedding params)
  #     hyperParamGridClusterOnly.df <- result$hyperParamGridClusterOnly.df (larger output with full PCA clustering gridsearch and results)
  #--------------------------------------------------------------------------------------------------------

  #-----------------------------------------REQUIRED PACKAGES----------------------------------------------
  # Seurat, parallel, doParallel, foreach, clusterSim, cluster, dplyr, scDEED, stringr, pracma
  #--------------------------------------------------------------------------------------------------------

  # STEP 1: grid search over FindNeighbors() and FindClusters() hyperparameters to find
  # the best combinations for each nClusters obtained

  # Extract embedding matrix outside of the loop(s)
  fullEmbeddingMatrix <- Embeddings(seuratObject, reduction = reduction)

  if (useParallel) {
    # Existing parallel path (kept for backwards compatibility); note: higher memory footprint
    nCores <- parallel::detectCores() - 2 # leave 2 cores unused to reduce risk of RAM overuse.
    myCluster <- parallel::makeCluster(nCores, type = "FORK")
    doParallel::registerDoParallel(cl = myCluster)
    clusteringMetrics <- foreach(
      findNeighborsDim = hyperParamGridPCA.df$findNeighborsDim,
      findClustersRes = hyperParamGridPCA.df$findClustersRes,
      .combine = "cbind"
    ) %dopar% {
      set.seed(6022)
      obj <- FindNeighbors(seuratObject, reduction = reduction, dims = 1:findNeighborsDim, verbose = FALSE)
      obj <- FindClusters(obj, resolution = findClustersRes, verbose = FALSE)
      metaData.df <- obj@meta.data
      clusterMat <- fullEmbeddingMatrix[, 1:findNeighborsDim, drop = FALSE]
      clusterVect <- as.numeric(metaData.df$seurat_clusters)
      nClusters <- dim(table(metaData.df$seurat_clusters))
      daviesBouldinIdx <- index.DB(clusterMat, clusterVect)
      calinskiHarabasz <- index.G1(clusterMat, clusterVect)
      medianSilScore <- NaN
      if (computeSilhouette) {
        silMat <- clusterMat
        if (!is.null(silhouetteSample) && silhouetteSample < nrow(silMat)) {
          sampleIdx <- sample.int(nrow(silMat), silhouetteSample)
          silMat <- silMat[sampleIdx, , drop = FALSE]
          silVec <- clusterVect[sampleIdx]
        } else {
          silVec <- clusterVect
        }
        silScore <- silhouette(silVec, dist(silMat))
        if (!anyNA(silScore)) {
          silScoreSummary <- summary(silScore)
          medianSilScore <- unname(silScoreSummary$si.summary["Median"])
        }
      }
      out <- c(
        nClusters = nClusters,
        DBIndex = daviesBouldinIdx$DB,
        CHMetric = calinskiHarabasz,
        medianSilScore = medianSilScore
      )
      rm(obj, metaData.df, clusterMat, clusterVect)
      gc(verbose = FALSE)
      out
    }
    clusteringMetrics <- t(clusteringMetrics)
    hyperParamGridPCA.df$DBIndex <- clusteringMetrics[, "DBIndex"]
    hyperParamGridPCA.df$CHMetric <- clusteringMetrics[, "CHMetric"]
    hyperParamGridPCA.df$medianSilScore <- clusteringMetrics[, "medianSilScore"]
    hyperParamGridPCA.df$nClusters <- clusteringMetrics[, "nClusters"]
    parallel::stopCluster(cl = myCluster)
  } else {
    # Memory-optimized sequential path with neighbor graph reuse
    nRows <- nrow(hyperParamGridPCA.df)
    DBIndexVec <- numeric(nRows)
    CHVec <- numeric(nRows)
    SilVec <- numeric(nRows)
    nClusVec <- numeric(nRows)
    lastDim <- NA_integer_
    for (rowIdx in seq_len(nRows)) {
      findNeighborsDim <- hyperParamGridPCA.df$findNeighborsDim[rowIdx]
      findClustersRes <- hyperParamGridPCA.df$findClustersRes[rowIdx]
      # Recompute neighbor graph only when dimension changes
      if (is.na(lastDim) || findNeighborsDim != lastDim) {
        seuratObject <- FindNeighbors(seuratObject, reduction = reduction, dims = 1:findNeighborsDim, verbose = FALSE)
        lastDim <- findNeighborsDim
      }
      seuratObject <- FindClusters(seuratObject, resolution = findClustersRes, verbose = FALSE)
      metaData.df <- seuratObject@meta.data
      clusterMat <- fullEmbeddingMatrix[, 1:findNeighborsDim, drop = FALSE]
      clusterVect <- as.numeric(metaData.df$seurat_clusters)
      nClusVec[rowIdx] <- dim(table(metaData.df$seurat_clusters))
      daviesBouldinIdx <- index.DB(clusterMat, clusterVect)
      calinskiHarabasz <- index.G1(clusterMat, clusterVect)
      DBIndexVec[rowIdx] <- daviesBouldinIdx$DB
      CHVec[rowIdx] <- calinskiHarabasz
      SilVec[rowIdx] <- NaN
      if (computeSilhouette) {
        silMat <- clusterMat
        if (!is.null(silhouetteSample) && silhouetteSample < nrow(silMat)) {
          sampleIdx <- sample.int(nrow(silMat), silhouetteSample)
          silMat <- silMat[sampleIdx, , drop = FALSE]
          silVec <- clusterVect[sampleIdx]
        } else {
          silVec <- clusterVect
        }
        silScore <- silhouette(silVec, dist(silMat))
        if (!anyNA(silScore)) {
          silScoreSummary <- summary(silScore)
          SilVec[rowIdx] <- unname(silScoreSummary$si.summary["Median"])
        }
      }
      gc(verbose = FALSE)
    }
    hyperParamGridPCA.df$DBIndex <- DBIndexVec
    hyperParamGridPCA.df$CHMetric <- CHVec
    hyperParamGridPCA.df$medianSilScore <- SilVec
    hyperParamGridPCA.df$nClusters <- nClusVec
  }

  # scale the clustering metrics to fall between 0 and 1
  hyperParamGridPCA.df$DBIndex <- (hyperParamGridPCA.df$DBIndex - min(hyperParamGridPCA.df$DBIndex, na.rm = TRUE)) /
    (max(hyperParamGridPCA.df$DBIndex, na.rm = TRUE) - min(hyperParamGridPCA.df$DBIndex, na.rm = TRUE))

  hyperParamGridPCA.df$CHMetric <- (hyperParamGridPCA.df$CHMetric - min(hyperParamGridPCA.df$CHMetric, na.rm = TRUE)) /
    (max(hyperParamGridPCA.df$CHMetric, na.rm = TRUE) - min(hyperParamGridPCA.df$CHMetric, na.rm = TRUE))

  hyperParamGridPCA.df$medianSilScore <- (hyperParamGridPCA.df$medianSilScore - min(hyperParamGridPCA.df$medianSilScore, na.rm = TRUE)) /
    (max(hyperParamGridPCA.df$medianSilScore, na.rm = TRUE) - min(hyperParamGridPCA.df$medianSilScore, na.rm = TRUE))

  # make an inverse DB score so that it is more easily compared to the other two metrics
  # (small DB is best, but large CH and silscore are best)
  hyperParamGridPCA.df$DBIndexInverse <- 1 - hyperParamGridPCA.df$DBIndex

  # create hyperParamGridPCA.df$meanEvalMetric, as the mean of $DBIndex, $CHMetric, $medianSilScore
  hyperParamGridPCA.df$meanEvalMetric <- rowMeans(hyperParamGridPCA.df[, c("DBIndexInverse", "CHMetric", "medianSilScore")], na.rm = TRUE)

  # sort and export (optional) the hyperparamgrid
  hyperParamGridPCA.df <- hyperParamGridPCA.df %>%
    arrange(nClusters, desc(meanEvalMetric))

  # create a subset of hyperParamGridPCA.df with the best option for each nClusters
  hyperParamGridOptSubset.df <- hyperParamGridPCA.df %>%
    group_by(nClusters) %>%
    filter(meanEvalMetric == max(meanEvalMetric)) %>%
    ungroup()

  # Pre-allocate min.dist and n_neighbors columns with NA values
  hyperParamGridOptSubset.df$min.dist <- NA
  hyperParamGridOptSubset.df$n_neighbors <- NA

  permutedObject <- NULL

  if (isTRUE(integrated)) {
    baseEmbedding <- Embeddings(seuratObject, reduction = reduction)
    if (is.null(baseEmbedding)) {
      stop(sprintf("Reduction '%s' not found on supplied Seurat object", reduction))
    }

    permutedEmbedding <- baseEmbedding
    set.seed(314)
    for (colIdx in seq_len(ncol(permutedEmbedding))) {
      rowOrder <- pracma::randperm(nrow(permutedEmbedding))
      permutedEmbedding[, colIdx] <- baseEmbedding[rowOrder, colIdx]
    }

    permutedObject <- seuratObject
    reductionKey <- tryCatch(Seurat::Key(seuratObject[[reduction]]), error = function(...) NULL)
    if (is.null(reductionKey)) {
      reductionKey <- paste0(reduction, "_")
    }
    permutedObject[[reduction]] <- CreateDimReducObject(
      embeddings = permutedEmbedding,
      key = reductionKey,
      assay = DefaultAssay(permutedObject)
    )
  }

  # STEP 2: use scDEED to optimize the UMAP hyperparameters
  # ---------------------------start of scDEED UMAP opt---------------------------------------
  # if nrow(hyperParamGridOptSubset.df) is smaller than maxNClusters, then set maxNClusters to nrow(hyperParamGridOptSubset.df)
  if (nrow(hyperParamGridOptSubset.df) < maxNClusters) {
    maxNClusters <- nrow(hyperParamGridOptSubset.df)
  }

  # Loop over the rows of hyperParamGridOptSubset.df, finding the best UMAP embedding for each nClusters
  for (i in minNClusters:maxNClusters) {
    # Extract necessary values from hyperParamGridOptSubset.df
    nClusters <- hyperParamGridOptSubset.df$nClusters[i]
    numPCs <- hyperParamGridOptSubset.df$findNeighborsDim[i]
    findClustersRes <- hyperParamGridOptSubset.df$findClustersRes[i]

    set.seed(6022)
    seuratObject <- FindNeighbors(seuratObject, reduction = reduction, dims = 1:numPCs)
    seuratObject <- FindClusters(seuratObject, resolution = findClustersRes)

    if (isTRUE(integrated)) {
      scDEEDResult <- scDEED(seuratObject, # input Seurat object (must have UMAP or t-SNE already run)
        K = numPCs, # number of PCs
        reduction.method = "umap", # 'umap' or 'tsne'
        min.dist =  min.dist.vector, # scDEED default is 0.1 & 0.4; Seurat default is 0.3
        n_neighbors = n_neighbors.vector, # scDEED defaults are c(5, 20, 30, 40, 50). Seurat default is 30.
        pre_embedding = reduction,
        permuted = permutedObject,
        similarity_percent = 0.5, # scDEED default
        dubious_cutoff = 0.05, # scDEED default
        trustworthy_cutoff = 0.95) # scDEED default
    } else {
      scDEEDResult <- scDEED(seuratObject, # input Seurat object (must have UMAP or t-SNE already run)
        K = numPCs, # number of PCs
        reduction.method = "umap", # 'umap' or 'tsne'
        min.dist =  min.dist.vector, # scDEED default is 0.1 & 0.4; Seurat default is 0.3
        n_neighbors = n_neighbors.vector, # scDEED defaults are c(5, 20, 30, 40, 50). Seurat default is 30.
        similarity_percent = 0.5, # scDEED default
        dubious_cutoff = 0.05, # scDEED default
        trustworthy_cutoff = 0.95) # scDEED default
    }

    # Extract optimal UMAP hyperparameters (minimizing number of dubious cells)
    if (!is.null(scDEEDResult) && !is.null(scDEEDResult$num_dubious)) {
      nd <- scDEEDResult$num_dubious
      if (all(c("number_dubious_cells", "min.dist", "n_neighbors") %in% colnames(nd))) {
        optIdx <- which(nd$number_dubious_cells == min(nd$number_dubious_cells, na.rm = TRUE))[1]
        if (!is.na(optIdx)) {
          hyperParamGridOptSubset.df$min.dist[i] <- nd$min.dist[optIdx]
          hyperParamGridOptSubset.df$n_neighbors[i] <- nd$n_neighbors[optIdx]
        }
      }
    }
  }

  # Reorder columns to place min.dist and n_neighbors after findClustersRes
  hyperParamGridOptSubset.df <- hyperParamGridOptSubset.df %>%
    dplyr::select(findNeighborsDim, findClustersRes, min.dist, n_neighbors, everything())

  # Return the hyperparam dataframes
  return(list(hyperParamGridOptClusterAndUMAP.df = hyperParamGridOptSubset.df, hyperParamGridClusterOnly.df = hyperParamGridPCA.df))
}

#' Run scDEED UMAP Optimization Benchmark
#'
#' Isolated scDEED call for benchmarking UMAP hyperparameter optimization.
#' Returns timing information and full scDEED results including dubious cell counts.
#'
#' @param seuratObject A Seurat object with PCA/harmony reduction already computed
#' @param numPCs Number of principal components to use for FindNeighbors
#' @param findClustersRes Resolution for FindClusters
#' @param min.dist.vector Vector of min.dist values to test
#' @param n_neighbors.vector Vector of n_neighbors values to test
#' @param reduction Name of the reduction to use (e.g., "pca", "harmony")
#' @param integrated Logical; if TRUE, creates permuted embedding for scDEED
#' @param seed Random seed for reproducibility
#' @return A list containing elapsed_seconds, optimal_min.dist, optimal_n_neighbors,
#'         and the full num_dubious table from scDEED
runScDEEDBenchmark <- function(
    seuratObject,
    numPCs,
    findClustersRes,
    min.dist.vector = c(0.1, 0.3, 0.5),
    n_neighbors.vector = c(10, 30),
    reduction = "harmony",
    integrated = TRUE,
    seed = 6022
    ) {
  # Prepare permuted object if using integrated data
  permutedObject <- NULL

  if (isTRUE(integrated)) {
    baseEmbedding <- Embeddings(seuratObject, reduction = reduction)
    if (is.null(baseEmbedding)) {
      stop(sprintf("Reduction '%s' not found on supplied Seurat object", reduction))
    }

    permutedEmbedding <- baseEmbedding
    set.seed(314)
    for (colIdx in seq_len(ncol(permutedEmbedding))) {
      rowOrder <- pracma::randperm(nrow(permutedEmbedding))
      permutedEmbedding[, colIdx] <- baseEmbedding[rowOrder, colIdx]
    }

    permutedObject <- seuratObject
    reductionKey <- tryCatch(Seurat::Key(seuratObject[[reduction]]), error = function(...) NULL)
    if (is.null(reductionKey)) {
      reductionKey <- paste0(reduction, "_")
    }
    permutedObject[[reduction]] <- CreateDimReducObject(
      embeddings = permutedEmbedding,
      key = reductionKey,
      assay = DefaultAssay(permutedObject)
    )
  }

  # Run FindNeighbors and FindClusters
  set.seed(seed)
  seuratObject <- FindNeighbors(seuratObject, reduction = reduction, dims = 1:numPCs, verbose = FALSE)
  seuratObject <- FindClusters(seuratObject, resolution = findClustersRes, verbose = FALSE)

  # Time the scDEED call
  startTime <- Sys.time()

  if (isTRUE(integrated)) {
    scDEEDResult <- scDEED(seuratObject,
      K = numPCs,
      reduction.method = "umap",
      min.dist = min.dist.vector,
      n_neighbors = n_neighbors.vector,
      pre_embedding = reduction,
      permuted = permutedObject,
      similarity_percent = 0.5,
      dubious_cutoff = 0.05,
      trustworthy_cutoff = 0.95)
  } else {
    scDEEDResult <- scDEED(seuratObject,
      K = numPCs,
      reduction.method = "umap",
      min.dist = min.dist.vector,
      n_neighbors = n_neighbors.vector,
      similarity_percent = 0.5,
      dubious_cutoff = 0.05,
      trustworthy_cutoff = 0.95)
  }

  endTime <- Sys.time()
  elapsedSeconds <- as.numeric(difftime(endTime, startTime, units = "secs"))

  # Extract optimal parameters
  optimalMinDist <- NA
  optimalNNeighbors <- NA
  numDubiousTable <- NULL

  if (!is.null(scDEEDResult) && !is.null(scDEEDResult$num_dubious)) {
    nd <- scDEEDResult$num_dubious
    numDubiousTable <- nd

    if (all(c("number_dubious_cells", "min.dist", "n_neighbors") %in% colnames(nd))) {
      optIdx <- which(nd$number_dubious_cells == min(nd$number_dubious_cells, na.rm = TRUE))[1]
      if (!is.na(optIdx)) {
        optimalMinDist <- nd$min.dist[optIdx]
        optimalNNeighbors <- nd$n_neighbors[optIdx]
      }
    }
  }

  return(list(
    elapsed_seconds = elapsedSeconds,
    optimal_min.dist = optimalMinDist,
    optimal_n_neighbors = optimalNNeighbors,
    num_dubious = numDubiousTable
  ))
}

# function for saving plots as both pdf and png
savePlot <- function(
    plot,
    plotDir,
    filename,
    height,
    width,
    units = "in",
    dpi = 600,
    formats = c("pdf", "png")
    ) {
  # Ensure plotDir exists
  if (!dir.exists(plotDir)) dir.create(plotDir, recursive = TRUE)

  # Save as PDF
  if ("pdf" %in% formats) {
    pdf(file.path(plotDir, paste0(filenameSuffix, "_", filename, ".pdf")),
      height = height,
      width = width)
    # The on.exit() calls prevent the graphics device from staying open if the plot fails to print
    on.exit(dev.off(), add = TRUE)
    print(plot)
    dev.off()
    on.exit(NULL)
  }

  # Save as PNG
  if ("png" %in% formats) {
    png(
      file.path(plotDir, paste0(filenameSuffix, "_", filename, ".png")),
      height = height,
      width = width,
      units = units,
      res = dpi
    )
    # The on.exit() calls prevent the graphics device from staying open if the plot fails to print
    on.exit(dev.off(), add = TRUE)
    print(plot)
    dev.off()
    on.exit(NULL)
  }
}

# Remove cells with: a) 3+ alphas, b) cells with 2 alphas and 2 betas, and c) 2 betas
# Remove iNKT and MAIT cells
callMultiplets <- function(tcrs,
                           nAlphaCut = 3,
                           nBetaCut = 2,
                           alphaAndBetaCut = 2,
                           callINKT = TRUE,
                           callMAIT = TRUE) {
  # Count chains
  chainCounts <- tcrs %>%
    dplyr::group_by(barcode) %>%
    dplyr::summarise(nAlpha = sum(chain ==  "TRA"),
      nBeta = sum(chain == "TRB"))

  multiplets <- chainCounts %>%
    dplyr::filter(nAlpha >= nAlphaCut |
      nBeta >= nBetaCut |
      (nAlpha >= alphaAndBetaCut & nBeta >= alphaAndBetaCut))

  tcrs$multiplet <- tcrs$barcode %in% multiplets$barcode

  if (callMAIT == TRUE) {
    tcrs <- tcrs %>%
      dplyr::mutate(isMAIT = (str_detect(v_gene, "TRAV1-2") &
        str_detect(j_gene, "TRAJ(33|12|20)")))
  }

  if (callINKT == TRUE) {
    tcrs <- tcrs %>%
      dplyr::mutate(isINKT = (str_detect(v_gene, "TRAV10") &
        str_detect(j_gene, "TRAJ18")))
  }
  return(tcrs)
}

# finds all pairs of alpha/beta TCRs
# TCRs is a dataframe that has a detected TCR chain for each row (alpha or beta or gamma, etc.) This is the output of mixcr, and maybe something analogous for 10x. But the columns might be named something different for 10x.
# adapted from A.Hu code
combinePairs <- function(tcrDf) {
  a <- tcrDf[tcrDf$chain %in% c("TRA", "a"), ]
  b <- tcrDf[tcrDf$chain %in% c("TRB", "b"), ]
  commonbarcodes <- intersect(a$barcode, b$barcode) # get barcodes that have both an alpha or beta chain
  a <- a[a$barcode %in% commonbarcodes, ]
  b <- b[b$barcode %in% commonbarcodes, ]

  # Compute how many alpha/beta pairs there should be
  atab <- table(a$barcode)
  btab <- table(b$barcode)
  nPairs <- sum(atab[commonbarcodes] * btab[commonbarcodes])

  # Initialize the pairs data frame
  pairs <- data.frame(CDR3b = rep("", nPairs),
    CDR3bnt = rep("", nPairs),
    TRBV = rep("", nPairs),
    TRBJ = rep("", nPairs),
    CDR3a = rep("", nPairs),
    CDR3ant = rep("", nPairs),
    TRAV = rep("", nPairs),
    TRAJ = rep("", nPairs),
    barcode = rep("", nPairs),
    fullLengthNTa = rep("", nPairs),
    fullLengthNTb = rep("", nPairs))

  # Iterate through the TCR alphas and betas and fill out the pair dataframe
  k <- 1
  for (barcode in commonbarcodes) {
    arows <- which(a$barcode == barcode)
    brows <- which(b$barcode == barcode)
    for (i in arows) {
      for (j in brows) {
        pairs[k, ] <- c(b$cdr3[j],
          b$cdr3_nt[j],
          b$v_gene[j],
          b$j_gene[j],
          a$cdr3[i],
          a$cdr3_nt[i],
          a$v_gene[i],
          a$j_gene[i],
          barcode,
          a$fullLengthNT[i],
          b$fullLengthNT[j])
        k <- k + 1
      }
    }
  }
  return(pairs)
}

# Create Sankey diagram showing hierarchical clustering splits
# With custom JS for full-lineage highlighting on hover
create_clustering_sankey <- function(cluster_tracking_df, plot_title = "Hierarchical Clustering",
                                     plotDir = NULL, filename = NULL,
                                     height = 8, width = 12) {
  # cluster_tracking_df should have columns: cell_id, and one column per nClusters value
  # e.g., clusters_n3, clusters_n4, clusters_n5, etc.

  if (!requireNamespace("networkD3", quietly = TRUE)) {
    warning("networkD3 package not available. Sankey plot will not be created.")
    return(NULL)
  }
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    warning("htmlwidgets package not available. Sankey plot will not be created.")
    return(NULL)
  }

  # Get cluster column names in order
  cluster_cols <- grep("^clusters_n", colnames(cluster_tracking_df), value = TRUE)
  if (length(cluster_cols) < 2) {
    warning("Need at least 2 clustering levels to create Sankey diagram")
    return(NULL)
  }

  # Sort by nClusters value
  cluster_nums <- as.integer(sub("clusters_n", "", cluster_cols))
  cluster_cols <- cluster_cols[order(cluster_nums)]

  # Build links dataframe for Sankey
  links_list <- list()

  for (i in 1:(length(cluster_cols) - 1)) {
    from_col <- cluster_cols[i]
    to_col <- cluster_cols[i + 1]

    # Create flow table
    flow_df <- cluster_tracking_df %>%
      dplyr::group_by(!!rlang::sym(from_col), !!rlang::sym(to_col)) %>%
      dplyr::summarise(value = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(
        source = paste0(from_col, "_", !!rlang::sym(from_col)),
        target = paste0(to_col, "_", !!rlang::sym(to_col))
      )

    links_list[[i]] <- flow_df %>%
      dplyr::select(source, target, value)
  }

  # Combine all links
  links <- dplyr::bind_rows(links_list)

  # Create nodes dataframe
  nodes <- data.frame(
    name = unique(c(links$source, links$target)),
    stringsAsFactors = FALSE
  )

  # Convert links to use node indices
  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1

  # Create Sankey diagram
  sankey_plot <- networkD3::sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "IDsource",
    Target = "IDtarget",
    Value = "value",
    NodeID = "name",
    fontSize = 12,
    nodeWidth = 30,
    nodePadding = 10
  )

  # Build cell-level lineage data for JavaScript

  # Each cell has a path through all clustering levels
  # Convert to JSON format: array of objects, each with the node names at each level
  cell_lineages <- lapply(1:nrow(cluster_tracking_df), function(i) {
    row_data <- cluster_tracking_df[i, cluster_cols, drop = FALSE]
    node_names <- sapply(cluster_cols, function(col) {
      paste0(col, "_", row_data[[col]])
    })
    as.list(node_names)
  })

  # Convert to JSON string
  cell_lineages_json <- jsonlite::toJSON(cell_lineages, auto_unbox = TRUE)

  # Custom JavaScript for full-lineage highlighting based on CELL-LEVEL data
  # When hovering over a node, find all cells that pass through that node,
  # then highlight only the other nodes those same cells pass through
  # AND rescale link widths to show lineage-specific cell counts
  full_lineage_js <- sprintf('
function(el, x) {
  setTimeout(function() {
    var svg = d3.select(el).select("svg");
    var link = svg.selectAll(".link");
    var node = svg.selectAll(".node");

    // Cell-level lineage data: each element is an object with node names for one cell
    var cellLineages = %s;

    // Build index: nodeName -> array of cell indices that pass through this node
    var nodeToCell = {};
    for (var cellIdx = 0; cellIdx < cellLineages.length; cellIdx++) {
      var cellPath = cellLineages[cellIdx];
      for (var key in cellPath) {
        var nodeName = cellPath[key];
        if (!nodeToCell[nodeName]) nodeToCell[nodeName] = [];
        nodeToCell[nodeName].push(cellIdx);
      }
    }

    // Store original link stroke-widths for restoration
    link.each(function(d) {
      d.originalStrokeWidth = d3.select(this).style("stroke-width");
      d.originalValue = d.value;
    });

    // Helper to make a link key from source->target names
    function linkKey(srcName, tgtName) {
      return srcName + "|||" + tgtName;
    }

    // Remove existing handlers and add new ones
    node.on("mouseover", null).on("mouseout", null);

    node.on("mouseover", function(arg1, arg2) {
      // D3 v4: function(d, i) - D3 v6: function(event, d)
      var nodeData;
      if (arg1 && typeof arg1.name === "string") {
        nodeData = arg1;
      } else if (arg2 && typeof arg2.name === "string") {
        nodeData = arg2;
      } else {
        nodeData = d3.select(this).datum();
      }

      if (!nodeData || !nodeData.name) {
        return;
      }

      var nodeName = nodeData.name;

      // Find all cells that pass through this node
      var cellsInNode = nodeToCell[nodeName] || [];
      var cellSet = {};
      for (var i = 0; i < cellsInNode.length; i++) {
        cellSet[cellsInNode[i]] = true;
      }

      // Find ALL nodes that these cells pass through AND count lineage-specific flows
      var lineageNodes = {};
      var lineageFlows = {};  // linkKey -> count of cells in lineage that use this link

      for (var i = 0; i < cellsInNode.length; i++) {
        var cellIdx = cellsInNode[i];
        var cellPath = cellLineages[cellIdx];

        // Sort keys by numeric value (clusters_n3, clusters_n4, ..., clusters_n17)
        // Extract the number from "clusters_nX" and sort numerically
        var keys = Object.keys(cellPath).sort(function(a, b) {
          var numA = parseInt(a.replace("clusters_n", ""), 10);
          var numB = parseInt(b.replace("clusters_n", ""), 10);
          return numA - numB;
        });

        // Mark all nodes in this cells path
        for (var j = 0; j < keys.length; j++) {
          lineageNodes[cellPath[keys[j]]] = true;
        }

        // Count flows between consecutive levels for this cell
        for (var j = 0; j < keys.length - 1; j++) {
          var srcNode = cellPath[keys[j]];
          var tgtNode = cellPath[keys[j + 1]];
          var lk = linkKey(srcNode, tgtNode);
          lineageFlows[lk] = (lineageFlows[lk] || 0) + 1;
        }
      }

      // Find max lineage flow for scaling
      var maxLineageFlow = 0;
      for (var lk in lineageFlows) {
        if (lineageFlows[lk] > maxLineageFlow) {
          maxLineageFlow = lineageFlows[lk];
        }
      }

      // Dim nodes not in lineage
      node.style("opacity", function(nd) {
        return (nd && nd.name && lineageNodes[nd.name]) ? 1.0 : 0.15;
      });

      // Update links: dim non-lineage links, rescale lineage links
      link.each(function(ld) {
        var srcName = ld.source.name;
        var tgtName = ld.target.name;
        var lk = linkKey(srcName, tgtName);
        var lineageCount = lineageFlows[lk] || 0;
        var inLineage = lineageCount > 0;

        var elem = d3.select(this);

        if (inLineage) {
          // Rescale stroke-width based on lineage-specific count
          // Use same scaling as original Sankey but with lineage count
          var originalWidth = parseFloat(ld.originalStrokeWidth) || ld.dy || 1;
          var originalValue = ld.originalValue || ld.value || 1;

          // Scale: newWidth = originalWidth * (lineageCount / originalValue)
          var scaledWidth = Math.max(1, originalWidth * (lineageCount / originalValue));

          elem.style("stroke-width", scaledWidth + "px");
          elem.style("stroke-opacity", 0.7);
        } else {
          elem.style("stroke-opacity", 0.02);
        }
      });
    });

    node.on("mouseout", function() {
      node.style("opacity", 1.0);
      // Restore original link widths and opacity
      link.each(function(ld) {
        var elem = d3.select(this);
        elem.style("stroke-width", ld.originalStrokeWidth);
        elem.style("stroke-opacity", 0.5);
      });
    });

  }, 100);
}
', cell_lineages_json)

  # Attach the custom JS to the widget
  sankey_plot <- htmlwidgets::onRender(sankey_plot, full_lineage_js)

  # Save if directory and filename provided
  if (!is.null(plotDir) && !is.null(filename)) {
    html_file <- file.path(plotDir, paste0(filename, ".html"))
    htmlwidgets::saveWidget(sankey_plot, html_file, selfcontained = TRUE)
    cat("Sankey plot saved to:", html_file, "\n")
  }

  return(sankey_plot)
}

# Create static Sankey diagram using ggsankey for clustering hierarchy
create_clustering_sankey_static <- function(cluster_tracking_df, plot_title = "Hierarchical Clustering",
                                            plotDir = NULL, filename = NULL,
                                            height = 8, width = 12) {
  # cluster_tracking_df should have columns: cell_id, and one column per nClusters value
  # e.g., clusters_n3, clusters_n4, clusters_n5, etc.

  if (!requireNamespace("ggsankey", quietly = TRUE)) {
    warning("ggsankey package not available. Installing...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("davidsjoberg/ggsankey")
  }
  library(ggsankey)

  # Get cluster column names in order
  cluster_cols <- grep("^clusters_n", colnames(cluster_tracking_df), value = TRUE)
  if (length(cluster_cols) < 2) {
    warning("Need at least 2 clustering levels to create Sankey diagram")
    return(NULL)
  }

  # Sort by nClusters value
  cluster_nums <- as.integer(sub("clusters_n", "", cluster_cols))
  cluster_cols <- cluster_cols[order(cluster_nums)]

  # Prepare data for ggsankey - convert to long format
  sankey_long <- cluster_tracking_df %>%
    ggsankey::make_long(!!!rlang::syms(cluster_cols))

  # Count cells at each node
  node_counts <- sankey_long %>%
    dplyr::filter(!is.na(node)) %>%
    dplyr::group_by(x, node) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      label = paste0(node, "\n(n=", format(count, big.mark = ","), ")")
    )

  # Join counts back to sankey data
  sankey_long <- sankey_long %>%
    dplyr::left_join(node_counts %>% dplyr::select(x, node, label), by = c("x", "node"))

  # Create color palette - use a colorblind-friendly palette
  unique_clusters <- unique(sankey_long$node)
  n_colors <- length(unique_clusters)
  cluster_colors <- setNames(
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(12, n_colors), "Set3"))(n_colors),
    unique_clusters
  )

  # Create axis labels from column names
  axis_labels <- paste0("n=", cluster_nums)

  # Create static Sankey plot
  p <- ggplot(sankey_long,
    aes(x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = factor(node),
      label = label)) +
    geom_sankey(flow.alpha = 0.5, node.color = "black", smooth = 6) +
    geom_sankey_label(size = 2.5, color = "black", fill = "white", alpha = 0.85) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    scale_x_discrete(labels = axis_labels) +
    labs(title = plot_title,
      subtitle = paste0("Showing how clusters split across ", length(cluster_cols), " granularity levels"),
      x = "Clustering Resolution (number of clusters)") +
    theme_void() +
    theme(
      axis.text.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
      plot.margin = margin(10, 10, 10, 10)
    )

  # Save if directory and filename provided
  if (!is.null(plotDir) && !is.null(filename)) {
    tryCatch(
      {
        savePlot(
          plot = p,
          plotDir = plotDir,
          filename = filename,
          height = height,
          width = width,
          units = "in",
          dpi = 600,
          formats = c("pdf", "png")
        )
        cat("Static Sankey plot saved to:", file.path(plotDir, paste0(filename, ".pdf/.png")), "\n")
      },
      error = function(e) {
        warning("Failed to save static Sankey plot: ", e$message)
      })
  }

  return(p)
}

# Extract UMAP coordinates and relevant metadata
create_faceted_umap <- function(seurat_obj,
                                reduction = "ref.umap_Treg",
                                group_by = "seurat_clusters_Treg",
                                facet_by = "stimulationFigures",
                                color_palette = palRNAClustersTreg,
                                pt_size = 1,
                                pt_alpha = 0.8) {
  # Extract UMAP coordinates
  umapCoords.tmp <- as.data.frame(Embeddings(seurat_obj, reduction = reduction))
  colnames(umapCoords.tmp) <- c("UMAP_1", "UMAP_2")

  # Add cell metadata
  umapCoords.tmp$cell_id <- rownames(umapCoords.tmp)
  metadata <- seurat_obj@meta.data
  metadata$cell_id <- rownames(metadata)

  # Combine coordinates and metadata
  plot_data <- merge(umapCoords.tmp, metadata, by = "cell_id")

  # Get the cluster column and convert to factor if needed
  plot_data[[group_by]] <- factor(plot_data[[group_by]])

  # Create the plot
  ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[group_by]])) +
    geom_point(size = pt_size, alpha = pt_alpha) +
    scale_color_manual(values = color_palette) +
    facet_wrap(reformulate(facet_by)) +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster", title = "") +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 8)
    )
}

# can be used for faceted expression UMAPs
create_feature_umap <- function(seurat_obj,
                                feature,
                                assay = "RNA", # or "FB"
                                reduction = "ref.umap",
                                facet_by = "stimulationFigures",
                                pt_size = 1,
                                pt_alpha = 0.8) {
  # Extract UMAP coordinates
  umapCoords.tmp <- as.data.frame(Embeddings(seurat_obj, reduction = reduction))
  colnames(umapCoords.tmp) <- c("UMAP_1", "UMAP_2")
  umapCoords.tmp$cell_id <- rownames(umapCoords.tmp)

  # Extract feature expression
  expr <- FetchData(seurat_obj, vars = feature, assay = assay)
  expr$cell_id <- rownames(expr)

  # Add metadata
  metadata <- seurat_obj@meta.data
  metadata$cell_id <- rownames(metadata)

  # Combine all
  plot_data <- merge(umapCoords.tmp, metadata, by = "cell_id")
  plot_data <- merge(plot_data, expr, by = "cell_id")

  # Plot
  ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[feature]])) +
    geom_point(size = pt_size, alpha = pt_alpha) +
    scale_color_viridis() +
    facet_wrap(reformulate(facet_by)) +
    labs(x = "UMAP 1", y = "UMAP 2", color = feature, title = paste(feature, "expression")) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 8)
    )
}
