# Cell cycle scoring function
runCellCycle <- function(sample, x) {
  print(sample)
  x@meta.data$sample <- sample
  # Subset cell cycle genes to just the ones in our data
  sGenes <- drS[drS$ortholog_gene %in% rownames(x),]$ortholog_gene
  g2mGenes <- drGM[drGM$ortholog_gene %in% rownames(x),]$ortholog_gene
  
  # Set default assay as RNA and normalize
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  
  # Cell Cycle scoring
  x <- CellCycleScoring(x,
                        g2m.features = g2mGenes,
                        s.features = sGenes,
                        set.ident = TRUE)
  results <- table(x$sample, x$Phase)
  print(results)
  
  return(x)
}


# Function to run clustering on ATAC portion
clusteringATAC <- function(sample, x, resolution) {
  print(sample)
  # Run clustering based on the ATAC data
  x <- FindNeighbors(x,
                     reduction = "atac_lsi",
                     dims = 2:30,
                     graph.name = c("atac_nn", "atac_snn"),
                     verbose = TRUE)
  
  # Find clusters
  x <- FindClusters(x,
                    resolution = resolution,
                    graph.name = c("atac_snn"),
                    verbose = TRUE)
  
  # Assign clusters to the name ATAC_clusters
  x$ATAC_clusters <- x$seurat_clusters
  x$seurat_clusters <- NULL
  
  # Number of colors
  nColors <- length(unique(x$ATAC_clusters))
  
  # Plot
  p <- DimPlot(x, 
               reduction = "atac_umap", 
               group.by = "ATAC_clusters") +
    scale_color_manual(values = col[1:nColors])
  fileName <- paste(sample, resolution, "ATAC_Clusters.png", sep = "_")
  
  p1 <- FeaturePlot(x,
                    reduction = "atac_umap",
                    features = "ClamScore") +
    scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))
  
  p2 <- DimPlot(x, 
               reduction = "atac_umap", 
               group.by = "ClamClass") 

  p4 <- DimPlot(x,
                reduction = "atac_umap",
                group.by = "Phase") #+
  scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))
  
  full <- cowplot::plot_grid(p, p1, p2, p4)
  ggsave(paste0(figDir, "ATAC_Clustering/", fileName), full, width = 12, height = 10)
  
  return(x)
}

# Function to cluster based on the RNA
clusteringRNA <- function(sample, x, resolution) {
  print(sample)
  # Run clustering based on the ATAC data
  x <- FindNeighbors(x,
                     reduction = "rna_pca",
                     dims = 2:30,
                     graph.name = c("rna_nn", "rna_snn"),
                     verbose = TRUE)
  
  # Find clusters
  x <- FindClusters(x,
                    resolution = resolution,
                    graph.name = c("rna_snn"),
                    verbose = TRUE)
  
  # Assign clusters to the name ATAC_clusters
  x$RNA_clusters <- x$seurat_clusters
  x$seurat_clusters <- NULL
  
  # Number of colors
  nColors <- length(unique(x$RNA_clusters))
  
  
  
  return(x)
}



# Get WNN

run_wnn <- function(sample, x, resolution) {
  print(sample)
  x <- FindMultiModalNeighbors(
    x,
    reduction.list = list("rna_pca", "atac_lsi"),
    dims.list = list(1:50, 2:30),
    modality.weight.name = c("RNA.weight", "ATAC.weight"),
    verbose = TRUE)
  
  # Run UMAP
  x <- RunUMAP(x, 
               nn.name = "weighted.nn", 
               reduction.name = "wnn.umap")
  

  # Find clusters
  x <- FindClusters(x,
                    resolution = resolution,
                    graph.name = c("wsnn"),
                    verbose = TRUE)
  
  x$WNN_clusters <- x$seurat_clusters
  x$seurat_clusters <- NULL
  
  # Number of colors
  nColors <- length(unique(x$WNN_clusters))
  
  # Plot
  p <- DimPlot(x, 
               reduction = "wnn.umap", 
               group.by = "WNN_clusters") +
    scale_color_manual(values = col[1:nColors])
  fileName <- paste(sample, resolution, "WNN_Clusters.png", sep = "_")

  p1 <- FeaturePlot(x,
                    reduction = "wnn.umap",
                    features = "ClamScore") +
    scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))
  
  p2 <- DimPlot(x, 
                reduction = "wnn.umap", 
                group.by = "ClamClass") 
  
  p4 <- DimPlot(x,
                reduction = "wnn.umap",
                group.by = "Phase") #+
    #scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))
  
  full <- cowplot::plot_grid(p, p1, p4)
  ggsave(paste0(figDir, "WNN_Clustering/", fileName), full, width = 12, height = 10)
  
  return(x)
}



# Function to run scDblFinder
run_scDblFinder <- function(sample, x) {
  print(sample)
  # Convert Seurat to singleCellExperiment
  sce <- as.SingleCellExperiment(x)
  
  # Run scDblFinder
  sce <- scDblFinder(sce)
  
  # Extract results
  result <- as.data.frame(sce@colData@listData)
  x$scDoubletClass <- result$scDblFinder.class
  x$scDoubletScore <- result$scDblFinder.score
  
  # Set number of colors
  nColors <- length(unique(x$RNA_clusters))
  
  # Plot UMAPs
  p <- DimPlot(x, 
               reduction = "rna_umap", 
               group.by = "RNA_clusters") +
    scale_color_manual(values = col[1:nColors])
  fileName <- paste(sample, "RNA_Clusters.png", sep = "_")
  
  p2 <- DimPlot(x,
                reduction = "rna_umap",
                group.by = "scDoubletClass") 
  
  p1 <- FeaturePlot(x,
                    reduction = "rna_umap",
                    features = "scDoubletScore") +
    scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))
  
  p4 <- DimPlot(x,
                reduction = "rna_umap",
                group.by = "Phase") #+
  scale_colour_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")))
  
  full <- cowplot::plot_grid(p, p1, p2, p4)
  ggsave(paste0(figDir, "RNA_Clustering/", fileName), full, width = 12, height = 10)
  
  return(x)
  
}

# Function to compare the doublet calls between scDoubletFinder and Clamulet
compare_doublet_calls <- function(x, sample) {
  message(sample)
  # Extract the metadata
  meta <- x@meta.data
  # Extract vectors of barcodes for singlet/doublet calls for both methods
  sing <- rownames(meta[meta$scDoubletClass == "singlet",])
  doub <- rownames(meta[meta$scDoubletClass == "doublet",])
  message(paste0("Number of scDoubletFinder doublets: ", length(doub)))
  clamSing <- rownames(meta[meta$ClamClass == "clam_Singlet",])
  clamDoub <- rownames(meta[meta$ClamClass == "clam_Doublet",])
  message(paste0("Number of Clamulet doublets: ", length(clamDoub)))
  
  
  # Find the shared doublet calls
  shared <- intersect(doub, clamDoub)
  numShared <- length(shared)
  message(paste0("Number of shared doublet calls: ", numShared))
  
  
  # Save shared calls as a csv
  shared <- as.data.frame(shared)
  shared$Sample <- sample
  colnames(shared) <- c("Common_Doublets", "Sample")  
  
  venn.plot1 <- venn.diagram(
    x = list("scDF" = doub, "Clamulet" = clamDoub),
    filename = paste0(figDir, sample, "_Venn_Diagram.png"),
    output = TRUE,
    category.names = c("scDF", "Clamulet"),
    col = "transparent",
    fill = c("red", "blue"),
    alpha = 0.5,
    label.col = "white",
    cex = 3,
    cat.cex = 3,
    main = c("Doublet Calls"),
    disable.logging = TRUE,
    resolution = 150,
    imagetype = "png"
  )
  
}



