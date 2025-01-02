# Function to read in filtered 10X datasets for each sample
read_10x_multi <- function(sample, ident) {
  print(sample)
  counts <- Read10X_h5(filename = paste0(sampleDir, sample, "/outs/filtered_feature_bc_matrix.h5"))
  x <- CreateSeuratObject(counts = counts$`Gene Expression`,
                          min.cells = 10,
                          min.features = 100,
                          assay = "RNA")
  x[["ident"]] <- paste0(ident)
  return(x)
}

# Function to preform basic basic filtering, normalization, PCA and UMAP
preprocess_RNA_basic <- function(x, ident) {
  print(ident)
  xSub <- subset(x,
                 subset = nCount_RNA >= 300)
  
  # Add mito qc column
  xSub$percentMT <- PercentageFeatureSet(object = xSub,
                                         pattern = "^mt-")
  
  # Add log10 gene/UMI (complexity measure of dataset)
  xSub$log10GenesPerUMI <- log10(xSub$nFeature_RNA) / log10(xSub$nCount_RNA)
  
  return(xSub)
}

# Function to filter the scRNA data
filter_rna <- function(x, ident) {
  print(ident)
  xSub <- subset(x,
                 subset = nFeature_RNA > 1000 &
                   nFeature_RNA < 7500 &
                   percentMT < 20 &
                   log10GenesPerUMI > 0.70)
  
  # Find variable features
  xSub <- FindVariableFeatures(xSub, nfeatures = 2000)
  
  # Scale data
  xSub <- NormalizeData(xSub)
  xSub <- ScaleData(xSub)
  
  # Run PCA
  xSub <- RunPCA(xSub)
  
  # Run UMAP
  xSub <- RunUMAP(xSub, dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")
  
  return(xSub)
}

