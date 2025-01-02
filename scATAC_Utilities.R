# Load libraries
library(scuttle)

run_atac_qc_filtered <- function(sample, files_dir, annotations){
  # print sample name 
  message(paste0("running sample ", sample))
  # read in counts 
  counts <- Read10X_h5(filename = paste0(files_dir, sample, "/outs/filtered_feature_bc_matrix.h5"))
  # extract gene annotations from EnsDb
  annotations <- annotations
  # change to UCSC style
  #seqlevelsStyle(annotations) <- 'UCSC'
  # create chromatin assay 
  chrom_assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = paste0(files_dir, sample, '/outs/atac_fragments.tsv.gz'),
    annotation = annotations,
    min.cells = 10,
    min.features = 200)
  # read in metadata 
  metadata <- read.csv(
    file = paste0(files_dir, sample, "/outs/per_barcode_metrics.csv"),
    header = TRUE)
  rownames(metadata) <- metadata[,1]
  # subset for cell barcodes in chromatin assay 
  any(is.na(rownames(metadata)))
  any(is.na(colnames(chrom_assay)))
  print(ncol(chrom_assay))
  print(nrow(metadata))
  metadata_sub <- metadata[rownames(metadata) %in% colnames(chrom_assay),]
  # create seurat object 
  seurat <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks", 
    meta.data = metadata_sub)
  seurat$sample <- paste0(sample)
  # add the gene information to the object
  Annotation(seurat) <- annotations
  # calc. reads in peaks, and balcklist ratio
  seurat$pct_reads_in_peaks <- seurat$atac_peak_region_fragments / seurat$atac_fragments * 100
  # calc. perc of read in blacklist regions 
  overlaps <- findOverlaps(query = seurat[['peaks']], subject = blacklist_Drerio)
  hit.regions <- queryHits(x = overlaps)
  data.matrix <- GetAssayData(object = seurat, assay = 'peaks', slot = "counts")[hit.regions, , drop = FALSE]
  seurat$blacklist_region_fragments <- colSums(data.matrix)
  seurat$blacklist_fraction <- seurat$blacklist_region_fragments / seurat$atac_peak_region_fragments
  # compute nucleosome signal score per cell
  seurat <- NucleosomeSignal(object = seurat)
  seurat$nucleosome_group <- ifelse(seurat$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  # compute TSS enrichment score per cell
  seurat <- TSSEnrichment(object = seurat, fast = FALSE)
  seurat$high.tss <- ifelse(seurat$TSS.enrichment > 2, 'High', 'Low')
  # return seurat object 
  seurat
}

# Function to calculate outliers
find_outliers <- function(x) {
  print(paste0("Number of total cells: ", ncol(x)))
  outliers <- isOutlier(x$atac_peak_region_fragments)
  outliersDF <- as.data.frame(outliers)
  outliersDF$barcode <- rownames(outliersDF)
  outliersDFtrue <- outliersDF[outliersDF$outliers == TRUE,]$barcode
  outliersDFfalse <- outliersDF[outliersDF$outliers == FALSE,]$barcode
  print(paste0("Number of outliers: ", length(outliersDFtrue)))
  cellsKeep <- colnames(x)[colnames(x) %in% outliersDFfalse]
  
  # Subset the outliers
  xSub <- subset(x, cells = cellsKeep)
  pre <- ncol(x)
  post <- ncol(xSub)
  removed <- pre - post
  print(paste0(removed, " cells removed as outliers."))
  print(paste0("New total number of cells: ", ncol(xSub)))
  return(xSub)
}

# Function to filter cells by QC metrics
## "find_outlier()" function was replaced with "isOutlier()" from library(scater)
filter_cells_by_qc <- function(xSub){
  xSub2 <- subset(xSub,
                  subset = atac_peak_region_fragments > 3000 &
                     pct_reads_in_peaks > 15 &
                     blacklist_fraction < 0.03 &
                     nucleosome_signal < 4 &
                     TSS.enrichment > 2 &
                     TSS.enrichment < 15)
  print(paste0("Post filtering cell count: ", ncol(xSub2)))
  cat("\n")
  return(xSub2)
}

# Create custom fragment histograms
custom_FragmentHistogram <- function(seurat_obj, sample){
  p <- FragmentHistogram(object = seurat_obj) + ggtitle(sample)
  save_plot(paste0(outDir, "fragment_distribution-", sample, ".png"), p, base_height = 5, base_width = 5)
  return(p)
}

# Call peaks with MACS3
custom_call_peaks <- function(seurat_obj, macs2.path){
  # call peaks using MACS3
  peaks <- CallPeaks(seurat_obj, macs2.path = macsPath)
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  peaks
}

# define function to quantify counts for each object with a combined peak set 
quantify_counts <- function(sample, seurat_obj){
  #sample <- samples[[1]]
  #seurat_obj <- seurat_list_sub_2[[1]]
  message(sample)
  # read in metadata 
  metadata <- read.csv(
    file = paste0(sampleDir, sample, "/outs/per_barcode_metrics.csv"),
    header = TRUE,
    row.names = 1)
  # subset for cell barcodes in chromatin assay 
  metadata_sub <- metadata[rownames(metadata) %in% colnames(seurat_obj),]
  
  # create fragment object
  frags <- CreateFragmentObject(
    path = paste0(sampleDir, sample, "/outs/atac_fragments.tsv.gz"),
    cells = rownames(metadata_sub)
  )
  # quantify counts 
  counts <- FeatureMatrix(
    fragments = frags,
    features = combinedPeaks,
    cells = rownames(metadata_sub)
  )
  # create seurat objects 
  assay <- CreateChromatinAssay(counts = counts, fragments = frags, annotation = anno)
  seurat <- CreateSeuratObject(assay, assay = "ATAC", meta.data = metadata_sub)
  
  # add sample identifier 
  seurat$sample <- sample
  # dim. red. 
  seurat <- RunTFIDF(seurat)
  seurat <- FindTopFeatures(seurat, min.cutoff = 20)
  seurat <- RunSVD(seurat)
  seurat <- RunUMAP(seurat, dims = 2:50, reduction = 'lsi')
  
  # rename cells 
  seurat_obj <- RenameCells(seurat, new.names = paste0(colnames(seurat), "_", unique(seurat$sample)))
  
  # return seurat object 
  seurat
}

AddMotifs.custom <- function(
    object,
    genome,
    pfm,
    verbose = TRUE,
    ...
) {
  if (!requireNamespace("motifmatchr", quietly = TRUE)) {
    stop("Please install motifmatchr.\n",
         "https://www.bioconductor.org/packages/motifmatchr/")
  }
  if (is.null(x = names(x = pfm))) {
    warning("No 'names' attribute found in PFMatrixList. ",
            "Extracting names from individual entries.", immediate. = TRUE)
    names(x = pfm) <- vapply(
      X = pfm, FUN = slot, FUN.VALUE = "character", "name"
    )
  }
  if (verbose) {
    message("Building motif matrix")
  }
  motif.matrix <- CreateMotifMatrix(
    features = object,
    pwm = pfm,
    genome = genome,
    use.counts = FALSE
  )
  if (verbose) {
    message("Finding motif positions")
  }
  
  # for positions, a list of granges is returned
  # each element of list is a PFM name
  # each entry in granges is the position within a feature that matches motif
  obj_keep <- as.character(seqnames(x = object)) %in% seqlevels(x = genome)
  motif.positions <- motifmatchr::matchMotifs(
    pwms = pfm,
    subject = object[obj_keep],
    out = 'positions',
    genome = genome, p.cutoff = 5e-10
  )
  if (verbose) {
    message("Creating Motif object")
  }
  motif <- CreateMotifObject(
    data = motif.matrix,
    positions = motif.positions,
    pwm = pfm
  )
  return(motif)
}
